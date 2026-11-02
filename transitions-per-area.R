#!/usr/bin/env Rscript
## transitions-per-area.R
##
## For every administrative area: what happened between year0 (2000) and
## year1 (2020) to each 30 m cell in the GLCLU data (Potapov et al. 2022)?
##
## Every cell is assigned a broad class in both years (forest, short
## vegetation, cropland, built-up, water, other, no data). The script counts
## the full 7 x 7 transition matrix per area and derives the shares behind
## figure 7 in the forest-change paper:
##
##   A  forest -> short vegetation
##   C  forest -> cropland
##   E  forest -> built-up
##   non-deforestation (everything else)
##
## Cells are weighted by the fraction of the cell that lies inside the area
## (exactextractr), exactly as in the published forest change data set, so
## forest_share_2000 and forest_change_2020 computed here should reproduce
## forestation_2000 and forestation_2020 in forest_change.csv. Use --check to
## verify that.
##
## Usage (from the shell):
##   Rscript transitions-per-area.R --areas=forest-change.gpkg \
##       --tiles_root="$HOME/mnt/windows/Deforestation and poverty - Forest coverage/5 degrees" \
##       --workers=14
##
##   Only the areas in a file, e.g. the DHS areas (names of an sp object):
##   Rscript transitions-per-area.R --areas=forest-change.gpkg --only_regions_file=dhs.grid.RDS
##
##   Quick test on one top-level region (Krâchéh in Cambodia):
##   Rscript transitions-per-area.R --areas=forest-change.gpkg --only_regions='^116\.Krâchéh'
##
##   Regenerate figure 7 (bug-fixed, works across tile borders):
##   Rscript transitions-per-area.R --areas=forest-change.gpkg \
##       --plot_region='116.Krâchéh.Chhloung.Kampong Damrei' --plot_file=fig7.png
##
## From R: source("transitions-per-area.R") defines the functions without
## running anything, so they can be used from an Rnw file, e.g.
##   res <- run_transitions(cfg = modifyList(default_cfg, list(areas = "...")))
##
## Expected layout of tiles_root (same as in deforestation.Rnw):
##   <tiles_root>/2000/15N_105E.tif, <tiles_root>/2020/15N_105E.tif, ...
## i.e. 5 x 5 degree tiles produced by split-tiffs.sh.

suppressPackageStartupMessages({
    library(sf)
    library(terra)
    library(exactextractr)
    library(data.table)
    library(future.apply)
})

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------

default_cfg <- list(
    areas        = "forest-change.gpkg",   # .gpkg/.shp (st_read) or .RDS with an sf object
    region_col   = "region",
    tiles_root   = "~/mnt/windows/Deforestation and poverty - Forest coverage/5 degrees",
    year0        = 2000L,
    year1        = 2020L,
    tile_deg     = 5,
    max_area_km2 = 1000,                   # split larger areas into pieces (as in deforestation.Rnw)
    chunk_size   = 25L,                    # max pieces per job; tiles are split into several jobs
    workers      = max(1L, parallelly::availableCores() - 2L),
    cache_dir    = "transitions-cache",    # one RDS per job; delete to force recomputation
    tmp_dir      = "transitions-tmp",      # temporary files, removed after each job; --tmp_dir=/tmp/...
                                           # is faster if /tmp is a RAM disk with enough room
    out_file     = "forest_transitions_2000_2020.csv",
    only_regions = "",                     # regex on region id; "" = all areas
    only_regions_file = "",                # file listing the areas to include (see region_ids_from_file)
    check        = "",                     # path to forest_change.csv for validation; "" = skip
    plot_region  = "",                     # region id to plot (figure 7); "" = no plot
    plot_file    = "transitions.png",
    log_file     = "transitions.log"
)

## ---------------------------------------------------------------------------
## Classes
## ---------------------------------------------------------------------------

## Class index is 0-based so that the transition code is from * 7 + to.
class_names <- c("forest", "shortveg", "cropland", "builtup", "water", "other", "nodata")
n_cls <- length(class_names)

## Lookup table, GLCLU map value (0-255) -> class index. Same code lists as
## in deforestation.Rnw.
glclu_lut <- local({
    lut <- rep(5L, 256L)                          # other: snow/ice (241) and unused codes
    lut[c(27:48, 127:148) + 1L] <- 0L             # forest: trees >= 5 m, terra firma + wetland
    lut[c(0:26, 100:126) + 1L]  <- 1L             # short vegetation (incl. trees < 5 m)
    lut[244L + 1L]              <- 2L             # cropland
    lut[250L + 1L]              <- 3L             # built-up
    lut[c(200:207, 254) + 1L]   <- 4L             # open surface water + ocean
    lut[255L + 1L]              <- 6L             # no data
    lut
})

to_class <- function(v) {
    out <- glclu_lut[as.integer(v) + 1L]          # NA and values > 255 give NA
    out[is.na(out)] <- 6L
    out
}

transition_names <- as.vector(t(outer(class_names, class_names, paste, sep = "__")))
## transition_names[k + 1] is the name of code k = from * 7 + to

## exactextractr summary function. With summarize_df = TRUE it gets a
## data.frame with one column per raster layer (y0, y1) and coverage_fraction.
## Returns the coverage-weighted cell count for each of the 49 transitions.
transition_summary <- function(df) {
    k <- to_class(df$y0) * n_cls + to_class(df$y1)
    w <- df$coverage_fraction
    full <- w >= 1
    out <- as.numeric(tabulate(k[full] + 1L, nbins = n_cls * n_cls))
    if (any(!full)) {
        s <- rowsum(w[!full], k[!full])
        idx <- as.integer(rownames(s)) + 1L
        out[idx] <- out[idx] + s[, 1]
    }
    out
}

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

log_msg <- function(cfg, ...) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
    message(line)
    cat(line, "\n", file = path.expand(cfg$log_file), append = TRUE, sep = "")
}

## Tile names (e.g. "15N_105E.tif") overlapping a bbox. Same naming scheme as
## glclu_tile_names_from_bbox() in deforestation.Rnw: name = north-west corner.
tile_names_from_bbox <- function(bb, tile_deg = 5, eps = 1e-10) {
    xmax2 <- max(bb[["xmin"]], bb[["xmax"]] - eps)
    ymax2 <- max(bb[["ymin"]], bb[["ymax"]] - eps)
    xs <- seq(floor(bb[["xmin"]] / tile_deg) * tile_deg, floor(xmax2 / tile_deg) * tile_deg, by = tile_deg)
    ys <- seq(floor(bb[["ymin"]] / tile_deg) * tile_deg, floor(ymax2 / tile_deg) * tile_deg, by = tile_deg)
    g <- expand.grid(x = xs, y = ys)
    lat <- as.integer(round(g$y + tile_deg))
    lon <- as.integer(round(g$x))
    paste0(sprintf("%02d%s", abs(lat), ifelse(lat >= 0, "N", "S")), "_",
           sprintf("%03d%s", abs(lon), ifelse(lon >= 0, "E", "W")), ".tif")
}

tile_bbox_sfc <- function(tile, tile_deg = 5) {
    p <- strsplit(sub("\\.tif$", "", tile), "_")[[1]]
    lat <- as.integer(sub("[NS]$", "", p[1])) * ifelse(grepl("S$", p[1]), -1L, 1L)
    lon <- as.integer(sub("[EW]$", "", p[2])) * ifelse(grepl("W$", p[2]), -1L, 1L)
    st_as_sfc(st_bbox(c(xmin = lon, xmax = lon + tile_deg,
                        ymin = lat - tile_deg, ymax = lat), crs = 4326))
}

tile_files <- function(cfg, tile) {
    c(y0 = file.path(cfg$tiles_root, cfg$year0, tile),
      y1 = file.path(cfg$tiles_root, cfg$year1, tile))
}

## Region ids from a file: .RDS with a character vector, a named object
## (e.g. sp::SpatialPolygons, whose names are the region ids), or a
## data.frame/sf with a "region" column; otherwise a text file, one id per line.
region_ids_from_file <- function(path) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        x <- readRDS(path)
        ids <- if (is.character(x)) x
               else if (is.data.frame(x) && "region" %in% names(x)) as.character(x$region)
               else if (inherits(x, "SpatialPolygons")) vapply(x@polygons, function(p) p@ID, "")
               else names(x)
    } else {
        ids <- readLines(path, encoding = "UTF-8")
    }
    if (is.null(ids) || !length(ids)) stop("Found no region ids in ", path)
    unique(enc2utf8(trimws(ids[nzchar(ids)])))
}

read_areas <- function(cfg) {
    a <- if (grepl("\\.rds$", cfg$areas, ignore.case = TRUE)) readRDS(cfg$areas) else st_read(cfg$areas, quiet = TRUE)
    a <- st_as_sf(a)
    a <- a[, cfg$region_col]
    names(a)[1] <- "region"
    a <- st_set_geometry(a, st_geometry(a))   # normalise geometry column name
    if (is.na(st_crs(a))) st_crs(a) <- 4326
    if (st_crs(a) != st_crs(4326)) a <- st_transform(a, 4326)
    a$region <- as.character(a$region)
    if (nzchar(cfg$only_regions_file)) {
        ids <- region_ids_from_file(path.expand(cfg$only_regions_file))
        a$region <- enc2utf8(a$region)
        hit <- ids %in% a$region
        message(sprintf("%s: %d region ids, %d found among the areas", cfg$only_regions_file, length(ids), sum(hit)))
        if (any(!hit)) message("  not found (first 10): ", paste(head(ids[!hit], 10), collapse = "; "))
        a <- a[a$region %in% ids, ]
    }
    if (nzchar(cfg$only_regions)) a <- a[grepl(cfg$only_regions, a$region), ]
    if (nrow(a) == 0) stop("No areas selected (check --areas and --only_regions)")
    a
}

## Split large areas into pieces of about max_area_km2 (as split_large_polygons()
## in deforestation.Rnw) so that exactextractr never holds huge windows in memory.
split_areas <- function(areas, max_area_km2) {
    old <- sf_use_s2(); on.exit(suppressMessages(sf_use_s2(old)), add = TRUE); suppressMessages(sf_use_s2(FALSE))
    geoms <- st_geometry(areas)
    areas_m2 <- suppressMessages(as.numeric(st_area(geoms)))
    big_idx <- which(areas_m2 > max_area_km2 * 1e6)

    split_one <- function(i) {
        suppressMessages(sf_use_s2(FALSE))
        one <- geoms[i]
        bb <- st_bbox(one)
        n <- ceiling(areas_m2[i] / (max_area_km2 * 1e6))
        cs <- sqrt(as.numeric((bb$xmax - bb$xmin) * (bb$ymax - bb$ymin)) / n)
        grid <- st_make_grid(one, cellsize = cs, square = TRUE)
        g <- suppressWarnings(suppressMessages(st_intersection(one, grid)))
        ## keep polygonal parts only (intersections can produce slivers as lines/points)
        if (any(st_geometry_type(g) == "GEOMETRYCOLLECTION"))
            g <- suppressWarnings(st_collection_extract(g, "POLYGON"))
        g <- g[st_geometry_type(g) %in% c("POLYGON", "MULTIPOLYGON") & !st_is_empty(g)]
        if (length(g)) st_cast(g, "MULTIPOLYGON") else st_cast(one, "MULTIPOLYGON")
    }
    ## runs in parallel if a future plan has been set (see run_transitions)
    big_pieces <- future_lapply(big_idx, split_one, future.seed = TRUE,
                                future.chunk.size = 10L)

    small <- setdiff(seq_along(geoms), big_idx)
    out <- st_sf(region = c(areas$region[small], rep(areas$region[big_idx], lengths(big_pieces))),
                 geometry = c(st_cast(geoms[small], "MULTIPOLYGON"), do.call(c, big_pieces)))
    out$piece_id <- seq_len(nrow(out))
    out
}

## Which pieces intersect which tile.
map_tiles <- function(pieces, tile_deg) {
    old <- sf_use_s2(); on.exit(suppressMessages(sf_use_s2(old)), add = TRUE); suppressMessages(sf_use_s2(FALSE))
    tiles <- unique(unlist(lapply(st_geometry(pieces), function(g) tile_names_from_bbox(st_bbox(g), tile_deg))))
    tiles_sf <- st_sf(tile = tiles, geometry = do.call(c, lapply(tiles, tile_bbox_sfc, tile_deg = tile_deg)))
    hits <- suppressMessages(st_intersects(tiles_sf, pieces))
    names(hits) <- tiles
    hits[lengths(hits) > 0]
}

## ---------------------------------------------------------------------------
## Work for one job: a chunk of the pieces that intersect one tile
## ---------------------------------------------------------------------------

## Split every tile's pieces into chunks of at most chunk_size pieces, and
## make the chunks small enough that there are at least as many jobs as
## workers (so that even a single tile uses all cores).
make_jobs <- function(tile_to_piece, chunk_size, workers) {
    total <- sum(lengths(tile_to_piece))
    size <- max(1L, min(as.integer(chunk_size), floor(total / workers)))
    jobs <- list()
    for (tile in names(tile_to_piece)) {
        idx <- as.integer(tile_to_piece[[tile]])
        n_chunks <- ceiling(length(idx) / size)
        grp <- rep(seq_len(n_chunks), each = size, length.out = length(idx))
        for (ch in seq_len(n_chunks))
            jobs[[length(jobs) + 1L]] <- list(tile = tile, idx = idx[grp == ch], chunk = ch, n_chunks = n_chunks)
    }
    jobs
}

## Remove temporary files written by this process (terra and raster put the
## process id in their temporary file names).
cleanup_tmp <- function(cfg) {
    dirs <- unique(c(tempdir(), cfg$tmp_dir))
    f <- list.files(dirs, recursive = TRUE, full.names = TRUE)
    f <- f[grepl(paste0("_", Sys.getpid(), "[_.]"), basename(f))]
    if (length(f)) unlink(f, force = TRUE)
    invisible(length(f))
}

tmp_size_gb <- function(cfg) {
    f <- list.files(unique(c(tempdir(), cfg$tmp_dir)), recursive = TRUE, full.names = TRUE)
    sum(file.size(f), na.rm = TRUE) / 2^30
}

process_job <- function(job, pieces, cfg) {
    on.exit(cleanup_tmp(cfg), add = TRUE)
    tile <- job$tile
    idx <- job$idx
    sub_p <- pieces[idx, ]
    cache <- file.path(cfg$cache_dir, sprintf("%s_%s_%s_%04d.RDS", sub("\\.tif$", "", tile),
                                              cfg$year0, cfg$year1, job$chunk))

    ## Cached result is valid if it covers all requested pieces (keyed on region
    ## + piece geometry hash so that it survives changes in --only_regions).
    piece_keys <- paste(sub_p$region, vapply(st_geometry(sub_p), digest_geom, ""), sep = "|")
    if (file.exists(cache)) {
        old <- tryCatch(readRDS(cache), error = function(e) NULL)
        if (!is.null(old) && all(piece_keys %in% old$key)) {
            res <- old[match(piece_keys, old$key)]
            res[, piece_id := sub_p$piece_id]
            return(res)
        }
    }

    f <- tile_files(cfg, tile)
    if (!all(file.exists(f))) {
        log_msg(cfg, "MISSING tile ", tile, " (", sum(!file.exists(f)), " of 2 years missing); ",
                length(idx), " pieces affected")
        return(data.table(piece_id = sub_p$piece_id, missing_tile = tile))
    }

    t0 <- Sys.time()
    r <- rast(unname(f))
    names(r) <- c("y0", "y1")

    if (any(st_geometry_type(sub_p) == "POLYGON")) sub_p <- st_cast(sub_p, "MULTIPOLYGON")

    vals <- exact_extract(r, sub_p, fun = transition_summary, summarize_df = TRUE, progress = FALSE)
    m <- if (is.list(vals)) do.call(rbind, vals) else t(vals)
    res <- as.data.table(m)
    setnames(res, transition_names)
    res[, key := piece_keys]

    dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(cache, ".", Sys.getpid(), ".tmp")
    saveRDS(res, tmp); file.rename(tmp, cache)

    res[, piece_id := sub_p$piece_id]
    log_msg(cfg, "tile ", tile, " chunk ", job$chunk, "/", job$n_chunks, ": ", nrow(sub_p), " pieces, ",
            format(round(sum(m) / 1e6, 1), nsmall = 1), " M cells in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min",
            { g <- tmp_size_gb(cfg); if (g > 0.1) sprintf(" [temporary files: %.1f GB]", g) else "" })
    res
}

## Small stable hash of a geometry (no extra package needed).
digest_geom <- function(g) {
    x <- unlist(g)
    sprintf("%d_%.8f_%.8f_%.8f", length(x), sum(x), sum(x * seq_along(x)) / length(x), x[1])
}

## ---------------------------------------------------------------------------
## Derived measures
## ---------------------------------------------------------------------------

add_derived <- function(dt) {
    col <- function(from, to) {
        nm <- as.vector(outer(from, to, paste, sep = "__"))
        rowSums(dt[, ..nm])
    }
    land <- c("forest", "shortveg", "cropland", "builtup", "other")   # not water, not no data

    dt[, valid_2000  := col(land, class_names)]
    dt[, valid_2020  := col(class_names, land)]
    dt[, forest_2000 := col("forest", class_names)]
    dt[, forest_2020 := col(class_names, "forest")]

    ## Same definitions as forest_change.csv (per cent and percentage points)
    dt[, forest_share_2000  := 100 * forest_2000 / valid_2000]
    dt[, forest_change_2020 := 100 * (forest_2020 / valid_2020 - forest_2000 / valid_2000)]

    ## The four colours of figure 7, as per cent of the land area in 2000.
    ## pct_forest_to_other (forest -> water/other) is shown separately so that
    ## the five shares sum to 100.
    dt[, pct_A_forest_to_shortveg := 100 * forest__shortveg / valid_2000]
    dt[, pct_C_forest_to_cropland := 100 * forest__cropland / valid_2000]
    dt[, pct_E_forest_to_builtup  := 100 * forest__builtup  / valid_2000]
    dt[, pct_forest_to_other      := 100 * col("forest", c("water", "other")) / valid_2000]
    dt[, pct_non_deforestation    := 100 - pct_A_forest_to_shortveg - pct_C_forest_to_cropland -
                                         pct_E_forest_to_builtup - pct_forest_to_other]

    ## What replaced the forest: per cent of the year-2000 forest (cf. table 5)
    dt[, pctforest_A := 100 * forest__shortveg / forest_2000]
    dt[, pctforest_C := 100 * forest__cropland / forest_2000]
    dt[, pctforest_E := 100 * forest__builtup  / forest_2000]
    dt[, pctforest_stable := 100 * forest__forest / forest_2000]

    ## Forest gain paths, per cent of land area in 2000
    dt[, pct_B_shortveg_to_forest := 100 * shortveg__forest / valid_2000]
    dt[, pct_D_cropland_to_forest := 100 * cropland__forest / valid_2000]
    dt[, pct_F_builtup_to_forest  := 100 * builtup__forest  / valid_2000]

    ## Initial land use mix (context for "what can replace the forest")
    dt[, pct_cropland_2000 := 100 * col("cropland", class_names) / valid_2000]
    dt[, pct_cropland_2020 := 100 * col(class_names, "cropland") / valid_2020]
    dt[]
}

## ---------------------------------------------------------------------------
## Main computation
## ---------------------------------------------------------------------------

run_transitions <- function(cfg = default_cfg) {
    t_start <- Sys.time()
    cfg$tiles_root <- path.expand(cfg$tiles_root)
    for (yr in c(cfg$year0, cfg$year1)) {
        d <- file.path(cfg$tiles_root, yr)
        if (!dir.exists(d)) stop("Tile directory not found: ", d, " (check --tiles_root)")
    }
    dir.create(path.expand(cfg$tmp_dir), showWarnings = FALSE, recursive = TRUE)
    cfg$tmp_dir <- normalizePath(path.expand(cfg$tmp_dir))
    terra::terraOptions(tempdir = cfg$tmp_dir)
    Sys.setenv(CPL_TMPDIR = cfg$tmp_dir)          # GDAL's own temporary files
    log_msg(cfg, "temporary files in ", cfg$tmp_dir, " and ", tempdir())
    areas <- read_areas(cfg)
    log_msg(cfg, "areas: ", nrow(areas))

    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    ## multicore forks the R process, so large globals (areas, pieces) are shared,
    ## not copied; the default 500 MiB export check does not apply.
    oopt <- options(future.globals.maxSize = +Inf); on.exit(options(oopt), add = TRUE)
    if (cfg$workers > 1) future::plan(future::multicore, workers = cfg$workers) else future::plan(future::sequential)
    ## Avoid oversubscription: each worker uses one GDAL thread
    terra::setGDALconfig("GDAL_NUM_THREADS", "1")

    n_big <- sum(suppressMessages(as.numeric(sf::st_area(sf::st_geometry(areas)))) > cfg$max_area_km2 * 1e6, na.rm = TRUE)
    log_msg(cfg, "splitting ", n_big, " areas larger than ", cfg$max_area_km2, " km2 into pieces (",
            cfg$workers, " workers) ...")
    pieces <- split_areas(areas, cfg$max_area_km2)
    log_msg(cfg, "finding tiles for ", nrow(pieces), " pieces ...")
    tile_to_piece <- map_tiles(pieces, cfg$tile_deg)
    log_msg(cfg, "pieces: ", nrow(pieces), ", tiles: ", length(tile_to_piece), ", workers: ", cfg$workers)

    jobs <- make_jobs(tile_to_piece, cfg$chunk_size, cfg$workers)
    log_msg(cfg, "jobs: ", length(jobs))

    res <- future_lapply(jobs, process_job, pieces = pieces,
                         cfg = cfg, future.seed = TRUE, future.scheduling = Inf)
    res <- rbindlist(res, fill = TRUE)
    ## make sure all transition columns exist, even if every tile was missing
    for (nm in setdiff(c(transition_names, "piece_id"), names(res))) res[, (nm) := numeric(0)]

    ## Pieces in a missing tile -> the area is incomplete
    missing <- if ("missing_tile" %in% names(res)) res[!is.na(missing_tile)] else res[0]
    res <- if ("missing_tile" %in% names(res)) res[is.na(missing_tile)] else res
    res[, region := pieces$region[piece_id]]

    out <- res[, lapply(.SD, sum), by = region, .SDcols = transition_names]
    out[, complete := !region %in% pieces$region[missing$piece_id]]
    ## Areas with no tile at all (should not happen) are kept with NA counts
    out <- merge(data.table(region = unique(areas$region)), out, by = "region", all.x = TRUE)
    out[is.na(complete), complete := FALSE]
    out <- add_derived(out)
    ## add_derived() labels the two years 2000 and 2020; use the actual years
    if (cfg$year0 != 2000L || cfg$year1 != 2020L)
        setnames(out, gsub("2020", "YEAR1", gsub("2000", "YEAR0", names(out))) |>
                      sub(pattern = "YEAR0", replacement = cfg$year0) |>
                      sub(pattern = "YEAR1", replacement = cfg$year1))

    fwrite(out, cfg$out_file)
    log_msg(cfg, "wrote ", nrow(out), " areas to ", cfg$out_file, " (",
            sum(!out$complete, na.rm = TRUE), " incomplete) in ",
            round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
    invisible(out)
}

## Compare with the published forest change data set.
check_against <- function(out, check_file, year0 = 2000L, year1 = 2020L) {
    ref <- fread(check_file, encoding = "UTF-8")
    setnames(ref, gsub("\\.", "_", names(ref)))
    ## local forest_change.csv: forestation_change_YYYY; Zenodo: forestation_YYYY (YYYY > 2000)
    old <- paste0("forestation_change_", year1); new <- paste0("forestation_", year1)
    if (!new %in% names(ref) && old %in% names(ref)) setnames(ref, old, new)
    cat(sprintf("Checking against %s\n", check_file))
    if (year0 == 2000L) {
        m <- merge(out[, .(region, a = get(paste0("forest_share_", year0)))],
                   ref[, .(region, b = forestation_2000)], by = "region")
        if (max(m$b, na.rm = TRUE) <= 1.0001) m[, b := 100 * b]
        cat(sprintf("Matched %d of %d areas\n", sum(complete.cases(m)), nrow(out)))
        cat(sprintf("forest share %d: max |diff| = %.4f pp, cor = %.6f\n", year0,
                    max(abs(m$a - m$b), na.rm = TRUE), cor(m$a, m$b, use = "complete")))
        if (new %in% names(ref)) {
            m <- merge(out[, .(region, a = get(paste0("forest_change_", year1)))],
                       ref[, .(region, b = get(new))], by = "region")
            cat(sprintf("forest change %d-%d: max |diff| = %.4f pp, cor = %.6f\n", year0, year1,
                        max(abs(m$a - m$b), na.rm = TRUE), cor(m$a, m$b, use = "complete")))
        } else cat("No column for ", year1, " in ", check_file, "\n")
    } else cat("Check only implemented for year0 = 2000\n")
    invisible(NULL)
}

## ---------------------------------------------------------------------------
## Figure 7, cell by cell (fixed version of the chunk focus-deforestation)
## ---------------------------------------------------------------------------

transition_raster <- function(region_id, cfg = default_cfg) {
    a <- read_areas(modifyList(cfg, list(only_regions = "", only_regions_file = "")))
    poly <- a[a$region == region_id, ]
    if (nrow(poly) != 1) stop("region not found: ", region_id)
    tiles <- tile_names_from_bbox(st_bbox(poly), cfg$tile_deg)
    cfg$tiles_root <- path.expand(cfg$tiles_root)
    f0 <- file.path(cfg$tiles_root, cfg$year0, tiles)
    f1 <- file.path(cfg$tiles_root, cfg$year1, tiles)
    v <- vect(poly)
    ## vrt() mosaics adjacent tiles, so areas crossing a tile border are complete
    r0 <- crop(if (length(f0) > 1) vrt(f0) else rast(f0), v)
    r1 <- crop(if (length(f1) > 1) vrt(f1) else rast(f1), v)
    c0 <- to_class(values(r0, mat = FALSE))
    c1 <- to_class(values(r1, mat = FALSE))
    ## 1 = A, 2 = C, 3 = E, 4 = everything else
    k <- rep(4L, length(c0))
    k[c0 == 0L & c1 == 1L] <- 1L
    k[c0 == 0L & c1 == 2L] <- 2L
    k[c0 == 0L & c1 == 3L] <- 3L
    out <- rast(r0); values(out) <- k
    out <- mask(out, v)
    levels(out) <- data.frame(id = 1:4, class = c("A forest → short vegetation",
                                                  "C forest → cropland",
                                                  "E forest → built-up",
                                                  "non-deforestation"))
    coltab(out) <- data.frame(value = 1:4, col = c(rgb(46, 168, 153, maxColorValue = 255),
                                                   rgb(230, 126, 34, maxColorValue = 255),
                                                   rgb(155, 89, 182, maxColorValue = 255),
                                                   "#F2F2F2"))
    attr(out, "polygon") <- poly
    out
}

plot_transitions <- function(region_id, cfg = default_cfg, file = NULL) {
    r <- transition_raster(region_id, cfg)
    if (!is.null(file)) png(file, width = 1800, height = 1100, res = 200)
    plot(r, mar = c(2, 2, 1, 12))
    plot(st_geometry(attr(r, "polygon")), add = TRUE, col = NA)
    if (!is.null(file)) dev.off()
    invisible(r)
}

## ---------------------------------------------------------------------------
## Command line
## ---------------------------------------------------------------------------

parse_args <- function(args, cfg = default_cfg) {
    for (a in args) {
        m <- regmatches(a, regexec("^--([A-Za-z0-9_]+)=(.*)$", a))[[1]]
        if (length(m) != 3) stop("Cannot parse argument: ", a, " (use --key=value)")
        key <- m[2]; val <- m[3]
        if (validUTF8(val)) Encoding(val) <- "UTF-8"   # region names with non-ASCII letters
        if (!key %in% names(cfg)) stop("Unknown option --", key, ". Options: ", paste(names(cfg), collapse = ", "))
        cfg[[key]] <- if (is.numeric(cfg[[key]])) as(val, class(cfg[[key]])) else val
    }
    cfg
}

if (sys.nframe() == 0L) {
    cfg <- parse_args(commandArgs(trailingOnly = TRUE))
    ## R's own tempdir() is fixed at start-up and defaults to /tmp, which can be
    ## a small RAM disk. Restart once with TMPDIR pointing to tmp_dir on disk.
    if (Sys.getenv("TRANSITIONS_TMPDIR_SET") != "1") {
        dir.create(path.expand(cfg$tmp_dir), showWarnings = FALSE, recursive = TRUE)
        tmp <- normalizePath(path.expand(cfg$tmp_dir))
        script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
        status <- system2(file.path(R.home("bin"), "Rscript"),
                          c(shQuote(script), shQuote(commandArgs(trailingOnly = TRUE))),
                          env = c(paste0("TMPDIR=", shQuote(tmp)), "TRANSITIONS_TMPDIR_SET=1"))
        quit(save = "no", status = status)
    }
    if (nzchar(cfg$plot_region)) {
        plot_transitions(cfg$plot_region, cfg, file = cfg$plot_file)
        message("wrote ", cfg$plot_file)
    } else {
        out <- run_transitions(cfg)
        if (nzchar(cfg$check)) check_against(out, cfg$check, cfg$year0, cfg$year1)
    }
}
