#!/usr/bin/env Rscript
## drivers-per-area.R
##
## For every administrative area, split the forest loss 2000-2020 in the GLCLU
## data (Potapov et al. 2022) by
##
##  1. the dominant driver of forest loss in the WRI / Google DeepMind map
##     "Global drivers of forest loss at 1 km resolution" v1.2 (2001-2024),
##     https://zenodo.org/records/15366671 (drivers_forest_loss_1km_2001_2024_v1_2.tif)
##  2. oil palm planted in (year0, year1], from "Global oil palm extent and
##     planting year 1990-2021" (Descals et al. 2024),
##     https://zenodo.org/records/13379129 (the planting-year tiles, 30 m)
##
## Every 30 m GLCLU cell that was forest in year0 and is not forest in year1
## gets the driver of the 1 km cell it lies in. Oil palm is counted from the
## planting-year layer, whatever GLCLU says in year1: young palms are short
## vegetation, older palms (>= 5 m) are "forest" in GLCLU.
##
## Cells are weighted by the fraction that lies inside the area (exactextractr),
## as in transitions-per-area.R, whose functions this script reuses.
##
## Usage:
##   ## 1. look at the input files first (class codes, CRS, year coding)
##   Rscript drivers-per-area.R --inspect=1 \
##       --drivers_file=~/data/drivers_forest_loss_1km_2001_2024_v1_2.tif \
##       --palm_dir=~/data/oilpalm_planting_year
##
##   ## 2. run, e.g. for the DHS areas
##   Rscript drivers-per-area.R --areas=forest-change.gpkg \
##       --tiles_root=~/annex/projekt/deforestation/Forest\ Coverage/5\ degrees/ \
##       --only_regions_file=dhs.grid.RDS \
##       --drivers_file=... --palm_dir=...
##
## Output (default forest_drivers_2000_2020.csv), one row per area, all in
## coverage-weighted 30 m cell counts and as per cent of the land area in year0:
##   loss_total, loss_<driver>        forest loss year0 -> year1 by driver
##   palm_from_forest                 forest in year0, oil palm planted after year0
##   palm_from_forest_<class>         ... by GLCLU class in year1 (forest/shortveg/cropland/other)
##   palm_from_nonforest              not forest in year0, oil palm planted after year0
##   palm_before                      oil palm planted in or before year0
##   loss_<path>_to_<class>           forest loss cross-classified, with no overlap: <path> is
##                                    "oil_palm" if oil palm was planted after year0, otherwise
##                                    the driver ("none" if no driver); <class> is the GLCLU class
##                                    in year1 (shortveg, cropland, builtup, other = water/other).
##                                    These sum to loss_total.
##   palm_from_forest_<type>_<period> oil palm planted on forest land by plantation type in the 2021
##                                    extent map (industrial, smallholder; untyped = planting year
##                                    but not oil palm in the extent map) and planting period
##                                    (early = year0+1..year0+10, late = the rest up to year1);
##                                    needs --palm_type_dir. These sum to palm_from_forest.
##   palm_from_nonforest_<type>       the same for land that was not forest in year0
##   pct_*                            the same, per cent of the land area in year0

## ---------------------------------------------------------------------------
## Reuse transitions-per-area.R (same directory as this script)
## ---------------------------------------------------------------------------

this_script <- local({
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) normalizePath(sub("^--file=", "", a[1])) else
        if (!is.null(sys.frame(1)$ofile)) normalizePath(sys.frame(1)$ofile) else "drivers-per-area.R"
})
source(file.path(dirname(this_script), "transitions-per-area.R"))

drivers_cfg <- modifyList(default_cfg, list(
    drivers_file     = "drivers_forest_loss_1km_2001_2024_v1_2.tif",
    drivers_band     = 1L,           # band 1 = most likely driver class
    ## labels for driver codes 1, 2, ...; check with --inspect=1 and the dataset README
    driver_labels    = "permanent_agriculture,hard_commodities,shifting_cultivation,logging,wildfire,settlements_infrastructure,other_natural_disturbances",
    palm_dir         = "",           # directory with the oil palm planting-year GeoTIFFs; "" = skip oil palm
    palm_pattern     = "\\.tif$",
    palm_year_offset = 0L,           # added to the raster value to get the calendar year (see --inspect)
    palm_type_dir    = "",           # directory with the oil palm extent 2021 GeoTIFFs (0 other,
                                     # 1 industrial, 2 smallholder); "" = no split by type
    palm_type_pattern = "\\.tif$",
    early_years      = 10L,          # planted in year0+1 .. year0+early_years counts as early
    gdal_cache_mb    = 256L,         # GDAL block cache per worker (the default is 5% of RAM *per process*)
    mem_share        = 0.6,          # share of RAM that all workers together may use for terra
    inspect          = 0L,
    cache_dir        = "drivers-cache-v3",   # v3: adds oil palm by type and planting period
    out_file         = "forest_drivers_2000_2020.csv",
    log_file         = "drivers.log"
))

driver_labels_of <- function(cfg) trimws(strsplit(cfg$driver_labels, ",")[[1]])

## ---------------------------------------------------------------------------
## Oil palm tile index (extent of every tile in EPSG:4326)
## ---------------------------------------------------------------------------

palm_index <- function(cfg, dir = cfg$palm_dir, pattern = cfg$palm_pattern, idx_name = "palm_index.RDS") {
    if (!nzchar(dir)) return(NULL)
    files <- list.files(path.expand(dir), pattern = pattern,
                        recursive = TRUE, full.names = TRUE)
    if (!length(files)) stop("No oil palm files matching ", pattern, " in ", dir)
    idx_file <- file.path(cfg$cache_dir, idx_name)
    if (file.exists(idx_file)) {
        idx <- readRDS(idx_file)
        if (setequal(idx$file, files)) return(idx)
    }
    ext4326 <- lapply(files, function(f) {
        r <- rast(f)
        p <- as.polygons(ext(r), crs = crs(r))
        as.vector(ext(project(p, "EPSG:4326")))
    })
    idx <- data.table(file = files, do.call(rbind, ext4326))
    dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(idx, idx_file)
    idx
}

palm_type_index <- function(cfg)
    palm_index(cfg, cfg$palm_type_dir, cfg$palm_type_pattern, "palm_type_index.RDS")

## ---------------------------------------------------------------------------
## Inspect the input files
## ---------------------------------------------------------------------------

inspect_inputs <- function(cfg) {
    cat("== Drivers map:", cfg$drivers_file, "\n")
    d <- rast(path.expand(cfg$drivers_file))
    print(d)
    s <- spatSample(d[[cfg$drivers_band]], 2e5, method = "regular", na.rm = FALSE, warn = FALSE)[[1]]
    cat("\nValue frequencies in band", cfg$drivers_band, "(regular sample of 200,000 cells):\n")
    print(table(s, useNA = "ifany"))
    lab <- driver_labels_of(cfg)
    cat("\nLabels the script will use: ", paste(seq_along(lab), lab, sep = " = ", collapse = "; "),
        "\nAny other value (incl. 0 and NA) is counted as 'none'.\n")

    if (nzchar(cfg$palm_dir)) {
        idx <- palm_index(cfg)
        cat("\n== Oil palm planting year:", nrow(idx), "files in", cfg$palm_dir, "\n")
        r <- rast(idx$file[1])
        print(r)
        cat("CRS of the files:\n"); print(table(vapply(idx$file, function(f) crs(rast(f), describe = TRUE)$code, "")))
        v <- values(r, mat = FALSE); v <- v[!is.na(v)]
        cat("\nValue frequencies in", basename(idx$file[1]), "(non-NA cells):\n"); print(table(v))
        cat("\nWith --palm_year_offset=", cfg$palm_year_offset, " these become the years ",
            paste(range(v[v > 0]) + cfg$palm_year_offset, collapse = "-"),
            ". Values <= 0 are treated as 'no oil palm'.\n", sep = "")
    }
    if (nzchar(cfg$palm_type_dir)) {
        tidx <- palm_type_index(cfg)
        cat("\n== Oil palm extent 2021:", nrow(tidx), "files in", cfg$palm_type_dir, "\n")
        r <- rast(tidx$file[1]); print(r)
        cat("\nValue frequencies in", basename(tidx$file[1]), "(expected 0 other, 1 industrial, 2 smallholder):\n")
        print(table(values(r, mat = FALSE), useNA = "ifany"))
    }
    invisible(NULL)
}

## ---------------------------------------------------------------------------
## Per piece
## ---------------------------------------------------------------------------

n_driver_slots <- function(cfg) length(driver_labels_of(cfg)) + 1L   # slot 1 = none

dest_classes <- c("shortveg", "cropland", "builtup", "other")

cross_names <- function(cfg) {
    paths <- c("none", driver_labels_of(cfg), "oil_palm")
    as.vector(t(outer(paths, dest_classes, function(p, d) paste0("loss_", p, "_to_", d))))
}

palm_types <- c("untyped", "industrial", "smallholder")      # extent map codes 0 (or missing), 1, 2

palm_type_names <- function() {
    c(as.vector(t(outer(palm_types, c("early", "late"), function(t, p) paste0("palm_from_forest_", t, "_", p)))),
      paste0("palm_from_nonforest_", palm_types))
}

summary_names <- function(cfg) {
    c("valid_2000", "loss_total", paste0("loss_", c("none", driver_labels_of(cfg))),
      "palm_from_forest", paste0("palm_from_forest_", c("forest", "shortveg", "cropland", "other")),
      "palm_from_nonforest", "palm_before", cross_names(cfg), palm_type_names())
}

make_summary_fun <- function(cfg) {
    nd <- n_driver_slots(cfg)
    y0 <- cfg$year0; y1 <- cfg$year1; off <- cfg$palm_year_offset
    function(df) {
        w  <- df$coverage_fraction
        c0 <- to_class(df$y0); c1 <- to_class(df$y1)
        land0  <- c0 <= 3L | c0 == 5L                      # not water (4), not no data (6)
        forest0 <- c0 == 0L
        loss <- forest0 & c1 != 0L & c1 != 6L
        d <- as.integer(df$drv); d[is.na(d) | d < 1L | d >= nd] <- 0L
        loss_by_d <- vapply(0:(nd - 1L), function(k) sum(w[loss & d == k]), 0)

        py <- df$palm + off
        has_palm <- !is.na(df$palm) & df$palm > 0
        planted  <- has_palm & py > y0 & py <= y1
        pf <- forest0 & planted
        c1g <- ifelse(c1 <= 2L, c1, 3L)                    # forest, shortveg, cropland, other
        ## cross-classified loss: path (driver 0..nd-1, oil palm = nd) x class in year1
        path <- ifelse(pf, nd, d)
        dest <- ifelse(c1 == 1L, 0L, ifelse(c1 == 2L, 1L, ifelse(c1 == 3L, 2L, 3L)))
        k <- path * 4L + dest
        cross <- numeric((nd + 1L) * 4L)
        if (any(loss)) {
            s <- rowsum(w[loss], k[loss])
            cross[as.integer(rownames(s)) + 1L] <- s[, 1]
        }
        ## oil palm by type (extent 2021) and planting period
        pt <- as.integer(df$ptype); pt[is.na(pt) | pt < 1L | pt > 2L] <- 0L
        early <- py <= y0 + cfg$early_years
        pnf <- !forest0 & land0 & planted
        ## order: untyped early, untyped late, industrial early, ... (as palm_type_names())
        by_type <- c(as.vector(vapply(0:2, function(k)
                         c(sum(w[pf & pt == k & early]), sum(w[pf & pt == k & !early])), c(0, 0))),
                     vapply(0:2, function(k) sum(w[pnf & pt == k]), 0))
        c(sum(w[land0]), sum(w[loss]), loss_by_d,
          sum(w[pf]), vapply(0:3, function(k) sum(w[pf & c1g == k]), 0),
          sum(w[pnf]),
          sum(w[has_palm & py <= y0 & land0]),
          cross, by_type)
    }
}

## Oil palm planting year on the grid of `template` (NA where no tile)
palm_on_grid <- function(template, bb, pidx) {
    if (is.null(pidx)) return(NULL)
    hit <- pidx[xmin < bb[["xmax"]] & xmax > bb[["xmin"]] & ymin < bb[["ymax"]] & ymax > bb[["ymin"]]]
    if (!nrow(hit)) return(NULL)
    out <- NULL
    for (f in hit$file) {
        p <- tryCatch(project(rast(f), template, method = "near"), error = function(e) NULL)
        if (is.null(p)) next
        out <- if (is.null(out)) p else cover(out, p)
    }
    out
}

process_piece <- function(g, files, drv, pidx, sfun, cfg) {
    nm <- summary_names(cfg)
    zero <- setNames(numeric(length(nm)), nm)
    v <- vect(g)
    r0 <- tryCatch(crop(rast(files[["y0"]]), v, snap = "out"), error = function(e) NULL)
    if (is.null(r0)) return(zero)                          # piece does not overlap this tile
    r1 <- crop(rast(files[["y1"]]), ext(r0))
    d  <- resample(drv, r0, method = "near")
    p  <- palm_on_grid(r0, as.vector(ext(r0)), pidx$year)
    if (is.null(p)) { p <- rast(r0); values(p) <- NA_integer_ }
    pt <- palm_on_grid(r0, as.vector(ext(r0)), pidx$type)
    if (is.null(pt)) { pt <- rast(r0); values(pt) <- NA_integer_ }
    s <- c(r0, r1, d, p, pt)
    names(s) <- c("y0", "y1", "drv", "palm", "ptype")
    res <- exact_extract(s, st_sf(geometry = st_sfc(g, crs = 4326)), fun = sfun,
                         summarize_df = TRUE, progress = FALSE)
    m <- if (is.list(res)) res[[1]] else as.numeric(res)
    setNames(as.numeric(m), nm)
}

process_job_drivers <- function(job, pieces, pidx, cfg) {
    on.exit(cleanup_tmp(cfg), add = TRUE)
    tile <- job$tile
    sub_p <- pieces[job$idx, ]
    cache <- file.path(cfg$cache_dir, sprintf("%s_%s_%s_%04d.RDS", sub("\\.tif$", "", tile),
                                              cfg$year0, cfg$year1, job$chunk))
    piece_keys <- paste(sub_p$region, vapply(st_geometry(sub_p), digest_geom, ""), sep = "|")
    if (file.exists(cache)) {
        old <- tryCatch(readRDS(cache), error = function(e) NULL)
        if (!is.null(old) && all(piece_keys %in% old$key)) {
            res <- old[match(piece_keys, old$key)]
            res[, piece_id := sub_p$piece_id]
            return(res)
        }
    }
    files <- tile_files(cfg, tile)
    if (!all(file.exists(files))) {
        log_msg(cfg, "MISSING tile ", tile, "; ", nrow(sub_p), " pieces affected")
        return(data.table(piece_id = sub_p$piece_id, missing_tile = tile))
    }
    t0 <- Sys.time()
    drv <- rast(path.expand(cfg$drivers_file))[[cfg$drivers_band]]
    sfun <- make_summary_fun(cfg)
    geoms <- st_geometry(sub_p)
    m <- do.call(rbind, lapply(seq_along(geoms), function(i)
        process_piece(geoms[[i]], files, drv, pidx, sfun, cfg)))
    res <- as.data.table(m)
    res[, key := piece_keys]
    dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(cache, ".", Sys.getpid(), ".tmp")
    saveRDS(res, tmp); file.rename(tmp, cache)
    res[, piece_id := sub_p$piece_id]
    log_msg(cfg, "tile ", tile, " chunk ", job$chunk, "/", job$n_chunks, ": ", nrow(sub_p), " pieces, ",
            format(round(sum(res$valid_2000) / 1e6, 1), nsmall = 1), " M land cells in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min",
            { g <- tmp_size_gb(cfg); if (g > 0.1) sprintf(" [temporary files: %.1f GB]", g) else "" })
    res
}

## ---------------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------------

run_drivers <- function(cfg = drivers_cfg) {
    t_start <- Sys.time()
    cfg$tiles_root <- path.expand(cfg$tiles_root)
    cfg$drivers_file <- path.expand(cfg$drivers_file)
    if (!file.exists(cfg$drivers_file)) stop("Drivers file not found: ", cfg$drivers_file)
    for (yr in c(cfg$year0, cfg$year1)) {
        d <- file.path(cfg$tiles_root, yr)
        if (!dir.exists(d)) stop("Tile directory not found: ", d, " (check --tiles_root)")
    }
    dir.create(path.expand(cfg$tmp_dir), showWarnings = FALSE, recursive = TRUE)
    cfg$tmp_dir <- normalizePath(path.expand(cfg$tmp_dir))
    terra::terraOptions(tempdir = cfg$tmp_dir)
    Sys.setenv(CPL_TMPDIR = cfg$tmp_dir)

    pidx <- list(year = palm_index(cfg), type = palm_type_index(cfg))
    if (is.null(pidx$year)) log_msg(cfg, "no --palm_dir given: oil palm columns will be 0")
    else log_msg(cfg, "oil palm planting year: ", nrow(pidx$year), " tiles")
    if (is.null(pidx$type)) log_msg(cfg, "no --palm_type_dir given: all oil palm is 'untyped'")
    else log_msg(cfg, "oil palm extent (type): ", nrow(pidx$type), " tiles")

    areas <- read_areas(cfg)
    log_msg(cfg, "areas: ", nrow(areas))

    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    oopt <- options(future.globals.maxSize = +Inf); on.exit(options(oopt), add = TRUE)
    ## Memory limits per worker, set before the workers are forked. By default every process may
    ## use 5% of RAM for GDAL's block cache and 60% of RAM for terra, which with many workers adds
    ## up to more than the machine has.
    terra::setGDALconfig("GDAL_CACHEMAX", as.character(cfg$gdal_cache_mb))
    ram.gb <- tryCatch(as.numeric(sub("\\D+(\\d+).*", "\\1",
                           grep("^MemTotal", readLines("/proc/meminfo"), value = TRUE))) / 2^20,
                       error = function(e) NA)
    if (!is.na(ram.gb)) terra::terraOptions(memmax = cfg$mem_share * ram.gb / max(1, cfg$workers))
    log_msg(cfg, "memory: ", if (is.na(ram.gb)) "unknown RAM" else sprintf("%.0f GB RAM", ram.gb),
            ", GDAL cache ", cfg$gdal_cache_mb, " MB and terra max ",
            if (is.na(ram.gb)) "default" else sprintf("%.1f GB", cfg$mem_share * ram.gb / max(1, cfg$workers)),
            " per worker")
    if (cfg$workers > 1) future::plan(future::multicore, workers = cfg$workers) else future::plan(future::sequential)
    terra::setGDALconfig("GDAL_NUM_THREADS", "1")

    pieces <- split_areas(areas, cfg$max_area_km2)
    tile_to_piece <- map_tiles(pieces, cfg$tile_deg)
    jobs <- make_jobs(tile_to_piece, cfg$chunk_size, cfg$workers)
    log_msg(cfg, "pieces: ", nrow(pieces), ", tiles: ", length(tile_to_piece), ", jobs: ", length(jobs),
            ", workers: ", cfg$workers)

    res <- future_lapply(jobs, process_job_drivers, pieces = pieces, pidx = pidx, cfg = cfg,
                         future.seed = TRUE, future.scheduling = Inf)
    res <- rbindlist(res, fill = TRUE)
    nm <- summary_names(cfg)
    for (x in setdiff(c(nm, "piece_id"), names(res))) res[, (x) := numeric(0)]
    missing <- if ("missing_tile" %in% names(res)) res[!is.na(missing_tile)] else res[0]
    res <- if ("missing_tile" %in% names(res)) res[is.na(missing_tile)] else res
    res[, region := pieces$region[piece_id]]

    out <- res[, lapply(.SD, sum), by = region, .SDcols = nm]
    out[, complete := !region %in% pieces$region[missing$piece_id]]
    out <- merge(data.table(region = unique(areas$region)), out, by = "region", all.x = TRUE)
    out[is.na(complete), complete := FALSE]

    ## per cent of the land area in year0
    for (x in setdiff(nm, "valid_2000")) out[, (paste0("pct_", x)) := 100 * get(x) / valid_2000]
    ## label the years
    if (cfg$year0 != 2000L) setnames(out, "valid_2000", paste0("valid_", cfg$year0))

    fwrite(out, cfg$out_file)
    log_msg(cfg, "wrote ", nrow(out), " areas to ", cfg$out_file, " (",
            sum(!out$complete), " incomplete) in ",
            round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
    invisible(out)
}

if (sys.nframe() == 0L) {
    cfg <- parse_args(commandArgs(trailingOnly = TRUE), drivers_cfg)
    if (cfg$inspect == 1L) {
        inspect_inputs(cfg)
        quit(save = "no", status = 0)
    }
    ## restart once with TMPDIR on disk (see transitions-per-area.R)
    if (Sys.getenv("TRANSITIONS_TMPDIR_SET") != "1") {
        dir.create(path.expand(cfg$tmp_dir), showWarnings = FALSE, recursive = TRUE)
        tmp <- normalizePath(path.expand(cfg$tmp_dir))
        status <- system2(file.path(R.home("bin"), "Rscript"),
                          c(shQuote(this_script), shQuote(commandArgs(trailingOnly = TRUE))),
                          env = c(paste0("TMPDIR=", shQuote(tmp)), "TRANSITIONS_TMPDIR_SET=1"))
        quit(save = "no", status = status)
    }
    run_drivers(cfg)
}
