#!/usr/bin/env Rscript
## iwi-per-area.R
##
## Population-weighted average IWI per administrative area and three-year
## period 1990-2019, from the satellite-based wealth maps of Africa by
## Pettersson et al. (2023):
##
##   "Temporal neighborhood-level material wealth maps of Africa (1990-2019)",
##   Harvard Dataverse, file wealth_map.tif (CC BY 4.0),
##   https://dataverse.harvard.edu/file.xhtml?fileId=7477600
##   10 bands: 1990-92, 1993-95, ..., 2017-19; IWI scaled 0-1; 6.72 km cells;
##   cells outside populated areas (GHSL) are 0 = no estimate.
##
## Population weights: GHS-POP R2023A, 30 arc seconds (~1 km), WGS84, epochs
## every 5 years, e.g.
##   https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/
##     GHS_POP_E2000_GLOBE_R2023A_4326_30ss/V1-0/GHS_POP_E2000_GLOBE_R2023A_4326_30ss_V1_0.zip
## Population for the middle year of each period is interpolated linearly
## between the two nearest epochs.
##
## Method: the IWI map is resampled (nearest neighbour) to the ~1 km population
## grid, so each 1 km cell carries the IWI of the 6.7 km cell it lies in. The
## IWI of an area is then
##     sum(pop * IWI * coverage) / sum(pop * coverage)
## over the 1 km cells with an IWI estimate, where coverage is the fraction of
## the cell inside the area (exactextractr). Only the population inside the
## area counts, which matters for areas smaller than a 6.7 km cell.
##
## Usage:
##   Rscript iwi-per-area.R --areas=forest-change.gpkg --only_regions_file=dhs.grid.RDS \
##       --wealth_file=~/data/wealth_map.tif --pop_dir=~/data/ghs_pop
##
## Output:
##   iwi_africa_long.csv   one row per area and period: iwi (0-100), pop,
##                         pop_share_with_iwi (share of the population living in
##                         cells with an IWI estimate)
##   iwi_africa_trend.csv  one row per area: IWI 1999-2001 and 2017-2019, OLS slopes
##                         (IWI units per year) for 2000-2018 and 1991-2018

this_script <- local({
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) normalizePath(sub("^--file=", "", a[1])) else
        if (!is.null(sys.frame(1)$ofile)) normalizePath(sys.frame(1)$ofile) else "iwi-per-area.R"
})
source(file.path(dirname(this_script), "transitions-per-area.R"))   # read_areas(), parse_args(), log_msg()

iwi_cfg <- modifyList(default_cfg, list(
    wealth_file = "wealth_map.tif",
    pop_dir     = "ghs_pop",
    pop_pattern = "GHS_POP_E([0-9]{4})_.*\\.tif$",
    first_year  = 1990L,          # first year of band 1
    period_len  = 3L,             # years per band
    workers     = 5L,             # periods processed in parallel (each needs a few GB of temporary disk)
    tmp_dir     = "iwi-tmp",
    out_file    = "iwi_africa_long.csv",
    trend_file  = "iwi_africa_trend.csv",
    log_file    = "iwi.log"
))

## Population epochs found in pop_dir: named vector year -> file
pop_epochs <- function(cfg) {
    f <- list.files(path.expand(cfg$pop_dir), pattern = cfg$pop_pattern, recursive = TRUE, full.names = TRUE)
    if (!length(f)) stop("No GHS-POP files matching ", cfg$pop_pattern, " in ", cfg$pop_dir,
                         " (unzip the downloads first)")
    yrs <- as.integer(sub(cfg$pop_pattern, "\\1", basename(f)))
    f <- f[order(yrs)]; yrs <- sort(yrs)
    setNames(f, yrs)
}

## Population raster for `year`, cropped to `e`, linear between the nearest epochs
pop_for_year <- function(year, epochs, e) {
    yrs <- as.integer(names(epochs))
    lo <- max(yrs[yrs <= year], min(yrs)); hi <- min(yrs[yrs >= year], max(yrs))
    p_lo <- crop(rast(epochs[[as.character(lo)]]), e)
    if (lo == hi) return(p_lo)
    p_hi <- crop(rast(epochs[[as.character(hi)]]), e)
    a <- (year - lo) / (hi - lo)
    (1 - a) * p_lo + a * p_hi
}

iwi_summary <- function(df) {
    w <- df$pop * df$coverage_fraction
    w[is.na(w)] <- 0
    has <- !is.na(df$iwi)
    c(iwi_x_pop = sum(df$iwi[has] * w[has]), pop_with_iwi = sum(w[has]), pop = sum(w))
}

one_period <- function(k, areas, epochs, cfg) {
    on.exit(cleanup_tmp(cfg), add = TRUE)
    t0 <- Sys.time()
    y_from <- cfg$first_year + cfg$period_len * (k - 1L)
    mid <- y_from + (cfg$period_len - 1L) / 2
    e <- ext(vect(areas)) + 0.1
    pop <- pop_for_year(mid, epochs, e)
    iwi <- crop(rast(path.expand(cfg$wealth_file))[[k]], e, snap = "out")
    iwi <- classify(iwi, cbind(0, NA)) * 100           # 0 = no estimate; scale to 0-100
    iwi <- resample(iwi, pop, method = "near")
    s <- c(iwi, pop); names(s) <- c("iwi", "pop")
    a <- st_cast(areas, "MULTIPOLYGON")
    v <- exact_extract(s, a, fun = iwi_summary, summarize_df = TRUE, progress = FALSE)
    m <- if (is.list(v)) do.call(rbind, v) else t(v)
    res <- data.table(region = areas$region, band = k,
                      period = sprintf("%d-%d", y_from, y_from + cfg$period_len - 1L), mid_year = mid,
                      iwi = m[, "iwi_x_pop"] / m[, "pop_with_iwi"],
                      pop = m[, "pop"], pop_share_with_iwi = m[, "pop_with_iwi"] / m[, "pop"])
    log_msg(cfg, "period ", res$period[1], " (population ", mid, "): ", nrow(res), " areas in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
    res
}

ols_slope <- function(x, y) {
    ok <- !is.na(y)
    if (sum(ok) < 2) return(NA_real_)
    x <- x[ok]; y <- y[ok]
    sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
}

run_iwi <- function(cfg = iwi_cfg) {
    t_start <- Sys.time()
    dir.create(path.expand(cfg$tmp_dir), showWarnings = FALSE, recursive = TRUE)
    cfg$tmp_dir <- normalizePath(path.expand(cfg$tmp_dir))
    terra::terraOptions(tempdir = cfg$tmp_dir)
    Sys.setenv(CPL_TMPDIR = cfg$tmp_dir)

    wealth <- rast(path.expand(cfg$wealth_file))
    n_bands <- nlyr(wealth)
    epochs <- pop_epochs(cfg)
    log_msg(cfg, "wealth map: ", n_bands, " bands; population epochs: ", paste(names(epochs), collapse = ", "))
    need <- cfg$first_year + cfg$period_len * (seq_len(n_bands) - 1L) + (cfg$period_len - 1L) / 2
    if (min(need) < min(as.integer(names(epochs))) || max(need) > max(as.integer(names(epochs))))
        log_msg(cfg, "WARNING: periods outside the population epochs use the nearest epoch")

    areas <- read_areas(cfg)
    ## keep areas that overlap the wealth map (Africa)
    wb <- st_as_sfc(st_bbox(as.vector(ext(wealth))[c(1, 3, 2, 4)] |> setNames(c("xmin", "ymin", "xmax", "ymax")),
                            crs = 4326))
    old <- suppressMessages(sf_use_s2(FALSE))
    inside <- lengths(suppressMessages(st_intersects(areas, wb))) > 0
    suppressMessages(sf_use_s2(old))
    areas <- areas[inside, ]
    log_msg(cfg, "areas overlapping the wealth map: ", nrow(areas))

    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    oopt <- options(future.globals.maxSize = +Inf); on.exit(options(oopt), add = TRUE)
    if (cfg$workers > 1) future::plan(future::multicore, workers = min(cfg$workers, n_bands)) else
        future::plan(future::sequential)

    res <- rbindlist(future_lapply(seq_len(n_bands), one_period, areas = areas, epochs = epochs,
                                   cfg = cfg, future.seed = TRUE, future.scheduling = Inf))
    setorder(res, region, band)
    fwrite(res, cfg$out_file)

    trend <- res[, list(
        iwi_1999_2001 = iwi[mid_year == 2000],
        iwi_2017_2019 = iwi[mid_year == 2018],
        iwi_slope_2000_2018 = ols_slope(mid_year[mid_year >= 2000], iwi[mid_year >= 2000]),
        iwi_slope_1991_2018 = ols_slope(mid_year, iwi),
        n_periods_with_iwi = sum(!is.na(iwi)),
        pop_2018 = pop[mid_year == 2018],
        min_pop_share_with_iwi = suppressWarnings(min(pop_share_with_iwi, na.rm = TRUE))
    ), by = region]
    fwrite(trend, cfg$trend_file)
    log_msg(cfg, "wrote ", cfg$out_file, " (", nrow(res), " rows) and ", cfg$trend_file, " (",
            nrow(trend), " areas, ", sum(trend$n_periods_with_iwi == 0), " without any IWI estimate) in ",
            round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
    invisible(list(long = res, trend = trend))
}

if (sys.nframe() == 0L) {
    cfg <- parse_args(commandArgs(trailingOnly = TRUE), iwi_cfg)
    if (Sys.getenv("TRANSITIONS_TMPDIR_SET") != "1") {
        dir.create(path.expand(cfg$tmp_dir), showWarnings = FALSE, recursive = TRUE)
        tmp <- normalizePath(path.expand(cfg$tmp_dir))
        status <- system2(file.path(R.home("bin"), "Rscript"),
                          c(shQuote(this_script), shQuote(commandArgs(trailingOnly = TRUE))),
                          env = c(paste0("TMPDIR=", shQuote(tmp)), "TRANSITIONS_TMPDIR_SET=1"))
        quit(save = "no", status = status)
    }
    run_iwi(cfg)
}
