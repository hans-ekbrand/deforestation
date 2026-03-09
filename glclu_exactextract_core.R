# Core workflow for GLCLU annual maps using exactextractr
# Forest share = forest pixel equivalents / valid land pixel equivalents
# Valid land excludes open surface water (200:207) and ocean (254)

library(sf)
library(raster)
library(exactextractr)
library(data.table)
library(future.apply)

forest_codes <- c(27:48, 127:148)   # >= 5m trees on terra firma + wetlands
water_codes  <- c(200:207, 254)     # open surface water + ocean

# Return GLCLU tile names overlapping an sf/sp bbox.
glclu_tile_names_from_bbox <- function(bbox, tile_deg = 10, eps = 1e-10) {
  # bbox can be st_bbox or sp::bbox; coerce to named vector with xmin/xmax/ymin/ymax
  if (is.matrix(bbox)) {
    xmin <- bbox["x", "min"]
    xmax <- bbox["x", "max"]
    ymin <- bbox["y", "min"]
    ymax <- bbox["y", "max"]
  } else {
    xmin <- bbox[["xmin"]]
    xmax <- bbox[["xmax"]]
    ymin <- bbox[["ymin"]]
    ymax <- bbox[["ymax"]]
  }

  xmax2 <- xmax - eps
  ymax2 <- ymax - eps
  if (xmax2 < xmin) xmax2 <- xmin
  if (ymax2 < ymin) ymax2 <- ymin

  x_starts <- seq(floor(xmin / tile_deg) * tile_deg,
                  floor(xmax2 / tile_deg) * tile_deg,
                  by = tile_deg)
  y_starts <- seq(floor(ymin / tile_deg) * tile_deg,
                  floor(ymax2 / tile_deg) * tile_deg,
                  by = tile_deg)

  tiles <- expand.grid(xmin_tile = x_starts, ymin_tile = y_starts)

  apply(tiles, 1, function(t) {
    lon_west  <- as.integer(t[["xmin_tile"]])
    lat_north <- as.integer(t[["ymin_tile"]] + tile_deg)

    lat_lab <- sprintf("%02d%s", abs(lat_north), ifelse(lat_north >= 0, "N", "S"))
    lon_lab <- sprintf("%03d%s", abs(lon_west),  ifelse(lon_west  >= 0, "E", "W"))
    paste0(lat_lab, "_", lon_lab, ".tif")
  })
}

ensure_glclu_tile <- function(year, tile, base_dir, overwrite = FALSE) {
  dir.create(file.path(base_dir, as.character(year)), recursive = TRUE, showWarnings = FALSE)
  
  local_file <- file.path(base_dir, as.character(year), paste0(tile, ".tif"))
  
  if (!file.exists(local_file) || overwrite) {
    url <- sprintf(
      "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/%s/%s",
      year, tile
    )
    
    tmp <- paste0(local_file, ".download")
    utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
    file.rename(tmp, local_file)
  }
  
  local_file
}

bboxes <- lapply(st_geometry(my.sf), st_bbox)
needed.tiles <- unique(unlist(sapply(bboxes, glclu_tile_names_from_bbox)))

files_2000 <- vapply(
  needed.tiles,
  ensure_glclu_tile,
  character(1),
  year = 2000,
  base_dir = data_dir
)

# bbox polygon for a tile name like 10N_000E.tif
st_bbox_from_tile_name <- function(tile_name, tile_deg = 10) {
  nm <- sub("\\.tif$", "", basename(tile_name))
  parts <- strsplit(nm, "_")[[1]]
  lat_north <- as.integer(sub("([0-9]+)[NS]", "\\1", parts[1]))
  if (grepl("S$", parts[1])) lat_north <- -lat_north
  lon_west <- as.integer(sub("([0-9]+)[EW]", "\\1", parts[2]))
  if (grepl("W$", parts[2])) lon_west <- -lon_west

  st_bbox(c(
    xmin = lon_west,
    xmax = lon_west + tile_deg,
    ymin = lat_north - tile_deg,
    ymax = lat_north
  ), crs = 4326)
}

# One exactextractr summary: returns forest and valid pixel equivalents per polygon.
forest_valid_summary_fun <- function(values, coverage_fraction) {
  is_water  <- values %in% water_codes
  is_forest <- values %in% forest_codes

  valid  <- sum(coverage_fraction[!is_water], na.rm = TRUE)
  forest <- sum(coverage_fraction[is_forest], na.rm = TRUE)
  c(forest = forest, valid = valid)
}

# Process one tile for one year over only intersecting polygons.
process_one_tile_one_year <- function(tile_file, polys_sf, poly_id_col = "poly_id") {
  bb <- st_bbox_from_tile_name(tile_file)
  tile_poly <- st_as_sfc(bb)
  hits <- lengths(st_intersects(polys_sf, tile_poly)) > 0
  if (!any(hits)) return(NULL)

  sub_polys <- polys_sf[hits, ]
  r <- raster(tile_file)

  vals <- exact_extract(r, sub_polys, forest_valid_summary_fun,
                        progress = FALSE)
  vals <- as.data.table(vals)
  vals[, (poly_id_col) := sub_polys[[poly_id_col]]]
  vals[, tile_file := basename(tile_file)]
  vals[]
}

# Process all relevant tiles for one year in parallel.
# yearly_dir is a directory like /path/to/GLCLU/2000/
compute_forest_share_for_year <- function(polys_sf,
                                          yearly_dir,
                                          poly_id_col = "poly_id",
                                          workers = max(1L, future::availableCores() - 1L),
                                          relevant_tile_names = NULL) {
  stopifnot(inherits(polys_sf, "sf"))
  if (is.null(polys_sf[[poly_id_col]])) {
    stop("poly_id_col must exist in polys_sf")
  }

  if (is.null(relevant_tile_names)) {
    # derive from full bbox of all polygons
    relevant_tile_names <- unique(glclu_tile_names_from_bbox(st_bbox(polys_sf)))
  }

  tile_files <- file.path(yearly_dir, relevant_tile_names)
  tile_files <- tile_files[file.exists(tile_files)]
  if (!length(tile_files)) stop("No matching tile files found in ", yearly_dir)

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  res_list <- future_lapply(tile_files, function(f) {
    process_one_tile_one_year(tile_file = f, polys_sf = polys_sf, poly_id_col = poly_id_col)
  }, future.seed = TRUE)

  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  if (!nrow(res)) {
    out <- data.table(poly_tmp = polys_sf[[poly_id_col]], forest = 0, valid = 0, share = NA_real_)
    data.table::setnames(out, "poly_tmp", poly_id_col)
    return(out)
  }

  out <- res[, .(forest = sum(forest, na.rm = TRUE),
                 valid  = sum(valid, na.rm = TRUE)), by = poly_id_col]
  out[, share := fifelse(valid > 0, forest / valid, NA_real_)]
  setorder(out, poly_id_col)
  out[]
}

ensure_glclu_tile <- function(year, tile, base_dir, overwrite = FALSE) {
  dir.create(file.path(base_dir, as.character(year)), recursive = TRUE, showWarnings = FALSE)
  
  local_file <- file.path(base_dir, as.character(year), paste0(tile, ".tif"))
  
  if (!file.exists(local_file) || overwrite) {
    url <- sprintf(
      "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/%s/%s.tif",
      year, tile
    )
    
    tmp <- paste0(local_file, ".download")
    utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
    file.rename(tmp, local_file)
  }
  
  local_file
}

# Convenience wrapper for two years and net change in percentage points.
compute_forest_share_change <- function(polys_sf,
                                        root_dir,
                                        year0,
                                        year1,
                                        poly_id_col = "poly_id",
                                        workers = max(1L, future::availableCores() - 1L),
                                        relevant_tile_names = NULL) {
  if (is.null(relevant_tile_names)) {
    # safer than just using the global bbox: union of per-polygon tile sets could be used,
    # but for one-country/admin run the global bbox is usually fine.
    relevant_tile_names <- unique(glclu_tile_names_from_bbox(st_bbox(polys_sf)))
  }

  dt0 <- compute_forest_share_for_year(polys_sf, file.path(root_dir, as.character(year0)),
                                       poly_id_col = poly_id_col,
                                       workers = workers,
                                       relevant_tile_names = relevant_tile_names)
  setnames(dt0, c("forest", "valid", "share"),
           c(sprintf("forest_%s", year0), sprintf("valid_%s", year0), sprintf("share_%s", year0)))

  dt1 <- compute_forest_share_for_year(polys_sf, file.path(root_dir, as.character(year1)),
                                       poly_id_col = poly_id_col,
                                       workers = workers,
                                       relevant_tile_names = relevant_tile_names)
  setnames(dt1, c("forest", "valid", "share"),
           c(sprintf("forest_%s", year1), sprintf("valid_%s", year1), sprintf("share_%s", year1)))

  out <- merge(dt0, dt1, by = poly_id_col, all = TRUE)
  out[, net_change_pp := 100 * (get(sprintf("share_%s", year1)) - get(sprintf("share_%s", year0)))]
  out[]
}

# Example usage:
# my.sf <- st_read("...")
# my.sf$poly_id <- seq_len(nrow(my.sf))
# out <- compute_forest_share_change(
#   polys_sf = my.sf,
#   root_dir = "/mnt/windows/Deforestation and poverty - Forest coverage",
#   year0 = 2000,
#   year1 = 2005,
#   poly_id_col = "poly_id",
#   workers = 8
# )
# my.sf <- merge(my.sf, out, by = "poly_id", all.x = TRUE)

my.sf <- readRDS("my.sf.RDS")
my.sf$poly_id <- seq_len(nrow(my.sf))
