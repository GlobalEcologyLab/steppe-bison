# Extinction time: 8734 (8810 - 8657)
# Extinction location: 135 E 67.6 N
# Peak abundance: 43,500 BP to 32,500 BP
# Occupancy pattern: 30 fossil sites

library(fs)
library(poems)
library(paleopop)
library(furrr)
library(purrr)
library(dplyr)
library(data.table)
library(readxl)
library(rgdal)
SOURCE_DIR <- getwd()
RESULTS_DIR <- file.path(SOURCE_DIR, "results")
data_dir <- file.path(SOURCE_DIR, "data")
projCRS <- "+proj=laea +lon_0=138.6132813 +lat_0=69.9047412 +datum=WGS84 +units=m +no_defs"
region <- readRDS(file.path(data_dir, "projected_paleo_region.RDS"))

files <- dir_ls(RESULTS_DIR, glob = "**.RData") %>% gtools::mixedsort()

parallel_cores <- 8
plan(multisession, workers = parallel_cores)

# Extinction time

extinction_time <- future_map_dbl(files, function(f) {
  extinction_index <- f %>% readRDS() %>% .$abundance %>% .[, 101:3851] %>%
    colSums(na.rm=T) %>% detect_index(~.==0)
  extinction_time <- -12*(3751 - extinction_index)-5001
  # if (extinction_time < -8810) {
  #   penalty <- abs(extinction_time) - 8810
  # } else if (extinction_time > -8657) {
  #   penalty <- 8657 - abs(extinction_time)
  # } else {
  #   penalty <- 0
  # }
  # return(penalty)
  return(extinction_time)
}, .progress = TRUE)

# Peak abundance

peak_abundance_penalty <- future_map_dbl(files, function(f) {
  peak_index <- f %>% readRDS() %>% .$abundance %>% .[, 101:3851] %>%
    colSums(na.rm=T) %>% zoo::rollapply(width = 10, FUN = mean, fill = NA) %>%
    which.max()
  peak_abundance <- -12*(3751 - peak_index)-5001
  if (peak_abundance < -43500) {
    penalty <- abs(peak_abundance) - 43500
  } else if (peak_abundance > -32500) {
    penalty <- 32500 - abs(peak_abundance)
  } else {
    penalty <- 0
  }
  return(penalty)
}, .progress = TRUE)

# Occupancy pattern

row_index <- c(2L, 7L, 8L, 11L, 12L, 13L, 64L, 66L, 68L, 69L, 70L, 71L, 76L, 77L, 89L, 96L, 133L, 140L, 141L, 142L, 147L,
               209L, 211L, 212L, 218L, 219L, 280L, 321L, 330L, 342L) %>% `+`(1) %>% prepend(1L)
fossil_data <- readxl::read_excel(file.path(data_dir, "bison_priscus_fossils.xlsx")) %>%
  dplyr::mutate(ID = 1:nrow(.)) %>% slice(row_index) %>%
  cbind(rgdal::project(as.matrix(select(., long, lat)), projCRS)) %>% select(-c(1:2))
indices <- raster::rasterize(x = dplyr::select(fossil_data, long, lat), y = region$region_raster,
                             field = fossil_data$ID) %>% raster::values() %>% purrr::map_lgl(~!is.na(.)) %>% which()
adjacent <- raster::rasterize(x = dplyr::select(fossil_data, long, lat), y = region$region_raster,
                              field = fossil_data$ID) %>%
  raster::adjacent(cells = indices, directions = 8, sorted=T, include=T, id=T)
fossil_raster <- raster::rasterize(x = dplyr::select(fossil_data, long, lat), y = region$region_raster,
                                   field = fossil_data$ID) %>%
  raster::mask(region$region_raster) %>% raster::as.data.frame(xy = T) %>%
  slice(adjacent[,3]) %>%
  cbind(adjacent) %>%
  mutate(layer = zoo::na.locf(if_else(from!=to, NA_real_, layer))) %>%
  select(x,y,layer)
fossil_data <- fossil_data %>% dplyr::left_join(fossil_raster, by = c("ID" = "layer")) %>%
  dplyr::mutate(max = round((50001 - (OxCal_age + OxCal_error))/12),
                min = round((50001 - (OxCal_age - OxCal_error))/12)) %>%
  dplyr::filter(max<3751, max>0, min<3751, min>0, !is.na(x))
id_coordinates <- data.table::as.data.table(cbind(pop_index = 1:3277, region$coordinates))
indexed_fossil_data <- as.data.frame(id_coordinates[fossil_data, on = c(x = "x", y = "y")])

occupancy_pattern <- future_map_dbl(files, function(f) {
  mat <- f %>% readRDS() %>% .$abundance %>% .[, 101:3851]
  presences <- purrr::map_lgl(
    purrr::map(1:length(indexed_fossil_data$pop_index),
               ~mat[indexed_fossil_data$pop_index[.],
                    indexed_fossil_data$max[.]:indexed_fossil_data$min[.]]),
    ~sum(., na.rm = T)>0
  )
  unique(dplyr::filter(
    dplyr::mutate(indexed_fossil_data, present = presences), present
  )$ID) -> sites
  return(length(sites))
}, .progress = TRUE)

# Extinction location

extinction_target <- indexed_fossil_data %>% slice(1:9) %>% select(x, y)

extinction_location <- future_map_dbl(files, function(f) {
  mat <- f %>% readRDS() %>% .$abundance %>% .[, 101:3851]
  extinction_index <- mat %>% colSums(na.rm=T) %>% detect_index(~.==0)
  if (extinction_index==0) {
    last_pop_indices <- which(as.logical(mat[, ncol(mat)]))
    if (length(last_pop_indices) > 1) {
      abundance_weights <- matrix(rep(mat[last_pop_indices, ncol(mat)], 2), ncol = 2)
      extinction_location <- .colSums(as.matrix(region$coordinates[last_pop_indices,])*abundance_weights, m = length(last_pop_indices), n = 2)/.colSums(abundance_weights, m = length(last_pop_indices), n = 2)
    } else {
      extinction_location <- as.numeric(region$coordinates[last_pop_indices,])
    }
  } else if (extinction_index > 1) {
    last_pop_indices <- which(as.logical(mat[, extinction_index - 1]))
    if (length(last_pop_indices) > 1) {
      abundance_weights <- matrix(rep(mat[last_pop_indices, extinction_index - 1], 2), ncol = 2)
      extinction_location <- .colSums(as.matrix(region$coordinates[last_pop_indices,])*abundance_weights, m = length(last_pop_indices), n = 2)/.colSums(abundance_weights, m = length(last_pop_indices), n = 2)
    } else {
      extinction_location <- as.numeric(region$coordinates[last_pop_indices,])
    }
  } else if (extinction_index == 1) {
    penalty <- 4214242
  }
  if (!exists("penalty")) {
    if (nrow(dplyr::filter(extinction_target, x==extinction_location[1], y==extinction_location[2]))>0) {
      penalty <- 0
    } else {
      penalty <- map_dbl(1:9, function(i) {
        sqrt((extinction_location[1]-extinction_target[i,1])^2 + (extinction_location[2]-extinction_target[i,2])^2)
      }) %>% min(na.rm=T)
    }
  }
  return(penalty)
}, .progress = TRUE)

summary_metrics <- cbind(extinction_penalty, peak_abundance_penalty, occupancy_pattern, extinction_location)

write.csv(summary_metrics, file.path(RESULTS_DIR, "summary_metrics.csv"))
