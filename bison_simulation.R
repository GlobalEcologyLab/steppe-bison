# Set directories for script and general and K Cuts and Human Density data files, and results
SOURCE_DIR <- getwd()
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")

# create results dir if needed
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Parallel cores available on machine
# 128 max
parallel_cores <- 95

# Install paleopop and dependencies
library(poems)
library(paleopop)
library(raster)
library(purrr)
library(stringr)

#### Step 1: Create a simulation model template with a study region ####
burn_in_steps <- 100
timesteps <- 3751 + burn_in_steps

# Step 1: Create a simulation model template with a study region

# Retrieve region, HS and human distributions
region <- readRDS(file.path(DATA_DIR, "projected_paleo_region.RDS"))
raster::plot(region$region_raster, main = "Bison region (cell indices)",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")
raster::plot(region$temporal_mask_raster()[[2626]], colNA = "blue")

# Re-usable distance matrix
distance_matrix <- DispersalGenerator$new(region = region, dispersal_max_distance = 500, # in km
                                          distance_scale = 1000)$calculate_distance_matrix()
summary(as.vector(distance_matrix))

# Distance-based environmental correlation (via a compacted Cholesky decomposition)
# mammoth method doc: b = 850 km
env_corr <- SpatialCorrelation$new(region = region, amplitude = 0.99, breadth = 850, distance_scale = 1000) # breadth in km
compact_decomposition <- env_corr$get_compact_decomposition(distance_matrix = distance_matrix) # default threshold

# Friction raster
friction <- readRDS(file.path(DATA_DIR, "bison_conductance_projected.RDS"))[region$region_indices,]
dim(friction)

# Human matrix
burn_in <- stack(replicate(100, raster(file.path(DATA_DIR, "projected_humans_mean.grd"))))
humans_mean_raster <- raster::stack(burn_in, raster::stack(file.path(DATA_DIR, "projected_humans_mean.grd")))
raster::plot(humans_mean_raster[[3850]], main = "Human mean final density",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")
raster::values(humans_mean_raster) <- raster::values(humans_mean_raster) # memory
humans_mean_matrix <- humans_mean_raster[region$region_indices]
humans_mean_matrix[!is.finite(humans_mean_matrix)] <- 0
burn_in <- stack(replicate(100, raster(file.path(DATA_DIR, "projected_humans_sd.grd"))))
humans_sd_raster <- raster::stack(burn_in, raster::stack(file.path(DATA_DIR, "projected_humans_sd.grd")))
raster::values(humans_sd_raster) <- raster::values(humans_sd_raster) # memory
humans_sd_matrix <- humans_sd_raster[region$region_indices]
humans_sd_matrix[!is.finite(humans_sd_matrix)] <- 0
rm(humans_mean_raster, humans_sd_raster); gc()

# Population (simulation) model template for fixed parameters
model_template <- PaleoPopModel$new(
  region = region,
  time_steps = timesteps, # include burn-in
  years_per_step = 12,
  populations = region$region_cells,
  # initial_abundance: generated
  transition_rate = 1.0,
  # standard_deviation: sampled
  compact_decomposition = compact_decomposition,
  # carrying_capacity: generated
  density_dependence = "logistic",
  # growth_rate_max: sampled
  harvest = TRUE,
  # harvest_max: sampled
  harvest_g = 0.4, # constant
  # harvest_z: sampled
  # harvest_max_n: sampled
  # human_density: generated
  dispersal_target_k = 10,
  # dispersal_data: generated
  # abundance_threshold: sampled,
  occupancy_threshold = 1,
  # niche_ref: sampled
  results_selection = c("abundance")
)

# Step 2: Create generators for initial abundance, carrying capacity, and dispersal

#### Step 2: Create generators for initial abundance, carrying capacity, and dispersal ####
##### Initial abundance and carrying capacity generated via example habitat suitability #####
# Here we set the fixed parameters that are template attached, so when we generate a clone of the generator, these will be parsed through
cuts <- list.files(K_CUTS_DIR) %>% stringr::str_split("_") %>% purrr::map(~.[3:4]) %>% purrr::map_chr(paste0, collapse = "_")

## Carrying-capacity generator
capacity_gen <- Generator$new(description = "capacity",
                              region = region,
                              generate_rasters = FALSE, # use but don't generate
                              burn_in_steps =  burn_in_steps,
                              generative_requirements = list(
                                hs_matrix = "file",
                                initial_abundance = "function",
                                carrying_capacity = "function"),
                              # these come from the latin hypercube sampler
                              inputs = c("density_max", "niche_ref"),
                              outputs = c("initial_abundance", "carrying_capacity"))
# Here we tell the generator to import the HS file and save it as "hs_matrix"
capacity_gen$add_file_template("hs_matrix",
                               path_template = file.path(NICHE_DIR, "clim_suit_%s_50k-5kBP_SCALED.RDS"),
                               path_params = "niche_ref",
                               file_type = "RDS")
# Here we subset the hs_matrix to have only the region cells, and we add the burn in.
# Also, we tell the generator to generate the carrying_capacity based on "density_max" and "hs_matrix".
capacity_gen$add_function_template("carrying_capacity",
                                   function_def = function(params) {
                                     hs_matrix <- params$hs_matrix[params$region$region_indices,]
                                     hs_matrix[!is.finite(hs_matrix)] <- 0
                                     # repeat the first timestep n times as a burn in
                                     hs_matrix <- cbind(replicate(params$burn_in_steps, hs_matrix[, 1]), hs_matrix)
                                     # round the density values
                                     round(params$density_max*hs_matrix)
                                   },
                                   call_params = c("density_max", "hs_matrix", "burn_in_steps", "region"))
# Here we tell the generator what function to use to generate initial_abundance
# based on the carrying capacity of the first time step
capacity_gen$add_function_template("initial_abundance",
                                   function_def = function(params) {
                                     params$carrying_capacity[, 1]
                                   },
                                   call_params = c("carrying_capacity"))
system.time({test_capacity <- capacity_gen$generate(input_values = list(density_max = 250, niche_ref = "0.9_4"))}) # ~4



# Distance-based dispersal generator: dispersal = p*exp(-1*distance/b) up to d_max (r)
#  known dispersal values
known_dispersals <- c(205, 190, 190, 190, 280, 280, 190, 155, 150, 125, 155, 125, 135, 135, 255, 260, 190)
summary(known_dispersals)
b_lookup <- data.frame(d_max = -Inf, b = 0:300)
for (i in 2:300) {
  b_lookup$d_max[i] <- which.max(exp(-1*(1:501)/b_lookup$b[i]) <= 0.19)
}
b_lookup$d_max[301] <- 501
dispersal_gen <- DispersalGenerator$new(region = region,
                                        dispersal_max_distance = 500, # km
                                        distance_classes = seq(10, 500, 10),
                                        distance_scale = 1000, # km
                                        dispersal_function_data = b_lookup,
                                        dispersal_friction = DispersalFriction$new(conductance = friction),
                                        inputs = c("dispersal_p",
                                                   "dispersal_r"),
                                        decimals = 3)
# dispersal_gen$calculate_distance_data(distance_matrix = distance_matrix) # pre-calculate
# saveRDS(dispersal_gen$distance_data, file.path(DATA_DIR, "dispersal_distance_data.RDS"))
dispersal_gen$distance_data <- readRDS(file = file.path(DATA_DIR, "dispersal_distance_data.RDS"))
system.time(test_dispersal <- dispersal_gen$generate(input_values =
                                                       list(dispersal_p = 0.5,
                                                            dispersal_r = 400))$dispersal_data)

##### Human Density #####
human_threshold_mean <- quantile(humans_mean_matrix[humans_mean_matrix > 0], 0.95, na.rm = FALSE)

# Build a human density generator
human_density_gen <- Generator$new(description = "Human Density Generator",
                                   humans_abundance = humans_mean_matrix,
                                   humans_var = humans_sd_matrix,
                                   spatial_correlation = env_corr,
                                   human_threshold = human_threshold_mean,
                                   generate_rasters = FALSE,
                                   generative_requirements = list(p_window = 'function',
                                                                  human_density = 'distribution'),
                                   inputs = c("p"),
                                   outputs = c("human_density"))

human_density_gen$add_function_template("p_window",
                                        function_def = function(params) {
                                          w <- params$p*10/100
                                          p_lower <- params$p - w
                                          p_upper <- params$p + w
                                          p_lower <- ifelse(p_lower < 0, 0, p_lower)
                                          p_upper <- ifelse(p_upper > 1, 1, p_upper)
                                          return(c(p_lower, p_upper))
                                        },
                                        call_params = c("p"))

human_density_gen$add_distribution_template("human_density",
                                            distr_type = "lognormal",
                                            distr_params = list(mean = "humans_abundance", sd = "humans_var"),
                                            sample = c("p_window"),
                                            normalize_threshold = "human_threshold")

## test the generator
system.time(test_humans <- human_density_gen$generate(input_values = list(p=0.1)))

rm(test_capacity, test_dispersal, test_humans); gc()

##### Simulation #####

# LHS samples

nsims <- 100
lhs_generator <- LatinHypercubeSampler$new()
lhs_generator$set_uniform_parameter("standard_deviation", lower = 0, upper = sqrt(0.06))
lhs_generator$set_uniform_parameter("growth_rate_max", lower = log(1.31), upper = log(2.84))
lhs_generator$set_uniform_parameter("abundance_threshold", lower = 0, upper = 500, decimals = 0)
lhs_generator$set_uniform_parameter("harvest_max", lower = 0, upper = 0.35)
lhs_generator$set_uniform_parameter("harvest_z", lower = 1, upper = 2)
lhs_generator$set_uniform_parameter("p", lower = 0, upper = 1) # ?
lhs_generator$set_uniform_parameter("density_max", lower = 500, upper = 3250) # also alias for harvest_max_n
lhs_generator$set_uniform_parameter("dispersal_p", lower = 0.05, upper = 0.25)
# lhs_generator$set_uniform_parameter("dispersal_b", lower = 65, upper = 145)
lhs_generator$set_uniform_parameter("dispersal_r", lower = 100, upper = 500)
lhs_generator$set_class_parameter("niche_ref", cuts)
sample_data <- lhs_generator$generate_samples(number = nsims, random_seed = 123)
sample_data$sample <- c(1:nsims)

#### Step 4: Build a simulation manager to run each simulation ####

# Create a simulation manager and run the sampled model simulations
sim_manager <- SimulationManager$new(sample_data = sample_data,
                                     model_template = model_template,
                                     generators = list(capacity_gen,
                                                       dispersal_gen,
                                                       human_density_gen),
                                     parallel_cores = parallel_cores,
                                     results_dir = RESULTS_DIR)

sim_manager$results_filename_attributes <- c("sample", "results")
run_output <- sim_manager$run()
run_output$summary
