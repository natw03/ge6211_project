# Load packages ----
library(MatchIt)
library(sf)
library(dplyr)

# Load datasets ----
kalimantan_buffer   <- st_read("Data/final_shapefiles/buffer_shp/buffer_5km_borneo.shp")
sumatra_buffer      <- st_read("Data/final_shapefiles/buffer_shp/buffer_5km_sumatra.shp")

kalimantan_674_shp  <- st_read("Data/final_shapefiles/study_area_shp/kalimantan_674_area.shp")
kalimantan_1477_shp <- st_read("Data/final_shapefiles/study_area_shp/kalimantan_1477_area.shp")
sumatra_1498_shp    <- st_read("Data/final_shapefiles/study_area_shp/sumatra_1498_area.shp")
sumatra_2403_shp    <- st_read("Data/final_shapefiles/study_area_shp/sumatra_2403_area.shp")

kalimantan_674_csv  <- read.csv("Data/combined_csv/modelling_table_674_clipped.csv")
kalimantan_1477_csv <- read.csv("Data/combined_csv/modelling_table_1477_clipped.csv")
sumatra_1498_csv    <- read.csv("Data/combined_csv/modelling_table_1498_clipped.csv")
sumatra_2403_csv    <- read.csv("Data/combined_csv/modelling_table_2403_clipped.csv")

# Functions ----
reproject_if_needed <- function(shp, target_epsg) {
  current_epsg <- st_crs(shp)$epsg
  if (is.na(current_epsg) || current_epsg != target_epsg) {
    message("Reprojecting from EPSG:", ifelse(is.na(current_epsg), "unknown", current_epsg),
            " to EPSG:", target_epsg)
    shp <- st_transform(shp, target_epsg)
  } else {
    message("Already in EPSG:", target_epsg, " — no reprojection needed.")
  }
  return(shp)
}

run_matching <- function(csv, site_id, project_start_year, n_pairs = 50, calipers = NULL) {
  
  pre_years <- (project_start_year - 5):(project_start_year - 1)
  cat("\n========================================\n")
  cat("Site:", site_id, "| Pre-treatment years:", min(pre_years), "-", max(pre_years), "\n")
  cat("========================================\n")
  
  ## Fill terrain variables ----
  # slope and aspect are static; -9999 sentinel replaced and value propagated across all years.
  csv <- csv %>%
    mutate(
      slope_fill  = ifelse(slope  == -9999, NA, slope),
      aspect_fill = ifelse(aspect == -9999, NA, aspect)
    ) %>%
    group_by(pixel_id) %>%
    mutate(
      slope_fill  = first(na.omit(slope_fill)),
      aspect_fill = first(na.omit(aspect_fill))
    ) %>%
    ungroup()
  
  n_no_slope <- csv %>%
    group_by(pixel_id) %>%
    summarise(all_na = all(is.na(slope_fill)), .groups = "drop") %>%
    filter(all_na) %>%
    nrow()
  cat("Pixels with no slope/aspect (imputed with median):", n_no_slope, "\n")
  
  ## Identify undisturbed pixels ----
  # Pixel is disturbed if burned == 1 OR fire_count > 0 in any pre-treatment year.
  disturbed <- csv %>%
    filter(year %in% pre_years) %>%
    filter(burned == 1 | fire_count > 0) %>%
    pull(pixel_id) %>%
    unique()
  cat("Ever-disturbed pixels excluded:", length(disturbed), "\n")
  
  pre_treatment <- csv %>%
    filter(year %in% pre_years) %>%
    filter(!pixel_id %in% disturbed)
  
  ## Create covariate table ----
  # Time-varying covariates averaged across pre-treatment window.
  covariate_table <- pre_treatment %>%
    group_by(pixel_id, in_study_area) %>%
    summarise(
      dem                          = mean(dem,                     na.rm = TRUE),
      slope                        = mean(slope_fill,              na.rm = TRUE),
      aspect                       = mean(aspect_fill,             na.rm = TRUE),
      dist_builtup                 = mean(dist_builtup,            na.rm = TRUE),
      dist_river                   = mean(dist_river,              na.rm = TRUE),
      mean_distance_to_forest_edge = mean(distance_to_forest_edge, na.rm = TRUE),
      mean_ndvi                    = mean(ndvi,                    na.rm = TRUE),
      mean_ndmi                    = mean(ndmi,                    na.rm = TRUE),
      mean_lst                     = mean(lst,                     na.rm = TRUE),
      mean_precip                  = mean(precip,                  na.rm = TRUE),
      lulc_forest                  = mean(lulc_forest,             na.rm = TRUE),
      lulc_oilpalm                 = mean(lulc_oilpalm,            na.rm = TRUE),
      lulc_pulpwood                = mean(lulc_pulpwood,           na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      !is.na(dem),
      dem > -100,
      !is.na(dist_builtup),
      !is.na(mean_distance_to_forest_edge),
      !is.na(mean_ndvi)
    ) %>%
    mutate(
      slope  = ifelse(is.nan(slope)  | is.na(slope),  median(slope,  na.rm = TRUE), slope),
      aspect = ifelse(is.nan(aspect) | is.na(aspect), median(aspect, na.rm = TRUE), aspect)
    )
  cat("Pixels available — treatment:", sum(covariate_table$in_study_area == 1),
      "| control:", sum(covariate_table$in_study_area == 0), "\n")
  
  ## Sample treatment pixels ----
  set.seed(42)
  treatment_pixels <- covariate_table %>% filter(in_study_area == 1)
  
  if (nrow(treatment_pixels) > n_pairs) {
    sampled_ids <- sample(treatment_pixels$pixel_id, n_pairs)
    match_input <- covariate_table %>%
      filter(in_study_area == 0 | pixel_id %in% sampled_ids)
  } else {
    match_input <- covariate_table
    message("Fewer than ", n_pairs, " treatment pixels available (n = ",
            nrow(treatment_pixels), "). Using all.")
  }
  
  ## Run matching ----
  # 1:1 nearest-neighbour without replacement, Mahalanobis distance (Guizar-Coutino et al. 2022).
  # calipers optionally passed per site to tighten balance on specific covariates (0.25 * SD).
  match_out <- matchit(
    in_study_area ~ dem + slope + dist_builtup +
      mean_distance_to_forest_edge + mean_ndvi,
    data        = match_input,
    method      = "nearest",
    distance    = "mahalanobis",
    caliper     = calipers,
    std.caliper = FALSE,
    replace     = FALSE,
    ratio       = 1
  )
  
  ## Balance check ----
  # Covariates with ASMD >= 0.25 indicate insufficient balance (Stuart 2010).
  match_summary <- summary(match_out, standardize = TRUE)
  std_diffs     <- match_summary$sum.matched[, "Std. Mean Diff."]
  
  cat("\n--- ASMD (matched) ---\n")
  print(round(std_diffs, 3))
  
  failed_covs <- names(std_diffs[abs(std_diffs) >= 0.25])
  if (length(failed_covs) == 0) {
    cat("✓ All covariates meet ASMD < 0.25 threshold.\n")
  } else {
    cat("⚠ Covariates exceeding threshold:", paste(failed_covs, collapse = ", "), "\n")
  }
  
  ## Extract matched pairs ----
  matched_data  <- match.data(match_out)
  treatment_ids <- matched_data %>% filter(in_study_area == 1) %>% pull(pixel_id)
  control_ids   <- matched_data %>% filter(in_study_area == 0) %>% pull(pixel_id)
  cat("Matched pairs:", length(treatment_ids), "\n")
  
  # Build paired table: one row per pair with treatment and control pixel_id side by side
  match_matrix  <- match_out$match.matrix
  matched_pairs <- data.frame(
    treatment_pixel_id = as.integer(rownames(match_matrix)),
    control_pixel_id   = as.integer(match_matrix[, 1])
  ) %>%
    filter(!is.na(control_pixel_id))
  
  ## Save output ----
  dir.create("outputs", showWarnings = FALSE)
  write.csv(matched_pairs, paste0("outputs/matched_pairs_", site_id, ".csv"), row.names = FALSE)
  cat("Matched pairs saved for site", site_id, "\n")
  
  return(list(
    match_out     = match_out,
    std_diffs     = std_diffs,
    failed_covs   = failed_covs,
    matched_pairs = matched_pairs
  ))
}

# Reproject ----
kalimantan_crs <- 23835
sumatra_crs    <- 23831

kalimantan_buffer   <- reproject_if_needed(kalimantan_buffer,   kalimantan_crs)
sumatra_buffer      <- reproject_if_needed(sumatra_buffer,      sumatra_crs)
kalimantan_674_shp  <- reproject_if_needed(kalimantan_674_shp,  kalimantan_crs)
kalimantan_1477_shp <- reproject_if_needed(kalimantan_1477_shp, kalimantan_crs)
sumatra_1498_shp    <- reproject_if_needed(sumatra_1498_shp,    sumatra_crs)
sumatra_2403_shp    <- reproject_if_needed(sumatra_2403_shp,    sumatra_crs)

# Run matching ----
# Calipers set at 0.25 * SD of the covariate across the full pixel pool, applied only
# to sites where those covariates failed the ASMD threshold without calipers.
results_674  <- run_matching(kalimantan_674_csv,  site_id = "674",  project_start_year = 2009, n_pairs = 50,
                             calipers = c(mean_ndvi = 0.0169, dist_builtup = 703.64))
results_1477 <- run_matching(kalimantan_1477_csv, site_id = "1477", project_start_year = 2010, n_pairs = 50)
results_1498 <- run_matching(sumatra_1498_csv,    site_id = "1498", project_start_year = 2014, n_pairs = 50,
                             calipers = c(dem = 4.50, slope = 0.73))
results_2403 <- run_matching(sumatra_2403_csv,    site_id = "2403", project_start_year = 2016, n_pairs = 50)

# Save combined ASMD table ----
asmd_all <- bind_rows(
  data.frame(site = "674",  covariate = names(results_674$std_diffs),  asmd = round(results_674$std_diffs,  3)),
  data.frame(site = "1477", covariate = names(results_1477$std_diffs), asmd = round(results_1477$std_diffs, 3)),
  data.frame(site = "1498", covariate = names(results_1498$std_diffs), asmd = round(results_1498$std_diffs, 3)),
  data.frame(site = "2403", covariate = names(results_2403$std_diffs), asmd = round(results_2403$std_diffs, 3))
) %>%
  mutate(passed = abs(asmd) < 0.25)

write.csv(asmd_all, "outputs/asmd_all_sites.csv", row.names = FALSE)
cat("ASMD table saved.\n")