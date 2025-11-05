# version 3 of code in this project 


Dataset_Luisa$duration_stimulants <- 
  Dataset_Luisa$last_prescription_before_17_stimulants - 
  Dataset_Luisa$first_prescription_stimulants


Dataset_Luisa$long_term_user_stimulants <- 
  ifelse(Dataset_Luisa$duration_stimulants > 12, 1, 0)



Dataset_Luisa$long_term_user_mpr <- 
  ifelse(Dataset_Luisa$mpr_before_17_stimulants >= 0.8, 1, 0)


Dataset_Luisa$early_starter <- 
  ifelse(Dataset_Luisa$first_prescription_stimulants < 12, 1, 0)



Dataset_Luisa$CI_and_stimulant_user <- ifelse(
  !is.na(Dataset_Luisa$first_prescription_CImeds) & 
    !is.na(Dataset_Luisa$first_prescription_stimulants), 1, 0)



Dataset_Luisa$gap_to_CImeds <- 
  Dataset_Luisa$first_prescription_CImeds - 
  Dataset_Luisa$first_prescription_stimulants




Dataset_Luisa$overlap_use <- ifelse(
  Dataset_Luisa$last_prescription_before_17_stimulants >= 
    Dataset_Luisa$first_prescription_CImeds, 1, 0)




Dataset_Luisa$treatment_trajectory <- with(Dataset_Luisa, ifelse(
  first_prescription_stimulants < 12 & duration_stimulants <= 12, "early_brief", ifelse(
    first_prescription_stimulants < 12 & duration_stimulants > 12, "early_sustained", ifelse(
      first_prescription_stimulants >= 12 & duration_stimulants <= 12, "late_brief", ifelse(
        first_prescription_stimulants >= 12 & duration_stimulants > 12, "late_sustained", NA)))))






# Standardize all PGS columns across traits and versions
pgs_vars <- grep("SCORESUM_PRS_", names(Dataset_Luisa), value = TRUE)

for (pgs in pgs_vars) {
  z_name <- paste0(pgs, "_z")
  Dataset_Luisa[[z_name]] <- scale(Dataset_Luisa[[pgs]])
}



for (pgs in pgs_vars) {
  int_name <- paste0(pgs, "_int")
  Dataset_Luisa[[int_name]] <- round(Dataset_Luisa[[pgs]] * 1000)
}




for (pgs in pgs_vars) {
  tertile_name <- paste0(pgs, "_tertile")
  Dataset_Luisa[[tertile_name]] <- cut(
    Dataset_Luisa[[pgs]],
    breaks = quantile(Dataset_Luisa[[pgs]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("Low", "Medium", "High")
  )
}




# Traits and versions
traits <- c("ADHD", "ASD", "MDD", "SCZ")
versions <- c("genr3", "genr4")

results_list <- list()

for (version in versions) {
  # Choose correct ancestry PCs for this version
  pc_vars <- paste0("C", 1:20, "_", version)
  
  for (trait in traits) {
    
    # Build variable name for this PGS
    pgs_var <- paste0("SCORESUM_PRS_", trait, "_", version)
    
    # MODEL 1: Logistic regression for treatment initiation
    formula_logit <- as.formula(
      paste("treatment_init ~", pgs_var, "+ GENDER + INCOME5 + EDUCM + EDUCP +", 
            paste(pc_vars, collapse = " + "))
    )
    
    model_logit <- glm(formula_logit, data = Dataset_Luisa, family = "binomial")
    
    # MODEL 2: Linear regression for symptom severity in treated children
    data_treated <- subset(Dataset_Luisa, treatment_init == 1)
    
    formula_lm <- as.formula(
      paste("Total_Problems_TScore_5 ~", pgs_var, "+ GENDER + INCOME5 + EDUCM + EDUCP +", 
            paste(pc_vars, collapse = " + "))
    )
    
    model_lm <- lm(formula_lm, data = data_treated)
    
    # Store results
    results_list[[paste(trait, version, "logit", sep = "_")]] <- summary(model_logit)
    results_list[[paste(trait, version, "lm", sep = "_")]] <- summary(model_lm)
  }
}




# Create binary variable for treatment initiation
Dataset_Luisa$treatment_init <- ifelse(
  Dataset_Luisa$first_prescription_stimulants == 1 | 
    Dataset_Luisa$first_prescription_CImeds == 1, 1, 0)




> Dataset_Luisa$treatment_init <- ifelse(
  +     !is.na(Dataset_Luisa$first_prescription_stimulants) |
    +         !is.na(Dataset_Luisa$first_prescription_CImeds), 1, 0)
> 
  > table(Dataset_Luisa$treatment_init, useNA = "ifany")




install.packages("nnet")  # if not already installed




library(nnet)

# Traits and versions
traits <- c("ADHD", "ASD", "MDD", "SCZ")
versions <- c("genr3", "genr4")

trajectory_results <- list()

for (version in versions) {
  pc_vars <- paste0("C", 1:20, "_", version)
  
  for (trait in traits) {
    pgs_var <- paste0("SCORESUM_PRS_", trait, "_", version, "_z")
    
    # Subset data with complete trajectory and PGS
    complete_data <- Dataset_Luisa[complete.cases(Dataset_Luisa[, c("treatment_trajectory", pgs_var, "GENDER", "INCOME5", "EDUCM", "EDUCP", pc_vars)]), ]
    
    if (nrow(complete_data) > 50) {
      formula_multinom <- as.formula(
        paste("treatment_trajectory ~", pgs_var, "+ GENDER + INCOME5 + EDUCM + EDUCP +", 
              paste(pc_vars, collapse = " + "))
      )
      
      multinom_model <- multinom(formula_multinom, data = complete_data)
      
      trajectory_results[[paste(trait, version, "multinom", sep = "_")]] <- summary(multinom_model)
    } else {
      warning(paste("Too few complete cases for", trait, version, "treatment trajectory model"))
    }
  }
}


13-07-25



# Function to extract and compute p-values from multinom summary
extract_multinom_pvals <- function(model_summary, trait, version) {
  coefs <- model_summary$coefficients
  ses   <- model_summary$standard.errors
  
  # Create dataframe of estimates, SEs, z-scores, p-values
  coef_df <- as.data.frame(coefs)
  se_df   <- as.data.frame(ses)
  
  output <- data.frame()
  
  for (row in 1:nrow(coef_df)) {
    z_vals <- coef_df[row, ] / se_df[row, ]
    p_vals <- 2 * (1 - pnorm(abs(z_vals)))
    
    row_df <- data.frame(
      trait     = trait,
      version   = version,
      category  = rownames(coef_df)[row],
      variable  = names(z_vals),
      estimate  = as.numeric(coef_df[row, ]),
      se        = as.numeric(se_df[row, ]),
      z         = as.numeric(z_vals),
      p         = as.numeric(p_vals)
    )
    
    output <- rbind(output, row_df)
  }
  
  return(output)
}


all_pvals <- data.frame()

for (name in names(trajectory_results)) {
  parts <- strsplit(name, "_")[[1]]
  trait <- parts[1]
  version <- parts[2]
  
  model_summary <- trajectory_results[[name]]
  
  extracted <- extract_multinom_pvals(model_summary, trait, version)
  
  if (!is.null(extracted)) {
    all_pvals <- rbind(all_pvals, extracted)
  }
}



extract_multinom_pvals <- function(model_summary, trait, version) {
  coefs <- model_summary$coefficients
  ses   <- model_summary$standard.errors
  
  # If empty or malformed, skip
  if (length(coefs) == 0 || length(ses) == 0) return(NULL)
  
  output <- data.frame()
  
  # Handle case when output is a single row (vector)
  if (is.null(dim(coefs))) {
    z_vals <- coefs / ses
    p_vals <- 2 * (1 - pnorm(abs(z_vals)))
    
    output <- data.frame(
      trait     = trait,
      version   = version,
      category  = "level1_vs_ref",
      variable  = names(coefs),
      estimate  = as.numeric(coefs),
      se        = as.numeric(ses),
      z         = as.numeric(z_vals),
      p         = as.numeric(p_vals)
    )
  } else {
    # For multi-row models
    coef_df <- as.data.frame(coefs)
    se_df   <- as.data.frame(ses)
    
    for (row in 1:nrow(coef_df)) {
      z_vals <- coef_df[row, ] / se_df[row, ]
      p_vals <- 2 * (1 - pnorm(abs(z_vals)))
      
      row_df <- data.frame(
        trait     = trait,
        version   = version,
        category  = rownames(coef_df)[row],
        variable  = names(z_vals),
        estimate  = as.numeric(coef_df[row, ]),
        se        = as.numeric(se_df[row, ]),
        z         = as.numeric(z_vals),
        p         = as.numeric(p_vals)
      )
      
      output <- rbind(output, row_df)
    }
  }
  
  return(output)
}



extract_multinom_pvals <- function(model_summary, trait, version) {
  coefs <- model_summary$coefficients
  ses   <- model_summary$standard.errors
  
  # Sanity check
  if (!all(names(coefs) == names(ses))) {
    stop("Mismatch in coefficient and SE names")
  }
  
  # Calculate z and p values
  z_vals <- coefs / ses
  p_vals <- 2 * (1 - pnorm(abs(z_vals)))
  p_vals[p_vals > 1] <- 1  # safety cap
  
  df <- data.frame(
    trait    = trait,
    version  = version,
    category = "level1_vs_ref",
    variable = names(coefs),
    estimate = coefs,
    se       = ses,
    z        = z_vals,
    p        = p_vals,
    row.names = NULL
  )
  
  return(df)
}




all_pvals <- data.frame()

for (name in names(trajectory_results)) {
  parts <- strsplit(name, "_")[[1]]
  trait <- parts[1]
  version <- parts[2]
  
  model_summary <- trajectory_results[[name]]
  extracted <- extract_multinom_pvals(model_summary, trait, version)
  
  if (!is.null(extracted)) {
    all_pvals <- rbind(all_pvals, extracted)
  }
}



extract_multinom_pvals <- function(model_summary, trait, version) {
  coefs <- model_summary$coefficients
  ses   <- model_summary$standard.errors
  
  if (is.null(coefs) || is.null(ses)) return(NULL)
  
  # Ensure both are numeric and same length
  stopifnot(length(coefs) == length(ses), all(names(coefs) == names(ses)))
  
  z_vals <- coefs / ses
  p_vals <- 2 * (1 - pnorm(abs(z_vals)))
  
  # Clean up invalid values
  p_vals[is.na(p_vals)] <- 1
  p_vals[p_vals > 1] <- 1
  p_vals[p_vals < 0] <- 0
  
  df <- data.frame(
    trait    = rep(trait, length(coefs)),
    version  = rep(version, length(coefs)),
    category = rep("level1_vs_ref", length(coefs)),
    variable = names(coefs),
    estimate = as.numeric(coefs),
    se       = as.numeric(ses),
    z        = as.numeric(z_vals),
    p        = as.numeric(p_vals),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  return(df)
}



all_pvals <- data.frame()

for (name in names(trajectory_results)) {
  parts <- strsplit(name, "_")[[1]]
  trait <- parts[1]
  version <- parts[2]
  
  model_summary <- trajectory_results[[name]]
  extracted <- extract_multinom_pvals(model_summary, trait, version)
  
  if (!is.null(extracted)) {
    all_pvals <- rbind(all_pvals, extracted)
  }
}




14-07-25


05-11-25

# Required packages
packages <- c("dplyr","broom","purrr","ggplot2","forcats","tidyr","car")
install_if_missing <- function(pk){ if(!requireNamespace(pk, quietly=TRUE)) install.packages(pk) }
invisible(lapply(packages, install_if_missing))
library(dplyr); library(broom); library(purrr); library(ggplot2); library(forcats); library(tidyr); library(car)

# --- 0) Quick sanity: check object exists
if(!exists("Dataset_Luisa")) stop("Dataset_Luisa not found in environment. Load your data first.")

# --- 1) Rebuild duration and trajectory variables (keep the 12 cutoffs as you confirmed)
# NOTE: these assume first_prescription_* and last_prescription_before_17_* are numeric (in months or years consistent across vars)
Dataset_Luisa <- Dataset_Luisa %>%
  mutate(
    duration_stimulants = last_prescription_before_17_stimulants - first_prescription_stimulants,
    long_term_user_stimulants = ifelse(!is.na(duration_stimulants) & duration_stimulants > 12, 1, 0),
    long_term_user_mpr = ifelse(!is.na(mpr_before_17_stimulants) & mpr_before_17_stimulants >= 0.8, 1, 0),
    early_starter = ifelse(!is.na(first_prescription_stimulants) & first_prescription_stimulants < 12, 1, 0),
    CI_and_stimulant_user = ifelse(!is.na(first_prescription_CImeds) & !is.na(first_prescription_stimulants), 1, 0),
    gap_to_CImeds = ifelse(!is.na(first_prescription_CImeds) & !is.na(first_prescription_stimulants),
                           first_prescription_CImeds - first_prescription_stimulants, NA_integer_),
    overlap_use = ifelse(!is.na(last_prescription_before_17_stimulants) & !is.na(first_prescription_CImeds) &
                         last_prescription_before_17_stimulants >= first_prescription_CImeds, 1, 0)
  )

# Create the categorical treatment_trajectory per your original logic
Dataset_Luisa <- Dataset_Luisa %>%
  mutate(treatment_trajectory = case_when(
    !is.na(first_prescription_stimulants) & first_prescription_stimulants < 12 & !is.na(duration_stimulants) & duration_stimulants <= 12 ~ "early_brief",
    !is.na(first_prescription_stimulants) & first_prescription_stimulants < 12 & !is.na(duration_stimulants) & duration_stimulants > 12  ~ "early_sustained",
    !is.na(first_prescription_stimulants) & first_prescription_stimulants >= 12 & !is.na(duration_stimulants) & duration_stimulants <= 12 ~ "late_brief",
    !is.na(first_prescription_stimulants) & first_prescription_stimulants >= 12 & !is.na(duration_stimulants) & duration_stimulants > 12  ~ "late_sustained",
    TRUE ~ NA_character_
  ))
Dataset_Luisa$treatment_trajectory <- factor(Dataset_Luisa$treatment_trajectory,
                                             levels = c("early_brief","early_sustained","late_brief","late_sustained"))

# --- 2) Treatment initiation factor per your choice (0 = no treatment, 1 = stimulant only, 2 = CI only, 3 = both)
Dataset_Luisa <- Dataset_Luisa %>%
  mutate(
    stimulant_init = ifelse(!is.na(first_prescription_stimulants), 1, 0),
    CI_init = ifelse(!is.na(first_prescription_CImeds), 1, 0),
    treatment_init_cat = case_when(
      stimulant_init == 0 & CI_init == 0 ~ 0L,
      stimulant_init == 1 & CI_init == 0 ~ 1L,
      stimulant_init == 0 & CI_init == 1 ~ 2L,
      stimulant_init == 1 & CI_init == 1 ~ 3L,
      TRUE ~ NA_integer_
    )
  )
Dataset_Luisa$treatment_init_cat <- factor(Dataset_Luisa$treatment_init_cat,
                                           levels = c(0,1,2,3),
                                           labels = c("no_treatment","stimulant_only","CI_only","both"))

# Also keep a binary for "any treatment" (useful in binary logistic models)
Dataset_Luisa <- Dataset_Luisa %>%
  mutate(treatment_any = ifelse(stimulant_init == 1 | CI_init == 1, 1, 0))

# --- 3) Standardize RAW numeric PRS only (safe approach)
# pattern to select raw PRS (genr3/genr4 raw)
pgs_raw <- grep("^SCORESUM_PRS_.*_genr[34]$", names(Dataset_Luisa), value = TRUE)
message("Found PRS variables: ", paste(head(pgs_raw,20), collapse = ", "))

for(pgs in pgs_raw){
  # coerce to numeric (in case of factors/characters)
  x <- suppressWarnings(as.numeric(Dataset_Luisa[[pgs]]))
  if(all(is.na(x))){
    warning(paste("PRS", pgs, "could not be coerced to numeric — check variable. Skipping."))
    next
  }
  zname <- paste0(pgs, "_z")
  Dataset_Luisa[[zname]] <- as.numeric(scale(x))
}

# --- 4) QC / validation checks before modelling
qc <- list()
qc$n_rows <- nrow(Dataset_Luisa)
qc$n_treated_any <- sum(Dataset_Luisa$treatment_any == 1, na.rm = TRUE)
qc$treatment_cat_counts <- table(Dataset_Luisa$treatment_init_cat, useNA = "ifany")
qc$trajectory_counts <- table(Dataset_Luisa$treatment_trajectory, useNA = "ifany")

print(qc)

# Missingness table for core covariates
core_covs <- c("GENDER","INCOME5","EDUCM","EDUCP")
missing_core <- sapply(core_covs, function(v) sum(is.na(Dataset_Luisa[[v]])))
message("Missingness in core covariates:")
print(missing_core)

# Show PRS summary (means, sd)
pgs_z <- grep("_genr[34]_z$", names(Dataset_Luisa), value = TRUE)
if(length(pgs_z)>0){
  print(as.data.frame(t(sapply(pgs_z, function(v) c(mean = mean(Dataset_Luisa[[v]], na.rm=TRUE),
                                                   sd = sd(Dataset_Luisa[[v]], na.rm=TRUE)))),
                      row.names = pgs_z))
}

# Optional VIF check example for one model (requires numeric covariates) - only run if car::vif available
# (we will run VIFs after model building below)

# --- 5) Binary logistic regression: one PRS per model per generation
traits <- c("ADHD","ASD","MDD","SCZ")
gens <- c("genr3","genr4")

# store tidy results
binary_results <- list()
binary_tidy <- list()

for(trait in traits){
  for(gen in gens){
    prs_var <- paste0("SCORESUM_PRS_", trait, "_", gen, "_z")
    # check var exists
    if(!prs_var %in% names(Dataset_Luisa)){
      warning(prs_var, " not found — skipping")
      next
    }

    # Subset to complete cases for the binary model (treatment_any + covariates + PRS)
    covs <- c("treatment_any", prs_var, "GENDER","INCOME5","EDUCM","EDUCP", paste0("C",1:20,"_",gen))
    covs_present <- covs[covs %in% names(Dataset_Luisa)]
    complete_data <- Dataset_Luisa %>% select(all_of(covs_present)) %>% filter(!is.na(treatment_any))

    # ensure at least some treated and untreated
    if(nrow(complete_data) < 50 | length(unique(complete_data$treatment_any)) < 2){
      warning("Too few complete cases or lacking outcome variability for ", trait, gen, " — skipping")
      next
    }

    # Build formula
    formula_str <- paste("treatment_any ~", prs_var, "+ GENDER + INCOME5 + EDUCM + EDUCP +",
                         paste(intersect(paste0("C",1:20,"_",gen), names(complete_data)), collapse = " + "))
    model <- glm(as.formula(formula_str), data = complete_data, family = binomial())
    s <- summary(model)

    # store
    binary_results[[paste(trait, gen, sep = "_")]] <- model

    # tidy the PRS row using broom::tidy
    tidy_df <- broom::tidy(model) %>%
      filter(term == prs_var) %>%
      mutate(
        trait = trait, generation = gen,
        OR = exp(estimate),
        CI_lower = exp(estimate - 1.96 * std.error),
        CI_upper = exp(estimate + 1.96 * std.error)
      ) %>%
      select(trait, generation, term, estimate, std.error, statistic, p.value, OR, CI_lower, CI_upper)

    binary_tidy[[paste(trait,gen,sep="_")]] <- tidy_df
  }
}

binary_summary <- bind_rows(binary_tidy)

# Print summary table
print(binary_summary)

# --- 6) Simple forest plot for ORs (binary models)
if(nrow(binary_summary) > 0){
  plot_df <- binary_summary %>%
    mutate(label = paste0(trait, " (", generation, ")"),
           label = fct_reorder(label, OR))

  gg_forest <- ggplot(plot_df, aes(x = OR, y = label)) +
    geom_point() +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Odds ratio (log scale)", y = "", title = "PRS -> Any treatment (binary logistic) - ORs with 95% CI") +
    theme_minimal()

  print(gg_forest)
  # Save plot (optional)
  ggsave("prs_binary_forest_plot.png", gg_forest, width = 8, height = 4, dpi = 300)
}

# --- 7) Save results to CSV for table insertion into paper/pptx
if(nrow(binary_summary)>0){
  write.csv(binary_summary, "binary_prs_results_summary.csv", row.names = FALSE)
  message("Saved binary_prs_results_summary.csv and prs_binary_forest_plot.png to working directory.")
}

# --- 8) Optional: VIF check for one representative model (ADHD_genr3 if present)
if("ADHD_genr3" %in% names(binary_results)){
  m <- binary_results[["ADHD_genr3"]]
  # VIF needs a model with no missing and only numeric/factor vars
  try({
    vif_vals <- car::vif(m)
    print(vif_vals)
  }, silent = TRUE)
}

# --- End of script
