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