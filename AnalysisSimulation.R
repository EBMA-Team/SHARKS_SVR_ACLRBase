To illustrate how the analysis will work - a pilot simulation is provided below.

The first thing is to retrieve some starting parameters from the SHARKS registry. A recent governance report \[link\] was modified to return descriptives for IKDC at baseline and some basic demographics from the Tibiofemoral Instability cohort.

```{r}

SHARKSIKDCMean <- 47.5
SHARKSIKDCSD <- 12.6
SHARKSSexProp <- round(358/519,digits = 2)
SHARKSAgeMean <- 27.3
SHARKSAgeSD <- 7.8

# BMI


```

Let's assume for now that the IKDCScore is correlated to age and sex - this may need to be tested in the SHARKS sample, or the correlation coefficient retrieved from the literature.

```{r}
#| label: simulation-dataset-sharks
#| echo: true

set.seed(2065)

SHARKSDef <- defData(
  varname = "age_offset",
  dist = "gamma",
  formula = "log(27.3 - 16)",  # Adjust mean for the shift
  variance = (SHARKSAgeSD^2) / (27.3 - 16)^2,  # Adjust variance
  link = "log"
)

SHARKSDef <- defData(
  SHARKSDef,
  varname = "Age",
  formula = "age_offset + 16",
  dist = "nonrandom"
)


SHARKSDef <- defData(
  SHARKSDef, 
  varname = "Male", 
  dist = "binary",
  formula = 0.69, 
  link = "logit"
  )

SHARKSDef <- defData(
  SHARKSDef,
  varname = "bmi_scale",
  formula = "log(26 + 0.15 * (Age - 25) + 3 * Male)",  # Base BMI relationship
  variance = 0.04,  # Controls the shape/spread
  dist = "gamma",
  link = "log"
)


SHARKSDef <- defData(
  SHARKSDef,
  varname = "BMI",
  formula = "pmax(pmin(bmi_scale, 40), 18)",
  dist = "nonrandom"
)

SHARKSDef <- defData(
  SHARKSDef,
  varname = "MeniscusCartilage",
  formula = "-1.5 + 0.02 * Age + 0.04 * Male + 0.06 * BMI",
  dist = "binary",
  link = "logit"
)

SHARKSDef <- defData(
  SHARKSDef,
  varname = "ikdc_scale",
  formula = "log(72 - 0.25 * Age + 4 * Male - 0.4 * BMI - 12 * MeniscusCartilage)",
  variance = 0.12,  # Controls spread
  dist = "gamma",
  link = "log"
)

SHARKSDef <- defData(
  SHARKSDef,
  varname = "IKDC",
  formula = "pmax(pmin(ikdc_scale, 100), 0)",
  dist = "nonrandom"
)


SHARKSData <- genData(
  n = 500,
  SHARKSDef
)


```

A quick tabulation

```{r}
#| label: fig-SHARKSIKDC
#| fig-cap: "Distribution of simulated IKDC at baseline assessment for SHARKS"

FigSHARKSIKDC <- SHARKSData |> ggplot(aes(x = IKDC)) +
  stat_halfeye()

knitr::knit_print(FigSHARKSIKDC)

```

Let's summarise the dataset so far

```{r}

TableLocal <- tbl_summary(
  SHARKSData |> dplyr::select(
    Age,
    BMI,
    Male,
    MeniscusCartilage,
    IKDC
  ),
  statistic = list(
    all_continuous() ~ "{mean} ({sd})"
  )
)

knitr::knit_print(TableLocal)


```

Now we have a local dataset - we can construct a broader dataset spanning other baseline reports

::: {#tbl-benchmarks}
  | Source | Sample | IKDC Summary | Age | BMI | Comments |
    |------------|------------|------------|------------|------------|------------|
    | SHARKS | 500 | 47.5 (12.6) | 27.3 (7.8) | 28.2 (5.8) | Simulated from SHARKS registry and Ting et al 2022 |
    | [@day2021] | 1126 | 47.9 (16.7) | 30.6 (12.6) | 25.0 (4.1) | Specialty MSK academic institution |
    |  |  |  |  |  | UK National Ligament Registry (KOOS) |
    |  |  |  |  |  | NZ ACL Registry (KOOS) |
    |  | 3687 | 51.1 (17.1) |  |  | MOON cohort |
    |  | 53 | 68.1 (11.4) |  |  | Delaware-Oslo cohort |
    
    Caption
  :::
    
    ## Simulation Study
    
    ```{r}
  #| label: simulation-setup-1
  
  
  # Case-Mix Adjustment Simulation Study - Tidyverse with Epoxy and gtsummary
  # Comparing local vs literature baseline patient-reported outcomes
  
  # Set seed for reproducibility
  set.seed(42)
  
  # =============================================================================
  # DEFINE DATA GENERATING MECHANISMS AND STUDY CHARACTERISTICS
  # =============================================================================
  
  # Define study types and their characteristics
  study_types <- tibble(
    study_type = c("RCT", "Registry", "Observational_Cohort"),
    selection_bias_effect = c(0, -2, -1),  # RCT least biased
    measurement_precision = c(0.8, 0.9, 0.7),  # Registry most precise
    missing_data_rate = c(0.05, 0.15, 0.20)
  )
  
  # Define data capture methods
  capture_methods <- tibble(
    capture_method = c("Electronic", "Paper"),
    measurement_error_sd = c(2, 4),  # Electronic more precise
    processing_delay_bias = c(0, 0.5)
  )
  
  # Define healthcare setting types
  setting_types <- tibble(
    setting_type = c("Academic_Medical_Center", "Community_Hospital", 
                     "Specialty_Clinic", "Primary_Care"),
    resource_level = c(4, 2, 3, 1),  # 1=low, 4=high resources
    case_complexity_mean = c(70, 50, 60, 40),  # baseline severity
    staff_training_effect = c(2, -1, 1, -2)
  )
  
  ```
  
  ```{r}
  #| label: simulation-setup-2
  
  # =============================================================================
  # DEFINE SIMULATION PARAMETERS
  # =============================================================================
  
  # Sample sizes for different scenarios
  n_local <- 500      # Your local dataset
  n_literature <- 2000 # Combined literature studies
  
  # Patient characteristics definitions using simstudy
  def_patient <- defData(varname = "age", dist = "normal", formula = 65, variance = 200) %>%
    defData(varname = "gender", dist = "binary", formula = 0.6) %>%  # 60% female
    defData(varname = "bmi", dist = "normal", formula = 28, variance = 25) %>%
    defData(varname = "comorbidity_count", dist = "poisson", formula = 2) %>%
    defData(varname = "ses_score", dist = "normal", formula = 50, variance = 100) %>%  # socioeconomic status
    defData(varname = "education_years", dist = "normal", formula = 14, variance = 9) %>%
    defData(varname = "disease_duration_months", dist = "gamma", formula = 24, variance = 300)
  
  # Healthcare setting effects
  def_setting <- defData(varname = "setting_id", dist = "categorical", 
                         formula = "0.3;0.3;0.2;0.2")  # distribution across settings
  
  ```
  
  ```{r}
  #| label: simulation-setup-3
  
  
  # =============================================================================
  # SIMULATION FUNCTIONS
  # =============================================================================
  
  # Function to generate patient-reported outcome based on characteristics
  generate_pro_score <- function(data, study_params, setting_params, capture_params) {
    
    # Base PRO score (e.g., 0-100 scale, higher = better quality of life)
    base_score <- 50
    
    data %>%
      mutate(
        # Patient characteristic effects
        age_effect = -0.2 * (age - 65),  # older patients score lower
        gender_effect = if_else(gender == 1, -3, 0),  # females report lower scores
        bmi_effect = -0.1 * pmax(0, bmi - 25),  # obesity effect
        comorbidity_effect = -4 * comorbidity_count,
        ses_effect = 0.1 * (ses_score - 50),
        education_effect = 0.5 * (education_years - 12),
        duration_effect = -0.05 * disease_duration_months,
        
        # Study design effects
        selection_bias = study_params$selection_bias_effect,
        
        # Healthcare setting effects
        setting_effect = setting_params$staff_training_effect,
        complexity_adjustment = -0.1 * (setting_params$case_complexity_mean - 50),
        
        # Data capture effects
        measurement_error = rnorm(n(), 0, capture_params$measurement_error_sd),
        processing_bias = capture_params$processing_delay_bias,
        
        # Calculate final PRO score
        pro_score_raw = base_score + age_effect + gender_effect + bmi_effect + 
          comorbidity_effect + ses_effect + education_effect + 
          duration_effect + selection_bias + setting_effect + 
          complexity_adjustment + measurement_error + processing_bias,
        
        # Bound scores between 0 and 100
        pro_score = pmax(0, pmin(100, pro_score_raw))
      ) %>%
      select(-ends_with("_effect"), -selection_bias, -setting_effect, 
             -complexity_adjustment, -measurement_error, -processing_bias, -pro_score_raw)
  }
  
  # Function to introduce missing data based on study characteristics
  introduce_missing_data <- function(data, missing_rate) {
    n_missing <- round(nrow(data) * missing_rate)
    
    if (n_missing > 0) {
      # Missing not at random - sicker patients more likely to have missing data
      data %>%
        mutate(
          missing_prob = plogis(-2 + 0.02 * (50 - pro_score)),
          row_id = row_number()
        ) %>%
        slice_sample(n = n_missing, weight_by = missing_prob) %>%
        pull(row_id) -> missing_indices
      
      data$pro_score[missing_indices] <- NA
    }
    
    return(data)
  }
  
  # =============================================================================
  # GENERATE SIMULATED DATASETS
  # =============================================================================
  
  generate_study_dataset <- function(n, study_type, setting_type, capture_method, 
                                     dataset_label) {
    
    # Get parameters for this study configuration
    study_params <- study_types %>% filter(study_type == !!study_type)
    setting_params <- setting_types %>% filter(setting_type == !!setting_type)
    capture_params <- capture_methods %>% filter(capture_method == !!capture_method)
    
    # Generate patient characteristics
    patients <- genData(n, def_patient) %>% as_tibble()
    
    # Adjust patient characteristics based on study type and setting
    if (study_type == "RCT") {
      # RCTs often have stricter inclusion criteria
      patients <- patients %>%
        filter(age >= 18, age <= 80, comorbidity_count <= 3)
    }
    
    if (setting_type == "Academic_Medical_Center") {
      # Academic centers often see more complex cases
      patients <- patients %>%
        mutate(
          comorbidity_count = comorbidity_count + rpois(n(), 1),
          disease_duration_months = disease_duration_months * 1.5
        )
    }
    
    # Generate PRO scores
    patients <- generate_pro_score(patients, study_params, setting_params, capture_params)
    
    # Introduce missing data
    patients <- introduce_missing_data(patients, study_params$missing_data_rate)
    
    # Add study metadata
    patients %>%
      mutate(
        study_type = study_type,
        setting_type = setting_type,
        capture_method = capture_method,
        dataset_label = dataset_label,
        study_id = paste(dataset_label, row_number(), sep = "_")
      )
  }
  
  # =============================================================================
  # SIMULATE LOCAL DATASET
  # =============================================================================
  
  # Your local dataset - assume single setting, single study design
  local_data <- generate_study_dataset(
    n = n_local,
    study_type = "Observational_Cohort",  # Typical for real-world data
    setting_type = "Community_Hospital",   # Specify your setting
    capture_method = "Electronic",         # Modern EMR capture
    dataset_label = "Local"
  )
  
  # Create summary using epoxy
  local_summary <- local_data %>%
    summarise(
      n = n(),
      mean_pro = round(mean(pro_score, na.rm = TRUE), 2),
      sd_pro = round(sd(pro_score, na.rm = TRUE), 2)
    )
  
  # =============================================================================
  # SIMULATE LITERATURE DATASETS
  # =============================================================================
  
  # Create multiple literature studies with different characteristics
  literature_studies <- tibble(
    n = c(800, 600, 400, 200),
    study_type = c("RCT", "Registry", "Observational_Cohort", "Observational_Cohort"),
    setting_type = c("Academic_Medical_Center", "Specialty_Clinic", 
                     "Community_Hospital", "Primary_Care"),
    capture_method = c("Electronic", "Electronic", "Paper", "Paper"),
    dataset_label = c("Literature_RCT_Academic", "Literature_Registry_Specialty",
                      "Literature_Obs_Community", "Literature_Obs_Primary")
  )
  
  # Generate all literature datasets using pmap
  literature_data <- literature_studies %>%
    pmap_dfr(generate_study_dataset)
  
  # Create literature summary
  literature_summary <- literature_data %>%
    summarise(
      n = n(),
      mean_pro = round(mean(pro_score, na.rm = TRUE), 2),
      sd_pro = round(sd(pro_score, na.rm = TRUE), 2)
    )
  
  # =============================================================================
  # COMBINE AND ANALYZE DATASETS
  # =============================================================================
  
  # Combine all data for analysis
  all_data <- bind_rows(local_data, literature_data)
  
  # Calculate summary statistics by dataset type using gtsummary
  summary_table <- all_data %>%
    select(dataset_label, study_type, setting_type, capture_method, 
           pro_score, age, gender, bmi, comorbidity_count, ses_score) %>%
    mutate(
      gender = factor(gender, levels = c(0, 1), labels = c("Male", "Female")),
      study_type = factor(study_type),
      setting_type = factor(setting_type),
      capture_method = factor(capture_method)
    ) %>%
    tbl_summary(
      by = dataset_label,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 2,
      label = list(
        pro_score ~ "PRO Score",
        age ~ "Age (years)",
        gender ~ "Gender",
        bmi ~ "BMI (kg/m²)",
        comorbidity_count ~ "Comorbidity Count",
        ses_score ~ "SES Score",
        study_type ~ "Study Type",
        setting_type ~ "Setting Type",
        capture_method ~ "Capture Method"
      )
    ) %>%
    add_overall() %>%
    modify_header(label ~ "**Variable**") %>%
    modify_caption("**Study Characteristics by Dataset**")
  
  # =============================================================================
  # CASE-MIX ADJUSTMENT ANALYSIS EXAMPLE
  # =============================================================================
  
  # Prepare data for case-mix adjustment
  analysis_data <- all_data %>%
    filter(!is.na(pro_score)) %>%
    mutate(is_local = if_else(dataset_label == "Local", 1, 0))
  
  # Unadjusted comparison
  unadjusted_comparison <- analysis_data %>%
    group_by(is_local) %>%
    summarise(mean_pro = mean(pro_score), .groups = 'drop')
  
  unadjusted_diff <- unadjusted_comparison %>%
    pivot_wider(names_from = is_local, values_from = mean_pro, names_prefix = "group_") %>%
    mutate(difference = group_1 - group_0) %>%
    pull(difference)
  
  # Case-mix adjusted comparison using linear regression
  casemix_model <- lm(pro_score ~ is_local + age + gender + bmi + 
                        comorbidity_count + ses_score + education_years + 
                        disease_duration_months, 
                      data = analysis_data)
  
  # Extract adjusted difference and confidence interval using broom
  model_results <- casemix_model %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "is_local")
  
  # Create model results table using gtsummary
  model_table <- casemix_model %>%
    tbl_regression(
      label = list(
        is_local ~ "Local vs Literature",
        age ~ "Age (years)",
        gender ~ "Gender (Female)",
        bmi ~ "BMI (kg/m²)",
        comorbidity_count ~ "Comorbidity Count",
        ses_score ~ "SES Score",
        education_years ~ "Education (years)",
        disease_duration_months ~ "Disease Duration (months)"
      )
    ) %>%
    add_glance_source_note(
      label = list(r.squared ~ "R²", AIC ~ "AIC"),
      fmt_fun = r.squared ~ function(x) style_number(x, digits = 3)
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_caption("**Case-Mix Adjustment Model Results**")
  
  adjusted_diff <- model_results$estimate
  adjusted_ci <- c(model_results$conf.low, model_results$conf.high)
  
  # =============================================================================
  # CASE-MIX CORRECTION MODELING AND APPLICATION
  # =============================================================================
  
  # Step 1: Build comprehensive case-mix adjustment model using literature data only
  literature_model_data <- analysis_data %>%
    filter(is_local == 0)  # Literature data only
  
  # Fit comprehensive case-mix model on literature data
  literature_casemix_model <- lm(
    pro_score ~ age + I(age^2) + gender + bmi + I(bmi^2) + 
      comorbidity_count + I(comorbidity_count^2) + 
      ses_score + education_years + disease_duration_months +
      # Interaction terms for key variables
      age:gender + bmi:comorbidity_count + ses_score:education_years +
      # Study design effects
      study_type + setting_type + capture_method,
    data = literature_model_data
  )
  
  # Model performance
  literature_model_summary <- literature_casemix_model %>%
    glance() %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))
  
  # Create literature model table
  literature_model_table <- literature_casemix_model %>%
    tbl_regression(
      label = list(
        age ~ "Age (years)",
        gender ~ "Gender (Female)",
        bmi ~ "BMI (kg/m²)",
        comorbidity_count ~ "Comorbidity Count",
        ses_score ~ "SES Score",
        education_years ~ "Education (years)",
        disease_duration_months ~ "Disease Duration (months)",
        study_type ~ "Study Type",
        setting_type ~ "Setting Type",
        capture_method ~ "Capture Method"
      )
    ) %>%
    add_glance_source_note(
      label = list(r.squared ~ "R²", AIC ~ "AIC"),
      fmt_fun = r.squared ~ function(x) style_number(x, digits = 3)
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_caption("**Literature-Based Case-Mix Adjustment Model**")
  
  # Step 2: Calculate expected scores for literature population
  # Get literature population means for standardization
  literature_population_means <- literature_model_data %>%
    summarise(
      age_pop = mean(age),
      gender_pop = mean(gender),
      bmi_pop = mean(bmi),
      comorbidity_count_pop = mean(comorbidity_count),
      ses_score_pop = mean(ses_score),
      education_years_pop = mean(education_years),
      disease_duration_months_pop = mean(disease_duration_months),
      .groups = 'drop'
    )
  
  # Create standardized reference population
  reference_population <- tibble(
    age = literature_population_means$age_pop,
    gender = literature_population_means$gender_pop,
    bmi = literature_population_means$bmi_pop,
    comorbidity_count = literature_population_means$comorbidity_count_pop,
    ses_score = literature_population_means$ses_score_pop,
    education_years = literature_population_means$education_years_pop,
    disease_duration_months = literature_population_means$disease_duration_months_pop,
    study_type = "Registry",  # Use most common/neutral study type
    setting_type = "Community_Hospital",  # Use most common setting
    capture_method = "Electronic"  # Use most precise method
  )
  
  # Calculate expected score for reference population
  reference_expected_score <- predict(literature_casemix_model, reference_population)
  
  # Step 3: Apply case-mix correction to local data
  local_corrected_data <- local_data %>%
    filter(!is.na(pro_score)) %>%
    mutate(
      # Calculate expected score for each local patient based on their characteristics
      # but using reference study design characteristics
      expected_score_reference = map_dbl(1:n(), function(i) {
        patient_data <- tibble(
          age = age[i],
          gender = gender[i], 
          bmi = bmi[i],
          comorbidity_count = comorbidity_count[i],
          ses_score = ses_score[i],
          education_years = education_years[i],
          disease_duration_months = disease_duration_months[i],
          study_type = "Registry",  # Reference study design
          setting_type = "Community_Hospital",  # Reference setting
          capture_method = "Electronic"  # Reference capture method
        )
        predict(literature_casemix_model, patient_data)
      }),
      
      # Calculate expected score for each local patient with their actual study characteristics
      expected_score_actual = predict(literature_casemix_model, 
                                      select(., age, gender, bmi, comorbidity_count, 
                                             ses_score, education_years, disease_duration_months,
                                             study_type, setting_type, capture_method)),
      
      # Case-mix adjustment: Remove the effect of patient characteristics differences
      # and study design differences
      case_mix_adjustment = expected_score_reference - expected_score_actual,
      
      # Apply correction to observed scores
      pro_score_corrected = pro_score + case_mix_adjustment,
      
      # Ensure corrected scores stay within valid range
      pro_score_corrected = pmax(0, pmin(100, pro_score_corrected))
    )
  
  # Step 4: Calculate correction impact
  correction_summary <- local_corrected_data %>%
    summarise(
      n = n(),
      original_mean = round(mean(pro_score), 2),
      original_sd = round(sd(pro_score), 2),
      corrected_mean = round(mean(pro_score_corrected), 2),
      corrected_sd = round(sd(pro_score_corrected), 2),
      mean_adjustment = round(mean(case_mix_adjustment), 2),
      median_adjustment = round(median(case_mix_adjustment), 2),
      adjustment_range_min = round(min(case_mix_adjustment), 2),
      adjustment_range_max = round(max(case_mix_adjustment), 2),
      literature_reference_mean = round(reference_expected_score, 2)
    )
  
  # Create correction summary table using gt
  correction_table <- tibble(
    Metric = c("Sample Size", "Original Mean (SD)", "Corrected Mean (SD)", 
               "Literature Reference", "Mean Adjustment", "Adjustment Range"),
    Value = c(
      as.character(correction_summary$n),
      glue::glue("{correction_summary$original_mean} ({correction_summary$original_sd})"),
      glue::glue("{correction_summary$corrected_mean} ({correction_summary$corrected_sd})"),
      as.character(correction_summary$literature_reference_mean),
      as.character(correction_summary$mean_adjustment),
      glue::glue("[{correction_summary$adjustment_range_min}, {correction_summary$adjustment_range_max}]")
    )
  ) %>%
    gt() %>%
    tab_header(
      title = "Case-Mix Correction Summary",
      subtitle = "Local Dataset Adjustment Results"
    ) %>%
    cols_label(
      Metric = "Metric",
      Value = "Value"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    )
  
  # =============================================================================
  # VISUALIZATION: ORIGINAL vs CORRECTED LOCAL DATA
  # =============================================================================
  
  # Prepare data for visualization
  plot_data <- local_corrected_data %>%
    select(study_id, pro_score, pro_score_corrected, case_mix_adjustment) %>%
    pivot_longer(cols = c(pro_score, pro_score_corrected), 
                 names_to = "score_type", 
                 values_to = "score") %>%
    mutate(
      score_type = case_when(
        score_type == "pro_score" ~ "Original",
        score_type == "pro_score_corrected" ~ "Case-Mix Corrected"
      ),
      score_type = factor(score_type, levels = c("Original", "Case-Mix Corrected"))
    )
  
  # Add literature reference line data
  literature_ref_mean <- reference_expected_score
  
  # Plot 1: Density comparison
  p3 <- plot_data %>%
    ggplot(aes(x = score, fill = score_type)) +
    geom_density(alpha = 0.7) +
    geom_vline(xintercept = literature_ref_mean, linetype = "dashed", 
               color = "black", size = 1) +
    annotate("text", x = literature_ref_mean + 5, y = 0.06, 
             label = paste("Literature Reference\n(", round(literature_ref_mean, 1), ")"),
             hjust = 0, size = 3) +
    scale_fill_manual(values = c("Original" = "red", "Case-Mix Corrected" = "darkgreen")) +
    labs(title = "Local PRO Scores: Original vs Case-Mix Corrected",
         subtitle = "Comparison with Literature Reference Population",
         x = "Patient-Reported Outcome Score",
         y = "Density",
         fill = "Score Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 2: Before/After scatter plot
  p4 <- local_corrected_data %>%
    ggplot(aes(x = pro_score, y = pro_score_corrected)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = literature_ref_mean, linetype = "dotted", color = "black") +
    geom_vline(xintercept = literature_ref_mean, linetype = "dotted", color = "black") +
    annotate("text", x = 20, y = 90, 
             label = "Identity Line\n(no change)", color = "red", size = 3) +
    annotate("text", x = 80, y = literature_ref_mean + 5, 
             label = "Literature Reference", color = "black", size = 3) +
    labs(title = "Individual Patient Score Changes",
         subtitle = "Original vs Case-Mix Corrected PRO Scores",
         x = "Original PRO Score",
         y = "Case-Mix Corrected PRO Score") +
    theme_minimal() +
    coord_fixed()
  
  # Plot 3: Adjustment magnitude by patient characteristics
  p5 <- local_corrected_data %>%
    ggplot(aes(x = age, y = case_mix_adjustment, color = factor(gender))) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("0" = "blue", "1" = "red"),
                       labels = c("Male", "Female")) +
    labs(title = "Case-Mix Adjustment by Patient Age and Gender",
         subtitle = "Positive values indicate upward correction",
         x = "Age (years)",
         y = "Case-Mix Adjustment (points)",
         color = "Gender") +
    theme_minimal()
  
  # Plot 4: Box plot comparison with literature
  comparison_plot_data <- bind_rows(
    local_corrected_data %>%
      select(pro_score, pro_score_corrected) %>%
      pivot_longer(cols = everything(), names_to = "score_type", values_to = "score") %>%
      mutate(
        dataset = case_when(
          score_type == "pro_score" ~ "Local (Original)",
          score_type == "pro_score_corrected" ~ "Local (Corrected)"
        )
      ),
    literature_model_data %>%
      select(pro_score) %>%
      mutate(dataset = "Literature", score_type = "literature", score = pro_score) %>%
      select(-pro_score)
  ) %>%
    mutate(dataset = factor(dataset, levels = c("Literature", "Local (Original)", "Local (Corrected)")))
  
  p6 <- comparison_plot_data %>%
    ggplot(aes(x = dataset, y = score, fill = dataset)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    scale_fill_manual(values = c("Literature" = "blue", 
                                 "Local (Original)" = "red", 
                                 "Local (Corrected)" = "darkgreen")) +
    labs(title = "PRO Score Distributions: Literature vs Local Data",
         subtitle = "Effect of Case-Mix Correction on Local Data",
         x = "Dataset",
         y = "PRO Score") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Final comparison statistics table
  final_comparison_table <- tibble(
    Dataset = c("Literature", "Local (Original)", "Local (Corrected)"),
    Mean = c(
      round(mean(literature_model_data$pro_score), 2),
      correction_summary$original_mean,
      correction_summary$corrected_mean
    ),
    SD = c(
      round(sd(literature_model_data$pro_score), 2),
      correction_summary$original_sd,
      correction_summary$corrected_sd
    ),
    `Difference from Literature` = c(
      0,
      round(correction_summary$original_mean - mean(literature_model_data$pro_score), 2),
      round(correction_summary$corrected_mean - mean(literature_model_data$pro_score), 2)
    )
  ) %>%
    gt() %>%
    tab_header(
      title = "Final Comparison: PRO Scores",
      subtitle = "Literature vs Local Data (Original and Corrected)"
    ) %>%
    cols_label(
      Dataset = "Dataset",
      Mean = "Mean",
      SD = "SD",
      `Difference from Literature` = "Difference from Literature"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    data_color(
      columns = `Difference from Literature`,
      domain = c(-10, 10),
      palette = c("red", "white", "blue"),
      na_color = "white"
    )
  
  # Save corrected dataset for further analysis
  local_data_with_correction <- local_corrected_data
  
  ```
  
  
  ```{r}
  
  # Create summary using epoxy
  epoxy_html(
    .data = list(
      local_n = nrow(local_data),
      lit_n = nrow(literature_data),
      orig_mean = correction_summary$original_mean,
      corr_mean = correction_summary$corrected_mean,
      lit_mean = round(reference_expected_score, 2),
      improvement = round(abs(correction_summary$original_mean - reference_expected_score) - 
                            abs(correction_summary$corrected_mean - reference_expected_score), 2)
    ),
    "
  ## Case-Mix Adjustment Simulation Results
  
  ### Dataset Generation Summary
  - **Local dataset:** {{local_n}} patients generated
  - **Literature datasets:** {{lit_n}} patients across multiple studies
  
  ### Key Findings
  - **Original local mean PRO score:** {{orig_mean}}
  - **Case-mix corrected local mean:** {{corr_mean}}  
  - **Literature reference mean:** {{lit_mean}}
  - **Improvement in alignment:** {{improvement}} points
  
  ### Available Objects
  - `local_data`: Original local dataset
  - `literature_data`: Combined literature studies  
  - `local_corrected_data`: Local data with case-mix corrections
  - `literature_casemix_model`: Case-mix adjustment model
  - `summary_table`: Study characteristics by dataset (gtsummary)
  - `model_table`: Case-mix adjustment model results (gtsummary)
  - `correction_table`: Correction summary (gt)
  - `final_comparison_table`: Final comparison (gt)
  "
  )
  
  # Print tables
  knitr::knit_print(summary_table)
  knitr::knit_print(model_table) 
  knitr::knit_print(literature_model_table)
  knitr::knit_print(correction_table)
  knitr::knit_print(final_comparison_table)
  
  # Print plots
  knitr::knit_print(p3)
  knitr::knit_print(p4)
  knitr::knit_print(p5)
  knitr::knit_print(p6)
  
  ```
  ## Simulation Study v2
  
  From simstudy help (clustered data)
  
  
  ```{r}
  #| label: simulationv2-prep-1
  
  # study type
  # capture method
  # setting type
  # country_region (?)
  
  gen.region <- defData(
    varname = "Region",
    dist = "normal",
    formula = 0,
    variance = 4,
    id = "idRegion"
  )
  
  gen.region <- defData(
    gen.region,
    varname = "nStudy",
    dist = "noZeroPoisson",
    formula = 4
  )
  
  set.seed(2065)
  
  dtRegion <- genData(8, gen.region)
  ```
  
  ```{r}
  #| label: simulationv2-prep-2
  
  gen.study <- defDataAdd(
    varname = "Study", 
    dist = "normal", 
    formula = 0, 
    variance = 4
  )
  
  gen.study <- defDataAdd(
    gen.study,
    varname = "nSample", 
    dist = "noZeroPoisson",
    formula = 200
  )
  
  
  dtStudy <- genCluster(dtRegion, "idRegion", numIndsVar = "nStudy", level1ID = "idStudy")
  dtStudy <- addColumns(gen.study, dtStudy)
  
  ```
  
  
  ```{r}
  #| label: simulationv2-prep-3
  
  
  gen.patient <- defDataAdd(
    varname = "Male",
    dist = "binary",
    formula = 0.5
  )
  
  gen.patient <- defDataAdd(
    varname = "Male", 
    dist = "binary",
    formula = 0.6
  )
  gen.patient <- defDataAdd(
    gen.patient, 
    varname = "age", 
    dist = "normal",
    formula = "27.5",
    variance = "12.5^2"
  )
  gen.patient <- defDataAdd(
    gen.patient, 
    varname = "IKDC", 
    dist = "normal",
    formula = "50 - 5*Male + Study + Region", 
    variance = 2
  )
  
  dtPatient <- genCluster(
    dtStudy, 
    cLevelVar = "idStudy",
    numIndsVar = "nSample",
    level1ID = "idPatient"
  )
  
  dtPatient <- addColumns(
    gen.patient, 
    dtPatient
  )
  
  ```
  
  
  ```{r}
  
  
  SimulatedTable <- gtsummary::tbl_summary(
    dtPatient |> dplyr::select(
      idRegion,
      Male,
      age,
      IKDC
    ),
    by = "idRegion",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})"
    )
  )
  
  knitr::knit_print(SimulatedTable)
  
  ```
  ```{r}
  #| label: fig-sim-IKDC
  #| fig-cap: "Distribution of simulated IKDC at baseline assessment across multiple regions"
  
  dtPatient1 <- dtPatient |> dplyr::mutate(
    idRegion = forcats::as_factor(idRegion)
  )
  
  FigSimIKDC <- dtPatient1 |> ggplot(aes(y = idRegion, x = IKDC)) +
    stat_halfeye()
  
  knitr::knit_print(FigSimIKDC)
  
  ```
  