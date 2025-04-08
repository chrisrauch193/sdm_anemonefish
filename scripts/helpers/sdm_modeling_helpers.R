# scripts/helpers/sdm_modeling_helpers.R

# ... (keep other helper functions like load_clean_individual_occ, thin_individual_occ, generate_sdm_background) ...

#' Run SDM Tuning using SDMtune::gridSearch with Train/Test Split
#'
#' Prepares data, splits into training and testing sets, trains a base model,
#' and executes SDMtune grid search evaluating hyperparameters on the test set.
#' NOTE: This uses a single train/test split, not k-fold CV for tuning,
#' as creating SDMmodelCV with Maxnet seems problematic based on errors.
#'
#' @param occs_thinned_sf sf object of thinned occurrence points.
#' @param predictor_stack SpatRaster object of predictors.
#' @param background_df Data frame with 'x', 'y' columns for background points.
#' @param config Configuration list (for SDM settings, test proportion).
#'
#' @return An SDMtune object containing tuning results, or NULL on error.
#' @export
run_sdm_SDMtune_grid <- function(occs_thinned_sf, predictor_stack, background_df, config) {
  cat("    Preparing data and running SDMtune::gridSearch with Train/Test split...\n") # Updated message
  
  if (is.null(occs_thinned_sf) || nrow(occs_thinned_sf) == 0 || is.null(predictor_stack) || is.null(background_df)) {
    warning("Invalid inputs for running SDMtune.", call. = FALSE); return(NULL)
  }
  
  occs_coords <- sf::st_coordinates(occs_thinned_sf)
  
  # 1. Prepare SWD object
  swd_data <- tryCatch({
    SDMtune::prepareSWD(
      species = "species",
      p = occs_coords,
      a = background_df,
      env = predictor_stack
    )
  }, error = function(e) {
    warning("Failed to prepare SWD object: ", e$message, call. = FALSE); return(NULL)
  })
  
  if (is.null(swd_data)) return(NULL)
  
  # 2. *** SPLIT Data into Training and Testing ***
  # Use a specified proportion for testing (e.g., 20% = 0.2)
  test_proportion <- 1 / config$sdm_n_folds # Use k from config to define test proportion
  cat("    Splitting data: ", (1-test_proportion)*100, "% Train,", test_proportion*100, "% Test (presence only).\n")
  datasets <- tryCatch({
    SDMtune::trainValTest(swd_data, test = test_proportion, only_presence = TRUE, seed = 123) # Added seed for reproducibility
  }, error = function(e) {
    warning("Failed to split data using trainValTest: ", e$message, call. = FALSE); return(NULL)
  })
  
  if (is.null(datasets)) return(NULL)
  train_swd <- datasets[[1]]
  test_swd <- datasets[[2]]
  
  # 3. Train an initial base model on the TRAINING data
  cat("    Training initial base model on training data...\n")
  initial_model <- tryCatch({
    SDMtune::train(method = config$sdm_method, # e.g., "Maxnet"
                   data = train_swd
                   # Use default hyperparameters for the initial model
    )
  }, error = function(e) {
    warning("Failed to train initial base model: ", e$message, call. = FALSE); return(NULL)
  })
  
  if (is.null(initial_model) || !inherits(initial_model, "SDMmodel")) {
    warning("Failed to create a valid SDMmodel object.", call. = FALSE); return(NULL)
  }
  
  # 4. Define the hyperparameter grid
  hyper_grid <- config$sdm_tune_grid
  
  # 5. Run gridSearch using the initial model and the TEST dataset
  tuned_results <- NULL
  tryCatch({
    cat("    Running gridSearch, evaluating on the test set...\n")
    tuned_results <- SDMtune::gridSearch(
      model = initial_model, # Pass the base SDMmodel trained on train_swd
      hypers = hyper_grid,
      metric = config$sdm_evaluation_metric, # Metric to evaluate on the test set
      test = test_swd,        # Provide the test SWD object
      save_models = TRUE      # Save models trained on training data with different hypers
    )
    cat("    SDMtune::gridSearch with train/test split completed.\n")
    # print(SDMtune::results(tuned_results)) # Optional: view results table
  }, error = function(e) {
    warning("SDMtune::gridSearch failed: ", e$message, call. = FALSE)
    tuned_results <- NULL
  })
  
  return(tuned_results)
}


#' Select Best Model and Predict SDM using SDMtune (Train/Test Version)
#'
#' Selects the best model based on test set performance from gridSearch results
#' and generates a prediction raster using a model retrained on ALL data
#' with the best hyperparameters.
#'
#' @param SDMtune_results The output object from SDMtune::gridSearch (train/test version).
#' @param predictor_stack SpatRaster object used for prediction.
#' @param config Configuration list.
#'
#' @return A SpatRaster object with the prediction, or NULL on error.
#' @export
predict_sdm_SDMtune <- function(SDMtune_results, predictor_stack, config) {
  cat("    Selecting best hyperparameters and predicting (retraining on full dataset)...\n")
  
  if (is.null(SDMtune_results) || !inherits(SDMtune_results, "SDMtune")) {
    warning("Invalid SDMtune results object provided for prediction.", call. = FALSE); return(NULL)
  }
  if (is.null(predictor_stack)) {
    warning("Predictor stack is required for prediction.", call. = FALSE); return(NULL)
  }
  
  best_hypers <- NULL
  best_model_from_tuning <- NULL
  tryCatch({
    res_df <- SDMtune::results(SDMtune_results)
    if(nrow(res_df) > 0 && length(SDMtune_results@models) > 0) {
      # Identify best hyperparameters based on the test metric (e.g., 'test_AUC')
      metric <- config$sdm_evaluation_metric
      metric_test_name <- paste0("test_", metric) # Name used in results() output
      
      if (!metric_test_name %in% colnames(res_df)) {
        # Fallback if exact name isn't present
        if ("test_AUC" %in% colnames(res_df)) metric_test_name <- "test_AUC"
        else if ("test_TSS" %in% colnames(res_df)) metric_test_name <- "test_TSS"
        else { metric_test_name <- NA }
        
        if(is.na(metric_test_name)){
          warning("Cannot find test metric '", paste0("test_", metric), "' or fallbacks in results. Using first hyperparameter set.", call.=FALSE)
          best_row_index <- 1
        } else {
          warning("Test metric '", paste0("test_", metric), "' not found, using '", metric_test_name, "' instead for selection.", call.=FALSE)
          best_row_index <- which.max(res_df[[metric_test_name]])
        }
      } else {
        best_row_index <- which.max(res_df[[metric_test_name]])
      }
      
      if(length(best_row_index) == 0 || best_row_index > nrow(res_df)) {
        warning("Could not determine best model index from results. Using first hyperparameter set.", call.=FALSE)
        best_row_index <- 1
      }
      
      # Extract the best hyperparameter combination
      best_hypers <- SDMtune_results@hypers[best_row_index, , drop = FALSE] # Keep as data frame row
      cat("    Best hyperparameters found (via", metric_test_name, "):", paste(names(best_hypers), best_hypers[1,], collapse=", "), "\n")
      # Get the corresponding model trained during gridSearch (it was trained on the training set)
      best_model_from_tuning <- SDMtune_results@models[[best_row_index]]
      
    } else {
      warning("No models or results found within the SDMtune results object.", call.=FALSE)
      return(NULL) # Cannot proceed without results
    }
  }, error = function(e) {
    warning("Error identifying best model hyperparameters: ", e$message, call. = FALSE)
    return(NULL) # Cannot proceed
  })
  
  if (is.null(best_hypers) || is.null(best_model_from_tuning)) {
    warning("Failed to retrieve best hyperparameters or model from tuning.", call.=FALSE)
    return(NULL)
  }
  
  # *** RETRAIN the model on the FULL dataset using the best hyperparameters ***
  # This is standard practice after tuning on a subset
  cat("    Retraining model on full dataset with best hyperparameters...\n")
  full_data <- SDMtune_results@data # Access the original full SWD object stored by gridSearch
  final_model <- tryCatch({
    SDMtune::train(method = config$sdm_method,
                   data = full_data,
                   # Apply the best hyperparameters found - need to pass them as named args
                   fc = best_hypers$fc, # Example, adjust if names differ
                   reg = best_hypers$reg   # Example, adjust if names differ
                   # Add other hypers as needed, extracting from best_hypers
    )
  }, error = function(e) {
    warning("Failed to retrain final model on full dataset: ", e$message, call.=FALSE); return(NULL)
  })
  
  if (is.null(final_model)) return(NULL)
  
  # Predict using the final model trained on all data
  prediction_raster <- NULL
  tryCatch({
    prediction_raster <- SDMtune::predict(
      object = final_model,    # Use the model retrained on full data
      data = predictor_stack,
      type = "cloglog"        # Or other relevant type
    )
    names(prediction_raster) <- "suitability"
    cat("    Prediction raster generated from final model.\n")
  }, error = function(e) {
    warning("SDMtune::predict failed on final model: ", e$message, call. = FALSE)
    prediction_raster <- NULL
  })
  
  return(prediction_raster)
}

#-------------------------------------------------------------------------------