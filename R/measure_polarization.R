#' Load the trained models and parameters
#'
#' Load the trained models and parameters for measuring the polarization in a cell type.
#'
#' @param celltype A string of cell type whose polarization will be measured.
#' Must be one of: 'B', 'NK', 'CD8T', 'CD4T', 'Treg', 'Tgd', 'pDC',
#' 'cDC1', 'cDC2', 'MigDC', 'LC', 'Macro', 'Mono', 'Neu'.
#' @param verbose Whether to message the loading progress.
#' Default: TRUE.
#'
#' @return A list containing 7 variables:
#' \describe{
#'   \item{mean}{The mean of Universial Cell Embeddings of all cells in the training
#'   dataset.}
#'   \item{sd}{The standard deviation of Universial Cell Embeddings of all cells
#'   in the training dataset.}
#'   \item{mean_unpolar}{The mean of Universial Cell Embeddings of unpolarized cells
#'   in the training dataset.}
#'   \item{pc_loadings}{The principal component loadings to transform Universial
#'   Cell Embeddings into principal component embeddings.}
#'   \item{models}{The trained machine learning models to predict the polarization
#'   scores based on principal component embeddings for each cell polarization
#'   state.}
#'   \item{unpolar_responses}{The polarization scores of unpolarized cells in the
#'   training dataset. The probability of input data being unpolarized is calculated
#'   by comparing to the distribution of these responses.}
#'   \item{calibration_quantiles}{The 1~99 quantiles of nonconformity scores in
#'   calibration data, used for conformal prediction.}
#' }
#'
#' @details This function will load a global variable named 'polar_params'
#' containing trained machine learning models and parameters.
#'
#' @importFrom utils data
#' @export
#'
#' @examples
#' params_cd8t <- LoadPolarParams('CD8T')
#'
LoadPolarParams <- function(celltype,
                            verbose = TRUE)
{
  if (!exists('polar_params'))
  {
    if (verbose)
    {
      message('No model and parameter found. Loading models and parameters from the saved data.')
    }
    data(polar_params)
  }
  if (!celltype %in% c(
    'B',
    'NK',
    'CD8T',
    'CD4T',
    'Treg',
    'Tgd',
    'pDC',
    'cDC1',
    'cDC2',
    'MigDC',
    'LC',
    'Macro',
    'Mono',
    'Neu'
  ))
  {
    stop(
      'celltype must be one of:
      B, NK, CD8T, CD4T, Treg, Tgd, pDC, cDC1, cDC2, MigDC, LC, Macro, Mono, Neu.'
    )
  }

  if (verbose)
  {
    message(paste('Loaded models and parameters for', celltype))
  }

  return(polar_params[[celltype]])
}

#' Measure polarization of immune cells
#'
#' This is the main function to measure the polarization of immune cells based on
#' their Universial Cell Embeddings. The function uses trained machine learning
#' models to predict the polarization of cells in the input data. Conformal
#' prediction is employed to make statistically valid predictions.
#'
#' @param object A Seurat object with an assay of UCE embedding for immune cell
#' polarization measurement. It is recommended that the object only contain one
#' cell type of interest.
#' @param celltype A string of cell type whose polarization will be measured.
#' Must be one of: B, NK, CD8T, CD4T, Treg, Tgd, pDC,
#' cDC1, cDC2, MigDC, LC, Macro, Mono, Neu.
#' @param embedding The name of dimension reduction or assay saving the cell
#' embeddings from the single-cell foundation model in the Seurat object.
#' Default: uce.
#' @param unpolarized_cell The name of cells that are considered unpolarized cells.
#' See details for selecting unpolarized cells. Default: NA.
#' @param return.df Whether to return the results as a dataframe. If TRUE,
#' return a dataframe. Else, return the Seurat object with updated meta.data.
#' Default: FALSE.
#' @param pc The index of principal component to be used by machine learning models.
#' Default: 1:20. DO NOT change this parameter if you uses the trained model provided
#' in the package! Only change it when you have trained custom models using
#' \code{CalculateParams}.
#' @param error_level The error level for conformal prediction. Increasing this
#' parameter will make the prediction less likely to be both classes ('Intermediate'),
#' but also increase the chance of empty class ('Uncertain'). Default: 0.05.
#' @param unpolarized_prob Whether to calculate the estimated probability of input
#' data being unpolarized. The values are calculated by comparing the predicted
#' polarization scores to unpolarized cells' distribution of polarization scores.
#' Default: FALSE.
#' @param verbose Whether to message the progress in different polarization states.
#' Default: TRUE.
#'
#' @return Return either a dataframe or a Seurat object with updated meta.data,
#' depending on param return.df.
#' By default, the dataframe includes 2 * #(polarization states) columns.
#' For each cell type, there are 4~6 polarization states.
#' For each polarization state, there is a column of polarization scores (range 0~1)
#' and a column of state predictions (‘Polarized’, ‘Intermediate’, ‘Unpolarized’,
#' or ‘Uncertain’) for all cells.
#' If parameter \code{unpolarized_prob} is TRUE, an additional column of probability
#' being unpolarized for each polarization state will be added to the dataframe.
#' For details of polarization states, please refer to \code{data(polar_states)}
#' and the reference (Cui et al. Nature. 2024).
#'
#' @details Providing the unpolarized cell names is optional, but highly recommended.
#' The function will correct cell embedding batch effect by mapping the center
#' of unpolarized cells in the input dataset to the center of unpolarized cells
#' in the training dataset, thereby improving the measurement of polarization scores.
#' Here are some criteria for selecting unpolarized cells:
#' 1. When there are untreated samples and samples treated with certain stimulation
#' like cytokines or drugs, use the cells in untreated samples.
#' 2. When there is a cell clusters with high expression of genes related to
#' steady state, like C1Q-high macrophages, use the cells in that cell cluster.
#'
#' If not provided, the cell embedding batch effect between the training dataset
#' and the input dataset will not be corrected. Despite that, we found that the
#' method could generate results similar as the results when providing unpolarized
#' cells in several independent datasets, potentially benefiting from UCE's
#' inherent ability in reducing technical variations.
#'
#' @references
#' Liu W, Zhao Z. Scupa: Single-cell unified polarization assessment of immune
#' cells using the single-cell foundation model. bioRxiv. 2024. doi:
#' 10.1101/2024.08.15.608093
#'
#' Cui A, Huang T, Li S, et al. Dictionary of immune responses to
#' cytokines at single-cell resolution. Nature. 2024;625(7994):377-384.
#' doi:10.1038/s41586-023-06816-9
#'
#' Rosen Y, Roohani Y, et al. Universal Cell Embeddings: A Foundation Model for
#' Cell Biology. BioRxiv. 2023. doi: 10.1101/2023.11.28.568918
#'
#' @importFrom stats glm p.adjust predict quantile sd
#' @importFrom e1071 svm
#' @importFrom Seurat AddMetaData Reductions Embeddings Assays
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load a Seurat object with UCE embeddings
#' ifnb <- readRDS('ifnb_treatment.rds')
#' # subset the object to get CD8+ T cells
#' cd8t <- subset(ifnb, seurat_annotations == 'CD8 T')
#' # Use the function to measure CD8+ T cell polarization
#' # Cells from the control samples are used as unpolarized cells.
#' cd8t <- MeasurePolar(cd8t, celltype='CD8T',
#' unpolarized_cell=WhichCells(cd8t, expression = stim == 'CTRL'))
#' # Plot the CD8+ T cell polarization scores
#' FeaturePlot(cd8t, c('T8.a_score','T8.b_score','T8.c_score','T8.e_score','T8.f_score'))
#' # Plot the CD8+ T cell polarization predictions
#' FeaturePlot(cd8t, c('T8.a_pred','T8.b_pred','T8.c_pred','T8.e_pred','T8.f_pred'))
#' }
#'
MeasurePolar <- function(object,
                         celltype,
                         embedding = 'uce',
                         unpolarized_cell = NA,
                         return.df = FALSE,
                         pc = 1:20,
                         error_level = 0.05,
                         unpolarized_prob = FALSE,
                         verbose = TRUE)
{
  # Validate celltype
  if (!celltype %in% c(
    'B',
    'NK',
    'CD8T',
    'CD4T',
    'Treg',
    'Tgd',
    'pDC',
    'cDC1',
    'cDC2',
    'MigDC',
    'LC',
    'Macro',
    'Mono',
    'Neu'
  ))
  {
    stop(
      'celltype must be one of:
      B, NK, CD8T, CD4T, Treg, Tgd, pDC, cDC1, cDC2, MigDC, LC, Macro, Mono, Neu.'
    )
  }
  # Load the saved machine learning models and parameters
  saved_params <- LoadPolarParams(celltype,
                                  verbose = verbose)
  if(embedding %in% Reductions(object))
  {
    input_emb <- t(Embeddings(object, reduction=embedding))
  }
  else if(embedding %in% Assays(object))
  {
    input_emb <- object[[embedding]]@data
  }
  else if('Xuce_' %in% Reductions(object))
  {
    message('Embedding', embedding, 'not found in Reductions and Assays. Using Xuce_ in Reductions!')
    input_emb <- t(Embeddings(object, reduction='Xuce_'))
  }
  else
    stop(
      paste('Embedding', embedding, 'not found in Reductions and Assays. Please input the correct embedding name!')
    )
  # If some unpolarized cells are provided, the center of input unpolarized cell
  # embeddings will be mapped to the the center of unpolarized cell embeddings
  # in the training dataset to correct for dataset difference.
  if (length(unpolarized_cell) == 1)
  {
    if (is.na(unpolarized_cell))
    {
      mapped_emb <- input_emb
    }
    else
    {
      mapped_emb <- input_emb - rowMeans(input_emb[, unpolarized_cell]) +
        saved_params$mean_unpolar
    }
  }
  else
  {
    mapped_emb <- input_emb - rowMeans(input_emb[, unpolarized_cell]) +
      saved_params$mean_unpolar
  }

  # Transform the UCEs to PC embeddings for dimension reduction
  scaled_emb <- (mapped_emb - saved_params$mean) / saved_params$sd
  pc_emb <- data.frame(t(as.matrix(scaled_emb)) %*%
                         saved_params$pc_loadings)[, paste('PC', pc, sep='_')]
  df_polar <- list()
  models <- saved_params$models
  unpolar_responses <- saved_params$unpolar_responses

  # Predict the polarization scores using trained models
  for (state in sort(names(models)))
  {
    if (verbose)
    {
      message(paste('Measure polarization state:', state))
    }
    model_class <- class(models[[1]])
    if (length(model_class) == 2 & model_class[1] == 'glm')
    {
      predictions_state <-
        predict(models[[state]], pc_emb, type = 'response')
    }
    else if (model_class == 'svm')
    {
      if (models[[1]]$type == 3)
      {
        predictions_state <-
          predict(models[[state]], as.matrix(pc_emb))
        predictions_state <- pmax(pmin(predictions_state, 1), 0)
      }
      else if (models[[1]]$type == 0)
      {
        pred <-
          predict(models[[state]], as.matrix(pc_emb), probability = TRUE)
        predictions_state <- attr(pred, 'probabilities')[, '1']
      }

    }
    else if (model_class == 'randomForest')
    {
      predictions_state <-
        predict(models[[state]], as.matrix(pc_emb), type = 'prob')[, '1']
    }
    predictions_state <- as.vector(predictions_state)
    df_polar[[paste(state, 'score', sep = '_')]] <- predictions_state

    # conformal prediction to get predicted polarization states
    cal_quantiles <- saved_params$calibration_quantiles[[state]]
    quantile_threshold <- cal_quantiles[as.integer(error_level*100)]
    # conformal_pred values:
    # 1: empty class, 2: unpolarized class,
    # 3: polarized class; 4: both classes
    conformal_pred <- 2*(predictions_state > 1 - quantile_threshold) +
                      (predictions_state < quantile_threshold) + 1
    conformal_pred <- factor(c('Uncertain','Unpolarized',
                        'Polarized', 'Intermediate')[conformal_pred],
               levels=c('Polarized', 'Intermediate','Unpolarized','Uncertain'))
    df_polar[[paste(state, 'pred', sep = '_')]] <- conformal_pred
    # calculate the probability of each cell being unpolarized by comparing with
    # the unpolarized cells' score distribution in the reference data
    if(unpolarized_prob)
    {
      pvals <- rep(1, length(predictions_state))
      for (i in 1:length(predictions_state))
      {
        pvals[i] <-
          1 - mean(predictions_state[i] >= unpolar_responses[[state]])
      }
      df_polar[[paste(state, 'p_unpolar', sep = '_')]] <- pvals
      #df_polar[[paste(state, 'padj', sep = '_')]] <-
      #  p.adjust(pvals, method = 'BH')
    }

  }
  df_polar <- data.frame(df_polar)
  rownames(df_polar) <- colnames(input_emb)

  # Return the results as a dataframe or update the Seurat object
  if (return.df)
    return(df_polar)

  object <- AddMetaData(object, df_polar)
  return(object)
}
