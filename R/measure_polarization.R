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
#' @return A list containing 6 variables:
#' \describe{
#'   \item{mean}{The mean of Universial Cell Embeddings of all cells in the training dataset.}
#'   \item{sd}{The standard deviation of Universial Cell Embeddings of all cells in the training dataset.}
#'   \item{mean_unpolar}{The mean of Universial Cell Embeddings of unpolarized cells in the training dataset.}
#'   \item{pc_loadings}{The principal component loadings to transform Universial Cell Embeddings into principal component embeddings.}
#'   \item{regression_models}{The trained support vector machine models to predict the polarization scores based on principal component embeddings for each cell polarization state.}
#'   \item{unpolar_responses}{The polarization scores of unpolarized cells in the training dataset. The polarization P values of input data are calculated by comparing to the distribution of these responses.}
#' }
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
      'celltype must be one of: B, NK, CD8T, CD4T, Treg, Tgd, pDC,
      cDC1, cDC2, MigDC, LC, Macro, Mono, Neu.'
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
#' models to predict the polarization of cells in the input data.
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
#' @param pc The index of principal component to be used by machine learning models.
#' Default: 1:20.
#' @param unpolarized_cell The name of cells that are considered unpolarized cells.
#' See details for selecting unpolarized cells. Default: NA.
#' @param return.df Whether to return the results as a dataframe. If TRUE,
#' return a dataframe. Else, return the Seurat object with updated meta.data.
#' Default: FALSE.
#' @param verbose Whether to message the progress in different polarization states.
#' Default: TRUE.
#'
#' @return Return either a dataframe or a Seurat object with updated meta.data,
#' depending on param return.df.
#' The dataframe includes 3 * N(polarization states) columns.
#' For each cell type, there are 4~6 polarization states.
#' For each polarization state, the polarization scores (range 0~1), P values,
#' and adjusted P values are calculated for all cells.
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
#' @references Cui A, Huang T, Li S, et al. Dictionary of immune responses to
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
#' # Plot the CD8+ T cell polarization P values
#' FeaturePlot(cd8t, c('T8.a_p','T8.b_p','T8.c_p','T8.e_p','T8.f_p'))
#' }
#'
MeasurePolar <- function(object,
                         celltype,
                         embedding = 'uce',
                         pc = 1:20,
                         unpolarized_cell = NA,
                         return.df = FALSE,
                         verbose = TRUE)
{
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
      'celltype must be one of: B, NK, CD8T, CD4T, Treg, Tgd, pDC,
         cDC1, cDC2, MigDC, LC, Macro, Mono, Neu.'
    )
  }
  saved_params <- LoadPolarParams(celltype,
                                  verbose = verbose)
  if(embedding %in% Reductions(object))
    input_emb <- t(Embeddings(object, reduction=embedding))
  else if(embedding %in% Assays(object))
    input_emb <- object[[embedding]]@data
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
  models <- saved_params$regression_models
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
        predictions_state[predictions_state < 0] <- 0
        predictions_state[predictions_state > 1] <- 1
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
      require('RandomForest')
      predictions_state <-
        predict(models[[state]], as.matrix(pc_emb), type = 'prob')[, '1']
    }
    # calculate the P value of each cell by comparing with the distribution of
    # scores of unpolarized cells in the training data
    pvals <- rep(1, length(predictions_state))
    for (i in 1:length(predictions_state))
    {
      pvals[i] <-
        1 - mean(predictions_state[i] >= unpolar_responses[[state]])
    }

    df_polar[[paste(state, 'score', sep = '_')]] <-
      as.vector(predictions_state)
    df_polar[[paste(state, 'p', sep = '_')]] <- pvals
    df_polar[[paste(state, 'padj', sep = '_')]] <-
      p.adjust(pvals, method = 'BH')
  }
  df_polar <- data.frame(df_polar)
  rownames(df_polar) <- colnames(input_emb)

  if (return.df)
    return(df_polar)

  object <- AddMetaData(object, df_polar)
  return(object)
}
