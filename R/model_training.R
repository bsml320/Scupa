# A helper function to calculate cosine similarity in a matrix
cosine_sim <- function(matrix)
{
  normalized_matrix <- sweep(matrix, 2, sqrt(colSums(matrix^2)), `/`)
  cosine_similarity_matrix <- t(normalized_matrix) %*% normalized_matrix
  return(cosine_similarity_matrix)
}

#' Assign the polarization states in Immune Dictionary
#'
#' Assign the polarization states to a cell type in Immune Dictionary based on
#' Universial Cell Embeddings and marker gene expression.
#'
#' Run this function only if you want to retrain the models using custom parameters!
#'
#' @param object A Seurat object of a cell type from the Immune Dictionary.
#' @param state_list A list of 4~6 elements of polarization states. Each element
#' should be a list of 2 variables: marker and cytokine. See the data(polar_states)
#' for a full list of polarization states.
#' @param include_other_cyt Logical. Whether to considering the non-driving cytokines
#' when assigning polarization states. If FALSE, only cells from samples treated
#' with driving cytokines may be assigned to the polarized states. If TRUE, all
#' cells treated by any cytokines may be assigned to the polarized states.
#' Default: FALSE.
#' @param dist_thres A numeric value in 0~1. Cells with cosine distance bigger
#' than this threshold will not be assigned to the polarized states.
#' Default: 0.2.
#' @param marker_quantile A numeric value in 0~1. Cells with the mean marker gene
#' expression smaller than this quantile will not be assigned into the polarized
#' states. Default: 0.9.
#' @param embedding The name of embedding saving the cell embeddings from the
#' single-cell foundation model in the Seurat object. Default: uce.
#' @param center_unpolar Logical. Whether to use unpolarized cell center as the
#' initial point for cell embedding vectors. The cosine distance is calculated
#' using the ell embedding vectors.If FALSE, the center of all cells is used as
#' the initial point for cell embedding vectors. Default: TRUE.
#'
#' @return A Seurat object with an additional variable 'polar' corresponding to
#' polarization states.
#'
#' @importFrom stats quantile
#' @importFrom Matrix colMeans
#' @importFrom Seurat Reductions Embeddings Assays
#' @export
#'
#' @examples
#' \dontrun{
#' # Please download the Immune Dictionary data at
#' cd8t <- AssignStates(cd8t, polar_states_cd8t, dist_thres=0.85)
#' }
#'
AssignStates <- function(object,
                         state_list,
                         include_other_cyt=FALSE,
                         dist_thres=0.8,
                         marker_quantile=0.9,
                         embedding='uce',
                         center_unpolar=TRUE)
{
  object$polar <- 'else'
  object$polar[object$cyt_display == 'PBS'] <- 'unpolarized'
  states <- names(state_list)
  min_dist_state <- dist_thres
  if(embedding %in% Reductions(object))
    emb_matrix <- t(Embeddings(object, reduction=embedding))
  else if(embedding %in% Assays(object))
    emb_matrix <- object[[embedding]]@data
  else
    stop(
      paste('Embedding', embedding, 'not found in Reductions and Assays.
            Please input the correct embedding name!')
    )
  unpolar_center <- apply(emb_matrix[,object$cyt_display == 'PBS'],1,mean)
  all_center <- apply(emb_matrix,1,mean)
  if(center_unpolar)
    using_center <- unpolar_center
  else
    using_center <- all_center

  for(state in states)
  {
    message(paste('Assigning state:', state))
    cell_state <- names(object$cyt_display)[object$cyt_display %in% state_list[[state]]$cytokine]
    mean_exp_marker <- colMeans(x = object[['RNA']]@data[state_list[[state]]$marker, ], na.rm = TRUE)

    if(include_other_cyt)
    {
      simlarity_cell <- cosine_sim(emb_matrix - using_center)[, cell_state]
    }
    else
    {
      simlarity_cell <- cosine_sim(emb_matrix[, cell_state] - using_center)
    }
    mean_dist_state <- 1 - rowMeans(simlarity_cell)
    cell_state_filtered <- names(mean_dist_state)[mean_dist_state < min_dist_state]
    if(include_other_cyt)
    {
      min_dist_state <- pmin(mean_dist_state, min_dist_state)
    }
    cell_state_final <- rownames(object@meta.data) %in% cell_state_filtered &
      mean_exp_marker > quantile(mean_exp_marker, marker_quantile)
    object$polar[cell_state_final] <- state
  }
  return(object)
}


#' Train machine learning models and calculating parameters
#'
#' Train machine learning models and calculating some necessary parameters to
#' measure immune cell polarization. The models are trained using the output from
#' \code{AssignStates}.
#'
#' Run this function only if you want to retrain the models using custom parameters!
#'
#' @param object A Seurat object of a cell type from the Immune Dictionary.
#' @param embedding The name of embedding saving the cell embeddings from the single-cell
#' foundation model in the Seurat object. Default: UCE.
#' @param pc The index of principal component as predictor variables for training
#' machine learning models. Default: 1:20.
#' @param classify_method A string from one of 'logistic','svc','svr','rf'. The
#' machine learning models to classify polarized cells and unpolarized cells.
#' Default: svr.
#' @param semi_supervised Logical. Whether to use the semi-supervised approach in
#' training machine learning models. If TRUE, the unlabeled cells from samples
#' treated with non-driving cytokines are used to train the final output machine
#' learning models. Default: FALSE.
#' @param verbose Whether to message the progress in different polarization states.
#' Default: TRUE.
#'
#' @return A list containing 6 variables:
#' \describe{
#'   \item{mean}{The mean of Universial Cell Embeddings of all cells in the training
#'   dataset.}
#'   \item{sd}{The standard deviation of Universial Cell Embeddings of all cells
#'   in the training dataset.}
#'   \item{mean_unpolar}{The mean of Universial Cell Embeddings of unpolarized cells
#'   in the training dataset.}
#'   \item{pc_loadings}{The principal component loadings to transform Universial
#'   Cell Embeddings into principal component embeddings.}
#'   \item{regression_models}{The trained machine learning models to predict the
#'   polarization scores based on principal component embeddings for each cell
#'   polarization state.}
#'   \item{unpolar_responses}{The polarization scores of unpolarized cells in the
#'   training dataset. The polarization P values of input data are calculated by
#'   comparing to the distribution of these responses.}
#' }
#'
#' @seealso Logistic regression model: \code{\link[stats]{glm}} with parameter: family = 'binomial'.
#'
#' SVC model: \code{\link[e1071]{svm}} with parameter: kernel='linear'.
#'
#' SVR model: \code{\link[e1071]{svm}} with parameters: kernel='linear', type='eps-regression'.
#'
#' Random forest model: \code{\link[randomForest]{randomForest}} with default parameters.
#'
#' @importFrom stats sd glm predict
#' @importFrom e1071 svm
#' @importFrom Seurat Reductions Embeddings Assays Loadings FetchData
#' @export
#'
#' @examples
#' \dontrun{
#' # Please download the Immune Dictionary data at
#' cd8t <- AssignStates(cd8t, polar_states_cd8t, dist_thres=0.85)
#' polar_params_cd8t <- CalculateParams(cd8t, classify_method='svr')
#' }
#'
CalculateParams <- function(object,
                            embedding='uce',
                            pc=1:20,
                            classify_method='svr',
                            semi_supervised=FALSE,
                            verbose=TRUE)
{
  if(embedding %in% Reductions(object))
    embs <- t(Embeddings(object, reduction=embedding))
  else if(embedding %in% Assays(object))
    embs <- object[[embedding]]@data
  else
    stop(
      paste('Embedding', embedding, 'not found in Reductions and Assays.
            Please input the correct embedding name!')
    )
  emb_means <- rowMeans(embs)
  emb_sds <- apply(embs, 1, sd)
  emb_means_unpolar <- rowMeans(embs[,object$polar == 'unpolarized'])
  emb_loadings <- Loadings(object, reduction = 'pca')[rownames(embs), ]
  state_train <- FetchData(object[,object$polar != 'else'],
                           c(paste('PC',pc,sep='_'),'polar'))
  regression_models <- list()
  unpolar_responses <- list()
  for(state in sort(unique(state_train$polar[state_train$polar != 'unpolarized'])))
  {
    if(verbose)
    {
      message(paste('Calculating parameters for state:', state))
    }
    state_train_i <- state_train[state_train$polar %in% c('unpolarized', state),]
    state_train_i$polar <- ifelse(state_train_i$polar=='unpolarized', 0, 1)
    emb_unpolar <- state_train_i[state_train_i$polar==0,paste('PC',pc,sep='_')]
    train_X <- as.matrix(state_train_i[,paste('PC',pc,sep='_')])
    train_Y <- state_train_i$polar
    sample_weights <- ifelse(train_Y==0, 1, sum(train_Y==0)/sum(train_Y==1))

    if(semi_supervised)
    {
      state_train_extended <- FetchData(object[,object$polar %in% c('unpolarized','else',state)],
                                        c(paste('PC',pc,sep='_'),'polar'))
      train_X_extended <- as.matrix(state_train_extended[,paste('PC',pc,sep='_')])
      train_Y_extended <- rep(NA, length(state_train_extended$polar))
      train_Y_extended[state_train_extended$polar=='unpolarized'] <- 0
      train_Y_extended[state_train_extended$polar==state] <- 1
    }

    if(classify_method == 'logistic')
    {
      df <- data.frame(train_X)
      df$Y <- train_Y
      model <- strip::strip(glm(Y ~ ., df, family = 'binomial'),
                     keep=c('predict','print'))
      if(semi_supervised)
      {
        prediction_extended <- predict(model, data.frame(train_X_extended), type='response')
        mean_prediction_state <- mean(prediction_extended[state_train_extended$polar==state])
        mean_prediction_unpolar <- mean(prediction_extended[state_train_extended$polar=='unpolarized'])
        # for unlabelded cells, assign to unpolarized if their prediction is up to mean prediction of unpolarized cells
        train_Y_extended[prediction_extended<=mean_prediction_unpolar] <- 0
        # for unlabelded cells, assign to polarized if their prediction is no less than mean prediction of polarized cells
        train_Y_extended[prediction_extended>=mean_prediction_state] <- 1
        cell_kept <- !is.na(train_Y_extended)
        train_X_extended <- train_X_extended[cell_kept,]
        train_Y_extended <- train_Y_extended[cell_kept]
        df <- data.frame(train_X_extended)
        df$Y <- train_Y_extended
        model <- strip(glm(Y ~ ., df, family = 'binomial'),
                       keep=c('predict','print'))
      }
      unpolar_prediction <- predict(model, emb_unpolar, type='response')
    }

    else if(classify_method == 'svc')
    {
      model <- svm(train_X,
                   factor(train_Y, levels=c(0,1)),
                   kernel='linear',
                   weights = sample_weights,
                   probability=TRUE)
      if(semi_supervised)
      {
        pred <- predict(model, train_X_extended, probability=TRUE)
        prediction_extended <- attr(pred, 'probabilities')[,'1']
        mean_prediction_state <- mean(prediction_extended[state_train_extended$polar==state])
        mean_prediction_unpolar <- mean(prediction_extended[state_train_extended$polar=='unpolarized'])
        # for unlabelded cells, assign to unpolarized if their prediction is up to mean prediction of unpolarized cells
        train_Y_extended[prediction_extended<=mean_prediction_unpolar] <- 0
        # for unlabelded cells, assign to polarized if their prediction is no less than mean prediction of polarized cells
        train_Y_extended[prediction_extended>=mean_prediction_state] <- 1
        cell_kept <- !is.na(train_Y_extended)
        train_X_extended <- train_X_extended[cell_kept,]
        train_Y_extended <- train_Y_extended[cell_kept]
        sample_weights <- ifelse(train_Y_extended==0, 1, sum(train_Y_extended==0)/sum(train_Y_extended==1))
        model <- svm(train_X_extended,
                     factor(train_Y_extended, levels=c(0,1)),
                     kernel='linear',
                     probability=TRUE,
                     weights = sample_weights)
      }
      pred <- predict(model, as.matrix(emb_unpolar), probability=TRUE)
      unpolar_prediction <- attr(pred, 'probabilities')[,'1']
    }

    else if(classify_method == 'svr')
    {
      model <- svm(train_X,
                   train_Y,
                   kernel='linear',
                   type='eps-regression',
                   weights = sample_weights)
      if(semi_supervised)
      {
        prediction_extended <- predict(model, train_X_extended)
        mean_prediction_state <- mean(prediction_extended[state_train_extended$polar==state])
        mean_prediction_unpolar <- mean(prediction_extended[state_train_extended$polar=='unpolarized'])
        # for unlabelded cells, assign to unpolarized if their prediction is up to mean prediction of unpolarized cells
        train_Y_extended[prediction_extended<=mean_prediction_unpolar] <- 0
        # for unlabelded cells, assign to polarized if their prediction is no less than mean prediction of polarized cells
        train_Y_extended[prediction_extended>=mean_prediction_state] <- 1
        cell_kept <- !is.na(train_Y_extended)
        train_X_extended <- train_X_extended[cell_kept,]
        train_Y_extended <- train_Y_extended[cell_kept]
        sample_weights <- ifelse(train_Y_extended==0,
                                 1, sum(train_Y_extended==0)/sum(train_Y_extended==1))
        model <- svm(train_X_extended,
                     train_Y_extended,
                     kernel='linear',
                     type='eps-regression',
                     weights = sample_weights)
      }
      unpolar_prediction <- predict(model, as.matrix(emb_unpolar))
      unpolar_prediction[unpolar_prediction<0] <- 0
      unpolar_prediction[unpolar_prediction>1] <- 1
    }

    else if(classify_method == 'rf')
    {
      model <- randomForest::randomForest(train_X,
                            factor(train_Y, levels=c(0,1)))
      if(semi_supervised)
      {
        prediction_extended <- predict(model, train_X_extended, type='prob')[,'1']
        mean_prediction_state <- mean(prediction_extended[state_train_extended$polar==state])
        mean_prediction_unpolar <- mean(prediction_extended[state_train_extended$polar=='unpolarized'])
        # for unlabelded cells, assign to unpolarized if their prediction is up to mean prediction of unpolarized cells
        train_Y_extended[prediction_extended<=mean_prediction_unpolar] <- 0
        # for unlabelded cells, assign to polarized if their prediction is no less than mean prediction of polarized cells
        train_Y_extended[prediction_extended>=mean_prediction_state] <- 1
        cell_kept <- !is.na(train_Y_extended)
        train_X_extended <- train_X_extended[cell_kept,]
        train_Y_extended <- train_Y_extended[cell_kept]
        sample_weights <- ifelse(train_Y_extended==0,
                                 1, sum(train_Y_extended==0)/sum(train_Y_extended==1))
        model <- randomForest::randomForest(train_X_extended,
                              factor(train_Y_extended), levels=c(0,1))
      }
      unpolar_prediction <- predict(model, as.matrix(emb_unpolar), type='prob')[,'1']
    }

    regression_models[[state]] <- model
    unpolar_responses[[state]] <- sort(unpolar_prediction)
  }

  object_saved_params <- list(mean=emb_means,
                              sd=emb_sds,
                              mean_unpolar=emb_means_unpolar,
                              pc_loadings=emb_loadings,
                              regression_models=regression_models,
                              unpolar_responses=unpolar_responses)
  return(object_saved_params)
}
