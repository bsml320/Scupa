#' Polarization states
#'
#' A list of immune cell polarization states with top marker genes and drivind cytokines.
#'
#' @format A list with 66 elements of polarization states. The polarization states
#' are for 12 cell types: B, NK, CD8T, CD4T, Treg, Tgd, pDC, cDC1, cDC2, MigDC, LC,
#' Macro, Mono, Neu. Each element is a list containing 2 variables:
#' \describe{
#'   \item{marker}{The 5 marker genes of the polarization state.}
#'   \item{cytokine}{The driving cytokines of the polarization state.}
#' }
#' @references Cui A, Huang T, Li S, et al. Dictionary of immune responses to
#' cytokines at single-cell resolution. Nature. 2024;625(7994):377-384.
#' doi:10.1038/s41586-023-06816-9
#'
"polar_states"

#' Trained models and parameters for polarization measurement
#'
#' A list of trained support vector machine models and parameters of training dataset for measuring polarization scores and P values in other dataset.
#'
#' @format A list with 14 elements. The element names are B, NK, CD8T, CD4T, Treg,
#' Tgd, pDC, cDC1, cDC2, MigDC, LC, Macro, Mono, Neu, corresponding to each cell
#' type. Each element is a list containing 6 variables:
#' \describe{
#'   \item{mean}{The mean of Universial Cell Embeddings of all cells in the training dataset.}
#'   \item{sd}{The standard deviation of Universial Cell Embeddings of all cells in the training dataset.}
#'   \item{mean_unpolar}{The mean of Universial Cell Embeddings of unpolarized cells in the training dataset.}
#'   \item{pc_loadings}{The principle component loadings to transform Universial Cell Embeddings into principle component embeddings.}
#'   \item{regression_models}{The trained support vector machine models to predict the polarization scores based on principle component embeddings for each cell polarization state.}
#'   \item{unpolar_responses}{The polarization scores of unpolarized cells in the training dataset. The polarization P values of input data are calculated by comparing to the distribution of these responses.}
#' }
#'
"polar_params"
