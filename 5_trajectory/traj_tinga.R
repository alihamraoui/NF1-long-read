#' @title Trajectory inference with TInGa
#' @description This function is a wrapper aroung dynverse functions and TInGa
#' @param sobj A Seurat object (no default)
#' @param expression_assay CHARACTER : in which assay is the normalized expression ? (default to 'RNA')
#' @param expression_slot CHARACTER : in which slot from the expression_assay is the normalized expression ? (default to 'data')
#' @param count_assay CHARACTER : in which assay are the raw counts ? (default to 'RNA')
#' @param count_slot CHARACTER : in which slot from the count_assay are the raw counts ? (default to 'counts')
#' @param dimred_name CHARACTER : which reduction to use as a TInGa prior ? (default to 'pca')
#' @param dimred_max_dim CHARACTER : how many dimensions to use from the reduction (default to NULL, i.e. all dimensions)
#' @param seed INTEGER : a seed to run TInGa (default to 1337L)
#' @param root_cell_id CHARACTER : who is the root cell ? (no default)
#' @param tinga_parameters LIST : a named list with TInGa parameters. Here are all possible values and default (min tested - max tested) :
#'              max_iter  :    10000      (10000)
#'              max_nodes :    10         (4 - 30)
#'              epsilon_b :    0.05       (0.005 - 1)
#'              epsilon_n :    0.001      (0.0001 - 1)
#'              age_max   :    200        (100 - 500)
#'              lambda    :    200        (100 - 500)
#'              alpha     :    0.5,       (0.1 - 0.9)
#'              beta      :    0.99
#'              Only max_nodes is not set to the default value of TInGa (default to 30 in TInGa).
#' @param verbose LOGICAL : whether to print verbose or not (default to TRUE)
#' @return The function returns a trajectory object. The trajectory is rooted and contains pseudotime if root_cell_id is set.
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix t
#' @importFrom dynwrap wrap_expression add_prior_information infer_trajectory add_root add_pseudotime
#' @importFrom TInGa gng_param2
#' @export
traj_tinga = function(sobj = NULL,
                      expression_assay = "RNA",
                      expression_slot = "data",
                      count_assay = "RNA",
                      count_slot = "counts",
                      dimred_name = "pca",
                      dimred_max_dim = NULL,
                      seed = 1337L,
                      root_cell_id = NULL,
                      tinga_parameters = list(max_iter = 10000,
                                              max_nodes = 10,
                                              epsilon_b = 0.05,
                                              epsilon_n = 0.001,
                                              age_max = 200,
                                              lambda = 200,
                                              alpha = 0.5,
                                              beta = 0.99),
                      verbose = TRUE) {
  if (is.null(sobj)) stop("No Seurat object provided !")
  if (!(dimred_name %in% names(sobj@reductions))) stop("Reduction not available !")

  ## Get normalized expression matrix
  if (verbose) {message(paste0("Expression matrix : ", expression_assay, "@", expression_slot))}
  expression = Seurat::GetAssayData(sobj, assay = expression_assay, slot = expression_slot)
  expression = Matrix::t(expression)
  expression = expression[, order(colnames(expression))]

  ## Get raw counts matrix only for genes present in normalized matrix (compatible with SCT)
  if (verbose) {message(paste0("Raw counts matrix : ", count_assay, "@", count_slot))}
  counts = Seurat::GetAssayData(sobj, assay = count_assay, slot = count_slot)
  counts = counts[colnames(expression), ]
  counts = Matrix::t(counts)

  # Build a dynwrap object
  dataset = dynwrap::wrap_expression(id = sobj@project.name,
                                     expression = expression,
                                     counts = counts)

  # Add prior information : a reduced space
  dimred = sobj[[dimred_name]]@cell.embeddings
  if (!is.null(dimred_max_dim)) {
    dimred = dimred[, c(1:dimred_max_dim)]
  } else {
    dimred_max_dim = ncol(dimred)
  }
  if (verbose) {message(paste0("Reduction : ", dimred_name, " with ", dimred_max_dim, " dimensions"))}
  dataset = dynwrap::add_prior_information(dataset = dataset,
                                           dimred = dimred)

  ## TInGa parameters
  default_params = list(max_iter = 10000,
                        max_nodes = 10,
                        epsilon_b = 0.05,
                        epsilon_n = 0.001,
                        age_max = 200,
                        lambda = 200,
                        alpha = 0.5,
                        beta = 0.99)
  for (param in names(tinga_parameters)) {
    if (param %in% names(default_params)) {
      default_params[[param]] = tinga_parameters[[param]]
    }
  }

  ## Running TInGa
  if (verbose) {message("Trajectory inference")}
  set.seed(seed)
  traj = dynwrap::infer_trajectory(dataset = dataset,
                                   method = TInGa::gng_param2(),
                                   seed = seed, give_priors = 'dimred',
                                   parameters = default_params)

  ## Eventually add a root
  if (!is.null(root_cell_id)) {
    if (verbose) {message("Set a root to trajectory")}
    if (root_cell_id %in% colnames(sobj)) {
      ## Add root cell
      traj = dynwrap::add_root(trajectory = traj, root_cell_id = root_cell_id)

      ## Add pseudotime
      traj = dynwrap::add_pseudotime(trajectory = traj)
    } else {
      warning(paste0(root_cell_id, " is not in Seurat object. Trajectory is not rooted."))
    }
  }

  ## Output
  return(traj)
}
