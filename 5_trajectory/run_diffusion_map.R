#' @title Generate a diffusion map
#' @description This function is a wrapper around destiny::DiffusionMap function to compute a diffusion map and add it in the Seurat object. Note that it take a lot of time to run on a full expression matrix. It may cause R session to abort... Moreover, diffusion map can be well ran on a reduced space, such as PCA.
#' @param sobj A Seurat object (no default)
#' @param input CHARACTER or MATRIX : if it is a character, the function will search of an assay or a reductions in the Seurat object. Otherwise, it is supposed to be a dataframe with cells in rows and dimensions in columns (default to NULL)
#' @param seed INTEGER : the seed to be used by destiny::DiffusionMap (default to 1337L)
#' @param verbose LOGICAL : whether to print messages or not (default to FALSE)
#' @param n_eigs INTEGER : number of eigenvectors/values to compute in the diffusion map (default to 50)
#' @param suppress_dpt LOGICAL : whether to skip calculation of necessary (but spacious) information for DPT in the returned object or not (default to TRUE)
#' @param return_dm LOGICAL : whether to return the dm object or not (default to FALSE)
#' @param ... other parameters passed to destiny::DiffusionMap
#' @return This function returns the input Seurat object with the new DimReducObject. If return_dm is set to TRUE, returns a list with $sobj containing the Seurat object and $dm containing the dm object.
#' To plot eigenvalues and dimensions, use :
#'        `ggplot2::qplot(y = destiny::eigenvalues(dm))`
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom destiny DiffusionMap
#' @export
run_diffusion_map = function(sobj,
                             input = NULL,
                             seed = 1337L,
                             verbose = FALSE,
                             n_eigs = 50,
                             suppress_dpt = TRUE,
                             return_dm = FALSE,
                             ...) {

  ## Identify input
  if (is.null(input)) {
    stop("input must be a character specifying the assay (@data) or the reduction to use")
  } else {
    if (is.character(input)) {
      dm_name = paste0(input, "_dm")
      if (input %in% names(sobj@assays)) {
        if (verbose) {
          message(paste0("Set input to the assay ", input, "@data"))
        }
        input = as.matrix(sobj[[input]]@data)
      } else {
        if (input %in% names(sobj@reductions)) {
          if (verbose) {
            message(paste0("Set input to the reduction ", input))
          }
          input = sobj[[input]]@cell.embeddings
        } else {
          stop("input not found in Seurat object")
        }
      }
    } else {
      dm_name = "dm"
    }
  }

  ## Compute diffusion map
  if (verbose) message("Run DiffusionMap")
  set.seed(seed) # diffusion map is random !
  dm = destiny::DiffusionMap(input,
                             n_eigs = n_eigs,
                             suppress_dpt = suppress_dpt,
                             verbose = verbose,
                             ...)

  ## Add dm to Seurat object
  if (verbose) message(paste0("Add diffusion map to sobj, with the name : ", dm_name))
  dm_embeddings = data.frame(dm@eigenvectors) # cells x eigenvectors matrix

  ## if duplicated rows (cells) in input, dm_embeddings does not contain theme
  if (nrow(dm_embeddings) != nrow(input)) {
    # missing_cells = missing cells names
    missing_cells = rownames(input)[!(rownames(input) %in% rownames(dm_embeddings))]

    # all_duplicated_rows = missing_cells and there matches
    all_duplicated_rows = names(which(duplicated(input) | duplicated(input, fromLast = TRUE)))

    # matching_cells = such that missing_cells U matching_cells = all_duplicated_rows
    matching_cells = all_duplicated_rows[!(all_duplicated_rows %in% missing_cells)]

    # find match : missing_cells is a named vector with 'missing_cell "matching_cell"' element
    names(missing_cells) = missing_cells
    missing_cells = lapply(missing_cells, FUN = function(one_cell) {
      its_match = lapply(matching_cells, FUN = function(one_candidate_match) {
        return(identical(input[one_cell, ], input[one_candidate_match, ]))
      }) %>% unlist() %>% which()
      its_match = matching_cells[its_match]
      return(its_match)
    }) %>% unlist()

    # restore dm_embeddings
    row_names = rownames(dm_embeddings)
    dm_embeddings = rbind(dm_embeddings,
                          dm_embeddings[missing_cells, ])
    rownames(dm_embeddings) = c(row_names, names(missing_cells))
    # good order
    dm_embeddings = dm_embeddings[rownames(input), ]
  }

  ## Add in Seurat
  colnames(dm_embeddings) = paste0("DM", 1:ncol(x = dm_embeddings)) # dimension names
  rownames(dm_embeddings) =  colnames(sobj) # cell names

  sobj[[dm_name]] = Seurat::CreateDimReducObject(embeddings = as.matrix(dm_embeddings),
                                                 key = "DM_",
                                                 assay = Seurat::DefaultAssay(sobj))
  rm(dm_embeddings)

  ## Output
  if (return_dm) {
    output = list(sobj = sobj,
                  dm = dm)
  } else {
    output = sobj
  }

  return(output)
}
