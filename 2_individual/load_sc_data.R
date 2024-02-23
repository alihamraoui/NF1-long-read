#' @title Load scRNA-Ses count matrix into a Seurat object
#' @description Wrapper to load a cellranger counts output matrix into a Seurat object
#' @param data_path CHARACTER : pathway where the three cellranger files are stored (matrix, features and barcodes)
#' @param sample_name CHARACTER : sample name
#' @param assay CHARACTER : name of the assay corresponding to the initial input data (default to 'RNA')
#' @param gene_column NUMERIC : 1 (gene ID) or 2 (gene symbol), to be passed to the gene.column parameter in Seurat::Read10X function (default to 2)
#' @param droplets_limit INTEGER : max number of droplets to activate empty droplets filtering (default to 1E+05)
#' @param emptydrops_fdr INTEGER : False Discovery Rate to filter out empty droplets candidates (default to 1E+03)
#' @param filter_replicates LOGICAL : whether you want or not to group columns by barcode, in case of presence of duplicated barcode (default to TRUE)
#' @param BPPARAM an object of class DoparParam to be passed to BPPARAM parameter of DropletUtils::emptyDrops function. Could be generate by our create.parallel.instance function
#' @param my_seed INTEGER : seed to be used by DropletUtils::emptyDrops (default to 1337L)
#' @param verbose LOGICAL : whether to print messages or not (default to TRUE)
#' @return A Seurat object containing count data and parameters used
#' @importFrom DelayedArray colSums
#' @importFrom dplyr %>% group_by slice
#' @importFrom DropletUtils barcodeRanks emptyDrops
#' @importFrom rlang .data
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom tibble tibble
#' @export
load_sc_data = function(data_path = NULL,
                        sample_name = NULL,
                        assay = 'RNA',
                        gene_column = 2,
                        droplets_limit = 1E+05,
                        emptydrops_fdr = 1E-03,
                        filter_replicates = TRUE,
                        BPPARAM = BiocParallel::SerialParam(),
                        my_seed = 1337L,
                        verbose = TRUE){
  if (!(file.exists(data_path))) stop("Data path does not exists.")
  if (is.null(sample_name)) stop("Please provide a sample_name, it will stored in @project.name.")

  # ======================================== Load count matrix =======================================
  if (verbose) message("Loading data ...")

  scmat = Seurat::Read10X(data_path,
                          gene.column = gene_column)
  scmat = scmat[, order(colnames(scmat))]

  # ======================================== Filter duplicated cell barcodes ======================
  # it keeps the most popular entry
  if (verbose) {
    message('Droplets matrix dimensions :')
    print(dim(scmat))
  }

  droplets.nb = ncol(scmat)

  if (filter_replicates) {
    dup_bc = unique(colnames(scmat)[duplicated(colnames(scmat))])

    if (length(dup_bc) > 0) {
      if (verbose) {
        message(paste0("Found ", length(dup_bc), ' (', sprintf("%.2f", length(dup_bc) / ncol(scmat) * 100),
                       "%) replicated cell barcodes ! Filtering ..."))
      }

      dedup_tbl = tibble(barcode = colnames(scmat),
                         ori.index = 1:ncol(scmat),
                         count = DelayedArray::colSums(scmat)) %>%
        dplyr::group_by(.data$barcode) %>%
        dplyr::slice(which.max(.data$count))
      scmat = scmat[, dedup_tbl$ori.index]
      rm(dedup_tbl)

      droplets.nb = ncol(scmat)

      if (verbose) {
        message('Droplets matrix dimensions (deduplicated) :')
        print(dim(scmat))
      }

    } else message('No replicated barcode found.')
  }

  # ======================================== Remove empty droplets =================================
  umi_total_nb = sum(scmat)

  if (verbose) {
    message('Total UMIs :')
    print(umi_total_nb)
  }

  if (ncol(scmat) > droplets_limit & !is.null(emptydrops_fdr)) {
    if (verbose) message("Removing empty droplets with emptyDrops")
    set.seed(my_seed)
    bc_rank = DropletUtils::barcodeRanks(scmat)

    set.seed(my_seed)
    bc_rank2 = DropletUtils::emptyDrops(scmat, BPPARAM = BPPARAM)
    scmat = scmat[, which(bc_rank2$FDR < emptydrops_fdr)]

    umi_kept_nb = sum(scmat)

    if (verbose) {
      message('Droplets matrix dimensions (filtered) :')
      print(dim(scmat))

      message('Total UMIs (filtered) :')
      print(umi_kept_nb)

      message('Fraction of UMIs in cells :')
      print(umi_kept_nb / umi_total_nb)
    }
    rm(bc_rank, bc_rank2)
  } else {
    umi_kept_nb = umi_total_nb
  }

  # ======================================== Create Seurat object =================================
  #sobj = Seurat::CreateSeuratObject(counts = scmat, project = sample_name, assay = assay)
  #rm(scmat)

  #sobj[[paste0('log_nCount_', assay)]] = log(sobj[[paste0('nCount_', assay)]])
  # sobj@misc$droplets = droplets.nb
  # sobj@misc$umi.raw = umi_total_nb
  # sobj@misc$umi.filtered = umi_kept_nb
  # sobj@misc$cells.ori = ncol(sobj)
  # sobj@misc$params = list(emptydrops_fdr = emptydrops_fdr, seed = my_seed)
  # sobj@misc$Rsession = utils::sessionInfo()

  ## Output
  #return(sobj)
  return(scmat)
}
