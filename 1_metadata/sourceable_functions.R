source("/Users/hamraoui/Analyses/script_audrey/cell_annot_custom.R")
source("/Users/hamraoui/Analyses/script_audrey/clusters_annot.R")

annotate.aud = function (ScObject, cell_markers) {
  if ("cell_type" %in% colnames(ScObject@meta.data)) {
    ScObject$cell_type = NULL
  }
  ScObject = cell_annot_custom(ScObject,
                               assay = "RNA",          
                               newname = "cell_type",  
                               markers = cell_markers)
}