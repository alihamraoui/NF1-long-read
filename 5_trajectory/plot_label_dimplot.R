#' @title Label some cells on a 2D plot
#' @description This function is similar to Seurat::DimPlot or Seurat::FeaturePlot. It is used to label particular cells, for instance root cell for trajectory inference.
#' NOTE : For the moment, there is a problem where stroke_color contains multiple colors. The multiple colors are not displayed and all strokes are grey...
#' @param sobj A Seurat object (no default)
#' @param reduction CHARACTER : the reduction to plot. It must be one of names(sobj@reductions). By default, takes the first one (default to NULL)
#' @param pt_size NUMERIC : point size (default to 0.5)
#' @param label_by CHARACTER : the column name where label criteria is. It must be something that can be found with Seurat::FetchData. By default, takes the cell barcode (default to NULL)
#' @param label_val VECTOR : a vector specifying the value label_by to highlight cells or not. in By default, takes the first cell barcode (default to NULL)
#' @param col_by CHARACTER : the column name to color cells by. It must be something that can be found with Seurat::FetchData. By default, takes Seurat::Idents (default to NULL)
#' @param col_color VECTOR : color vector to color cells with (default to NULL, i.e. ggplot default palette)
#' @param add_label LOGICAL : whether to label cells by the label_val or not (default to TRUE)
#' @param add_stroke LOGICAL : whether to add a circle around labelled cells (default to TRUE)
#' @param stroke_color CHARACTER : a column name to highlight cells by (found by Seurat::FetchData). If not found, it must be a color name to color stroke by (default to "yellow")
#' @return This function returns a ggplot object.
#' @importFrom Seurat Idents FetchData
#' @importFrom ggplot2 ggplot aes geom_point theme_classic labs scale_color_gradientn scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' @export
plot_label_dimplot = function(sobj,
                              reduction = NULL,
                              pt_size = 0.5,
                              label_by = NULL,
                              label_val = NULL,
                              col_by = NULL,
                              col_color = NULL,
                              add_label = TRUE,
                              add_stroke = TRUE,
                              stroke_color = "yellow") {
  ## Check parameters
  if (is.null(sobj)) stop("No Seurat object provided")
  if (!is.logical(add_label)) stop("add_label must be logical")
  if (!is.logical(add_stroke)) stop("add_stroke must be logical")

  # Reduction
  if (is.null(reduction)) {
    reduction = names(sobj@reductions)[1]
  } else if (!(reduction %in% names(sobj@reductions))) {
    stop("Reduction not is Seurat object. Please take a value in names(sobj@reductions).")
  }
  df_coordinates = as.data.frame(sobj@reductions[[reduction]]@cell.embeddings)[, c(1,2)]
  reduction_keys = colnames(sobj@reductions[[reduction]])[c(1,2)] # x and y labs

  # Label
  if (is.null(label_by)) {
    df_coordinates$label = colnames(sobj) # default to cell name
    if (is.null(label_val)) {
      label_val = colnames(sobj)[1] # default to first cell
    }
  } else {
    df_coordinates$label = Seurat::FetchData(sobj, label_by)[, 1]
  }

  # Color
  if (is.null(col_by)) {
    df_coordinates$groupby = Seurat::Idents(sobj) # default to active ident
  } else {
    df_coordinates$groupby = Seurat::FetchData(sobj, col_by)[, 1]
  }

  # Stroke
  if (add_stroke) {
    output = tryCatch(Seurat::FetchData(sobj, stroke_color)[, 1],
                      error = function(e) return("not found"))
    if (output[1] != "not found") {
      df_coordinates$strokeby = output
      df_coordinates$strokeby = factor(df_coordinates$strokeby,
                                       levels = aquarius::gg_color_hue(n = length(unique(df_coordinates$strokeby))))
    } else {
      df_coordinates$strokeby = stroke_color
    }
    colnames(df_coordinates) = c("Dim1", "Dim2", "label", "groupby", "strokeby")
  } else {
    colnames(df_coordinates) = c("Dim1", "Dim2", "label", "groupby")
  }

  ## Figure
  # Base : all cells
  myplot = ggplot2::ggplot(data = df_coordinates,
                           mapping = aes(x = .data$Dim1, y = .data$Dim2, fill = .data$groupby)) +
    ggplot2::geom_point(pch = 21, size = 2*pt_size, stroke = 0) +
    ggplot2::theme_classic() +
    ggplot2::xlab(reduction_keys[1]) + ggplot2::ylab(reduction_keys[2])

  # Color
  if (is.numeric(df_coordinates$groupby)) {
    if (is.null(col_color)) {
      col_color = aquarius::gg_color_hue(n = 2)
    }
    myplot = myplot +
      ggplot2::scale_fill_gradientn(name = col_by,
                                    colors = col_color)
  } else {
    if (is.null(col_color)) {
      col_color = aquarius::gg_color_hue(n = length(unique(df_coordinates$groupby)))
    }
    if (is.null(names(col_color))) {
      names(col_color) = unique(df_coordinates$groupby)
    }
    myplot = myplot +
      ggplot2::scale_fill_manual(name = col_by,
                                 breaks = names(col_color),
                                 values = col_color)
  }

  # Add label
  if (add_label) {
    myplot = myplot +
      # Label some points
      ggrepel::geom_text_repel(data = subset(df_coordinates, label %in% label_val),
                               aes(x = .data$Dim1, y = .data$Dim2, label = .data$label),
                               box.padding   = 0.35, seed = 1337L,
                               point.padding = 0.5,
                               segment.color = 'grey50',
                               show.legend = FALSE)
  }

  # Change stroke
  if (add_stroke) {
    myplot = myplot +
      # Add stroke to point
      ggplot2::geom_point(data = subset(df_coordinates, label %in% label_val),
                          mapping = aes(x = .data$Dim1, y = .data$Dim2,
                                        fill = .data$groupby, col = .data$strokeby),
                          inherit.aes = TRUE, pch = 21, stroke = 2,
                          show.legend = FALSE) +
      ggplot2::scale_color_manual(values = df_coordinates[df_coordinates$label %in% label_val, "strokeby"])
  }

  # Output
  return(myplot)
}
