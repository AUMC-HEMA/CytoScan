#' Plot heatmap
#'
#' @param CS CytoScan object
#' @param featMethod Which features to use
#' @param plotData Which data to plot ("test", "reference" or "all")
#'
#' @return Heatmap
#' @export
plotHeatmap <- function(CS, featMethod, plotData = "test") {
  
  # --- Step 1: Select input features ---
  test_data <- CS$features$test[[featMethod]]
  n_test <- nrow(test_data)
  
  ref_data <- CS$features$reference[[featMethod]]
  if (is.null(ref_data) || nrow(ref_data) == 0) {
    ref_data <- NULL
    n_ref <- 0
  } else {
    n_ref <- nrow(ref_data)
  }
  
  if (plotData == "test") {
    inputData <- test_data
    outlier_flags <- if (!is.null(CS$outliers[[featMethod]])) {
      factor(as.logical(CS$outliers[[featMethod]][rownames(test_data)]), levels = c(FALSE, TRUE))
    } else {
      factor(rep(FALSE, n_test), levels = c(FALSE, TRUE))
    }
  } else if (plotData == "reference") {
    if (is.null(ref_data)) stop("Reference data not available")
    inputData <- ref_data
    outlier_flags <- factor(rep(FALSE, n_ref), levels = c(FALSE, TRUE))
  } else if (plotData == "all") {
    inputData <- rbind(test_data, ref_data)
    outlier_flags <- factor(c(
      if (!is.null(CS$outliers[[featMethod]])) as.logical(CS$outliers[[featMethod]][rownames(test_data)]) else rep(FALSE, n_test),
      rep(FALSE, n_ref)
    ), levels = c(FALSE, TRUE))
  } else {
    stop("plotData must be one of 'test', 'reference', 'all'")
  }
  
  # --- Step 2: Standardize per channel ---
  channels <- CS$metadata[[featMethod]]$channels
  for (channel in channels) {
    cols <- grep(channel, colnames(inputData), value = TRUE)
    mat <- as.matrix(inputData[, cols])
    inputData[, cols] <- (mat - mean(mat)) / sd(mat)
  }
  
  # --- Step 3: Split columns by channel ---
  n_cols_per_channel <- ncol(inputData) / length(channels)
  split <- rep(1:length(channels), each = n_cols_per_channel)
  
  # --- Step 4: Row annotation for outliers ---
  row_anno <- ComplexHeatmap::rowAnnotation(
    Outlier = outlier_flags,
    col = list(Outlier = c("FALSE" = "gray", "TRUE" = "red")),
    show_annotation_name = TRUE
  )
  
  # --- Step 5: Build heatmap ---
  g <- ComplexHeatmap::Heatmap(
    as.matrix(inputData),
    cluster_rows = as.dendrogram(hclust(dist(inputData), method = "ward.D")),
    cluster_columns = FALSE,
    column_split = split,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_title = channels,
    column_title_side = "top",
    column_title_rot = 90,
    column_title_gp = grid::gpar(fontsize = 10),
    left_annotation = row_anno,
    heatmap_legend_param = list(title = "Z-score"),
    row_dend_width = grid::unit(3, "cm")
  )
  
  return(g)
}



#' Plot PCA
#'
#' @param CS CytoScan object
#' @param featMethod Which features to use
#' @param fitData Which data to fit PCA on ("test", "reference" or "all")
#' @param plotData Which data to plot based on PCA fit ("test", "reference" or "all")
#' @param PCx Which PC to plot along x-axis
#' @param PCy Which PC to plot along y-axis
#' @param color Category used for point color ("outlier", "novelty", "labels")
#' @param shape Category used for point shape ("outlier", "novelty", "labels")
#' @param plotBars Whether to plot loading coefficients as bars along
#' @param plotArrows Whether to plot loadings as arrows
#' @param nLoadings Number of loadings to plot. Selects most influential
#'
#' @return PCA plot
#' @export
plotPCA <- function(CS, featMethod, fitData = "test", plotData = "test",
                    PCx = 1, PCy = 2, color = NULL, shape = NULL,
                    plotBars = FALSE, plotArrows = FALSE, nLoadings = NULL) {
  
  # 1. Select features for PCA fit
  inputData <- switch(fitData,
                      test = CS$features$test[[featMethod]],
                      reference = CS$features$reference[[featMethod]],
                      all = rbind(CS$features$test[[featMethod]], CS$features$reference[[featMethod]]))
  
  # 2. Handle zero/constant columns
  zero_var <- which(apply(inputData, 2, function(x) sd(x) == 0))
  if (length(zero_var) > 0) inputData[, zero_var] <- 1e-6
  
  # 3. Fit PCA
  pca <- prcomp(inputData, center = TRUE, scale. = TRUE)
  x_var <- round((pca$sdev[PCx]^2 / sum(pca$sdev^2)) * 100, 1)
  y_var <- round((pca$sdev[PCy]^2 / sum(pca$sdev^2)) * 100, 1)
  
  features <- colnames(inputData)
  
  # 4. Prepare plotting data
  get_plot_df <- function(data, slot_name) {
    df <- data
    df$slot <- slot_name
    
    # Safe color assignment
    if (is.null(color)) {
      df$color <- slot_name
    } else if (color == "outlier") {
      df$color <- CS$outliers[[featMethod]]
    } else if (color == "novelty") {
      df$color <- CS$novelties[[featMethod]]
    } else if (color == "labels") {
      df$color <- CS$labels[[slot_name]]
    } else {
      df$color <- slot_name
    }
    
    # Safe shape assignment
    if (is.null(shape)) {
      df$shape <- slot_name
    } else if (shape == "outlier") {
      df$shape <- CS$outliers[[featMethod]]
    } else if (shape == "novelty") {
      df$shape <- CS$novelties[[featMethod]]
    } else if (shape == "labels") {
      df$shape <- CS$labels[[slot_name]]
    } else {
      df$shape <- slot_name
    }
    
    df
  }
  
  plot_df <- switch(plotData,
                    test = get_plot_df(CS$features$test[[featMethod]], "test"),
                    reference = get_plot_df(CS$features$reference[[featMethod]], "reference"),
                    all = rbind(get_plot_df(CS$features$test[[featMethod]], "test"),
                                get_plot_df(CS$features$reference[[featMethod]], "reference")))
  
  # 5. Project PCA
  scores <- predict(pca, newdata = plot_df[, features])
  plot_df$x <- scores[, PCx]
  plot_df$y <- scores[, PCy]
  
  # 6. Base PCA plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = color, shape = shape)) +
    ggplot2::geom_point() +
    ggplot2::xlab(paste0("PC", PCx, " (", x_var, "%)")) +
    ggplot2::ylab(paste0("PC", PCy, " (", y_var, "%)")) +
    ggplot2::theme_bw()
  
  # 7. Optional arrows
  if (plotArrows) {
    loadings <- as.data.frame(pca$rotation[, c(PCx, PCy)])
    colnames(loadings) <- c("PCx", "PCy")
    loadings$feature <- rownames(loadings)
    if (!is.null(nLoadings)) {
      loadings$magnitude <- sqrt(loadings$PCx^2 + loadings$PCy^2)
      loadings <- loadings[order(-loadings$magnitude), ][1:nLoadings, ]
    }
    scale_factor <- max(abs(plot_df$x), abs(plot_df$y)) * 1.5
    loadings$PCx <- loadings$PCx * scale_factor
    loadings$PCy <- loadings$PCy * scale_factor
    p <- p +
      ggplot2::geom_segment(data = loadings,
                            ggplot2::aes(x = 0, y = 0, xend = PCx * 0.2, yend = PCy * 0.2),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                            color = "black") +
      ggrepel::geom_label_repel(data = loadings,
                                ggplot2::aes(x = PCx * 0.18, y = PCy * 0.18, label = feature),
                                fill = scales::alpha("white", 0.9),
                                color = "black",
                                box.padding = 0.1, point.padding = 0.2)
  }
  
  # 8. Optional loading bars
  if (plotBars) {
    loadings <- as.data.frame(pca$rotation)
    loadings$magnitude <- sqrt(loadings[, PCx]^2 + loadings[, PCy]^2)
    if (!is.null(nLoadings)) loadings <- loadings[order(-loadings$magnitude), ][1:nLoadings, ]
    xBar <- ggplot2::ggplot(loadings, ggplot2::aes(x = rownames(loadings), y = .data[[paste0("PC", PCx)]])) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                     panel.background = ggplot2::element_blank())
    yBar <- ggplot2::ggplot(loadings, ggplot2::aes(x = rownames(loadings), y = .data[[paste0("PC", PCy)]])) +
      ggplot2::geom_bar(stat = "identity") + ggplot2::coord_flip() +
      ggplot2::theme(panel.background = ggplot2::element_blank())
    return(gridExtra::grid.arrange(yBar, p, ggplot2::ggplot() + ggplot2::theme_void(), xBar,
                                   ncol = 2, heights = c(1, 0.35), widths = c(0.35, 1)))
  }
  
  return(p)
}


#' @noRd
aggregatePlotdata <- function(CS, plotData, channel, n=1000, color=NULL,
                              featMethod = NULL, nFiles = NULL){
  if (plotData == "all"){
    files <- c(CS$paths$test, CS$paths$reference)
  } else {
    files <- CS$paths[[plotData]]
  }
  
  if (!is.null(nFiles)){
    files <- sample(files, nFiles)
  }
  
  data <- list()
  for (file in files){
    ff <- readInput(CS, file, n=n)
    df <- data.frame(ff@exprs[,channel], check.names=FALSE)
    colnames(df) <- "value"
    df$index <- file
    data[[file]] <- df
  }
  data <- data.frame(dplyr::bind_rows(data), check.names = FALSE)
  data$short_name <- basename(as.character(data$index))
  
  if (is.null(color)){
    data$group <- plotData
  } else if (color == "input"){
    data$group <- ifelse(data$index %in% CS$paths$test, "test",
                         ifelse(data$index %in% CS$paths$reference, "reference", NA))
  } else if (color == "outlier"){
    data$group <- ifelse(data$index %in% CS$paths$reference, "reference",
                         ifelse(data$index %in% names(CS$outliers[[featMethod]][CS$outliers[[featMethod]] == TRUE]), 
                                "outlier", "inlier"))
  } else if (color == "novelty"){
    data$group <- ifelse(data$index %in% CS$paths$reference, "reference",
                         ifelse(data$index %in% names(CS$novelties[[featMethod]][CS$novelties[[featMethod]] == TRUE]), 
                                "novelty", "non-novelty"))
  }
  return(data)
}


#' Plot Boxplot
#'
#' @param CS CytoScan object
#' @param plotData Which data to plot based on PCA fit ("test", "reference" or "all")
#' @param n Number of cells to sample from each file
#' @param color Category used for color ("input", "outlier", "novelty")
#' @param featMethod Which features to use when plotting anomalies
#' @param nFiles Number of files to plot
#'
#' @return Boxplot
#' @import dplyr
#' @export
plotBoxplot <- function(CS, plotData, channel, n=1000, color=NULL,
                        featMethod = NULL, nFiles = NULL) {
  
  data <- aggregatePlotdata(CS, plotData=plotData, channel=channel, color=color,
                            featMethod=featMethod, nFiles=nFiles) %>%
    dplyr::mutate(short_name = basename(index))
  
  data$short_name <- factor(data$short_name,
                            levels = data %>%
                              dplyr::group_by(short_name) %>%
                              dplyr::summarise(median_value = median(value), .groups="drop") %>%
                              dplyr::arrange(median_value) %>%
                              dplyr::pull(short_name))
  
  med  <- median(data$value, na.rm = TRUE)
  iqr  <- IQR(data$value, na.rm = TRUE)
  ymin <- med - 3 * iqr
  ymax <- med + 3 * iqr
  
  ggplot2::ggplot(data, ggplot2::aes(x = short_name, y = value, fill = group)) +
    ggplot2::geom_boxplot(color = "black", outlier.shape = NA) +
    ggplot2::coord_cartesian(ylim = c(ymin, ymax)) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste("Expression for", channel), x = "", y = "") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, size = 8, hjust = 1),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 80, l = 5)
    )
}