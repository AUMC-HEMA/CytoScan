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
plotPCA <- function(CS, featMethod, fitData, plotData, PCx = 1, PCy = 2,
                    color = NULL, shape = NULL, plotBars = FALSE, 
                    plotArrows = FALSE, nLoadings = NULL){
  # Fit PCA on input data
  if (fitData == "test"){
    inputData <- CS$features$test[[featMethod]]
  } else if (fitData == "reference"){
    inputData <- CS$features$reference[[featMethod]]
  } else if (fitData == "all"){
    inputData <- rbind(CS$features$test[[featMethod]], CS$features$ref[[featMethod]])
  }
  pca <- stats::prcomp(inputData, center = TRUE, scale = TRUE)
  x_var <- pca$sdev[PCx]^2
  y_var <- pca$sdev[PCy]^2
  
  features <- colnames(inputData)
  
  if (plotData == "test" | plotData == "all"){
    testData <- CS$features$test[[featMethod]]
    testData$slot <- "test"
    if (is.null(color)){
      testData$color <- "test"
    } else if (color == "outlier"){
      testData$color <- CS$outliers[[featMethod]]
    } else if (color == "novelty"){
      testData$color <- CS$novelties[[featMethod]]
    } else if (color == "labels"){
      testData$color <- CS$labels$test
    }
    if (is.null(shape)){
      testData$shape <- "test"
    } else if (shape == "outlier"){
      testData$shape <- CS$outliers[[featMethod]]
    } else if (shape == "novelty"){
      testData$shape <- CS$novelties[[featMethod]]
    } else if (shape == "labels"){
      testData$shape <- CS$labels$test
    }
  }
  if (plotData == "reference" | plotData == "all"){
    referenceData <- CS$features$reference[[featMethod]]
    referenceData$slot <- "reference"
    if (is.null(color)){
      referenceData$color <- "reference"
    } else if (color == "outlier"){
      referenceData$color <- "reference"
    } else if (color == "novelty"){
      referenceData$color <- "reference"
    } else if (color == "labels"){
      referenceData$color <- CS$labels$reference
    }
    if (is.null(shape)){
      referenceData$shape <- "reference"
    } else if (shape == "outlier"){
      referenceData$shape <- "reference"
    } else if (shape == "novelty"){
      referenceData$shape <- "reference"
    } else if (shape == "labels"){
      referenceData$shape <- CS$labels$reference
    }
  }
  if (plotData == "test"){
    plotData <- testData
  } else if (plotData == "reference"){
    plotData <- referenceData
  } else if (plotData == "all"){
    plotData <- rbind(testData, referenceData)
  }
  output <- stats::predict(pca, newdata = plotData[, features])
  plotData$x <- output[, paste0("PC", as.character(PCx))]
  plotData$y <- output[, paste0("PC", as.character(PCy))]
  plotData <- data.frame(plotData, check.names = FALSE)
  
  PCAPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y, color = color, 
                                                    shape = shape)) +
    ggplot2::geom_point() +
    ggplot2::xlab(paste0("PC", PCx, " (", round(x_var, 1), "%)")) +
    ggplot2::ylab(paste0("PC", PCy, " (", round(y_var, 1), "%)")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(colour = "black", 
                                                        fill = NA))
  
  if (plotArrows) {
    loadings <- as.data.frame(pca$rotation[, c(PCx, PCy)])
    colnames(loadings) <- c("PCx", "PCy")
    loadings$feature <- rownames(loadings)
    if (!is.null(nLoadings)) {
      loadings$magnitude <- sqrt(loadings$PCx^2 + loadings$PCy^2)
      loadings <- loadings[order(-loadings$magnitude), ][1:nLoadings, ]
    }
    # Normalize and scale loadings for better visualization
    scale_factor <- min(max(abs(plotData$x)), max(abs(plotData$y))) * 1.5
    loadings$PCx <- loadings$PCx * scale_factor
    loadings$PCy <- loadings$PCy * scale_factor
    
    # Update geom_label_repel to position labels closer to the end of the arrows
    PCAPlot <- PCAPlot +
      ggplot2::geom_segment(data = loadings, 
                            ggplot2::aes(x = 0, y = 0, 
                                         xend = PCx * max(plotData$x) * 0.2, 
                                         yend = PCy * max(plotData$y) * 0.2),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), 
                            color = "black", size = 0.5) +
      ggrepel::geom_label_repel(data = loadings, 
                                ggplot2::aes(x = PCx * max(plotData$x) * 0.2 * 0.9, 
                                             y = PCy * max(plotData$y) * 0.2 * 0.9, 
                                             label = feature),
                                fill = scales::alpha("white", 1), 
                                color = "black", 
                                box.padding = 0.1,
                                point.padding = 0.2)
  }
  if (plotBars) {
    loadings <- as.data.frame(pca$rotation)
    loadings$magnitude <- sqrt(loadings[, paste0("PC", as.character(PCx))]^2 + 
                                 loadings[, paste0("PC", as.character(PCy))]^2)
    if (!is.null(nLoadings)) {
      loadings <- loadings[order(-loadings$magnitude), ][1:nLoadings, ]
    } else {
      loadings <- loadings[order(-loadings[, paste0("PC", as.character(PCx))]), ]
    }
    sorted_PCx <- loadings[order(abs(loadings[, paste0("PC", as.character(PCx))]),
                                 decreasing = TRUE), ]
    sorted_PCy <- loadings[order(abs(loadings[, paste0("PC", as.character(PCy))]),
                                 decreasing = FALSE), ]
    loadings_x <- data.frame(Variable = rownames(sorted_PCx), 
                             Loading = sorted_PCx[, paste0("PC", as.character(PCx))])
    loadings_x$Variable <- factor(loadings_x$Variable, levels = loadings_x$Variable)
    loadings_y <- data.frame(Variable = rownames(sorted_PCy), 
                             Loading = sorted_PCy[, paste0("PC", as.character(PCy))])
    loadings_y$Variable <- factor(loadings_y$Variable, levels = loadings_y$Variable)
    
    xBar <- ggplot2::ggplot(loadings_x, ggplot2::aes(x = Variable, y = Loading)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, 
                                                         hjust = 1),
                     axis.title.x = ggplot2::element_blank(), 
                     axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour = "black", 
                                                          fill = NA))
    yBar <- ggplot2::ggplot(loadings_y, ggplot2::aes(x = Variable, y = Loading)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                     axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour = "black", 
                                                          fill = NA))
    pGrid <- gridExtra::grid.arrange(yBar, PCAPlot, 
                                     ggplot2::ggplot() + ggplot2::theme_void(), 
                                     xBar,
                                     ncol = 2,
                                     heights = c(1, 0.35),
                                     widths = c(0.35, 1))
    return(pGrid)
  } else {
    return(PCAPlot)
  }
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