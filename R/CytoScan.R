#' Initialize CytoScan object
#' 
#' Function to initialize CytoScan workflow. Returns CytoScan object
#' 
#' @return Default initialized CytoScan object
#' 
#' @export
CytoScan <- function(){
  CS <- list()
  class(CS) <- "CytoScan"
  CS[["preprocessFunction"]] <- processInput
  CS[["parallel"]] <- list("parallelVars" = c("channels", "CS", "readInput"),
                           "parallelPackages" = c("flowCore"))
  return(CS)
}

 
#' Add additional annotations for test data 
#' 
#' @param CS CytoScan object
#' @param labels Vector of numeric or string labels
#'
#' @return CytoScan object
#'
#' @seealso \code{\link{addReferencelabels}}
#'
#' @export
addTestlabels <- function(CS, labels){
  CS[["labels"]][["test"]] <- labels
  return(CS)
}


#' Add additional annotations for reference data 
#' 
#' @param CS CytoScan object
#' @param labels Vector of numeric or string labels
#'
#' @return CytoScan object
#'
#' @seealso \code{\link{addTestlabels}}
#'
#' @export
addReferencelabels <- function(CS, labels){
  CS[["labels"]][["reference"]] <- labels
  return(CS)
}


#' @noRd
processInput <- function(ff){
  spill <- ff@description$SPILL
  ff <- flowCore::compensate(ff, spill)
  ff <- flowCore::transform(ff, flowCore::transformList(colnames(spill), 
                                            flowCore::arcsinhTransform(a = 0, 
                                                                       b = 1/150, 
                                                                       c = 0)))
  return(ff)
}


#' @export
readInput <- function(CS, path, n = NULL){
  set.seed(42)
  ff <- flowCore::read.FCS(path, which.lines = n)
  ff <- CS[["preprocessFunction"]](ff)
  return(ff)
}


#' @noRd
addData <- function(CS, input, type, read, reload, aggSize){
  # Check if the paths are already stored in the CytoScan object
  if (type %in% names(CS$paths) & !reload){
    for (path in input){
      if (!path %in% CS$paths[[type]]){
        # message(paste("Concatenating additional file path", path))
        CS$paths[[type]] <- c(CS$paths[[type]], path)
        if (read == TRUE){
          ff <- readInput(CS, path, CS$nAgg)
          CS$data[[type]][[path]] <- data.frame(ff@exprs, check.names = FALSE)
          }
        }
      }
    } else {
      CS$paths[[type]] <- input
      if (read == TRUE){
        CS$aggSize <- aggSize
        CS$nAgg <- ceiling(aggSize / length(input))
        for (path in input){
          ff <- readInput(CS, path, CS$nAgg)
          CS$data[[type]][[path]] <- data.frame(ff@exprs, check.names = FALSE)
        }
      }
    }
  return(CS)
}


#' Add test data to CytoScan object
#'
#' @param CS CytoScan object
#' @param input List of FCS file paths
#' @param read Whether load subsampled data already (default = FALSE)
#' @param reload Reload all data (default = FALSE)
#' @param aggSize Size of aggregate sample (default = 10000 cells)
#'
#' @return CytoScan object
#'
#' @seealso \code{\link{addReferencedata}}
#' 
#' @export
addTestdata <- function(CS, input, read = FALSE, reload = FALSE,
                        aggSize = 10000){
  CS <- addData(CS, input, "test", read, reload, aggSize)
  return(CS)
}


#' Add reference data to CytoScan object
#'
#' @param CS CytoScan object
#' @param input List of FCS file paths
#' @param read Whether load subsampled data already (default = FALSE)
#' @param reload Reload all data (default = FALSE)
#' @param aggSize Size of aggregate sample (default = 10000 cells)
#'
#' @return CytoScan object
#'
#' @seealso \code{\link{addTestdata}}
#' 
#' @export
addReferencedata <- function(CS, input, read = FALSE, reload = FALSE,
                             aggSize = 10000){
  CS <- addData(CS, input, "reference", read, reload, aggSize)
  return(CS)
}


#' Flag anomalies based on generated features 
#'
#' @param CS CytoScan object
#' @param featMethod Which features to use for anomaly detection
#' @param flagMethod Which type of anomaly detection to use ("outlier" or "novelty")
#' @param outlier_threshold Isolation score threshold for outlier detection.
#' @param novelty_threshold GMM quantile threshold for novelty detection.
#'
#' @return CytoScan object
#' 
#' @seealso \code{\link{generateFeatures}}
#'
#' @export
#' @import mclust
Flag <- function(CS, featMethod, flagMethod, outlier_threshold=0.5,
                 novelty_threshold=0.05){
  testFeatures <- CS$features$test[[featMethod]]
  if (flagMethod == "outlier" | flagMethod == "outliers"){
    forest <- isotree::isolation.forest(testFeatures, sample_size = 1,
                                        ntrees = 1000,
                                        ndim = 1, seed = 42)
    scores <- stats::predict(forest, testFeatures)
    # Scores > 0.5 are most likely outliers according to the original paper
    outliers <- as.factor(ifelse(scores >= outlier_threshold, TRUE, FALSE))
    CS$outliers[[featMethod]] <- outliers
  }
  if (flagMethod == "novelty" | flagMethod == "novelties"){
    refFeatures <- CS$features$reference[[featMethod]]
    # Fit GMM on reference features
    gmm <- mclust::Mclust(refFeatures, verbose=FALSE)
    llrTrain <- log(mclust::dens(refFeatures, gmm$modelName, parameters = gmm$parameters))
    threshold <- quantile(llrTrain, novelty_threshold)
    llrTest <- log(mclust::dens(testFeatures, gmm$modelName, parameters = gmm$parameters))
    novelties <- llrTest <= threshold
    names(novelties) <- CS$paths$test
    CS$novelties[[featMethod]] <- novelties
  }
  return(CS)
}
