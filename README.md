<img src="man/figures/logo.png" align="center" height="200"/>

`CytoScan` is an R package for detecting "flagging" of anomalous flow cytometry files (".FCS"). 

## What is an anomaly?

In CytoScan, we consider two types of anomalies:
1) Outliers
2) Novelties

Outlier files that are defined as anomalies with respect to all other files *within* the same dataset.

Novelties are files that are defined as anomalies with respect to another dataset. We advice using a high-quality, validated dataset.

## How to use CytoScan

You can install `CytoScan` directly from Github using the `devtools` library in R. 

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("AUMC-HEMA/CytoScan")
```

After installation, we refer to the documentation for further use.
