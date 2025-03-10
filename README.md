
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RRExtraction

<!-- badges: start -->
<!-- badges: end -->

Create Relative Relief surface for distinct scales and multiple scales
including their average. To date the profile-based approaches have not
been tested. Refer to the examples in the readme for working with DEM
surfaces.

## Installation

You can install the development version of RRExtraction from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jacoblehner/RRExtraction")
```

## Example usage:

``` r
# Load libraries
library(RRExtraction)
library(terra)
library(viridis)

# Load your DEM using terra package
data <- rast("/path/to/DEM")

# Extract RR for given scale (e.g., 7)
rr7 <- relRelief(inRas = data, sc = 7) # Scale = 7

# Extract RR for scales between 2 scales with a sequence step and the average RR
rr.Avg5.19.2 <- rr.Average(inRas = data, inMn = 5, inMx = 19, stp=2)

# Example of writing features to file
writeRaster(rr7, filename="/path/to/dir/RR_7.tif", overwrite = T)
writeRaster(rs.Rasts, filename="/path/to/dir/avgRR_5_19_2.tif", overwrite = T)
```

## Example 1: Relative relief at distinct scales

``` r
library(RRExtraction)
library(terra)
#> terra 1.8.29
library(viridis)
#> Loading required package: viridisLite
f <- system.file("ex/elev_vinschgau.tif", package="terra")
r <- rast(f) 
a <- disagg(r, 2.5)
data <- resample(r, a, "bilinear")
plot(data, col=turbo(200), legend = TRUE)
```

<img src="man/figures/README-example1-1.png" width="100%" />

``` r

# Scale of X * raster resolution
rr7 <- relRelief(inRas = data, sc = 7) # Scale = 7
rr17 <- relRelief(inRas = data, sc = 17) # Scale = 17
rr51 <- relRelief(inRas = data, sc = 51) # Scale = 51

rrSamples <- c(data, rr7, rr17, rr51)
names(rrSamples) <- c('DEM', '7', '17', '51')
plot(rrSamples, col=gray(1:200/200))
```

<img src="man/figures/README-example1-2.png" width="100%" />

## Example 2: Relative relief at various scales and the average

``` r
library(RRExtraction)
library(terra)
library(viridis)
f <- system.file("ex/elev_vinschgau.tif", package="terra")
r <- rast(f) 
a <- disagg(r, 2.5)
data <- resample(r, a, "bilinear")
plot(data, col=turbo(200), legend = TRUE)
```

<img src="man/figures/README-example2-1.png" width="100%" />

``` r

# Get RR for scales between a 5 and 19 scale stepping by 2 and the average RR
rr.Avg5.19.2 <- rr.Average(inRas = data, inMn = 5, inMx = 19, stp=2)
#> Scale: 5 Scale: 7 Scale: 9 Scale: 11 Scale: 13 Scale: 15 Scale: 17 Scale: 19 
names(rr.Avg5.19.2) <- c('5', '7', '9', '11', '13', '15', '17', '19', 'Avg')
plot(rr.Avg5.19.2, col=gray(1:200/200))
```

<img src="man/figures/README-example2-2.png" width="100%" />

``` r

# Get RR for scales between a 5 and 33 scale stepping by 4 and the average RR
rr.Avg5.33.4 <- rr.Average(inRas = data, inMn = 5, inMx = 33, stp=4)
#> Scale: 5 Scale: 9 Scale: 13 Scale: 17 Scale: 21 Scale: 25 Scale: 29 Scale: 33 
names(rr.Avg5.33.4) <- c('5', '9', '13', '17', '21', '25', '29', '33', 'Avg')
plot(rr.Avg5.33.4, col=gray(1:200/200))
```

<img src="man/figures/README-example2-3.png" width="100%" />
