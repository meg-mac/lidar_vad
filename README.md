# lidar_vad

Lidar VAD Raw Data Conversion to Wind Speed and Direction

Coded for Halo Stream Line XR Doppler LiDAR system output but can be modified for other data structures


To determine wind speed and direction, the lidar must conduct measurements from three different directions. To do this, the azimuth angle (Θ) must change continuously at a fixed elevation angle (ϕ) with constant angular speed (ω). In the case of horizontal homogeneous wind, the mean radial wind velocity will have a sine wave dependence on the azimuth angle. However, turbulence and other factors will alter this relationship and thus it was necessary to minimize errors through a least-squares technique (Banakh and Smalikho, 2013)

Banakh, V., & Smalikho, I. (2013). Coherent Doppler wind lidars in a turbulent atmosphere. Artech House.
## Installation

You can install the development version from
[GitHub](https://github.com/)
with:

<!-- the released version of rvad from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rvad")
```

And -->

``` r
#install.packages("devtools")
#devtools::install_github("meg-mac/lidar_vad")
```

## Example

