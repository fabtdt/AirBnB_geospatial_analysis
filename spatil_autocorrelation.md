Analyze and Test Spatial Autocorrelation of AirBnB Prices
================
Fabio Taddei Dalla Torre - Mat.: 214924

It is plausible that the prices of AirBnB are somehow spatially
autocorrelated. The analysis, to test whether this spatial
autocorrelation exist or not, will be carried out using the data of the
AirBnB that has been used in the Jupyter notebook. In the notebook is
possible to see how this dataset has been created.

The analysis will be carried on in the following way:

-   Loading the data and computing centroids
-   Build spatial weight matrices given different definition of
    neighborhood for consistency
-   Calculate Global Moran’s I index
-   Build Moran Scatterplot to look for local auto correlation
-   Calculate local Moran’s I index

**Loading the data**

``` r
# Loading the Cambrisge shapefile that contains all the data
cambridge <- readOGR("Cambridge", "Cambridge")
```

    ## OGR data source with driver: ESRI Shapefile 
    ## Source: "C:\Users\fabio\Documents\UNI\Data_Science\Geospatial_analysis\project_exam\Cambridge_Fabio_TaddeiDallaTorre\Cambridge", layer: "Cambridge"
    ## with 13 features
    ## It has 3 fields

``` r
# Check the data
#View(cambridge@data)
```

First of all I need to compute the centroids of the ares in order to
have a representative point for each districts.

``` r
centroids <- coordinates(cambridge)
```

## 1. Creating the spatial weight matrix

The first step in testing if there is spatial auto correlation between
the different prices is to define a spatial weigh matrix.

### 1.1 Defining spatial neighbors

Various definitions of neighborhood are possible, hence we need to
perform different solutions:

-   k-Nearest neighbors
-   Critical cut-off distance
-   Contiguity-based neighborhood
-   Free-Form spatial weight matrices

Different solution will be performed in order to test for consistency.

``` r
#To visualize the different neighborhood amd centroids
plot(cambridge, border = "purple")
points(centroids, cex=0.8)
```

![](spatil_autocorrelation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

**1.1.1 k-Nearest neighbors**

This definition of neighborhood ensures that each spatial unit has
exactly the same number k of neighbors. For example for k = 1 the
neighbor is the closest units.

I will use k = 1 and k = 4.

``` r
# longlat = F it will use the euclidean distnace, longlat = T it will take into consideraiton for the curvature of the Earth
knn1cb <- knn2nb(knearneigh(centroids,k=1,longlat=T)) 
knn4cb <- knn2nb(knearneigh(centroids,k=4,longlat=T))
```

``` r
#Cheking visually the 1-nn

plot(cambridge, border="grey")
plot(knn1cb, centroids, add=TRUE)
```

![](spatil_autocorrelation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

**1.1.2 Critical cut-off distance**

Two spatial units are considered as neighbors if their distance is
equal, or less than equal, to a certain fixed distance which represents
a critical cut-off. For this method first we need to compute the maximum
distance between two centroids, by doing this we assure that all regions
have at least one neighbor. This value will the minimum cut-off distance

Computing the minimum cut-off distance

``` r
# we need to compute the distance among neigbour and find the maximum
min_cut_off <- max(unlist(nbdists(knn1cb, centroids, longlat=T)))
min_cut_off
```

    ## [1] 1.280868

The minimum cut-off distance has been resulted 1.280868. Moreover I will
create different concept of neighborhood by using different cut-off for
consistency.

``` r
dnb1.3 <- dnearneigh(centroids, 0, 1.3, longlat=TRUE); dnb1.3
```

    ## Neighbour list object:
    ## Number of regions: 13 
    ## Number of nonzero links: 28 
    ## Percentage nonzero weights: 16.56805 
    ## Average number of links: 2.153846

``` r
dnb2 <- dnearneigh(centroids, 0, 2, longlat=TRUE); dnb2
```

    ## Neighbour list object:
    ## Number of regions: 13 
    ## Number of nonzero links: 54 
    ## Percentage nonzero weights: 31.95266 
    ## Average number of links: 4.153846

``` r
dnb5 <- dnearneigh(centroids, 0, 5, longlat=TRUE); dnb5
```

    ## Neighbour list object:
    ## Number of regions: 13 
    ## Number of nonzero links: 138 
    ## Percentage nonzero weights: 81.6568 
    ## Average number of links: 10.61538

**1.1.3 Contiguity-based neighborhood**

In this case two units are considered neighbor if they share a common
boundary.

``` r
contnb <- poly2nb(cambridge, queen=T)
# Plotting the graph pf centroids in order to visualise the results
plot(cambridge, border="grey")
plot(contnb, centroids, add=TRUE)
```

![](spatil_autocorrelation_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

**1.1.4 Free-Form spatial weight matrices** In this case the spatial
weight matrix will be calculated as an invers function of the distance.

``` r
distM <- as.matrix(dist(centroids)) #distance matrix
# Calcualting weights
W1 <- 1/(1+(distM)); diag(W1) <- 0

#Row-standardize
W1s <- W1/rowSums(W1) 
```

### 1.2 Defining spatial weights Matrix

Given the previous definition it is possible to compute different
version of spatial wight matrix.

``` r
# SWM made by knn
knn1cb.listw <- nb2listw(knn1cb,style="W")
knn4cb.listw <- nb2listw(knn4cb,style="W")

# SWM made by critical cut-off distance
dnb1.3.listw <- nb2listw(dnb1.3,style="W")
dnb2.listw <- nb2listw(dnb2,style="W")
dnb5.listw <- nb2listw(dnb5,style="W")

# SWM made by contiguity base approach
contnb.listw <- nb2listw(contnb,style="W")

# SWM made by Free-Form 
listW1s <- mat2listw(W1s)
```

## 2. Test for spatial autocorrelation

### 2.1 Global spatial autocorrelation

The existence of a global spatial autocorrelation among prices in the
different neighborhood will be carried on by calculating the Global
Moran’s I index given the previously computed spatial weight matrix.

In particular the Moran’s I test:

-   H0: no spatial autocorrelation
-   H1: spatial autocorrelation

``` r
# With the option *randomisation = FALSE* the test is performing assuming under normality assumption
# K-nn
moran.test(cambridge$mean_price, knn1cb.listw, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: knn1cb.listw    
    ## 
    ## Moran I statistic standard deviate = -1.1532, p-value = 0.8756
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.46148339       -0.08333333        0.10752442

``` r
moran.test(cambridge$mean_price, knn4cb.listw, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: knn4cb.listw    
    ## 
    ## Moran I statistic standard deviate = -0.42795, p-value = 0.6657
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.14830089       -0.08333333        0.02304640

``` r
# Critical cut-off distance
moran.test(cambridge$mean_price, dnb1.3.listw, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: dnb1.3.listw    
    ## 
    ## Moran I statistic standard deviate = -1.8912, p-value = 0.9707
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.56329557       -0.08333333        0.06440781

``` r
moran.test(cambridge$mean_price, dnb2.listw, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: dnb2.listw    
    ## 
    ## Moran I statistic standard deviate = -0.67626, p-value = 0.7506
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.19307854       -0.08333333        0.02633572

``` r
moran.test(cambridge$mean_price, dnb5.listw, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: dnb5.listw    
    ## 
    ## Moran I statistic standard deviate = 0.15096, p-value = 0.44
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.07725456       -0.08333333        0.00162155

``` r
# Contiguity based approach
moran.test(cambridge$mean_price, contnb.listw, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: contnb.listw    
    ## 
    ## Moran I statistic standard deviate = -1.1089, p-value = 0.8663
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.26037467       -0.08333333        0.02549196

``` r
#Free-Form
moran.test(cambridge$mean_price, listW1s, randomisation=FALSE)
```

    ## 
    ##  Moran I test under normality
    ## 
    ## data:  cambridge$mean_price  
    ## weights: listW1s    
    ## 
    ## Moran I statistic standard deviate = -0.093837, p-value = 0.5374
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##     -8.348794e-02     -8.333333e-02      2.714528e-06

In every scenario the P-Value is quite high, that means the null
hypothesis cannot be rejected and so that there are evidence of no
global spatial autocorrelation.

It is possible to perform the test under the randomization assumption.
In this case the observed values of x are randomly permuted despite the
underlying distribution in the population.

``` r
# K-nn
moran.test(cambridge$mean_price, knn1cb.listw, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: knn1cb.listw    
    ## 
    ## Moran I statistic standard deviate = -1.2176, p-value = 0.8883
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.46148339       -0.08333333        0.09645247

``` r
moran.test(cambridge$mean_price, knn4cb.listw, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: knn4cb.listw    
    ## 
    ## Moran I statistic standard deviate = -0.45249, p-value = 0.6745
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.14830089       -0.08333333        0.02061433

``` r
# Critical cut-off distance
moran.test(cambridge$mean_price, dnb1.3.listw, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: dnb1.3.listw    
    ## 
    ## Moran I statistic standard deviate = -1.9977, p-value = 0.9771
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.56329557       -0.08333333        0.05772075

``` r
moran.test(cambridge$mean_price, dnb2.listw, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: dnb2.listw    
    ## 
    ## Moran I statistic standard deviate = -0.71501, p-value = 0.7627
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.19307854       -0.08333333        0.02355866

``` r
moran.test(cambridge$mean_price, dnb5.listw, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: dnb5.listw    
    ## 
    ## Moran I statistic standard deviate = 0.15868, p-value = 0.437
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##      -0.077254560      -0.083333333       0.001467546

``` r
# Contiguity based approach
moran.test(cambridge$mean_price, contnb.listw, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: contnb.listw    
    ## 
    ## Moran I statistic standard deviate = -1.172, p-value = 0.8794
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##       -0.26037467       -0.08333333        0.02282006

``` r
#Free-Form
moran.test(cambridge$mean_price, listW1s, randomisation=TRUE)
```

    ## 
    ##  Moran I test under randomisation
    ## 
    ## data:  cambridge$mean_price  
    ## weights: listW1s    
    ## 
    ## Moran I statistic standard deviate = -0.098714, p-value = 0.5393
    ## alternative hypothesis: greater
    ## sample estimates:
    ## Moran I statistic       Expectation          Variance 
    ##     -8.348794e-02     -8.333333e-02      2.452943e-06

Even in this case the P-values are quite high compared to the usual
threshold value.

### 2.2 Local Spatial Autocorrelation

**Moran’s scatterplot**

Local spatial autocorrelation can be investigated by using the Moran’s
Scatterplot.

``` r
# Producing the Moran scatterpot by using the different spatial weith matrices calcluated at point 1.2
par(mfrow=c(2,2))
# K.nn
  g1 <- moran.plot(cambridge$mean_price, listw=knn1cb.listw, main="Moran scatterplot 1-nn", return_df=F)
  g2 <- moran.plot(cambridge$mean_price, listw=knn1cb.listw, main="Moran scatterplot 4-nn", return_df=F)
# Critical cut-off distance
  g3 <- moran.plot(cambridge$mean_price, listw=dnb1.3.listw, main="Moran scatterplot 1.3 km", return_df=F)
  g4 <- moran.plot(cambridge$mean_price, listw=dnb2.listw, main="Moran scatterplot 2 Km", return_df=F)
```

![](spatil_autocorrelation_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
  g5 <- moran.plot(cambridge$mean_price, listw=dnb5.listw, main="Moran scatterplot 5 Km", return_df=F)
# Contiguity based approach
  g6 <- moran.plot(cambridge$mean_price, listw=contnb.listw, main="Moran scatterplot Contiguity-Based", return_df=F)
# Free-Form
  g7 <- moran.plot(cambridge$mean_price, listW1s, main="Moran scatterplot Free-Form", return_df=F)
```

![](spatil_autocorrelation_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

On the horizontal axis of the scatterplot there is the variable of
interest (*mean\_price*) while on the vertical axis there is the spatial
lag of this variable. In this case all the scatterplots present a
negative inclined line, that means a possible negative correlation. In
all the cases the inclination of the line cannot suggest any clear sign
of local correlation.

**Local Moran’s I**

With the Local Moran’s I is possible to prove the statistical
significance of the insight given by the Moran Scatterplot

``` r
# K-nn
lmI_1nn <- localmoran(cambridge$mean_price, knn1cb.listw)
head(lmI_1nn)
```

    ##           Ii        E.Ii    Var.Ii       Z.Ii Pr(z > 0)
    ## 1  0.3838011 -0.08333333 0.7852445  0.5271563 0.2990425
    ## 2 -0.1769969 -0.08333333 0.7852445 -0.1056983 0.5420892
    ## 3 -0.2393658 -0.08333333 0.7852445 -0.1760810 0.5698849
    ## 4  0.0704648 -0.08333333 0.7852445  0.1735596 0.4311058
    ## 5 -0.1769969 -0.08333333 0.7852445 -0.1056983 0.5420892
    ## 6  0.3838011 -0.08333333 0.7852445  0.5271563 0.2990425

``` r
lmI_4nn <- localmoran(cambridge$mean_price, knn4cb.listw)
head(lmI_4nn)
```

    ##              Ii        E.Ii    Var.Ii        Z.Ii Pr(z > 0)
    ## 1 -0.0003359197 -0.08333333 0.1569407  0.20950617 0.4170266
    ## 2  0.0104137072 -0.08333333 0.1569407  0.23664091 0.4064677
    ## 3  0.1091419241 -0.08333333 0.1569407  0.48585555 0.3135348
    ## 4 -0.0933237773 -0.08333333 0.1569407 -0.02521837 0.5100596
    ## 5 -0.4574341342 -0.08333333 0.1569407 -0.94432372 0.8274979
    ## 6 -1.2485563572 -0.08333333 0.1569407 -2.94131352 0.9983659

``` r
# Critical cut-off distance
lmI_1.3 <- localmoran(cambridge$mean_price, dnb1.3.listw)
head(lmI_1.3)
```

    ##              Ii        E.Ii    Var.Ii        Z.Ii Pr(z > 0)
    ## 1 -0.0003359197 -0.08333333 0.1569407  0.20950617 0.4170266
    ## 2 -0.0532660440 -0.08333333 0.3663753  0.04967424 0.4801910
    ## 3 -0.2393658173 -0.08333333 0.7852445 -0.17608103 0.5698849
    ## 4  0.0834445515 -0.08333333 0.3663753  0.27553412 0.3914530
    ## 5 -0.4574341342 -0.08333333 0.1569407 -0.94432372 0.8274979
    ## 6 -1.3085008475 -0.08333333 0.2267522 -2.57288113 0.9949572

``` r
lmI_2 <- localmoran(cambridge$mean_price, dnb2.listw)
head(lmI_2)
```

    ##            Ii        E.Ii     Var.Ii        Z.Ii Pr(z > 0)
    ## 1 -0.11024874 -0.08333333 0.08712913 -0.09118412 0.5363269
    ## 2  0.06442567 -0.08333333 0.11505374  0.43561576 0.3315578
    ## 3  0.10914192 -0.08333333 0.15694066  0.48585555 0.3135348
    ## 4  0.01144422 -0.08333333 0.11505374  0.27941847 0.3899618
    ## 5 -0.42028401 -0.08333333 0.11505374 -0.99338125 0.8397379
    ## 6 -0.71387262 -0.08333333 0.08712913 -2.13614374 0.9836662

``` r
lmI_5 <- localmoran(cambridge$mean_price, dnb5.listw)
head(lmI_5)
```

    ##             Ii        E.Ii     Var.Ii       Z.Ii Pr(z > 0)
    ## 1  0.019523958 -0.08333333 0.02366409  0.6686366 0.2518637
    ## 2 -0.003825545 -0.08333333 0.01731759  0.6041795 0.2728622
    ## 3 -0.017726581 -0.08333333 0.04058810  0.3256486 0.3723451
    ## 4 -0.009013399 -0.08333333 0.01731759  0.5647570 0.2861195
    ## 5 -0.056868978 -0.08333333 0.01731759  0.2011026 0.4203092
    ## 6 -0.212846209 -0.08333333 0.03127990 -0.7322853 0.7680028

``` r
# Contiguity based approach
lmI_cb <- localmoran(cambridge$mean_price, contnb.listw)
head(lmI_cb)
```

    ##            Ii        E.Ii     Var.Ii        Z.Ii Pr(z > 0)
    ## 0 -0.11024874 -0.08333333 0.08712913 -0.09118412 0.5363269
    ## 1  0.02406978 -0.08333333 0.15694066  0.27111225 0.3931524
    ## 2  0.10914192 -0.08333333 0.15694066  0.48585555 0.3135348
    ## 3 -0.09332378 -0.08333333 0.15694066 -0.02521837 0.5100596
    ## 4 -0.41155414 -0.08333333 0.06718297 -1.26629968 0.8972971
    ## 5 -1.24855636 -0.08333333 0.15694066 -2.94131352 0.9983659

``` r
# Free - Form
lmI_cb <- localmoran(cambridge$mean_price, listW1s)
head(lmI_cb)
```

    ##             Ii        E.Ii     Var.Ii       Z.Ii Pr(z > 0)
    ## 0 -0.008136800 -0.08333333 0.01733918  0.5710624 0.2839787
    ## 1 -0.003462896 -0.08333333 0.01732587  0.6067901 0.2719951
    ## 2 -0.052959655 -0.08333333 0.01734378  0.2306352 0.4087991
    ## 3 -0.009250426 -0.08333333 0.01733090  0.5627396 0.2868061
    ## 4 -0.058656453 -0.08333333 0.01732763  0.1874652 0.4256480
    ## 5 -0.151030601 -0.08333333 0.01733906 -0.5141127 0.6964134

From the output of all the Local Moran’s I test it can be seen that
there are no evidence against the null assumption of no spatial
autocorrelation hence there are no evidence that could suggest local
spatial autocorrelation.

As already stated in the Jupyter notebook this phenomena of no spatial
autocorrelation could be the consequence of the proximity of the city of
Cambridge to the Boston one. Moreover by looking at the map of the are
it is possible to see that, the East-Cambridge district that is one of
the expensive, is the closest to the city of Boston. Because of that it
is possible that the prices of the AirBnB are influenced by the one of
the Boston, hence in order to investigate the presence of spatial
autocorrelation of prices should be better to do this by taking into
account the prices of Boston as well.
