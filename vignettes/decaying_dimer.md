Decaying-Dimerization Reaction Set
================

<!-- github markdown built using 
rmarkdown::render("vignettes/decaying_dimer.Rmd", output_format = "github_document")
-->

The Decaying-Dimerization Reaction Set consists of three species and
four reaction channels,

``` 
      S1 --c1--> 0
 S1 + S1 --c2--> S2
      S2 --c3--> S1 + S1
      S2 --c4--> S3
```

Load package

``` r
library(GillespieSSA)
```

Define parameters

``` r
parms <- c(c1=1.0, c2=0.002, c3=0.5, c4=0.04)
```

Define initial state vector

``` r
x0 <- c(s1=10000, s2=0, s3=0)
```

Define state-change matrix

``` r
nu <- matrix(c(-1, -2, +2,  0,
                0, +1, -1, -1,
                0,  0,  0, +1),
                nrow=3,byrow=TRUE)
```

Define propensity functions

``` r
a  <- c("c1*s1", "c2*s1*s1", "c3*s2", "c4*s2")
```

Final time

``` r
tf <- 10
```

Simulation name

``` r
simName <- "Decaying-Dimerizing Reaction Set"
```

Run simulations with the Direct method

``` r
set.seed(1)
out <- ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = "D",
  simName = simName,
  verbose = TRUE,
  consoleInterval = 1
) 
```

    ## Running D method with console output every 1 time step
    ## Start wall time: 2019-07-13 21:29:48...
    ## t=0 : 10000,0,0
    ## (0.593s) t=1.000055 : 838,3776,138
    ## (0.916s) t=2.00009 : 819,3231,264
    ## (1.262s) t=3.000308 : 683,2789,388
    ## (1.513s) t=4.000468 : 679,2356,511
    ## (1.737s) t=5.00046 : 597,2008,596
    ## (1.94s) t=6.000619 : 548,1669,679
    ## (2.117s) t=7.000096 : 503,1378,735
    ## (2.284s) t=8.000416 : 444,1110,797
    ## (2.433s) t=9.002152 : 377,886,854
    ## t=10.00217 : 314,697,892
    ## tf: 10.00217
    ## TerminationStatus: finalTime
    ## Duration: 2.545 seconds
    ## Method: D
    ## Nr of steps: 30439
    ## Mean step size: 0.0003285972+/-0.0004577009
    ## End wall time: 2019-07-13 21:29:50
    ## --------------------

``` r
ssa.plot(out, show.title = TRUE, show.legend = FALSE)
```

![](decaying_dimer_files/figure-gfm/direct-1.png)<!-- -->

Run simulations with the Explict tau-leap method

``` r
set.seed(1)
out <- ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = "ETL",
  tau = 0.003,
  simName = simName,
  verbose = TRUE,
  consoleInterval = 1
) 
```

    ## Running ETL method with console output every 1 time step
    ## Start wall time: 2019-07-13 21:29:51...
    ## t=0 : 10000,0,0
    ## (0.025s) t=1.002 : 862,3808,137
    ## (0.05s) t=2.001 : 778,3307,261
    ## (0.278s) t=3 : 776,2806,391
    ## (0.304s) t=4.002 : 654,2418,491
    ## (0.329s) t=5.001 : 604,2008,570
    ## (0.355s) t=6 : 554,1678,648
    ## (0.375s) t=7.002 : 495,1396,714
    ## (0.396s) t=8.001 : 456,1146,763
    ## (0.421s) t=9 : 386,935,808
    ## t=10.002 : 351,733,846
    ## tf: 10.002
    ## TerminationStatus: finalTime
    ## Duration: 0.441 seconds
    ## Method: ETL
    ## Nr of steps: 3334
    ## Mean step size: 0.003+/-0
    ## End wall time: 2019-07-13 21:29:51
    ## --------------------

``` r
ssa.plot(out, show.title = TRUE, show.legend = FALSE)
```

![](decaying_dimer_files/figure-gfm/etl-1.png)<!-- -->

Run simulations with the Binomial tau-leap method

``` r
set.seed(1)
out <- ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = "BTL",
  simName = simName,
  verbose = TRUE,
  consoleInterval = 1
) 
```

    ## Running BTL method with console output every 1 time step
    ## Start wall time: 2019-07-13 21:29:51...
    ## t=0 : 10000,0,0
    ## (0.207s) t=1.000269 : 800,3766,169
    ## (0.254s) t=2.001535 : 811,3203,319
    ## (0.294s) t=3.001799 : 771,2714,429
    ## (0.319s) t=4.001483 : 656,2335,530
    ## (0.345s) t=5.002348 : 584,1959,624
    ## (0.363s) t=6.003555 : 502,1650,697
    ## (0.379s) t=7.004771 : 459,1358,765
    ## (0.396s) t=8.002029 : 400,1113,815
    ## (0.407s) t=9.002283 : 353,891,850
    ## t=10.00815 : 301,713,887
    ## tf: 10.00815
    ## TerminationStatus: finalTime
    ## Duration: 0.416 seconds
    ## Method: BTL
    ## Nr of steps: 3014
    ## Mean step size: 0.003320553+/-0.002322824
    ## End wall time: 2019-07-13 21:29:52
    ## --------------------

``` r
ssa.plot(out, show.title = TRUE, show.legend = FALSE)
```

![](decaying_dimer_files/figure-gfm/btl-1.png)<!-- -->

Run simulations with the Optimized tau-leap method

``` r
set.seed(1)
out <- ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = "OTL",
  simName = simName,
  verbose = TRUE,
  consoleInterval = 1
) 
```

    ## Running OTL method with console output every 1 time step
    ## Start wall time: 2019-07-13 21:29:52...
    ## t=0 : 10000,0,0
    ## (0.011s) t=1.007136 : 838,3785,146
    ## (0.021s) t=2.016863 : 812,3233,281
    ## (0.031s) t=3.002185 : 695,2795,390
    ## (0.039s) t=4.012203 : 664,2384,482
    ## (0.097s) t=5.016351 : 568,2046,578
    ## (0.152s) t=6.001846 : 547,1712,646
    ## (0.162s) t=7.018118 : 455,1453,702
    ## (0.17s) t=8.025641 : 457,1163,755
    ## (0.177s) t=9.006531 : 396,923,799
    ## t=10.01378 : 364,716,841
    ## tf: 10.01378
    ## TerminationStatus: finalTime
    ## Duration: 0.184 seconds
    ## Method: OTL
    ## Nr of steps: 406
    ## Mean step size: 0.02466448+/-0.003695297
    ## Nr suspended tau leaps: 0(0%)
    ## End wall time: 2019-07-13 21:29:52
    ## --------------------

``` r
ssa.plot(out, show.title = TRUE, show.legend = FALSE)
```

![](decaying_dimer_files/figure-gfm/otl-1.png)<!-- -->
