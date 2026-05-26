Class 12 Genome Informatics
================
Melissa (PID: A19149673)
2026-05-08

## Section 4: Population Scale Analysis

> Q13. Read this file into R and determine the sample size for each
> genotype and their corresponding median expression levels for each of
> these genotypes.

``` r
expr <- read.table("https://bioboot.github.io/bggn213_W19/class-material/rs8067378_ENSG00000172057.6.txt")
head(expr)
```

    ##    sample geno      exp
    ## 1 HG00367  A/G 28.96038
    ## 2 NA20768  A/G 20.24449
    ## 3 HG00361  A/A 31.32628
    ## 4 HG00135  A/A 34.11169
    ## 5 NA18870  G/G 18.25141
    ## 6 NA11993  A/A 32.89721

``` r
nrow(expr)
```

    ## [1] 462

``` r
table(expr$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
tapply(expr$exp, expr$geno, median)
```

    ##      A/A      A/G      G/G 
    ## 31.24847 25.06486 20.07363

The sample sizes are 108 for A/A, 233 for A/G and 121 for G/G. The
median expression levels are about 31.25 for A/A, 25.06 for A/G and
20.07 for G/G.

> Q14.Generate a boxplot with a box per genotype, what could you infer
> from the relative expression value between A/A and G/G displayed in
> this plot? Does the SNP effect the expression of ORMDL3?

``` r
library(ggplot2)
```

Let’s make a boxplot

``` r
ggplot(expr, aes(geno, exp, fill=geno)) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.4)
```

![](Class-12-HW---Melissa-Bui_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The ORMDL3 expression is highest in A/A, intermediate in A/G and lowest
in G/G. Hence, this suggests the rs8067378 genotype is associated with
ORMDL3 expression, with the G/G genotype showing lower expression.
