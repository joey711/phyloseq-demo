
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# Fast Parallel UniFrac (in R)

## Background
A common tool in microbial ecology studies involving many samples is to calculate the "UniFrac" distance between all pairs of samples, and then perform various analyses on the resulting distance matrix.

[The phyloseq package](http://joey711.github.com/phyloseq/) includes a native `R` implementation of the better, faster, cleaner [Fast UniFrac algorithm](http://www.nature.com/ismej/journal/v4/n1/abs/ismej200997a.html). For legacy reasons and comparison, it also includes the original UniFrac algorithm, although there is unlikely to be a reason to use the original implementation because it is slower, and both approaches arrive at the same result. There are also two very different types of the standard UniFrac calculation: 

[Weighted UniFrac](http://aem.asm.org/content/73/5/1576.abstract) - which does take into account differences in abundance of species between samples, but takes longer to calculate; and

[Unweighted UniFrac](http://aem.asm.org/content/71/12/8228.short) - which only considers the presence/absence of species between sample pairs. 

Both can be useful, and share slightly different insight. Both weighted and unweighted UniFrac are included, and all UniFrac calculations have the option of running "in parallel" for faster results on computers that have multiple cores/processors available.

## Load phyloseq
Check your version number.

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.3.20'
```


## Examples
In phyloseq, all variants of the UniFrac distance can be called with the same generic function, `UniFrac`.
```R
UniFrac(physeq, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
```
where `physeq` is your `phyloseq-class` experiment-level object, imported or constructed by the phyloseq package.

Here is an example taken from the `UniFrac` documentation showing how to calculate it for the 3-sample "esophagus" dataset:

Perform UniFrac on esophagus data

```r
data("esophagus")
(y <- UniFrac(esophagus, TRUE))
```

```
##        B      C
## C 0.2035       
## D 0.2603 0.2477
```

```r
UniFrac(esophagus, TRUE, FALSE)
```

```
##        B      C
## C 0.1050       
## D 0.1401 0.1422
```

```r
UniFrac(esophagus, FALSE)
```

```
##        B      C
## C 0.5176       
## D 0.5182 0.5422
```


If you want to compare with the unweighted-unifrac available in picante

```r
picante::unifrac(as(t(otuTable(esophagus)), "matrix"), phy_tree(esophagus))
```

```
##        B      C
## C 0.5176       
## D 0.5182 0.5422
```



## Performance Comparison

The following graphic shows the computation time required to calculate the distance matrix for a 21-sample experiment at varying numbers of species. As you can see, for datasets with less than 1000 species (or more accurately, "OTUs"), the choice of algorithm has only inconsequential effect on computation time, with all approaches returning a result in under 10 seconds. However as more and more species are included (more complex phylogenetic tree) in the dataset, the computation time increases exponentially; and at a much faster rate for variants of the original UniFrac algorithm. As pointed out in the article describing Fast UniFrac, there is much less difference in performance between unweighted/weighted Fast UniFrac, making weighted UniFrac distances much more accessible (than previously) for large datasets. Also note that the computation time for all UniFrac methods will increase combinatorically as the number of samples increases. You can quickly note the number of pairwise UniFrac calculations necessary with the following:


```r
data("GlobalPatterns")
dim(combn(nsamples(GlobalPatterns), 2))[2]
```

```
## [1] 325
```


If you want to test the speed on your system/data, here is some example code for randomly subsampling your experiment:

```R
# Define the number of species/samples in your randomly subsampled experiment
N <- 300 # species
M <- 20  # samples

# Randomly subset N species from the original
ex   <- prune_species(sample(species.names(physeq), N), physeq)

# Randomly subset M samples from the original
ex   <- prune_species(sample(sample.names(ex), M), ex)

# Pick a random root from the sub-sampled dataset
tree <- root(tre(ex), sample(species.names(ex), 1), resolve.root=TRUE)
OTU  <- otuTable(ex)

# Rebuild to make sure you didn't screw up something.
ex   <- phyloseq(OTU, tree)

# Calculate the UniFrac dist matrix and print the elapsed system time
system.time(UniFrac(ex, weighted=TRUE))["elapsed"]
```

![UniFrac time trial](https://github.com/downloads/joey711/phyloseq/pretty_time_trial_v01.png)
Computation time required to calculate the distance matrix for a 21-sample experiment at varying numbers of species. For comparison, the performance of the original unweighted UniFrac algorithm as implemented in the "picante" package is also included.
