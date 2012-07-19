Demo: phyloseq – A Bioconductor package for handling and analysis of high-throughput phylogenetic sequence data 
========================================================

Paul J. McMurdie and Susan Holmes
Statistics Department, Stanford University,
Stanford, CA 94305, USA
    
E-mail: mcmurdie@stanford.edu
https://github.com/joey711/phyloseq
susan@stat.stanford.edu\
http://www-stat.stanford.edu/~susan/

# Summary

This is a demonstration manual for the phyloseq package. It is an R-markdown
source-file with example code-chunks that can be reused if you have
phyloseq installed.

# phyloseq Documentation

Vignettes are included in phyloseq. A quick way to load them from within `R` is:


```r
vignette("phyloseqbasics")
vignette("phyloseqanalysis")
```


This demonstration document supports a live demonstration of tools in
phyloseq, and supplements other documentation resources
available for [the phyloseq package](https://github.com/joey711/phyloseq) (e.g. wiki, vignettes, publications, function-level documentation, etc.).


# Installation
For the foreseeable near future, phyloseq is under active development.
Users are encouraged to consistently update their version from 
the development website on GitHub. The following code should
install the newest "bleeding edge" version of the phyloseq package
onto your system, including dependencies.
Further instructions are available at [the installation wiki](https://github.com/joey711/phyloseq/wiki/Installation)

```r
source("http://bioconductor.org/biocLite.R")
biocLite("multtest")
biocLite("genefilter")
Make sure you have devtools installed
install.packages("devtools")
Load the devtools package
library("devtools")
Build and install phyloseq
install_github("phyloseq", "joey711")
```


#  Load phyloseq, and import data.
Of course we need to start this tutorial by loading the phyloseq package
(assuming that it is installed.)

```r
library("phyloseq")
```


There is package-level documentation available. Note the following difference:

```r
`?`("phyloseq-package")
`?`(phyloseq  # this is a function)
```

The latter loads instead the documentation for the constructor function named `phyloseq()`


## Basic Data Import


### Importing the Output from [QIIME](http://qiime.org/)

```r
otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", 
    package = "phyloseq")
mapfile <- system.file("extdata", "master_map.txt", package = "phyloseq")
trefile <- system.file("extdata", "GP_tree_rand_short.newick.gz", 
    package = "phyloseq")
qiimex <- import_qiime(otufile, mapfile, trefile, showProgress = FALSE)
print(qiimex)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [500 species and 26 samples]
##                      species are rows
## Sample Data:         [26 samples by 7 sample variables]:
## Taxonomy Table:     [500 species by 7 taxonomic ranks]:
## Phylogenetic Tree:  [500 tips and 499 internal nodes]
##                      rooted
```



### Importing the Output from [mothur](http://www.mothur.org/)
Example: Manually re-import [the "esophagus dataset"](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC384727), which is already inclduded in the phyloseq package.

```r
mothlist <- system.file("extdata", "esophagus.fn.list.gz", package = "phyloseq")
mothgroup <- system.file("extdata", "esophagus.good.groups.gz", package = "phyloseq")
mothtree <- system.file("extdata", "esophagus.tree.gz", package = "phyloseq")
cutoff <- "0.10"
esophman <- import_mothur(mothlist, mothgroup, mothtree, cutoff)
print(esophman)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [58 species and 3 samples]
##                      species are rows
## Phylogenetic Tree:  [58 tips and 57 internal nodes]
##                      rooted
```


Let's test if they are identical as expected

```r
identical(esophagus, esophman)
```

```
## Error: object 'esophagus' not found
```



### Importing [biom-format](http://biom-format.org/) files
The biom-format is intended to be a complete representation of the OTU-clustering results, and so import can be performed with just one connection/file path.


```r
rich_sparse_biom <- system.file("extdata", "rich_sparse_otu_table.biom", 
    package = "phyloseq")
rich_sparse <- import_biom(rich_sparse_biom, taxaPrefix = "greengenes")
print(rich_sparse)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [5 species and 6 samples]
##                      species are rows
## Sample Data:         [6 samples by 4 sample variables]:
## Taxonomy Table:     [5 species by 7 taxonomic ranks]:
```


Note: the current biom-format definition lacks an phylogenetic tree. I am working with the biom-format team on including a phylogenetic tree in the next version of the format, and have contributed a [preliminary R package for biom-format I/O support](https://github.com/biom-format/biom-format/pull/27) as a pull-request to [the biom-format development page on GitHub](https://github.com/biom-format/biom-format). 

The biom-format definition allows for both sparse and dense representations of the abundance data, and is also flexible enough to allow a "minimal" (abundance table onle) and "rich" forms (includes sample and taxonomy data). *All of these forms are supported and automatically recognized/interpreted in phyloseq* through the `import_biom` function.


```r
# Define file path to all four format combinations
rich_dense_biom <- system.file("extdata", "rich_dense_otu_table.biom", 
    package = "phyloseq")
rich_sparse_biom <- system.file("extdata", "rich_sparse_otu_table.biom", 
    package = "phyloseq")
min_dense_biom <- system.file("extdata", "min_dense_otu_table.biom", 
    package = "phyloseq")
min_sparse_biom <- system.file("extdata", "min_sparse_otu_table.biom", 
    package = "phyloseq")
# Each import is only one line, we will import four different example
# 'datasets'
rich_dense <- import_biom(rich_dense_biom, taxaPrefix = "greengenes")
rich_sparse <- import_biom(rich_sparse_biom, taxaPrefix = "greengenes")
min_dense <- import_biom(min_dense_biom, taxaPrefix = "greengenes")
min_sparse <- import_biom(min_sparse_biom, taxaPrefix = "greengenes")
# print summary if phyloseq, the component class otherwise.
biom_ex_print <- function(i) {
    if (class(i) != "phyloseq") {
        class(i)
    } else {
        i
    }
}
sapply(list(rich_dense, rich_sparse, min_dense, min_sparse), biom_ex_print)
```

```
## [[1]]
## phyloseq-class experiment-level object
## OTU Table:          [5 species and 6 samples]
##                      species are rows
## Sample Data:         [6 samples by 4 sample variables]:
## Taxonomy Table:     [5 species by 7 taxonomic ranks]:
## 
## [[2]]
## phyloseq-class experiment-level object
## OTU Table:          [5 species and 6 samples]
##                      species are rows
## Sample Data:         [6 samples by 4 sample variables]:
## Taxonomy Table:     [5 species by 7 taxonomic ranks]:
## 
## [[3]]
## [1] "otuTable"
## 
## [[4]]
## [1] "otuTable"
## 
```


## A More Complicated Import Example
One of the example datasets included in the phyloseq package is derived from [the study first describing human microbiome "Enterotypes"](http://www.nature.com/nature/journal/v473/n7346/full/nature09944.html), and that dataset is called simply `enterotype`. It will be called in later examples using the `data` command.

A more recent study investigating human microbiome "Enterotypes" is titled [Linking Long-Term Dietary Patterns with Gut Microbial Enterotypes](http://www.sciencemag.org/content/334/6052/105.short) by Wu et al., Science, 334 (6052), 105–108. One of the three corresponding authors has the last name "Bushman", which also happens to be the title of the QIIME-processed version of this dataset at [the microbio.me/qiime database](http://www.microbio.me/qiime/).

We will import this data to illustrate a more complicated situation in which we need to import 3 different data components from two different file types (one is [biom-format](http://biom-format.org/), the other is  [sample data contained in a tab-delimited "Mapping File" format produced by QIIME](http://qiime.org/tutorials/tutorial.html).

For convenience and stability, these "Bushman" data files have been saved onto the AMI for the conference, and will be imported from their local system location. An even more complicated "direct import" example is provided in the next subsection, but produces the same result and is not run by the embedded code.


```r
biom_file <- "/R-pkgs/BioC_phyloseq_materials/Bushman/study_1011_closed_reference_otu_table.biom"
map_file <- "/R-pkgs/BioC_phyloseq_materials/Bushman/study_1011_mapping_file.txt"
# Now import the .biom-formatted otuTable-taxonomyTable file.
biom_otu_tax <- import_biom(biom_file, "greengenes")
# Add sample data to the dataset using merge
bmsd <- import_qiime_sampleData(map_file)
class(bmsd)
```

```
## [1] "sampleData"
## attr(,"package")
## [1] "phyloseq"
```

```r
dim(bmsd)
```

```
## [1] 102 225
```


## Direct ftp download, unzip, and import
The `.biom` and sample data files are also [provided online (ftp)](ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_1011_split_library_seqs_and_mapping.zip), and a useful way to download and import into phyloseq directly from the ftp address in the following example code. This is an example in which we download a zip file with both biom- and qiime-formatted data, unzip it in a temporary directory from with in R, import the relavant files using phyloseq importers, and then delete the temporary files. This code *should* be platform independent, but occasionally there are finicky Windows issues that arise.

```r
# This is not actually run
zipftp <- "ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_1011_split_library_seqs_and_mapping.zip"
# First create a temporary directory in which to store the unpacked
# file(s) from the .zip
tmpdir <- tempdir()
# Second create a temp file where you will put the .zip-file itself
temp <- tempfile()
# Now download the file and unzip to tmpdir directory
download.file(zipftp, temp)
unzip(temp, exdir = tmpdir)
# Define the biom file-path
biom_file <- file.path(tmpdir, list.files(tmpdir, pattern = ".biom"))
# Define the mapping file-path
map_file <- file.path(tmpdir, list.files(tmpdir, pattern = "mapping"))
# Now import the .biom-formatted otuTable/taxonomyTable file.
biom_otu_tax <- import_biom(biom_file, "greengenes")
# Add sample data to the dataset using merge
bmsd <- import_qiime_sampleData(map_file)
# Remove the temperorary file and directory where you unpacked the zip
# files
unlink(temp)
unlink(tmpdir)
```


## Merging datasets or components
We need to merge these two separate Bushman dataset objects into one "phyloseq" object.
Presently, the two data objects contain the `otuTable`, `taxonomyTable`, and `sampleData` components, respectively. If we had three objects that were all components (think single tables, or a tree), then we would use the constructor function, `phyloseq`. However, because the `.biom` file contained two tables (including an `otuTable`), the `import_biom` function returned a valid `"phyloseq-class"` instance instead that contained both components. Whenever you need to add or merge data componentes from one (or more) phyloseq-class objects, the merging function, `merge_phyloseq`, is recommended, rather than the constructor (`phyloseq`).

```r
Bushman <- merge_phyloseq(biom_otu_tax, bmsd)
```



#  Basic interaction with phyloseq data

Let's look at some basic print and accessor functions/methods provided in phyloseq.


## Print method

```r
Bushman
```

```
## phyloseq-class experiment-level object
## OTU Table:          [1873 species and 100 samples]
##                      species are rows
## Sample Data:         [100 samples by 225 sample variables]:
## Taxonomy Table:     [1873 species by 7 taxonomic ranks]:
```



## Convenience accessors

```r
nspecies(Bushman)
```

```
## [1] 1873
```

```r
nsamples(Bushman)
```

```
## [1] 100
```

```r
sample.names(Bushman)[1:10]
```

```
##  [1] "C.3075.01.S1.405156" "C.3068.01.S1.405137" "C.3003.01.P1.405142"
##  [4] "C.4007.01.P1.405187" "C.4001.01.P1.405188" "C.3067.01.S1.405183"
##  [7] "C.3043.01.S1.405217" "C.3006.01.P1.405164" "C.3078.01.S1.405181"
## [10] "C.3086.01.S1.405130"
```

```r
species.names(Bushman)[1:10]
```

```
##  [1] "248563" "110059" "223351" "358030" "367581" "16076"  "313844"
##  [8] "49837"  "517282" "296166"
```



## Interacting with the sample variables
This is useful later in plotting

```r
sample.variables(Bushman)[1:10]
```

```
##  [1] "X.SampleID"                  "BarcodeSequence"            
##  [3] "LinkerPrimerSequence"        "ASSIGNED_FROM_GEO"          
##  [5] "ASPARTAME_MG_AVE"            "TOT_CONJUGLINOLEICA_G_AVE"  
##  [7] "AGE"                         "VITC_ASCORBIC_ACID_MG_AVE"  
##  [9] "OXALIC_ACID_MG_AVE"          "PUFA_EICOSAPENTAENOIC_G_AVE"
```

```r
length(sample.variables(Bushman))
```

```
## [1] 225
```

```r
# How can we look at the values for a particular variable
getVariable(Bushman, sample.variables(Bushman)[5])
```

```
##   [1] 267.759   0.000   0.000   0.000   0.000   0.000   0.000   0.000
##   [9]  67.961   1.167 253.080   0.000   0.000  10.263   0.000 472.708
##  [17]  60.739   0.000 178.793 282.437   0.000  94.628   0.000   0.000
##  [25]  62.989   0.000   0.000  16.990   0.000   0.000   0.000   0.000
##  [33] 113.269  60.739   0.000  33.981   0.000   0.000   0.000   0.000
##  [41]   0.000 107.208   0.000   0.000   0.000   0.000   0.000  58.875
##  [49]   0.000   0.000   0.000 257.124   0.000   1.771   0.000   0.000
##  [57]   0.000   0.000   0.000   0.000   0.000   0.000 217.392   0.000
##  [65]   0.000  31.671   0.000 161.971   0.000   0.000   0.000   0.000
##  [73]  33.981   0.000  60.739   0.000   0.000   0.000 176.812  51.820
##  [81]   0.000 110.667 157.867   0.000  34.100   0.000   0.000   0.000
##  [89]   1.167   0.000 830.102   0.000  60.739   0.000   0.000   0.000
##  [97]   0.000   0.000   0.000  70.496
```



## Interacting with the taxonomic ranks
This is useful later in plotting

```r
rank.names(Bushman)
```

```
## [1] "ta1" "ta2" "ta3" "ta4" "ta5" "ta6" "ta7"
```

```r
getTaxa(Bushman, "ta1")
```

```
## [1] "Bacteria"
```


The `rank.names` returned in the previous chunk are a little weird. Let's assign more meaningful taxonomic rank names to the Bushman dataset, since we know (in this case) that the "greengenes" taxonomy tools were used for the assignment we actually then know the names of the seven taxonomic ranks used.

```r
colnames(taxTab(Bushman)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
    o = "Order", f = "Family", g = "Genus", s = "Species")
getTaxa(Bushman, "Kingdom")
```

```
## [1] "Bacteria"
```

```r
getTaxa(Bushman, "Phylum")
```

```
##  [1] "Firmicutes"      "Bacteroidetes"   "Tenericutes"    
##  [4] "Proteobacteria"  "Cyanobacteria"   "Actinobacteria" 
##  [7] "Lentisphaerae"   "Fusobacteria"    "TM7"            
## [10] "Verrucomicrobia" "Synergistetes"  
```



## Abundance Accessors
The purposes of the `sampleSums` and `speciesSums` function are pretty straightforward, but the `getSamples` and `getSpecies` functions can bet a bit confusing.

```r
sampleSums(Bushman)[1:10]
```

```
## C.3075.01.S1.405156 C.3068.01.S1.405137 C.3003.01.P1.405142 
##               13919               10488                6776 
## C.4007.01.P1.405187 C.4001.01.P1.405188 C.3067.01.S1.405183 
##                2098                3911               11690 
## C.3043.01.S1.405217 C.3006.01.P1.405164 C.3078.01.S1.405181 
##               11825                6748               14934 
## C.3086.01.S1.405130 
##               13554 
```

```r
speciesSums(Bushman)[1:10]
```

```
## 248563 110059 223351 358030 367581  16076 313844  49837 517282 296166 
##     53     71    692     90  15559   1219      3     18      4      1 
```

```r
getSpecies(Bushman, sample.names(Bushman)[5])[1:10]
```

```
## 248563 110059 223351 358030 367581  16076 313844  49837 517282 296166 
##      0      0      0      0      0      0      0      0      0      0 
```

```r
getSamples(Bushman, species.names(Bushman)[5])[1:10]
```

```
## C.3075.01.S1.405156 C.3068.01.S1.405137 C.3003.01.P1.405142 
##                 100                 337                 104 
## C.4007.01.P1.405187 C.4001.01.P1.405188 C.3067.01.S1.405183 
##                  64                   0                 374 
## C.3043.01.S1.405217 C.3006.01.P1.405164 C.3078.01.S1.405181 
##                   0                   0                 123 
## C.3086.01.S1.405130 
##                   9 
```


Note how a *sample name* is required by `getSpecies`, and vice versa. This might seem confusing at first, but `getSpecies` is returning all the OTU abundances from *one* sample, while `getSamples` is returning the abundances from all samples for *one* OTU.


# Simple Summary Graphics
Load additional graphics-related packages

```r
library("ggplot2")
library("scales")
library("grid")
```


##  Some examples for plotting richness estimates from un-trimmed data.

```r
plot_richness_estimates(Bushman)  #, 'sample.names', 'SampleType')
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 


```r
(p <- plot_richness_estimates(Bushman, x = "SEX"))
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 


```r
p + geom_boxplot(data = p$data, aes(x = SEX, y = value, color = NULL), 
    alpha = 0.1)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 

Some others you might try (not run in demo)

```r
plot_richness_estimates(Bushman, x = "INSOLUBLE_DIETARY_FIBER_G_AVE")
plot_richness_estimates(Bushman, x = "AGE_IN_YEARS")
plot_richness_estimates(Bushman, x = "VEGETABLE_PROTEIN_G_AVE")
```



## Plotting an Annotated Phylogenetic Tree
Trying to plot too many taxa (tree tips) at once obscures meaning. Let's look at just the Chlamydiae phylum in the incldued `GlobalPatterns` dataset. Note that this also requires subsetting the `GlobalPatterns` dataset using the `subset_species` function, part of the "preprocessing" tools described in the following section.

```r
data(GlobalPatterns)
GlobalPatterns
```

```
## phyloseq-class experiment-level object
## OTU Table:          [19216 species and 26 samples]
##                      species are rows
## Sample Data:         [26 samples by 7 sample variables]:
## Taxonomy Table:     [19216 species by 7 taxonomic ranks]:
## Phylogenetic Tree:  [19216 tips and 19215 internal nodes]
##                      rooted
```

```r
GP.chl <- subset_species(GlobalPatterns, Phylum == "Chlamydiae")
```


Map the sample environment ("SampleType") to point color, and taxonomic Family to point shape. Additionally, label the tips with the Genus name and scale the point size by abundance

```r
plot_tree(GP.chl, color = "SampleType", shape = "Family", label.tips = "Genus", 
    size = "abundance")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 

Not as informative :)

```r
plot_tree(GP.chl, "treeonly")
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 



## Abundance bar plots
For direct quantitative observation/comparison of abundances.

```r
data(enterotype)
TopNOTUs <- names(sort(speciesSums(enterotype), TRUE)[1:10])
ent10 <- prune_species(TopNOTUs, enterotype)
plot_taxa_bar(ent10, "Genus", x = "SeqTech", fill = "TaxaGroup")
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 


```r
plot_taxa_bar(ent10, "Genus", x = "SeqTech", fill = "TaxaGroup") + 
    facet_wrap(~Enterotype)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 


This one takes a little while to calculate (not run).

```r
plot_taxa_bar(Bushman, "Phylum", NULL, 0.9, "SEX", "INSOLUBLE_DIETARY_FIBER_G_AVE")
```



# Preprocessing Abundance Data 
This section includes examples preprocessing (filtering, trimming, subsetting, etc) phyloseq data. Let's start by resetting the `GlobalPatterns` example data and adding a human category.

```r
data(GlobalPatterns)
GlobalPatterns
```

```
## phyloseq-class experiment-level object
## OTU Table:          [19216 species and 26 samples]
##                      species are rows
## Sample Data:         [26 samples by 7 sample variables]:
## Taxonomy Table:     [19216 species by 7 taxonomic ranks]:
## Phylogenetic Tree:  [19216 tips and 19215 internal nodes]
##                      rooted
```

```r
# prune OTUs that are not present in at least one sample
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable:
sampleData(GP)$human <- factor(getVariable(GP, "SampleType") %in% 
    c("Feces", "Mock", "Skin", "Tongue"))
```



## rarefy abundances to even depth
Although of perhaps dubious necessity, it is common for OTU abundance tables
to be randomly subsampled to even sequencing depth prior various analyses, especially UniFrac.
Here is an example comparing the UniFrac-PCoA results with and without
"rarefying" the abundances (Requires phyloseq v1.1.28+)

Test with esophagus dataset

```r
data(esophagus)
eso <- rarefy_even_depth(esophagus)
# plot(as(otuTable(eso), 'vector'), as(otuTable(esophagus), 'vector'))
```


```r
UniFrac(eso)
```

```
##        B      C
## C 0.6078       
## D 0.5807 0.5278
```

```r
UniFrac(esophagus)
```

```
##        B      C
## C 0.5176       
## D 0.5182 0.5422
```


Test with GlobalPatterns dataset

```r
data(GlobalPatterns)
GP.chl <- subset_species(GlobalPatterns, Phylum == "Chlamydiae")
# remove the samples that have less than 20 total reads from Chlamydiae
GP.chl <- prune_samples(names(which(sampleSums(GP.chl) >= 20)), GP.chl)
# (p <- plot_tree(GP.chl, color='SampleType', shape='Family',
# label.tips='Genus', size='abundance'))
GP.chl.r <- rarefy_even_depth(GP.chl)
# plot(as(otuTable(GP.chl.r), 'vector'), as(otuTable(GP.chl), 'vector'))
```


To compare MDS of unweighted-UniFrac for GP.chl and GP.chl.r (default distance is unweighted UniFrac) 

```r
plot_ordination(GP.chl, ordinate(GP.chl, "MDS"), color = "SampleType") + 
    geom_point(size = 5)
```

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-30.png) 


```r
plot_ordination(GP.chl.r, ordinate(GP.chl.r, "MDS"), color = "SampleType") + 
    geom_point(size = 5)
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31.png) 


How does rarefying affect the larger untrimmed dataset? (not run)

```r
GP.r <- rarefy_even_depth(GP)
plot_ordination(GP, ordinate(GP), color = "SampleType") + geom_point(size = 5)
```


```r
plot_ordination(GP.r, ordinate(GP.r), color = "SampleType") + geom_point(size = 5)
```



## `prune_species` vs. `subset_species`
These are two very different methods for subsetting OTUs in a dataset.


```r
topN <- 20
most_abundant_taxa <- sort(speciesSums(GP), TRUE)[1:topN]
print(most_abundant_taxa)
```

```
##  549656  331820  279599  360229  317182   94166  158660  329744  550960 
## 2481214 1001492  927850  556441  528629  444142  443938  436539  384205 
##  189047  326977  317658  244304  171551  263681   98605   12812  536311 
##  312161  291577  234024  222365  215470  212089  203583  196550  196355 
##  192573  298875 
##  179312  177879 
```

```r
GP20 <- prune_species(names(most_abundant_taxa), GP)
# Exploratory tree #1 (Note, too ma) plot_tree(ex2, color='SampleType',
# label.tips='Family', size='abundance')
nspecies(GP20)
```

```
## [1] 20
```

```r
length(getTaxa(GP20, "Phylum"))
```

```
## [1] 5
```


That was `prune_species`, when we know the OTU-IDs of the ones we want, or a logical from a test
Alternatively, you can subset based on taxonomic rank, using `subset_species` and an expression.

```r
GP.chl <- subset_species(GP, Phylum == "Chlamydiae")
# Exploratory tree #2 plot_tree(GP.chl, color='SampleType',
# shape='Family', label.tips='Genus', size='abundance')
nspecies(GP.chl)
```

```
## [1] 21
```

```r
length(getTaxa(GP.chl, "Phylum"))
```

```
## [1] 1
```



## `filterfunSample` and `genefilterSample`
This tool takes after the `genefilter` function from the genefilter package, but emphasizes within-microbiome conditions. The following code illustrates . The function `topp` is a filter function that returns the most abundant `p` fraction of taxa. The `filterfunSample` function takes one or or more functions like `topp` and binds them in order to define a filtering protocol function, in this case called: `f1`. This function, `f1`, is then passed to `genefilterSample` along with the dataset that is going to be pruned as well as a value for `A`, the number of samples in which an OTU must pass the filtering conditions.

```r
topp(0.1)
```

```
## function (x) 
## {
##     if (na.rm) {
##         x = x[!is.na(x)]
##     }
##     x >= sort(x, decreasing = TRUE)[ceiling(length(x) * p)]
## }
## <environment: 0x1088f7298>
```

```r
f1 <- filterfunSample(topp(0.1))
print(f1)
```

```
## function (x) 
## {
##     fun = flist[[1]]
##     fval = fun(x)
##     for (fun in flist[-1]) {
##         fval = fval & fun(x)
##     }
##     return(fval)
## }
## <environment: 0x109e5c788>
## attr(,"class")
## [1] "filterfun"
```

```r
wh1 <- genefilterSample(GP, f1, A = round(0.5 * nsamples(GP)))
sum(wh1)
```

```
## [1] 793
```

```r
ex2 <- prune_species(wh1, GP)
print(GP)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [18988 species and 26 samples]
##                      species are rows
## Sample Data:         [26 samples by 8 sample variables]:
## Taxonomy Table:     [18988 species by 7 taxonomic ranks]:
## Phylogenetic Tree:  [18988 tips and 18987 internal nodes]
##                      rooted
```

```r
print(ex2)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [793 species and 26 samples]
##                      species are rows
## Sample Data:         [26 samples by 8 sample variables]:
## Taxonomy Table:     [793 species by 7 taxonomic ranks]:
## Phylogenetic Tree:  [793 tips and 792 internal nodes]
##                      rooted
```



### Filtering low-variance OTUs
Suppose we wanted to use the variance of OTUs across samples as a condition for filtering. For example, to remove OTUs that do not change much across all (or most) samples. We start this example by subsetting just the *Crenarchaeota* from all samples:

```r
gpac <- subset_species(GP, Phylum == "Crenarchaeota")
```


Now define a function to get the across-sample variance of each OTU/taxa/species.

```r
specvar <- sapply(species.names(gpac), function(i, physeq) {
    var(getSamples(physeq, i))
}, gpac)
p1 <- qplot(x = log10(variance), data = data.frame(variance = specvar), 
    binwidth = abs(do.call("-", as.list(range(log10(specvar))))/20))
print(p1)
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-38.png) 


However, the value of the variance is highly-dependent on the sequencing effort of each sample (the total number of reads sequenced from a particular sample). Thus we segway to transformations (e.g. convert to fractional abundance prior to filtering)


## Transformations
Useful for: Standardization / Normalization / Smoothing / Shrinking. Second-order function: `transformSampleCounts`.

Example normalizing to sample fraction before estimating variance for filtering.

```r
gpacf <- transformSampleCounts(gpac, function(x) {
    x/sum(x)
})
specvar <- sapply(species.names(gpacf), function(i, physeq) {
    var(getSamples(physeq, i))
}, gpacf)
qplot(x = log10(variance), data = data.frame(variance = specvar))
```

```
## stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
## this.
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-39.png) 


```r
p2 <- qplot(x = log10(variance), data = data.frame(variance = specvar), 
    binwidth = abs(do.call("-", as.list(range(log10(specvar))))/20))
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
```

![plot of chunk unnamed-chunk-40](figure/unnamed-chunk-40.png) 


Now how would we filter the taxa with variance smaller than 0.001?

```r
gpac_filt <- prune_species(specvar > 0.001, gpac)
```

Show results with a heatmp

```r
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-42](figure/unnamed-chunk-42.png) 

Note that the code in this section is not an endorsement of this particular threshold (0.001),
just a demonstration for how you might remove OTUs/taxa/species that do not change much across the samples in your experiment.

For normalization, also consider features in the "edgeR" package, and for standardization the `decostand` function in vegan-package


# Graphics for Inference and Exploration
In the following section(s) we will illustrate using graphical tools provided
by the phyloseq package. These are meant to be flexible ways to explore and summarize
the data.

For sake of time of calculation/rendering in this section, let's subset the data further to most abundant 5 phyla.

```r
data(GlobalPatterns)  # Reload GlobalPatterns
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
sampleData(GP)$human <- factor(getVariable(GP, "SampleType") %in% 
    c("Feces", "Mock", "Skin", "Tongue"))
top5ph <- sort(tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum), 
    decreasing = TRUE)[1:5]
GP1 <- subset_species(GP, Phylum %in% names(top5ph))
```


##  Distance Functions
Before we get to network plots, let's discuss distances. Many tools use distances to perform their calculations. In phyloseq, ordinations, heatmap organization, and network plots all use the `distance` function internally for calculating OTU or Sample distance matrices (actually represented as a "dist" object).


```r
help("distance")  # Same as '?distance'
`?`(distance)
```



```r
data(esophagus)
distance(esophagus)  # Unweighted UniFrac
```

```
##        B      C
## C 0.5176       
## D 0.5182 0.5422
```

Here are some other examples. There are some 45 or so methods. 

```r
distance(esophagus, weighted = TRUE)  # weighted UniFrac
distance(esophagus, "jaccard")  # vegdist jaccard
distance(esophagus, "bray")  # vegdist bray-curtis
distance(esophagus, "gower")  # vegdist option 'gower'
distance(esophagus, "g")  # designdist method option 'g'
distance(esophagus, "minkowski")  # invokes a method from the base dist() function.
distance(esophagus, "(A+B-2*J)/(A+B)")  # designdist custom distance
distance("help")
distance("list")
```


## `plot_network`
The following code illustrates use the `make_network` and `plot_network` commands from phyloseq. The `GP` variable is an only-slightly-modified version of the `GlobalPatterns` dataset. The threshold was determined empirically to show something interesting for demonstration. In practice, this value has a huge effect on the resulting network, and its usefulness, and it is highly recommended that you investigate the results from multiple values.

Also it is recommended that you take a look at [the plot_network wiki page](https://github.com/joey711/phyloseq/wiki/plot_network) 


```r
ig <- make_network(GP, "samples", "bray", threshold = 0.9)
(p3 <- plot_network(ig, GP, color = "SampleType", shape = "human", 
    line_weight = 0.4, label = NULL))
```

![plot of chunk plotnetwork-GPsamples](figure/plotnetwork-GPsamples.png) 


A similar network representation of samples from the "Enterotypes" dataset.

```r
data(enterotype)
ig <- make_network(enterotype, max.dist = 0.3)
plot_network(ig, enterotype, color = "SeqTech", shape = "Enterotype", 
    line_weight = 0.4, label = NULL)
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

![plot of chunk plotnetwork-entsamples](figure/plotnetwork-entsamples.png) 


An example showing a network representation of OTUs, representing communities of bacteria that occurr in similar profiles of samples.

```r
data(GlobalPatterns)
# prune to just the top 100 most abundant OTUs across all samples (crude).
GP100 <- prune_species(names(sort(speciesSums(GlobalPatterns), TRUE))[1:100], 
    GlobalPatterns)
jg <- make_network(GP100, "species", "jaccard", 0.3)
plot_network(jg, GP100, "species", color = "Phylum", line_weight = 0.4, 
    label = NULL)
```

![plot of chunk plotnetwork-otus](figure/plotnetwork-otus.png) 



## The `ordinate` function
Using `GP100` from the previous section, let's calculate the unweighted-UniFrac distance for each sample pair in the dataset, and then perform Multidimensional Scaling (aka Principle Coordinates Analysis) on the resulting distance. For details about calculating the UniFrac distance on larger datasets using parallel-computing options in supported by phyloseq, see [the wiki page on Fast Parallel UniFrac in phyloseq](https://github.com/joey711/phyloseq/wiki/Fast-Parallel-UniFrac)

```r
GP.MDS <- ordinate(GP100, "MDS", "unifrac")
```

Here are just a few examples of other supported combinations.

```r
GP.NMDS <- ordinate(GP, "NMDS", "gower")
GP.NMDS <- ordinate(GP, "NMDS", "bray")  # perform NMDS on bray-curtis distance
GP.NMDS.UF.ord <- ordinate(GP, "NMDS")  # UniFrac. Takes a while.
GP.NMDS.wUF.ord <- ordinate(GP, "NMDS", "unifrac", weighted = TRUE)  # weighted-UniFrac
GP.NMDS.gower <- ordinate(GP, "NMDS", "gower")
```



## `plot_ordination` function
The `plot_ordination` function has many options, and supports many combinations of ordinations, and sample or OTU co-variables mapped to color and shape. Many additional examples (with results) are included on [the plot_ordination wiki page](https://github.com/joey711/phyloseq/wiki/plot_ordination). For quicker reference, some example "1-liners" are also included at bottom of this section.

This combination of MDS/PCoA ordination of the UniFrac distance is recently very popular in microbiome analyses. 

```r
plot_ordination(GP100, GP.MDS, "samples", color = "SampleType") + 
    geom_line() + geom_point(size = 5)
```

![plot of chunk unnamed-chunk-49](figure/unnamed-chunk-49.png) 


Get the names of the most-abundant phyla, and use for subsetting.

```r
top.TaxaGroup <- sort(tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], 
    sum, na.rm = TRUE), decreasing = TRUE)
top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1 * 10^6]
# Now prune further, to just the most-abundant phyla
GP2 <- subset_species(GP, Phylum %in% names(top.TaxaGroup))
topsp <- names(sort(speciesSums(GP2), TRUE)[1:200])
GP2 <- prune_species(topsp, GP2)
```


Let's try DCA in one-liner syntax

```r
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "samples", color = "SampleType") + 
    geom_point(size = 4)
```

![plot of chunk unnamed-chunk-51](figure/unnamed-chunk-51.png) 


```r
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "species", color = "Phylum") + 
    geom_point(size = 4)
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "split")
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "split", color = "SampleType")
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "biplot", shape = "Phylum")
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "split", color = "Phylum", 
    label = "SampleType")
plot_ordination(GP2, ordinate(GP2, "DCA"), type = "split", color = "SampleType", 
    shape = "Phylum", label = "SampleType")
```



##   Mapping continuous variables to color (can't do it to shape)

```r
Bushman.ord <- ordinate(Bushman, "CCA")
plot_ordination(Bushman, Bushman.ord, "samples", color = "OMEGA3_FATTY_ACIDS_G_AVE")
```

![plot of chunk unnamed-chunk-53](figure/unnamed-chunk-53.png) 

```r
# plot_ordination(Bushman, Bushman.ord, 'samples', color='VIT_B6_MG_AVE')
```



##  `plot_heatmap` Function
In a [2010 article in BMC Genomics](http://www.biomedcentral.com/1471-2105/11/45), Rajaram and Oono describe an approach to creating a heatmap using ordination methods (namely, NMDS and PCA) to organize the rows and columns instead of (hierarchical) cluster analysis. In many cases the ordination-based ordering does a much better job than h-clustering at providing an order of elements that is easily interpretable. The authors provided an immediately useful example of their approach as [the NeatMap package for R](http://cran.r-project.org/web/packages/NeatMap/index.html). The NeatMap package can be used directly on the abundance table (`"otuTable"`-class) of phylogenetic-sequencing data, but the NMDS or PCA ordination options that it supports are not based on ecological distances. To fill this void, and because phyloseq already provides support for a large number of [ecological distances](https://github.com/joey711/phyloseq/wiki/distance) and [ordination methods](https://github.com/joey711/phyloseq/wiki/ordinate), phyloseq now includes the `plot_heatmap()` function: an ecology-oriented variant of the NeatMap approach to organizing a heatmap and build it using ggplot2 graphics tools. The [distance](https://github.com/joey711/phyloseq/wiki/distance) and [method](https://github.com/joey711/phyloseq/wiki/ordinate) arguments are the same as for the [plot_ordination](https://github.com/joey711/phyloseq/wiki/plot_ordination) function, and support large number of distances and ordination methods, respectively, with a strong leaning toward ecology. This function also provides the options to re-label the OTU and sample axis-ticks with a taxonomic name and/or sample variable, respectively, in the hope that this might hasten your interpretation of the patterns (See the documentation for the `sample.label` and `species.label` arguments, and the examples below). Note that this function makes no attempt to overlay dendrograms from hierarchical clustering next to the axes, as hierarchical clustering is not used to organize these plots. Also note that each re-ordered axis repeats at the edge, and so apparent clusters at the far right/left or top/bottom of the heat-map may actually be the same. For now, the placement of this edge can be considered arbitrary, so beware of this artifact of the graphic and visually check if there are two "mergeable" clusters at the edges of a particular axis. If you benefit from this phyloseq-specific implementation of [the NeatMap approach](http://cran.r-project.org/web/packages/NeatMap/index.html), please cite [the NeatMap article](http://www.biomedcentral.com/1471-2105/11/45), as well as phyloseq.

Further examples are provided at [the plot_heatmap wiki page](https://github.com/joey711/phyloseq/wiki/plot_heatmap)


```r
plot_heatmap(GP2, "MDS", "bray", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-54](figure/unnamed-chunk-54.png) 


Some alternative `plot_heatmap` transformations.

```r
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans = log_trans(10))
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans = identity_trans())
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans = boxcox_trans(0.15))
```



#  Validation
In this section, we will look at examples for using R to validate/test hypotheses
we may have generated through some of the previous exploration.


##  Multiple Testing
In this example we will perform testing on fractional abundances to remove effect of differences in total sequencing across samples for same taxa. For the sake of time, let's subset our testing to the first 200 most abundant OTUs, and transform these counts to their relative abundance in their source sample (note that this introduce bias in real-life scenarios).


```r
topsp <- names(sort(speciesSums(GP), TRUE)[1:200])
GP3f <- transformSampleCounts(GP, function(x) {
    x/sum(x)
})
GP3f <- prune_species(topsp, GP3f)
```


We are going to use the multtest wrapper included in phyloseq, `mt`. Try `?mt` for help on this function. To use this wrapper to calculate the multiple-inference-adjusted P-values, using the "human" sample variable:

```r
GP.fwer.table <- mt(GP3f, "human")
```

```
## B=10000
## b=100	b=200	b=300	b=400	b=500	b=600	b=700	b=800	b=900	b=1000	
## b=1100	b=1200	b=1300	b=1400	b=1500	b=1600	b=1700	b=1800	b=1900	b=2000	
## b=2100	b=2200	b=2300	b=2400	b=2500	b=2600	b=2700	b=2800	b=2900	b=3000	
## b=3100	b=3200	b=3300	b=3400	b=3500	b=3600	b=3700	b=3800	b=3900	b=4000	
## b=4100	b=4200	b=4300	b=4400	b=4500	b=4600	b=4700	b=4800	b=4900	b=5000	
## b=5100	b=5200	b=5300	b=5400	b=5500	b=5600	b=5700	b=5800	b=5900	b=6000	
## b=6100	b=6200	b=6300	b=6400	b=6500	b=6600	b=6700	b=6800	b=6900	b=7000	
## b=7100	b=7200	b=7300	b=7400	b=7500	b=7600	b=7700	b=7800	b=7900	b=8000	
## b=8100	b=8200	b=8300	b=8400	b=8500	b=8600	b=8700	b=8800	b=8900	b=9000	
## b=9100	b=9200	b=9300	b=9400	b=9500	b=9600	b=9700	b=9800	b=9900	b=10000	
## r=2	r=4	r=6	r=8	r=10	r=12	r=14	r=16	r=18	r=20	
## r=22	r=24	r=26	r=28	r=30	r=32	r=34	r=36	r=38	r=40	
## r=42	r=44	r=46	r=48	r=50	r=52	r=54	r=56	r=58	r=60	
## r=62	r=64	r=66	r=68	r=70	r=72	r=74	r=76	r=78	r=80	
## r=82	r=84	r=86	r=88	r=90	r=92	r=94	r=96	r=98	r=100	
## r=102	r=104	r=106	r=108	r=110	r=112	r=114	r=116	r=118	r=120	
## r=122	r=124	r=126	r=128	r=130	r=132	r=134	r=136	r=138	r=140	
## r=142	r=144	r=146	r=148	r=150	r=152	r=154	r=156	r=158	r=160	
## r=162	r=164	r=166	r=168	r=170	r=172	r=174	r=176	r=178	r=180	
## r=182	r=184	r=186	r=188	r=190	r=192	r=194	r=196	r=198	r=200	
```


No add some taxonomic columns to the result for interpretation (rows are OTUs)

```r
jranks <- c("Phylum", "Family", "Genus")
GP.fwer.table <- data.frame(GP.fwer.table, taxTab(GP3f)[rownames(GP.fwer.table), 
    jranks])
subset(GP.fwer.table, adjp < 0.05)
```

```
##        index teststat  rawp   adjp plower        Phylum           Family
## 108747   150    2.574 4e-04 0.0302 0.0243    Firmicutes Streptococcaceae
## 348374    52    1.263 4e-04 0.0302 0.0243 Bacteroidetes   Bacteroidaceae
## 158660     7    2.313 5e-04 0.0356 0.0298 Bacteroidetes   Bacteroidaceae
##                Genus
## 108747 Streptococcus
## 348374   Bacteroides
## 158660   Bacteroides
```


##  What if we want FDR instead of FWER?
Or to use other tools in multtest-package?

```r
library("multtest")
```

```
## Loading required package: Biobase
```

```
## Loading required package: BiocGenerics
```

```
## Attaching package: 'BiocGenerics'
```

```
## The following object(s) are masked from 'package:stats':
## 
## xtabs
```

```
## The following object(s) are masked from 'package:base':
## 
## Filter, Find, Map, Position, Reduce, anyDuplicated, cbind, colnames,
## duplicated, eval, get, intersect, lapply, mapply, mget, order, paste,
## pmax, pmax.int, pmin, pmin.int, rbind, rep.int, rownames, sapply, setdiff,
## table, tapply, union, unique
```

```
## Welcome to Bioconductor
## 
## Vignettes contain introductory material; view with 'browseVignettes()'. To
## cite Bioconductor, see 'citation("Biobase")', and for packages
## 'citation("pkgname")'.
```

```
## Attaching package: 'Biobase'
```

```
## The following object(s) are masked from 'package:phyloseq':
## 
## sampleNames
```

```r
mtm <- mt(GP3f, "human")
```

```
## B=10000
## b=100	b=200	b=300	b=400	b=500	b=600	b=700	b=800	b=900	b=1000	
## b=1100	b=1200	b=1300	b=1400	b=1500	b=1600	b=1700	b=1800	b=1900	b=2000	
## b=2100	b=2200	b=2300	b=2400	b=2500	b=2600	b=2700	b=2800	b=2900	b=3000	
## b=3100	b=3200	b=3300	b=3400	b=3500	b=3600	b=3700	b=3800	b=3900	b=4000	
## b=4100	b=4200	b=4300	b=4400	b=4500	b=4600	b=4700	b=4800	b=4900	b=5000	
## b=5100	b=5200	b=5300	b=5400	b=5500	b=5600	b=5700	b=5800	b=5900	b=6000	
## b=6100	b=6200	b=6300	b=6400	b=6500	b=6600	b=6700	b=6800	b=6900	b=7000	
## b=7100	b=7200	b=7300	b=7400	b=7500	b=7600	b=7700	b=7800	b=7900	b=8000	
## b=8100	b=8200	b=8300	b=8400	b=8500	b=8600	b=8700	b=8800	b=8900	b=9000	
## b=9100	b=9200	b=9300	b=9400	b=9500	b=9600	b=9700	b=9800	b=9900	b=10000	
## r=2	r=4	r=6	r=8	r=10	r=12	r=14	r=16	r=18	r=20	
## r=22	r=24	r=26	r=28	r=30	r=32	r=34	r=36	r=38	r=40	
## r=42	r=44	r=46	r=48	r=50	r=52	r=54	r=56	r=58	r=60	
## r=62	r=64	r=66	r=68	r=70	r=72	r=74	r=76	r=78	r=80	
## r=82	r=84	r=86	r=88	r=90	r=92	r=94	r=96	r=98	r=100	
## r=102	r=104	r=106	r=108	r=110	r=112	r=114	r=116	r=118	r=120	
## r=122	r=124	r=126	r=128	r=130	r=132	r=134	r=136	r=138	r=140	
## r=142	r=144	r=146	r=148	r=150	r=152	r=154	r=156	r=158	r=160	
## r=162	r=164	r=166	r=168	r=170	r=172	r=174	r=176	r=178	r=180	
## r=182	r=184	r=186	r=188	r=190	r=192	r=194	r=196	r=198	r=200	
```


Re-order to original, and use raw p-values for adjustment via `mt.rawp2adjp()`

```r
procedure <- c("Bonferroni", "Hochberg", "BH")
p.mtm <- mt.rawp2adjp(mtm[order(mtm[, "index"]), "rawp"], procedure)
# Re-order so that you can return original table (ordered p-value table)
p.adjp.ord <- p.mtm$adjp[order(p.mtm$index), ]
# Give it the original row names from m
rownames(p.adjp.ord) <- species.names(GP3f)
# Return the table of adjusted p-values for each hypothesis.
GP3f.mt.table <- data.frame(p.adjp.ord, taxTab(GP3f)[rownames(p.adjp.ord), 
    jranks])
# Re-rorder based on BH
GP3f.mt.table <- GP3f.mt.table[order(GP3f.mt.table[, "BH"]), ]
subset(GP3f.mt.table, BH < 0.05)
```

```
##          rawp Bonferroni Hochberg      BH        Phylum             Family
## 158660 0.0005       0.10   0.0990 0.03333 Bacteroidetes     Bacteroidaceae
## 348374 0.0004       0.08   0.0796 0.03333 Bacteroidetes     Bacteroidaceae
## 108747 0.0004       0.08   0.0796 0.03333    Firmicutes   Streptococcaceae
## 194053 0.0008       0.16   0.1576 0.04000    Firmicutes    Lachnospiraceae
## 331820 0.0021       0.42   0.4032 0.04462 Bacteroidetes     Bacteroidaceae
## 322235 0.0017       0.34   0.3315 0.04462 Bacteroidetes     Bacteroidaceae
## 259569 0.0017       0.34   0.3315 0.04462 Bacteroidetes      Rikenellaceae
## 291090 0.0018       0.36   0.3492 0.04462 Bacteroidetes Porphyromonadaceae
## 561077 0.0019       0.38   0.3667 0.04462 Cyanobacteria                   
## 357795 0.0029       0.58   0.5452 0.04462    Firmicutes    Lachnospiraceae
## 518438 0.0026       0.52   0.4966 0.04462    Firmicutes    Lachnospiraceae
## 261912 0.0027       0.54   0.5130 0.04462    Firmicutes    Lachnospiraceae
## 542011 0.0028       0.56   0.5292 0.04462 Cyanobacteria                   
## 317182 0.0037       0.74   0.6882 0.04933 Cyanobacteria                   
## 470172 0.0037       0.74   0.6882 0.04933    Firmicutes    Lachnospiraceae
##                  Genus
## 158660     Bacteroides
## 348374     Bacteroides
## 108747   Streptococcus
## 194053       Roseburia
## 331820     Bacteroides
## 322235     Bacteroides
## 259569       Alistipes
## 291090 Parabacteroides
## 561077                
## 357795       Roseburia
## 518438       Roseburia
## 261912           Dorea
## 542011                
## 317182                
## 470172     Coprococcus
```


Some Alternative packages for multiple inference correction: "qvalue" package, "multcomp" package


#  Getting phyloseq data into other R tools

A big concern for many users is how they can easily get phyloseq-formatted data into other R tools. The following examples are meant to illustrate doing that with some commonly-requested tasks.


##  Porting Data to vegan Functions

For OTU abundance tables, vegan expects samples as rows, and OTUs/species/taxa as columns (so does the picante package). The following is an example function, called `veganotu`, for extracting the OTU table from a phyloseq data object, and converting it to a properly oriented standard matrix recognized by the vegan package.


```r
veganotu <- function(physeq) {
    OTU <- otuTable(physeq)
    if (speciesAreRows(OTU)) {
        OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
}
```


Now use this function for data input to vegan functions

```r
library("vegan")
```


Try the `bioenv` function from vegan to test for sample variables that correlate well
with the microbial community distances.
Note from the `bioenv` documentation:
There are 2^p-1 subsets of p variables, and an exhaustive search may take
a very, very, very long time (parameter upto offers a partial relief).
The formula interface can also help, by specifying variables in the correlation.

Let's try this with the Bushman data, since it has ample continuous variables with which to correlate microbiom distances. 


```r
keepvariables <- which(sapply(sampleData(Bushman), is.numeric))
bushsd <- data.frame(sampleData(Bushman))[keepvariables]
```

There are so many sample variables measured, this could take a long time!
(`1.225996e+55` possible subsets)

```r
bioenv(veganotu(Bushman), bushsd)
```


Instead let's focus on a few that we are interested in comparing,
using the formula interface again. We'll look at the variable
names again to remind us what is available

```r
names(bushsd)[1:10]
```

```
##  [1] "ASPARTAME_MG_AVE"            "TOT_CONJUGLINOLEICA_G_AVE"  
##  [3] "AGE"                         "VITC_ASCORBIC_ACID_MG_AVE"  
##  [5] "OXALIC_ACID_MG_AVE"          "PUFA_EICOSAPENTAENOIC_G_AVE"
##  [7] "PHOSPHORUS_G_AVE"            "SOLUBLE_DIETARY_FIBER_G_AVE"
##  [9] "SFA_MARGARIC_ACID_G_AVE"     "CLA_TRANS10_CIS12_G_AVE"    
```



```r
bioenv(veganotu(Bushman) ~ DEPTH + AGE + TOTAL_FAT_G_AVE + INSOLUBLE_DIETARY_FIBER_G_AVE, 
    bushsd)
```

```
## 
## Call:
## bioenv(formula = veganotu(Bushman) ~ DEPTH + AGE + TOTAL_FAT_G_AVE +      INSOLUBLE_DIETARY_FIBER_G_AVE, data = bushsd) 
## 
## Subset of environmental variables with best correlation to community data.
## 
## Correlations:      spearman 
## Dissimilarities:   bray 
## 
## Best model has 2 parameters (max. 4 allowed):
## AGE INSOLUBLE_DIETARY_FIBER_G_AVE
## with correlation  0.1537 
## 
```


##  Ordination Example on the Gap Statistic
Here is an example performing the gap statistic on ordination results
calculated using phyloseq tools, followed by an example of how a 
ggplot-based wrapper for this example might be included in the phyloseq
package.


```r
library("cluster")
# Load data
data(enterotype)
# ordination
ent.dca <- ordinate(enterotype)
```


Gap Statistic: How many clusters are there?

```r
pam1 <- function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
x <- scores(ent.dca)
gskmn <- clusGap(x[, 1:2], FUN = kmeans, nstart = 20, K.max = 6, 
    B = 500)
gskmn <- clusGap(x[, 1:2], FUN = pam1, K.max = 6, B = 500)
plot(gskmn, main = "Gap statistic for the 'Enterotypes' data")
mtext("k = 2 is best ... but  k = 3  pretty close")
```

![plot of chunk unnamed-chunk-67](figure/unnamed-chunk-67.png) 


That's nice. Just in case it is useful, let's look at what the wrapper-function might look like.

```r
gap_statistic_ordination <- function(ord, FUNcluster, K.max = 6, 
    axes = c(1:2), B = 500, verbose = interactive(), ...) {
    require("cluster")
    # If 'pam1' was chosen, use this internally defined call to pam
    if (FUNcluster == "pam1") {
        FUNcluster <- function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
    }
    
    # Use the scores function to get the ordination coordinates
    x <- scores(ord)
    
    # If axes not explicitly defined (NULL), then use all of them
    if (is.null(axes)) {
        axes <- 1:ncol(x)
    }
    
    # Finally, perform, and return, the gap statistic calculation using
    # cluster::clusGap
    clusGap(x[, axes], FUN = FUNcluster, K.max = K.max, B = B, verbose = verbose, 
        ...)
}
# Define a plot method for results...
plot_clusgap <- function(clusgap, title = "Gap Statistic calculation results") {
    require("ggplot2")
    gstab <- data.frame(gs$Tab, k = 1:nrow(gs$Tab))
    p <- ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
    p <- p + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim))
    p <- p + opts(title = title)
    return(p)
}
```


Now test the result. Should work on ordination classes recognized by `scores` function, and provide a ggplot2 graphic instead of a base graphic.

```r
gs <- gap_statistic_ordination(ent.dca, "pam1", B = 50, verbose = FALSE)
print(gs, method = "Tibs2001SEmax")
```

```
## Clustering Gap statistic ["clusGap"].
## B=50 simulated reference sets, k = 1..6
##  --> Number of clusters (method 'Tibs2001SEmax', SE.factor=1): 2
##       logW E.logW    gap  SE.sim
## [1,] 4.430  4.898 0.4680 0.01834
## [2,] 3.614  4.567 0.9524 0.02240
## [3,] 3.437  4.376 0.9393 0.02669
## [4,] 3.294  4.204 0.9094 0.02740
## [5,] 3.215  4.074 0.8588 0.01730
## [6,] 3.104  3.970 0.8663 0.02258
```

```r
plot_clusgap(gs)
```

![plot of chunk unnamed-chunk-69](figure/unnamed-chunk-691.png) 

```r
# Base graphics means of plotting.
plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("k = 2 is best ... but  k = 3  pretty close")
```

![plot of chunk unnamed-chunk-69](figure/unnamed-chunk-692.png) 
