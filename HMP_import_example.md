<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

Example importing into [phyloseq](http://joey711.github.com/phyloseq/) the files produced by [Qiime](http://qiime.org/) run on [the HMPv35 dataset](http://hmpdacc.org/micro_analysis/microbiome_analyses.php), which is avilable from [HMP-DACC](http://hmpdacc.org/).
========================================================

Note that this is a document produced by [the markdown package](http://cran.r-project.org/web/packages/markdown/index.html) from an [R Markdown file](http://rstudio.org/docs/authoring/using_markdown). Markdown is a simple formatting syntax for authoring web pages, and in this case the structure includes reproducible code.

## Load [the phyloseq package](http://joey711.github.com/phyloseq/)

```r
library("phyloseq")
```


## Download the Files to a Temporary Location 
Note that `R`'s file handling is sophisticated enough that you do not need to unpack/unzip the `gz`- or `bz2`-compressed files. They are understood automatically.

### Download the OTU Abundance File.
Create a temporary local file where there large otu-tax file will be downloaded and read-from.


```r
temp_otutax <- tempfile()
otu_url <- "http://downloads.hmpdacc.org/data/HMQCP/otu_table_psn_v35.txt.gz"
download.file(otu_url, temp_otutax)
```

### Download the Sample Map File


```r
temp_map <- tempfile()
map_url <- "http://downloads.hmpdacc.org/data/HMQCP/v35_map_uniquebyPSN.txt.bz2"
download.file(map_url, temp_map)
```


### Download the Temporary Tree.

```r
temp_tree <- tempfile()
tree_url <- "http://downloads.hmpdacc.org/data/HMQCP/rep_set_v35.tre.gz"
download.file(tree_url, temp_tree)
```


It turns out that the tree has weird quotes added around tip names, for example

```r
# Read .gz tree directly (no extra uncompression step needed)
HMPv35.tree <- readTree(temp_tree)
head(species.names(HMPv35.tree))
```

```
## [1] "'OTU_97.15099'" "'OTU_97.13686'" "'OTU_97.30326'" "'OTU_97.26112'"
## [5] "'OTU_97.34719'" "'OTU_97.12776'"
```

These extra quotes ( `"'"` ) must be removed/replaced. This can be accomplished using the following line

```r
HMPv35.tree$tip.label <- substr(species.names(HMPv35.tree), 2, nchar(species.names(HMPv35.tree)) - 
    1)
```


## Import the Data and Create phyloseq Object.
I have included a timing function call so that you can see how long this large file took to import using my laptop. Also, the number of "chunks" is shown by the number of "dots" printed after the function call. Note that the phylogenetic tree was already imported, slightly modified in the previous code chunk, and we want this to provide this `HMPv35.tree` as the tree argument, rather than the raw tree file, which won't exactly match the `species.names` in the OTU-Tax file.


```r
system.time(HMPv35 <- import_qiime(temp_otutax, temp_map, HMPv35.tree, 
    chunk.size = 2000L))
```

```
## Processing map file...
## Processing otu/tax file...
## 
## Reading and parsing file in chunks ... Could take some time. Please be patient...
## 
## Building OTU Table in chunks. Each chunk is one dot.
## .......................Building Taxonomy Table...
## Processing phylogenetic tree file...
```

```
##    user  system elapsed 
##   663.4   311.6  2139.8 
```


## Investigate Some Basic Features of HMPv35
For more advanced tools for exploring the data, see [additional phyloseq documentation](http://joey711.github.com/phyloseq/)


```r
nspecies(HMPv35)
```

```
## [1] 45336
```

```r
nsamples(HMPv35)
```

```
## [1] 4743
```

```r
sample.variables(HMPv35)
```

```
## [1] "X.SampleID"     "RSID"           "visitno"        "sex"           
## [5] "RUNCENTER"      "HMPbodysubsite" "Mislabeled"     "Contaminated"  
## [9] "Description"   
```


## Save HMPv35 to an RData File
Save `HMPv35` to an `".RData"`" file so that you don't have to run this import ever again.


```r
save(HMPv35, file = "~/Dropbox/R/import_HMP_examples/HMPv35.RData")
```


## Alternatively, Download the Processed File from Us
We have already run this import code several times during testing and timing. This [already-imported `HMPv35` from QIIME output dataset](http://cloud.github.com/downloads/joey711/phyloseq/HMPv35.RData) is available from [the phyloseq downloads page](https://github.com/joey711/phyloseq/downloads).  Happily, this file will load in very quickly compared with importing it again from the data files.

