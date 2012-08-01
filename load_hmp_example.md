<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

## Load the pre-imported HMPv35 data from the RData file

```r
library("phyloseq")
```

```
## Warning: A specification for S3 class "connection" in package
## 'BiocGenerics' seems equivalent to one from package 'RJSONIO' and is not
## turning on duplicate class definitions for this class
```

```
## Warning: A specification for S3 class "file" in package 'BiocGenerics'
## seems equivalent to one from package 'RJSONIO' and is not turning on
## duplicate class definitions for this class
```

```
## Warning: A specification for S3 class "pipe" in package 'BiocGenerics'
## seems equivalent to one from package 'RJSONIO' and is not turning on
## duplicate class definitions for this class
```

```
## Warning: A specification for S3 class "textConnection" in package
## 'BiocGenerics' seems equivalent to one from package 'RJSONIO' and is not
## turning on duplicate class definitions for this class
```


```r
temp_hmp_file <- tempfile()
hmp_url <- "http://cloud.github.com/downloads/joey711/phyloseq/HMPv35.RData"
download.file(hmp_url, temp_hmp_file)
load(temp_hmp_file)
```


### Clean up.

```r
unlink(temp_hmp_file)
rm("hmp_url", "temp_hmp_file")
```


### Investigate the phyloseq-class representation of the data

```r
show(HMPv35)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [45336 species and 4743 samples]
##                      species are rows
## Sample Data:         [4743 samples by 9 sample variables]:
## Taxonomy Table:     [45336 species by 6 taxonomic ranks]:
## Phylogenetic Tree:  [45336 tips and 45099 internal nodes]
##                      unrooted
```


### For example, what if we wanted to root the tree?
These rank names need to be fixed. Must not have been clear during the import. This will be updated soon. The result is the default generic names added when none are provided during the import. In this case, they should be a character vector with elements like "Root", "Phylum", "Class", etc.


```r
rank.names(HMPv35)
```

```
## [1] "ta1" "ta2" "ta3" "ta4" "ta5" "ta6"
```


Here is a way to check for the presence of particular taxonomic elements in the dataset.

```r
all(as(taxTab(HMPv35), "matrix")[1:1000, 1] == "Root")
```

```
## [1] TRUE
```

```r
any("p__Crenarchaeota" %in% as(taxTab(HMPv35), "matrix")[1:1000, 
    2])
```

```
## [1] FALSE
```

```r
any("p__Bacteroidetes" %in% as(taxTab(HMPv35), "matrix")[1:1000, 
    2])
```

```
## [1] TRUE
```


Is the tree rooted?

```r
is.rooted(tre(HMPv35))
```

```
## [1] FALSE
```


This is how you would randomly select an OTU as outgroup to root the tree:

```r
root(tre(HMPv35), sample(species.names(HMPv35), 1), resolve.root = TRUE)
```

```
## 
## Phylogenetic tree with 45336 tips and 45100 internal nodes.
## 
## Tip labels:
## 	OTU_97.15099, OTU_97.13686, OTU_97.30326, OTU_97.26112, OTU_97.34719, OTU_97.12776, ...
## 	Node labels:
## 	NA, 0.697, 0.902, 0.996, 0.768, 0.971, ...
## 
## Rooted; includes branch lengths.
```

