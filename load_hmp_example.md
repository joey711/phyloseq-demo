
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

## Load phyloseq package

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.3.21'
```

```r
date()
```

```
## [1] "Fri Mar  8 22:54:31 2013"
```


## Load pre-imported HMPv35 data
The data was previously imported when [the HMP import demo](HMP_import_example.html) was interpreted and built. At the end of that demo is a `save` command that stores the `HMPv35` R data object in a compact R binary form that is easy to re-load back into R. That same data file is posted on [the phyloseq downloads page]("http://cloud.github.com/downloads/joey711/phyloseq/").

If you've already downloaded [the HMPv35.RData binary file](https://raw.github.com/joey711/phyloseq-demo/gh-pages/HMPv35.RData) to the working directory of your current R session, then you can easily `load` the `HMPv35` data into your workspace roughly 100X faster than to re-import the data from the raw files:

```r
system.time(load("HMPv35.RData"))
```

```
##    user  system elapsed 
##  11.295   0.495  11.795
```



### Investigate the phyloseq-class representation of the data

```r
HMPv35
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45336 taxa and 4743 samples ]
## sample_data() Sample Data:       [ 4743 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 45336 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 45336 tips and 45099 internal nodes ]
## refseq()      DNAStringSet:      [ 45336 reference sequences ]
```


```r
rank_names(HMPv35)
```

```
## [1] "Rank1"  "Phylum" "Class"  "Order"  "Family" "Genus"
```


Here is a way to check for the presence of particular taxonomic elements in the dataset.

```r
all(tax_table(HMPv35)[, "Rank1"] == "Root")
```

```
## [1] TRUE
```

```r
any("Crenarchaeota" %in% tax_table(HMPv35)[, "Phylum"])
```

```
## [1] FALSE
```

```r
any("Bacteroidetes" %in% tax_table(HMPv35)[, "Phylum"])
```

```
## [1] TRUE
```

```r
any("Firmicutes" %in% tax_table(HMPv35)[, "Phylum"])
```

```
## [1] TRUE
```



### Root the tree
Is the tree rooted?

```r
is.rooted(phy_tree(HMPv35))
```

```
## [1] FALSE
```


This is how you could root (verb) the tree in `HMPv35` by randomly selecting an OTU as outgroup to root the tree

```r
phy_tree(HMPv35) = root(phy_tree(HMPv35), sample(taxa_names(HMPv35), 1), resolve.root = TRUE)
phy_tree(HMPv35)
```

```
## 
## Phylogenetic tree with 45336 tips and 45100 internal nodes.
## 
## Tip labels:
## 	OTU_97.15099, OTU_97.13686, OTU_97.30326, OTU_97.26112, OTU_97.34719, OTU_97.12776, ...
## Node labels:
## 	NA, 0.941, 0.914, 0.918, 0.769, 0.954, ...
## 
## Rooted; includes branch lengths.
```



### Richness

```r
sample_variables(HMPv35)
```

```
## [1] "X.SampleID"     "RSID"           "visitno"        "sex"           
## [5] "RUNCENTER"      "HMPbodysubsite" "Mislabeled"     "Contaminated"  
## [9] "Description"
```

```r
levels(sample_data(HMPv35)$HMPbodysubsite)
```

```
##  [1] "Anterior_nares"               "Attached_Keratinized_gingiva"
##  [3] "Buccal_mucosa"                "Hard_palate"                 
##  [5] "Left_Antecubital_fossa"       "Left_Retroauricular_crease"  
##  [7] "Mid_vagina"                   "Palatine_Tonsils"            
##  [9] "Posterior_fornix"             "Right_Antecubital_fossa"     
## [11] "Right_Retroauricular_crease"  "Saliva"                      
## [13] "Stool"                        "Subgingival_plaque"          
## [15] "Supragingival_plaque"         "Throat"                      
## [17] "Tongue_dorsum"                "Vaginal_introitus"
```

```r
HMPv35saliva = subset_samples(HMPv35, HMPbodysubsite == "Saliva")
HMPv35saliva
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45336 taxa and 290 samples ]
## sample_data() Sample Data:       [ 290 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 45336 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 45336 tips and 45100 internal nodes ]
## refseq()      DNAStringSet:      [ 45336 reference sequences ]
```

```r
# plot_richness(HMPv35, x='HMPbodysubsite', color='sex')
```


