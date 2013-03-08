
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

## Load phyloseq package

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.3.20'
```

```r
date()
```

```
## [1] "Fri Mar  8 11:59:20 2013"
```


## Load pre-imported HMPv35 data
The data was previously imported when [the HMP import demo](HMP_import_example.html) was interpreted and built. At the end of that demo is a `save` command that stores the `HMPv35` R data object in a compact R binary form that is easy to re-load back into R. That same data file is posted on [the phyloseq downloads page]("http://cloud.github.com/downloads/joey711/phyloseq/").

```r
temp_hmp_file <- tempfile()
hmp_url <- "HMPv35.RData"
download.file(hmp_url, temp_hmp_file)
```

```
## Error: unsupported URL scheme
```

```r
load(temp_hmp_file)
```

```
## Warning: cannot open compressed file
## '/var/folders/pc/j6k8xlt13kdg_y8755vzgprw0000gn/T//RtmpAAP2FM/file17cce2f5fcbb',
## probable reason 'No such file or directory'
```

```
## Error: cannot open the connection
```


### Clean up.

```r
unlink(temp_hmp_file)
rm("hmp_url", "temp_hmp_file")
```

```
## Warning: object 'hmp_url' not found
```

```
## Warning: object 'temp_hmp_file' not found
```


### Investigate the phyloseq-class representation of the data

```r
show(HMPv35)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for
## function 'show': Error: object 'HMPv35' not found
```


### For example, what if we wanted to root the tree?
These rank names need to be fixed. Must not have been clear during the import. This will be updated soon. The result is the default generic names added when none are provided during the import. In this case, they should be a character vector with elements like "Root", "Phylum", "Class", etc.


```r
rank.names(HMPv35)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for
## function 'tax_table': Error: object 'HMPv35' not found
```


Here is a way to check for the presence of particular taxonomic elements in the dataset.

```r
all(as(tax_table(HMPv35), "matrix")[1:1000, 1] == "Root")
```

```
## Error: error in evaluating the argument 'object' in selecting a method for
## function 'tax_table': Error: object 'HMPv35' not found
```

```r
any("p__Crenarchaeota" %in% as(tax_table(HMPv35), "matrix")[1:1000, 2])
```

```
## Error: error in evaluating the argument 'object' in selecting a method for
## function 'tax_table': Error: object 'HMPv35' not found
```

```r
any("p__Bacteroidetes" %in% as(tax_table(HMPv35), "matrix")[1:1000, 2])
```

```
## Error: error in evaluating the argument 'object' in selecting a method for
## function 'tax_table': Error: object 'HMPv35' not found
```


Is the tree rooted?

```r
is.rooted(phy_tree(HMPv35))
```

```
## Error: error in evaluating the argument 'physeq' in selecting a method for
## function 'phy_tree': Error: object 'HMPv35' not found
```


This is how you would randomly select an OTU as outgroup to root the tree:

```r
root(phy_tree(HMPv35), sample(taxa_names(HMPv35), 1), resolve.root = TRUE)
```

```
## Error: error in evaluating the argument 'physeq' in selecting a method for
## function 'phy_tree': Error: object 'HMPv35' not found
```

