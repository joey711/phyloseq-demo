## Load the pre-imported HMPv35 data from an RData file
library("phyloseq")
temp_hmp_file <- tempfile()
hmp_url <- "http://cloud.github.com/downloads/joey711/phyloseq/HMPv35.RData"
download.file(hmp_url, temp_hmp_file)
load(temp_hmp_file)
### Clean up.
unlink(temp_hmp_file)
rm("hmp_url", "temp_hmp_file")
### Investigate the phyloseq-class representation of the data
show(HMPv35)
### For example, what if we wanted to root the tree?
rank.names(HMPv35)
all(as(taxTab(HMPv35), "matrix")[1:1000, 1] == "Root")
any( "p__Crenarchaeota" %in% as(taxTab(HMPv35), "matrix")[1:1000, 2] )
any( "p__Bacteroidetes" %in% as(taxTab(HMPv35), "matrix")[1:1000, 2] )
is.rooted(tre(HMPv35))
### This is how you would randomly select an OTU as outgroup to root the tree:
root(tre(HMPv35), sample(species.names(HMPv35), 1), resolve.root=TRUE)
