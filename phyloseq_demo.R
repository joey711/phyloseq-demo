## Section: Installation
# For the foreseeable near future, phyloseq is under active development.
# Users are encouraged to consistently update their version from 
# the development website on GitHub. The following code should
# install the newest "bleeding edge" version of the phyloseq package
# onto your system, including dependencies.
# Further instructions are available at [the installation wiki](https://github.com/joey711/phyloseq/wiki/Installation)
source("http://bioconductor.org/biocLite.R")
biocLite("multtest")
biocLite("genefilter")
# Make sure you have devtools installed
install.packages("devtools")
# Load the devtools package
library("devtools")
# Build and install phyloseq
install_github("phyloseq", "joey711")


## Section: Load phyloseq, and import data.
# Of course we need to start this tutorial by loading the phyloseq package
# (assuming that it is installed.)
library("phyloseq")

# For looking at package-level documentation, try 
# `?"phyloseq-package"`
# which is distinct from
# `?phyloseq`
# which is actually a function in the phyloseq-package.

### Subsection: Importing data example: Bushman et al., combined enterotypes dataset
# Here is an example in which we download a zip file with both biom- and qiime-formatted
# data, unzip it in a temporary directory from with in R,
# import the relavant files using phyloseq importers, and then delete the temporary files.
# This code *should* be platform independent, but occasionally there are finicky 
# Windows issues that arise. 

# The "genus and phylum" abundance tables are provided online (ftp) here:
zipftp <- "ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_1011_split_library_seqs_and_mapping.zip"
# First create a temporary directory in which to store the unpacked file(s) from the .zip
tmpdir <- tempdir()
# Second create a temp file where you will put the .zip-file itself
temp <- tempfile()
# Now download the file and unzip to tmpdir directory
download.file(zipftp, temp)
unzip(temp, exdir = tmpdir )

# Define the biom file-path
biom_file <- file.path(tmpdir, list.files(tmpdir, pattern=".biom"))
# Define the mapping file-path
map_file <- file.path(tmpdir, list.files(tmpdir, pattern="mapping"))
# Now import the .biom-formatted otuTable/taxonomyTable file.
biom_otu_tax <- import_biom(biom_file, "greengenes")
# Add sample data to the dataset using merge
bmsd <- import_qiime_sampleData(map_file)
Bushman <- merge_phyloseq(biom_otu_tax, bmsd)

# Remove the temperorary file and directory where you unpacked the zip files
unlink(temp)
unlink(tmpdir, TRUE)


## Section: Basic interaction with phyloseq data
# Let's look at some basic print and accessor functions/methods provided in phyloseq

# Print method
Bushman

# Convenience accessors
nspecies(Bushman)
nspecies(Bushman)
sample.names(Bushman)[1:10]
species.names(Bushman)[1:10]

# Interacting with the sample variables (useful later in plotting)
sample.variables(Bushman)[1:10]
length(sample.variables(Bushman))
# How can we look at the values for a particular variable
getVariable(Bushman, sample.variables(Bushman)[5])

# Interacting with the taxonomic ranks (useful later in plotting)
rank.names(Bushman)
getTaxa(Bushman, "ta1")

# Let's assign more meaningful taxonomic rank names to Bushman
colnames(taxTab(Bushman)) <- c(k="Kingdom", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species")
getTaxa(Bushman, "Kingdom")
getTaxa(Bushman, "Phylum")

# Abundance Sums, accessors
sampleSums(Bushman)[1:10]
speciesSums(Bushman)[1:10]
getSpecies(Bushman, sample.names(Bushman)[5])[1:10]
getSamples(Bushman, species.names(Bushman)[5])[1:10]


## Section: simple summary graphics
# Load additional graphics-related packages
library("ggplot2"); library("scales"); library("grid")

### Subsection: Some examples for plotting richness estimates from un-trimmed data.
plot_richness_estimates(Bushman)#, "sample.names", "SampleType")
(p <- plot_richness_estimates(Bushman, x="SEX"))
p + geom_boxplot(data=p$data, aes(x=SEX, y=value, color=NULL), alpha=0.1)
plot_richness_estimates(Bushman, x="AGE_IN_YEARS")
plot_richness_estimates(Bushman, x="INSOLUBLE_DIETARY_FIBER_G_AVE")
plot_richness_estimates(Bushman, x="VEGETABLE_PROTEIN_G_AVE")

### Subsection: Some examples plotting an annotated tree
# Trying to plot too many taxa (tree tips) at once obscures any meaning.
# Let's look at just the Chlamydiae phylum in the GlobalPatterns dataset
# (This also requires subsetting the GlobalPatterns dataset using subset_species())
data(GlobalPatterns)
GlobalPatterns
GP.chl <- subset_species(GlobalPatterns, Phylum=="Chlamydiae")

# Map the sample environment ("SampleType") to point color,
# and taxonomic Family to point shape. 
# Additionally, label the tips with the Genus name
# and scale the point size by abundance
plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance")
# Not as informative :)
plot_tree(GP.chl, "treeonly")

### Subsection: Abundance bar plots for direct observation/comparison of abundances
data(enterotype)
TopNOTUs <- names(sort(speciesSums(enterotype), TRUE)[1:10])
ent10   <- prune_species(TopNOTUs, enterotype)
plot_taxa_bar(ent10, "Genus", x="SeqTech", fill="TaxaGroup")
plot_taxa_bar(ent10, "Genus", x="SeqTech", fill="TaxaGroup") + facet_wrap(~Enterotype)
# This one takes a little while to calculate...
# plot_taxa_bar(Bushman, "Phylum", NULL, 0.9, "SEX", "INSOLUBLE_DIETARY_FIBER_G_AVE") 




## Section:  Preprocessing Abundance Data 
# Examples preprocessing (filtering, trimming, subsetting, etc) phyloseq data.
# Let's start by Resetting the example data and adding a human category.
data(GlobalPatterns)
GlobalPatterns
# prune OTUs that are not present in at least one sample
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable:
sampleData(GP)$human <- factor(getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))


### Subsection: rarefy abundances to even depth
# Although of perhaps dubious necessity, it is common for OTU abundance tables
# to be randomly subsampled to even sequencing depth prior various analyses, especially UniFrac.
# Here is an example comparing the UniFrac-PCoA results with and without
# "rarefying" the abundances (Requires phyloseq v1.1.28+)

# Test with esophagus dataset
data(esophagus)
eso <- rarefy_even_depth(esophagus)
plot(as(otuTable(eso), "vector"), as(otuTable(esophagus), "vector"))
UniFrac(eso); UniFrac(esophagus)
# Test with GlobalPatterns dataset
data(GlobalPatterns)
GP.chl <- subset_species(GlobalPatterns, Phylum=="Chlamydiae")
# remove the samples that have less than 20 total reads from Chlamydiae
GP.chl <- prune_samples(names(which(sampleSums(GP.chl)>=20)), GP.chl)
# # (p <- plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance"))
GP.chl.r <- rarefy_even_depth(GP.chl)
plot(as(otuTable(GP.chl.r), "vector"), as(otuTable(GP.chl), "vector"))
# Try ordination of GP.chl and GP.chl.r (default distance is unweighted UniFrac)
plot_ordination(GP.chl, ordinate(GP.chl, "MDS"), color="SampleType") + geom_point(size=5)
plot_ordination(GP.chl.r, ordinate(GP.chl.r, "MDS"), color="SampleType") + geom_point(size=5)

# How does rarefying affect larger datasets?
GP.r <- rarefy_even_depth(GP)
plot_ordination(GP, ordinate(GP), color="SampleType") + geom_point(size=5)
plot_ordination(GP.r, ordinate(GP.r), color="SampleType") + geom_point(size=5)


### Subsection: prune_species() VERSUS subset_species()
# Using plot_tree to show the results, we'll try two different 
# methods for subsetting OTUs in a dataset.

topN <- 20
most_abundant_taxa <- sort(speciesSums(GP), TRUE)[1:topN]
print(most_abundant_taxa)
ex2 <- prune_species(names(most_abundant_taxa), GP)
# Exploratory tree #1 (Note, too ma)
plot_tree(ex2, color="SampleType", label.tips="Family", size="abundance")

# That was prune_species(), when we know the OTU-IDs of the ones we want, or a logical from a test
# Alternatively, you can subset based on taxonomic rank, using subset_species() and an expression.
GP.chl <- subset_species(GP, Phylum=="Chlamydiae")
# Exploratory tree #2
plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance")


### Subsection: filterfunSample() and genefilterSample()
f1  <- filterfunSample(topp(0.1))
print(f1)
wh1 <- genefilterSample(GP, f1, A=round(0.5*nsamples(GP)))
sum(wh1)
ex2 <- prune_species(wh1, GP)


### (Optional) Subsection: Variance-based filtering

# Subset just the Crenarchaeota from all samples:

gpac <- subset_species(GP, Phylum=="Crenarchaeota")
#
#
specvar <- sapply(species.names(gpac), function(i, physeq){
	var(getSamples(physeq, i))
}, gpac)
p1 <- qplot(x=log10(variance), data=data.frame(variance=specvar),
		binwidth=abs(do.call("-", as.list(range(log10(specvar))))/20)
		)

# But:		
# The value of the variance is highly-dependent on the sequencing effort
# of each sample
# (the total number of reads sequenced from a particular sample)
# Segway to transformations (e.g. fractional abundance prior to filtering)




## Section: Transformations
# Useful for: Standardization / Normalization / Smoothing / Shrinking
# second-order function: `transformSampleCounts`

# Should normalize to sample fraction before estimating variance for filtering.
gpacf   <- transformSampleCounts(gpac, function(x){x/sum(x)})
specvar <- sapply(species.names(gpacf), function(i, physeq){var(getSamples(physeq, i))}, gpacf)
qplot(x=log10(variance), data=data.frame(variance=specvar))

p2 <- qplot(x=log10(variance), data=data.frame(variance=specvar),
	binwidth=abs(do.call("-", as.list(range(log10(specvar))))/20)
	)


grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))

# # # # # That's interesting. 
# # # # # Now how would we filter the taxa with variance smaller than 0.001?
gpac_filt <- prune_species(specvar > 0.001, gpac)
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family")
# Note that this is not an endorsement of this particular threshold (0.001)
# Just a demonstration

# For normalization, Susan says to look at the edgeR package.

# Also for standardization, decostand() function in vegan-package



## Section: Graphical Exploration of data
# In the following section(s) we will illustrate using graphical tools provided
# by the phyloseq package. These are meant to be flexible ways to explore and summarize
# the data.

### Subsection: Abundance values graphical inspection ("ground-truthing")
# using `plot_taxa_bar`


# For sake of time, let's subset further to most abundant 5 phyla
top5ph <- sort(tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
GP1    <- subset_species(GP, Phylum %in% names(top5ph))
plot_taxa_bar(GP1, "Phylum", NULL, threshold=0.9, "human", "SampleType") 
plot_taxa_bar(GP1, "Phylum", NULL, threshold=0.9, "human", "SampleType", facet_formula= TaxaGroup ~ .) 


### Subsection: distance functions
# For ordinations, network plots

# ```{r, eval=FALSE}
# help("distance") # Same as "?distance"
# ?distance
# ```

data(esophagus)
distance(esophagus) # Unweighted UniFrac
distance(esophagus, weighted=TRUE) # weighted UniFrac
distance(esophagus, "jaccard") # vegdist jaccard
distance(esophagus, "bray")    # vegdist bray-curtis
distance(esophagus, "gower")   # vegdist option "gower"
distance(esophagus, "g") # designdist method option "g"
distance(esophagus, "minkowski") # invokes a method from the base dist() function.
distance(esophagus, "(A+B-2*J)/(A+B)") # designdist custom distance
distance("help")
distance("list")


### Subsection: plot_sample_network
ig <- make_sample_network(GP, "bray", 0.9)
(p3  <- plot_sample_network(ig, GP, color="SampleType", shape="human", line_weight=0.4, label=NULL))

# Try enterotype example for network
data(enterotype)
ig <- make_sample_network(enterotype, max.dist=0.3)
(p  <- plot_sample_network(ig, enterotype, color="SeqTech",
    shape="Enterotype", line_weight=0.4, label=NULL))


### Subsection: ordinate function
GP.NMDS <- ordinate(GP, "NMDS", "bray") # perform NMDS on bray-curtis distance
# GP.NMDS.UF.ord   <- ordinate(GP, "NMDS") # UniFrac. Takes a while.
# GP.NMDS.wUF.ord  <- ordinate(GP, "NMDS", "unifrac", weighted=TRUE) # weighted-UniFrac
GP.NMDS.gower <- ordinate(GP, "NMDS", "gower")



### Subsection: plot_ordination function
plot_ordination(GP, GP.NMDS, "samples", color="SampleType") + geom_line() + geom_point(size=5)
# ?plot_ordination # Example "1-liners" are at bottom.

# Get the names of the most-abundant
top.TaxaGroup <- sort(
	tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum, na.rm = TRUE),
	decreasing = TRUE)
top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
# Now prune further, to just the most-abundant phyla
GP2 <- subset_species(GP, Phylum %in% names(top.TaxaGroup))
topsp <- names(sort(speciesSums(GP2), TRUE)[1:200])
GP2   <- prune_species(topsp, GP2)
GP.dpcoa <- ordinate(GP2, "DPCoA")
plot_ordination(GP2, GP.dpcoa, type="taxa", color="Phylum")
# Customize with ggplot2 layers added directly to output
library("ggplot2")
plot_ordination(GP2, GP.dpcoa, type="samples", color="SampleType") + geom_line() + geom_point(size=5)
p <- plot_ordination(GP2, GP.dpcoa, type="samples", color="SampleType", shape=human)
print(p)
p + geom_line() + geom_point(size=5)
plot_ordination(GP2, GP.dpcoa, type="species", color="Phylum") + geom_line() + geom_point(size=5)
plot_ordination(GP2, GP.dpcoa, type="biplot", shape="Phylum", label="SampleType")
plot_ordination(GP2, GP.dpcoa, type="biplot", shape="Phylum")
plot_ordination(GP2, GP.dpcoa, type="biplot", color="Phylum")
plot_ordination(GP2, GP.dpcoa, type="biplot", label="Phylum")
plot_ordination(GP2, GP.dpcoa, type="split", color="Phylum", label="SampleType")
plot_ordination(GP2, GP.dpcoa, type="split", color="SampleType", shape="Phylum", label="SampleType")

# # # # # dpcoa kinda ugly. Let's try DCA in one-liner syntax
# Try one-liner
plot_ordination(GP2, ordinate(GP2, "DCA"), type="samples", color="SampleType") + geom_point(size=4)
plot_ordination(GP2, ordinate(GP2, "DCA"), type="species", color="Phylum") + geom_point(size=4)
plot_ordination(GP2, ordinate(GP2, "DCA"), type="split")
plot_ordination(GP2, ordinate(GP2, "DCA"), type="split", color="SampleType")
plot_ordination(GP2, ordinate(GP2, "DCA"), type="biplot", shape="Phylum")
plot_ordination(GP2, ordinate(GP2, "DCA"), type="split", color="Phylum", label="SampleType")
plot_ordination(GP2, ordinate(GP2, "DCA"), type="split", color="SampleType",
	shape="Phylum", label="SampleType")


### Subsection:  Mapping continuous variables to color (can't do it to shape)
Bushman.ord <- ordinate(Bushman, "DCA")
plot_ordination(Bushman, Bushman.ord, "samples", color="OMEGA3_FATTY_ACIDS_G_AVE")
plot_ordination(Bushman, Bushman.ord, "samples", color="VIT_B6_MG_AVE")


### Subsection: plot_heatmap function
# Describe NeatMap approach, then show examples.
plot_heatmap(GP2, "NMDS", "bray", "SampleType", "Family")
plot_heatmap(GP2, "NMDS", "jaccard", "SampleType", "Family")

# plot_heatmap transformations.
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=log_trans(10))
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=identity_trans())
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=boxcox_trans(0.15))



## Section: Validation
# In this section, we will look at examples for using R to validate/test hypotheses
# we may have generated through some of the previous exploration.


### Subsection: multiple multiple testing
# Do testing on fractional abundance to remove 
# effect of differences in total sequencing across samples for same taxa.
# First, for the sake of time, let's subset our testing to the first 200 
# most abundant OTUs, and transform these counts to their relative abundance
# in their source sample.
topsp <- names(sort(speciesSums(GP), TRUE)[1:200])
GP3f  <- transformSampleCounts(GP, function(x){x/sum(x)})
GP3f  <- prune_species(topsp, GP3f)
head(otuTable(GP))
head(otuTable(GP2))
head(otuTable(GP3))

# We are going to use the multtest wrapper included in phyloseq, `mt`.
# Try `?mt` for help on this function.
# To use this wrapper to calculate the multiple-inference-adjusted P-values,
# using the "human" sample variable:
GP.fwer.table <- mt(GP3f, "human")
jranks <- c("Phylum", "Family", "Genus")
GP.fwer.table <- data.frame(GP.fwer.table, taxTab(GP3f)[rownames(GP.fwer.table), jranks])
subset(GP.fwer.table, adjp < 0.05)

### Subsection: What if we want FDR instead of FWER?
# Or to use other tools in multtest-package?
library("multtest")
mtm   <- mt(GP3f, "human")
# Re-order to original, and use raw p-values for adjustment via mt.rawp2adjp()
procedure <- c("Bonferroni", "Hochberg", "BH")
p.mtm <- mt.rawp2adjp(mtm[order(mtm[,"index"]), "rawp"], procedure) 
# Re-order so that you can return original table
# (ordered p-value table)
p.adjp.ord <- p.mtm$adjp[order(p.mtm$index), ]
# Give it the original row names from m
rownames(p.adjp.ord) <- species.names(GP3f)
# Return the table of adjusted p-values for each hypothesis.
GP3f.mt.table <- data.frame(p.adjp.ord, taxTab(GP3f)[rownames(p.adjp.ord), jranks])
# Re-rorder based on BH
GP3f.mt.table <- GP3f.mt.table[order(GP3f.mt.table[, "BH"]), ]
subset(GP3f.mt.table, BH < 0.05)
# Some Alternative packages for multiple inference correction:
# "qvalue" package, "multcomp" package


## Section: Getting phyloseq data into other R tools outside of phyloseq
# A big concern for many users is how they can easily get phyloseq-formatted
# data into other R tools. The following examples are meant to illustrate
# doing that with some commonly-requested tasks.

### Subsection: Porting data to vegan functions

# For OTU abundance tables, vegan expects samples as rows, and OTUs/species/taxa as columns
# (so does the picante package).
# The following is an example function, called `veganotu`, for 
# extracting the OTU table from a phyloseq data object, 
# and converting it to a properly oriented standard matrix recognized by 
# the vegan package.
veganotu <- function(physeq){
	OTU <- otuTable(physeq)
	if( speciesAreRows(OTU) ){ OTU <- t(OTU) }
	return( as(OTU, "matrix") )
}

# Now use this function for data input to vegan functions
library("vegan")

# Try the `bioenv` function from vegan to test for sample variables that correlate well
# with the microbial community distances.
# Note from the `bioenv` documentation:
# There are 2^p-1 subsets of p variables, and an exhaustive search may take
# a very, very, very long time (parameter upto offers a partial relief).
# The formula interface can also help, by specifying variables in the correlation.

# Add dummy continuous variable, `contsim`, to the `GP` dataset.

sampleData(GP)$contsim <- sample(1:100, nsamples(GP), TRUE)
sampleData(GP)$contsim1 <- sample(100:1, nsamples(GP), TRUE)
bioenv(veganotu(GP) ~ contsim + contsim1, data.frame(sampleData(GP)))

# Let's try this with the Bushman data. 
# There are so many sample variables measured, this could take a long time!
# (`1.225996e+55` possible subsets)

keepvariables <- which(sapply(sampleData(Bushman), is.numeric))
bushsd <- data.frame(sampleData(Bushman))[keepvariables]
bioenv(veganotu(Bushman), bushsd)

# Instead let's focus on a few that we are interested in comparing,
# using the formula interface again. We'll look at the variable
# names again to remind us what is available
names(bushsd)
bioenv(veganotu(Bushman) ~ DEPTH + AGE + TOTAL_FAT_G_AVE + INSOLUBLE_DIETARY_FIBER_G_AVE, bushsd)


### Subsection: Ordination Example on the Gap Statistic
# Here is an example performing the gap statistic on ordination results
# calculated using phyloseq tools, followed by an example of how a 
# ggplot-based wrapper for this example might be included in the phyloseq
# package.

library("cluster")
# Load data
data(enterotype)
# ordination
ent.dca <- ordinate(enterotype)
# Plot
plot_ordination(enterotype, ent.dca, color="Enterotype")
plot_ordination(enterotype, ent.dca, color="SeqTech")

# Gap Statistic: How many clusters are there?
pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
x <- scores(ent.dca)
# gskmn <- clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn <- clusGap(x[, 1:2], FUN=pam1, K.max = 6, B = 500)
plot(gskmn, main = "Gap statistic for the 'Enterotypes' data")
mtext("k = 2 is best ... but  k = 3  pretty close")

# What would the wrapper-function look like?
gap_statistic_ordination <- function(ord, FUNcluster, K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
	require("cluster")
	# If "pam1" was chosen, use this internally defined call to pam
	if(FUNcluster == "pam1"){
		FUNcluster <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))		
	}
	
	# Use the scores function to get the ordination coordinates
	x <- scores(ord)
	
	# If axes not explicitly defined (NULL), then use all of them
	if(is.null(axes)){axes <- 1:ncol(x)}
	
	# Finally, perform, and return, the gap statistic calculation using cluster::clusGap	
	clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}
# Define a plot method for results...
plot_clusgap <- function(clusgap, title="Gap Statistic calculation results"){
	require("ggplot2")
	gstab <- data.frame(gs$Tab, k=1:nrow(gs$Tab))
	p <- ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
	p <- p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
	p <- p + opts(title=title)
	return(p)
}

# Now test the result. Should work on ordination classes recognized by `scores` function
gs <- gap_statistic_ordination(ent.dca, "pam1", B=50, verbose=FALSE)
print(gs, method="Tibs2001SEmax")
plot_clusgap(gs)
# Base graphics means of plotting.
plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("k = 2 is best ... but  k = 3  pretty close")	


