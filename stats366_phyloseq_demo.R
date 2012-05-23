###################################################
# Load packages
###################################################
library("phyloseq"); library("ggplot2"); library("scales"); library("grid")

###################################################
# Examples for looking at package-level documentation
###################################################
?"phyloseq-package"
?phyloseq # this is a function

###################################################
### code chunk number 9: phyloseq_basics.Rnw:234-236
###################################################
  data(GlobalPatterns)
  GlobalPatterns
  

###################################################
### code chunk number 4: phyloseq_analysis.Rnw:103-109
###################################################
# prune OTUs that are not present in at least one sample
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable:
human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
# Add new human variable to sample data:
sampleData(GP)$human <- factor(human)


# Some examples for plotting richness estimates from un-trimmed data.
###################################################
### code chunk number 5: phyloseq_analysis.Rnw:121-124
###################################################
p <- plot_richness_estimates(GP, "human", "SampleType")
(p <- p + geom_boxplot(data=p$data, aes(x=human, y=value, color=NULL), alpha=0.1))
print(p)
p
###################################################
### code chunk number 6: phyloseq_analysis.Rnw:126-127
###################################################
# ggsave("phyloseq_analysis-richness_estimates.pdf", p, width=11, height=7) # If you wanted to save a pdf of the graphic


###################################################
###################################################
# Filtering example code
###################################################
###################################################

# # # # # prune_species() VERSUS subset_species()
###################################################
### code chunk number 14: phyloseq_basics.Rnw:360-363
###################################################
topN <- 20
most_abundant_taxa <- sort(speciesSums(GP), TRUE)[1:topN]
print(most_abundant_taxa)
ex2 <- prune_species(names(most_abundant_taxa), GP)
# Exploratory tree #1 (Note, too ma)
(p <- plot_tree(ex2, color="SampleType", label.tips="Family", size="abundance"))

# # # # # That was prune_species(), when we know the OTU-IDs of the ones we want, or a logical from a test
# # # # # Alternatively, you can subset based on taxonomic rank, using subset_species() and an expression.
###################################################
### code chunk number 7: phyloseq_analysis.Rnw:144-145
###################################################
GP.chl <- subset_species(GP, Phylum=="Chlamydiae")
# Exploratory tree #2
(p <- plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance"))


# # # # # Use filterfunSample() and genefilterSample()
###################################################
### code chunk number 17: phyloseq_basics.Rnw:391-396
###################################################
f1  <- filterfunSample(topp(0.1))
print(f1)
wh1 <- genefilterSample(GP, f1, A=round(0.5*nsamples(GP)))
sum(wh1)
ex2 <- prune_species(wh1, GP)


# # # # # Variance-based filtering
###################################################
# Subset just the Crenarchaeota from all samples:
###################################################
gpac <- subset_species(GP, Phylum=="Crenarchaeota")
#
#
specvar <- sapply(species.names(gpac), function(i, physeq){
	var(getSamples(physeq, i))
}, gpac)
p1 <- qplot(x=log10(variance), data=data.frame(variance=specvar),
		binwidth=abs(do.call("-", as.list(range(log10(specvar))))/20)
		)
# # # # # But:		
# # # # # The value of the variance is highly-dependent on the sequencing effort
# # # # # of each sample
# # # # # (the total number ofreads sequenced from a particular sample)
# # # # # Segway to transformations (e.g. fractional abundance prior to filtering)
###################################################
###################################################
# Transformations
# Useful for:
# Standardization / Normalization / Smoothing / Shrinking
###################################################
###################################################
# # # # # second-order function:
# # # # # transformSampleCounts()
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

# # # # # That's interesting. Now how would we filter the taxa with variance smaller than 0.001?
gpac_filt <- prune_species(specvar > 0.001, gpac)
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family")
#
# # # # # # For normalization, Susan says to look at the edgeR package.
#
# # # # # Also for standardization, decostand() function in vegan-package

###################################################
###################################################
# Graphical Exploration of data
#
###################################################
###################################################
# heatmap
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=identity_trans())
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=boxcox_trans(0.15))
