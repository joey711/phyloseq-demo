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

# # # # # That's interesting. 
# # # # # Now how would we filter the taxa with variance smaller than 0.001?
gpac_filt <- prune_species(specvar > 0.001, gpac)
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family")
#
# # # # # Note that this is not an endorsement of this particular threshold (0.001)
# # # # # Just a demonstration
#
# # # # # For normalization, Susan says to look at the edgeR package.
#
# # # # # Also for standardization, decostand() function in vegan-package

###################################################
###################################################
# Graphical Exploration of data
#
###################################################
###################################################

# # # # # Abundance values graphical inspection ("ground-truthing")
###################################################
# plot_taxa_bar
###################################################
# Subset further to top 5 phyla
top5ph <- sort(tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
GP1    <- subset_species(GP, Phylum %in% names(top5ph))
plot_taxa_bar(GP1, "Phylum", NULL, threshold=0.9, "human", "SampleType") 
plot_taxa_bar(GP1, "Phylum", NULL, threshold=0.9, "human", "SampleType", 
							facet_formula= TaxaGroup ~ .) 


###################################################
# distance function
###################################################
?distance
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
help("distance") # Same as "?distance"


# # # # # # plot_sample_network
ig <- make_sample_network(GP, "bray", 0.9)
(p3  <- plot_sample_network(ig, GP, color="SampleType", shape="human", line_weight=0.4, label=NULL))


# Try enterotype example
data(enterotype)
ig <- make_sample_network(enterotype, max.dist=0.3)
(p  <- plot_sample_network(ig, enterotype, color="SeqTech",
    shape="Enterotype", line_weight=0.4, label=NULL))


###################################################
# ordinate function
###################################################
GP.NMDS <- ordinate(GP, "NMDS", "bray") # perform NMDS on bray-curtis distance
# GP.NMDS.UF.ord   <- ordinate(GP, "NMDS") # UniFrac. Takes a while.
# GP.NMDS.wUF.ord  <- ordinate(GP, "NMDS", "unifrac", weighted=TRUE) # weighted-UniFrac
GP.NMDS.gower <- ordinate(GP, "NMDS", "gower")


###################################################
# plot_ordination function
###################################################
(p1 <- plot_ordination(GP, GP.NMDS, "samples", color="SampleType") +
  geom_line() + geom_point(size=5) )

?plot_ordination # Example "1-liners" are at bottom.


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


# # # # # Describe NeatMap approach, then show examples.
###################################################
# plot_heatmap function
###################################################
plot_heatmap(GP2, "NMDS", "bray", "SampleType", "Family")
plot_heatmap(GP2, "NMDS", "jaccard", "SampleType", "Family")

# plot_heatmap transformations.
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=log_trans(10))
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=identity_trans())
plot_heatmap(gpac_filt, "NMDS", "bray", "SampleType", "Family", trans=boxcox_trans(0.15))


###################################################
# Validation
###################################################
# Do testing on fractional abundance to remove 
# effect of differences in total sequencing across samples for same taxa
# Recall that:
# # # # GP2 <- subset_species(GP, Phylum %in% names(top.TaxaGroup))
# # # # topsp <- names(sort(speciesSums(GP2), TRUE)[1:200])
# # # # GP2   <- prune_species(topsp, GP2)
topsp <- names(sort(speciesSums(GP), TRUE)[1:200])
GP3f  <- transformSampleCounts(GP, function(x){x/sum(x)})
# GP3f <- subset_species(GP3f, Phylum %in% names(top.TaxaGroup))
GP3f  <- prune_species(topsp, GP3f)
head(otuTable(GP))
head(otuTable(GP2))
head(otuTable(GP3))

# # # # # Example: through multiple testing
?mt

# Calculate the multiple-inference-adjusted P-values
GP.fwer.table <- mt(GP3f, "human")
# getTaxa(GP3f, "Family")
jranks <- c("Phylum", "Family", "Genus")
GP.fwer.table <- data.frame(GP.fwer.table, taxTab(GP3f)[rownames(GP.fwer.table), jranks])
subset(GP.fwer.table, adjp < 0.05)


# # # # # What if we want FDR?
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



# # # # # ALTERNATIVES: "qvalue" package, "multcomp" package


# Susan likes the hypergeometric:
?fisher.test

