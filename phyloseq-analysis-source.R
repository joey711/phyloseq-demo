### R code from vignette source 'phyloseq_analysis.Rnw'

###################################################
### code chunk number 1: phyloseq_analysis.Rnw:76-77 (eval = FALSE)
###################################################
## vignette("phyloseq_basics")


###################################################
### code chunk number 2: phyloseq_analysis.Rnw:84-86
###################################################
library("phyloseq")
library("ggplot2")


###################################################
### code chunk number 3: phyloseq_analysis.Rnw:97-98
###################################################
data(GlobalPatterns)


###################################################
### code chunk number 4: phyloseq_analysis.Rnw:103-109
###################################################
# prune OTUs that are not present in at least one sample
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable:
human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
# Add new human variable to sample data:
sampleData(GP)$human <- factor(human)


###################################################
### code chunk number 5: phyloseq_analysis.Rnw:121-124
###################################################
p <- plot_richness_estimates(GP, "human", "SampleType")
(p <- p + geom_boxplot(data=p$data, 
          aes(x=human, y=value, color=NULL), alpha=0.1))


###################################################
### code chunk number 6: phyloseq_analysis.Rnw:126-127
###################################################
ggsave("phyloseq_analysis-richness_estimates.pdf", p, width=11, height=7)


###################################################
### code chunk number 7: phyloseq_analysis.Rnw:144-145
###################################################
GP.chl <- subset_species(GP, Phylum=="Chlamydiae")


###################################################
### code chunk number 8: phyloseq_analysis.Rnw:149-151
###################################################
(p <- plot_tree(GP.chl, color="SampleType", shape="Family",
					label.tips="Genus", size="abundance"))


###################################################
### code chunk number 9: phyloseq_analysis.Rnw:153-154
###################################################
ggsave("phyloseq_analysis-plot_tree_example.pdf", p, width=7, height=7)


###################################################
### code chunk number 10: phyloseq_analysis.Rnw:170-171
###################################################
data(enterotype)


###################################################
### code chunk number 11: EntAbundPlot
###################################################
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 30
barplot(sort(speciesSums(enterotype), TRUE)[1:N]/nsamples(enterotype), las=2)


###################################################
### code chunk number 12: phyloseq_analysis.Rnw:190-191
###################################################
rank.names(enterotype)


###################################################
### code chunk number 13: phyloseq_analysis.Rnw:197-200
###################################################
TopNOTUs <- names(sort(speciesSums(enterotype), TRUE)[1:10]) 
ent10   <- prune_species(TopNOTUs, enterotype)
print(ent10)


###################################################
### code chunk number 14: phyloseq_analysis.Rnw:205-206
###################################################
sample.variables(ent10)


###################################################
### code chunk number 15: entbarplot0
###################################################
(p <- plot_taxa_bar(ent10, "Genus", x="SeqTech", fill="TaxaGroup") +
			facet_wrap(~Enterotype) )


###################################################
### code chunk number 16: phyloseq_analysis.Rnw:215-216
###################################################
ggsave("phyloseq_analysis-entbarplot.pdf", p, width=8, height=6)


###################################################
### code chunk number 17: GPheatmap
###################################################
data("GlobalPatterns")
gpac <- subset_species(GlobalPatterns, Phylum=="Crenarchaeota")
(p <- plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family"))


###################################################
### code chunk number 18: phyloseq_analysis.Rnw:246-247
###################################################
ggsave("phyloseq_analysis-plot_heatmap01.pdf", p, width=7, height=10)


###################################################
### code chunk number 19: ent-network
###################################################
data(enterotype)
ig <- make_sample_network(enterotype, max.dist=0.3)
(p  <- plot_sample_network(ig, enterotype, color="SeqTech",
    shape="Enterotype", line_weight=0.4, label=NULL))


###################################################
### code chunk number 20: phyloseq_analysis.Rnw:274-275
###################################################
ggsave("phyloseq_analysis-plot_sample_network.pdf", p, width=7, height=7)


###################################################
### code chunk number 21: phyloseq_analysis.Rnw:295-298 (eval = FALSE)
###################################################
## my.physeq <- import("Biom", BIOMfilename="myBiomFile.biom")
## my.ord    <- ordinate(my.physeq)
## plot_ordination(my.physeq, my.ord, color="myFavoriteVarible")


###################################################
### code chunk number 22: phyloseq_analysis.Rnw:302-306 (eval = FALSE)
###################################################
## ?import
## ?ordinate
## ?distance
## ?plot_ordination


###################################################
### code chunk number 23: phyloseq_analysis.Rnw:316-317
###################################################
data(GlobalPatterns)


###################################################
### code chunk number 24: phyloseq_analysis.Rnw:319-320 (eval = FALSE)
###################################################
## GPUF <- UniFrac(GlobalPatterns)


###################################################
### code chunk number 25: phyloseq_analysis.Rnw:322-324
###################################################
# Loads the pre-computed distance matrix, GPUF
load("Unweighted_UniFrac.RData")


###################################################
### code chunk number 26: phyloseq_analysis.Rnw:326-327
###################################################
GloPa.pcoa <- pcoa(GPUF)


###################################################
### code chunk number 27: PCoAScree
###################################################
barplot(GloPa.pcoa$values$Relative_eig)


###################################################
### code chunk number 28: GPfig5ax1213
###################################################
(p12 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", color="SampleType") +
  geom_line() + geom_point(size=5) + scale_colour_hue(legend = FALSE) )
(p13 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", axes=c(1, 3),
  color="SampleType") + geom_line() + geom_point(size=5) )


###################################################
### code chunk number 29: phyloseq_analysis.Rnw:351-353
###################################################
ggsave("phyloseq_analysis-GPfig5ax12.pdf", p12, width=7, height=7)
ggsave("phyloseq_analysis-GPfig5ax13.pdf", p13, width=9, height=7)


###################################################
### code chunk number 30: GP_UF_NMDS0
###################################################
# (Re)load UniFrac distance matrix and GlobalPatterns data
data(GlobalPatterns)
load("Unweighted_UniFrac.RData") # reloads GPUF variable
GP.NMDS <- metaMDS(GPUF, k=2) # perform NMDS, set to 2 axes
(p <- plot_ordination(GlobalPatterns, GP.NMDS, "samples", color="SampleType") +
  geom_line() + geom_point(size=5) )


###################################################
### code chunk number 31: GP_UF_NMDS1
###################################################
ggsave("phyloseq_analysis-GP_UF_NMDS.pdf", p, width=9, height=7)


###################################################
### code chunk number 32: GPCAscree0
###################################################
data(GlobalPatterns)
# Take a subset of the GP dataset, top 200 species
topsp <- names(sort(speciesSums(GlobalPatterns), TRUE)[1:200])
GP    <- prune_species(topsp, GlobalPatterns)
# Subset further to top 5 phyla, among the top 200 OTUs.
top5ph <- sort(tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
GP     <- subset_species(GP, Phylum %in% names(top5ph))
# Re-add human variable to sample data:
sampleData(GP)$human <- factor(human)


###################################################
### code chunk number 33: GPCAscree
###################################################
# Now perform a unconstrained correspondence analysis
gpca  <- ordinate(GP, "CCA")
barplot(gpca$CA$eig/sum(gpca$CA$eig), las=2)


###################################################
### code chunk number 34: GPCA1234
###################################################
(p12 <- plot_ordination(GP, gpca, "samples", color="SampleType") + 
  geom_line() + geom_point(size=5) )
(p34 <- plot_ordination(GP, gpca, "samples", axes=c(3, 4), color="SampleType") + 
  geom_line() + geom_point(size=5) )


###################################################
### code chunk number 35: GPCA1234s
###################################################
ggsave("phyloseq_analysis-GPCA12.pdf", p12, width=9, height=7)
ggsave("phyloseq_analysis-GPCA34.pdf", p34, width=9, height=7)


###################################################
### code chunk number 36: GPCAspecplot0
###################################################
p1  <- plot_ordination(GP, gpca, "species", color="Phylum")
(p1 <- ggplot(p1$data, p1$mapping) + geom_point(size=5, alpha=0.5) + 
  facet_wrap(~Phylum) + scale_colour_hue(legend = FALSE) )


###################################################
### code chunk number 37: GPCAspecplot1
###################################################
# Save as raster to control file size.
ggsave("phyloseq_analysis-GPCAspecplot.pdf", p1, width=10, height=7)


###################################################
### code chunk number 38: GPCAspecplotTopo0
###################################################
(p3 <- ggplot(p1$data, p1$mapping) + geom_density2d() +
	facet_wrap(~Phylum) + scale_colour_hue(legend = FALSE) )


###################################################
### code chunk number 39: GPCAspecplotTopo1
###################################################
# Do this one PDF, because it shouldn't be large file.
ggsave("phyloseq_analysis-GPCAspecplotTopo.pdf", p3, width=10, height=7)


###################################################
### code chunk number 40: GPCAjitter0
###################################################
library("reshape")
# Melt the species-data.frame, DF, to facet each CA axis separately
mdf <- melt(p1$data[, c("CA1", "CA2", "Phylum", "Family", "Genus")], 
            id=c("Phylum", "Family", "Genus") )
# Select some special outliers for labelling
LF <- subset(mdf, variable=="CA2" & value < -1.0)
# build plot: boxplot summaries of each CA-axis, with labels
p <- ggplot(mdf, aes(Phylum, value, color=Phylum)) + geom_boxplot() + 
  facet_wrap(~variable, 2) + scale_colour_hue(legend = FALSE) +
  theme_bw() + opts( axis.text.x = theme_text(angle = -90, hjust = 0) )
# Add the text label layer, and render ggplot graphic
(p <- p + geom_text(aes(Phylum, value+0.1, color=Phylum, label=Family), 
                    data=LF, vjust=0, size=2) )


###################################################
### code chunk number 41: GPCAjitter
###################################################
# Save as raster to control file size.
#ggsave("phyloseq_analysis-GPCAjitter.png", p, width=6, height=6, dpi=75)
ggsave("phyloseq_analysis-GPCAjitter.pdf", p, width=7, height=7)


###################################################
### code chunk number 42: GPtaxaplot0
###################################################
(p <- plot_taxa_bar(GP, "Phylum", NULL, threshold=0.9, "human", "SampleType", 
							facet_formula= TaxaGroup ~ .) )


###################################################
### code chunk number 43: GPtaxaplot
###################################################
ggsave("phyloseq_analysis-GPtaxaplot.pdf", p, width=8, height=10)


###################################################
### code chunk number 44: GPdpcoa01
###################################################
GP.dpcoa <- DPCoA(GP)
# GP.dpcoa <- ordinate(GP, "DPCoA") # Alternative; ordinate() function
pdpcoa <- plot_ordination(GP, GP.dpcoa, type="biplot",
     color="SampleType", shape="Phylum")
shape.fac <- pdpcoa$data[, deparse(pdpcoa$mapping$shape)]
man.shapes <- c(19, 21:25)
names(man.shapes) <- c("samples", levels(shape.fac)[levels(shape.fac)!="samples"])
p2dpcoa <- pdpcoa + scale_shape_manual(values=man.shapes)


###################################################
### code chunk number 45: GPdpcoa02
###################################################
ggsave("phyloseq_analysis-GPdpcoaBiplot.pdf", p2dpcoa, width=9, height=7)


###################################################
### code chunk number 46: distancefun
###################################################
data(esophagus)
distance(esophagus) # Unweighted UniFrac
distance(esophagus, weighted=TRUE) # weighted UniFrac
distance(esophagus, "jaccard") # vegdist jaccard
distance(esophagus, "g") # betadiver method option "g"
distance("help")


###################################################
### code chunk number 47: phyloseq_analysis.Rnw:595-600 (eval = FALSE)
###################################################
## data(esophagus)
## UniFrac(esophagus, weighted=TRUE)
## # distance(esophagus, weighted=TRUE) # Alternative using the distance() function
## UniFrac(esophagus, weighted=FALSE)
## # distance(esophagus) # Alternative using the distance() function


###################################################
### code chunk number 48: phyloseq_analysis.Rnw:604-606
###################################################
round( UniFrac(esophagus, weighted=TRUE), 3)
round( UniFrac(esophagus, weighted=FALSE), 3)


###################################################
### code chunk number 49: phyloseq_analysis.Rnw:617-625
###################################################
# (Re)load UniFrac distance matrix and GlobalPatterns data
data(GlobalPatterns)
load("Unweighted_UniFrac.RData") # reloads GPUF variable
# Manually define color-shading vector based on sample type.
colorScale    <- rainbow(length(levels(getVariable(GlobalPatterns, "SampleType"))))
cols          <- colorScale[getVariable(GlobalPatterns, "SampleType")] 
GP.tip.labels <- as(getVariable(GlobalPatterns, "SampleType"), "character")
GP.hclust     <- hclust(GPUF, method="average")


###################################################
### code chunk number 50: GPfig4
###################################################
plot(as.phylo(GP.hclust), show.tip.label=TRUE, tip.color="white")
tiplabels(GP.tip.labels, col=cols, frame="none", adj=-0.05, cex=0.7)


###################################################
### code chunk number 51: GPfig4jaccCLC
###################################################
jaccCLC <- hclust(distance(GlobalPatterns, "jaccard"))
plot( as.phylo(jaccCLC), show.tip.label=TRUE, tip.color="white" )
tiplabels(GP.tip.labels, col=cols, frame="none", adj=-0.05, cex=0.7)


###################################################
### code chunk number 52: phyloseq_analysis.Rnw:670-677 (eval = FALSE)
###################################################
## data(enterotype)
## # Filter samples that don't have Enterotype classification.
## x <- subset_samples(enterotype, !is.na(Enterotype))
## 
## # Calculate the multiple-inference-adjusted P-values
## ent.p.table <- mt(x, "Enterotype", test="f")
## print(head(ent.p.table, 10))

