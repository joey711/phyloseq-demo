### R code from vignette source 'phyloseq_basics.Rnw'

###################################################
### code chunk number 1: phyloseq_basics.Rnw:72-73 (eval = FALSE)
###################################################
## vignette("phyloseq_analysis")


###################################################
### code chunk number 2: phyloseq_basics.Rnw:86-87
###################################################
  library("phyloseq")


###################################################
### code chunk number 3: phyloseq_basics.Rnw:117-121 (eval = FALSE)
###################################################
##   otufilename <- "../data/ex1_otutable.txt"
##   mapfilename <- "../data/ex1_sampleData.txt"
##   trefilename <- "../data/ex1_tree.tre"
##   MyExpmt1 <- import_qiime(otufilename, mapfilename, trefilename)


###################################################
### code chunk number 4: phyloseq_basics.Rnw:148-154 (eval = FALSE)
###################################################
## mothur_list_file  <- "~/Downloads/mothur/Esophagus/esophagus.an.list"
## mothur_group_file <- "~/Downloads/mothur/Esophagus/esophagus.good.groups"
## mothur_tree_file  <- "~/Downloads/mothur/Esophagus/esophagus.tree"
## # # Actual examples follow:
## show_mothur_list_cutoffs(mothur_list_file)
## MyExpmt1 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file)


###################################################
### code chunk number 5: phyloseq_basics.Rnw:175-178 (eval = FALSE)
###################################################
##   pyrotagger_tab_file <- "path/to/my/filename.txt"
##   myData1 <- import_pyrotagger_tab(pyrotagger_tab_file,
##     strict_taxonomy=FALSE, keep_potential_chimeras=FALSE)


###################################################
### code chunk number 6: phyloseq_basics.Rnw:182-183 (eval = FALSE)
###################################################
##   myData1 <- import_pyrotagger_tab(pyrotagger_tab_file)  


###################################################
### code chunk number 7: phyloseq_basics.Rnw:194-195 (eval = FALSE)
###################################################
##   myOTU1 <- import_RDP_cluster("path/to/my/filename.clust")


###################################################
### code chunk number 8: phyloseq_basics.Rnw:217-221 (eval = FALSE)
###################################################
##   data(GlobalPatterns)
##   data(esophagus)
##   data(enterotype)
##   data(soilrep) 


###################################################
### code chunk number 9: phyloseq_basics.Rnw:234-236
###################################################
  data(GlobalPatterns)
  GlobalPatterns


###################################################
### code chunk number 10: phyloseq_basics.Rnw:279-283 (eval = FALSE)
###################################################
## otu1 <- otuTable(raw_abundance_matrix, speciesAreRows=FALSE)
## sam1 <- sampleData(raw_sample_data.frame) 
## tax1 <- taxTab(raw_taxonomy_matrix)
## tre1 <- read.nexus(my_nexus_file)


###################################################
### code chunk number 11: phyloseq_basics.Rnw:288-289 (eval = FALSE)
###################################################
## ex1b <- phyloseq(my_otuTable, my_sampleData, my_taxonomyTable, my_tree)


###################################################
### code chunk number 12: phyloseq_basics.Rnw:293-294 (eval = FALSE)
###################################################
## ex1c <- phyloseq(my_otuTable, my_sampleData)


###################################################
### code chunk number 13: phyloseq_basics.Rnw:356-357
###################################################
topN <- 20


###################################################
### code chunk number 14: phyloseq_basics.Rnw:360-363
###################################################
data(GlobalPatterns)
most_abundant_taxa <- sort(speciesSums(GlobalPatterns), TRUE)[1:topN]
ex2 <- prune_species(names(most_abundant_taxa), GlobalPatterns)


###################################################
### code chunk number 15: phyloseq_basics.Rnw:369-371
###################################################
topFamilies <- taxTab(ex2)[, "Family"]
as(topFamilies, "vector")


###################################################
### code chunk number 16: phyloseq_basics.Rnw:381-387 (eval = FALSE)
###################################################
## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
## f1  <- filterfunSample(topk(2))
## wh1 <- genefilterSample(testOTU, f1, A=2)
## wh2 <- c(T, T, T, F, F)
## prune_species(wh1, testOTU)
## prune_species(wh2, testOTU)


###################################################
### code chunk number 17: phyloseq_basics.Rnw:391-396
###################################################
data(GlobalPatterns)
f1  <- filterfunSample(topp(0.1))
wh1 <- genefilterSample(GlobalPatterns, f1, A=(1/2*nsamples(GlobalPatterns)))
sum(wh1)
ex2 <- prune_species(wh1, GlobalPatterns)


###################################################
### code chunk number 18: phyloseq_basics.Rnw:399-400
###################################################
  print(ex2)


###################################################
### code chunk number 19: phyloseq_basics.Rnw:407-412 (eval = FALSE)
###################################################
## data(GlobalPatterns)
## f1  <- filterfunSample(topf(0.9))
## wh1 <- genefilterSample(GlobalPatterns, f1, A=(1/3*nsamples(GlobalPatterns)))
## sum(wh1)
## prune_species(wh1, GlobalPatterns)


###################################################
### code chunk number 20: phyloseq_basics.Rnw:418-419
###################################################
ex3 <- subset_samples(GlobalPatterns, SampleType%in%c("Freshwater", "Ocean", "Freshwater (creek)"))


###################################################
### code chunk number 21: phyloseq_basics.Rnw:422-423
###################################################
ex3


###################################################
### code chunk number 22: phyloseq_basics.Rnw:429-430
###################################################
subset(sampleData(GlobalPatterns), SampleType%in%c("Freshwater", "Ocean", "Freshwater (creek)"))


###################################################
### code chunk number 23: phyloseq_basics.Rnw:436-437
###################################################
ex4 <- subset_species(GlobalPatterns, Phylum=="Firmicutes")


###################################################
### code chunk number 24: phyloseq_basics.Rnw:440-441
###################################################
ex4


###################################################
### code chunk number 25: phyloseq_basics.Rnw:447-449
###################################################
randomSpecies100 <- sample(species.names(GlobalPatterns), 100, replace=FALSE)
ex5 <- prune_species(randomSpecies100, GlobalPatterns)


###################################################
### code chunk number 26: phyloseq_basics.Rnw:459-461 (eval = FALSE)
###################################################
## data(GlobalPatterns)
## ex2 <- transformsamplecounts(GlobalPatterns, I)


###################################################
### code chunk number 27: phyloseq_basics.Rnw:468-469
###################################################
ex4  <- transformsamplecounts(GlobalPatterns, threshrankfun(500))


###################################################
### code chunk number 28: phyloseq_basics.Rnw:480-481 (eval = FALSE)
###################################################
## ex6 <- taxglom(GlobalPatterns, taxlevel="Genus")


###################################################
### code chunk number 29: phyloseq_basics.Rnw:487-488 (eval = FALSE)
###################################################
## ex7 <- tipglom(GlobalPatterns, speciationMinLength = 0.05)


###################################################
### code chunk number 30: phyloseq_basics.Rnw:527-529 (eval = FALSE)
###################################################
## 	source("http://bioconductor.org/biocLite.R")
## 	biocLite("phyloseq")


###################################################
### code chunk number 31: phyloseq_basics.Rnw:539-542 (eval = FALSE)
###################################################
## 	source("http://bioconductor.org/biocLite.R")
## 	biocLite("multtest")
## 	biocLite("genefilter")


###################################################
### code chunk number 32: phyloseq_basics.Rnw:547-556 (eval = FALSE)
###################################################
##   install.packages("ape")
##   install.packages("doParallel")
##   install.packages("foreach")
##   install.packages("ggplot2")
##   install.packages("igraph")
##   install.packages("picante")
##   install.packages("vegan")
##   install.packages("RJSONIO")
##   install.packages("plyr")


###################################################
### code chunk number 33: phyloseq_basics.Rnw:561-567 (eval = FALSE)
###################################################
##   # Make sure you have devtools installed
##   install.packages("devtools")
##   # Load the devtools package
##   library("devtools")
##   # Build and install phyloseq
##   install_github("phyloseq", "joey711")


###################################################
### code chunk number 34: phyloseq_basics.Rnw:573-577 (eval = FALSE)
###################################################
##   install.packages("doParallel")
##   install.packages("doMC")
##   install.packages("doSNOW")
##   install.packages("doMPI")

