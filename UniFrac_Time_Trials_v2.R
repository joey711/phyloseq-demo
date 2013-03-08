################################################################################
# Accuracy and time trials for different UniFrac implementations
################################################################################
rm(list=ls())
library("phyloseq"); library("foreach"); library("doParallel")
save.data.name <- "~/Dropbox/R/unifrac_time_trials_2_2012_02_07_6pm.RData"
# save.data.name <- "~/Dropbox/R/unifrac_time_trials_2_2012_02_03_11pm.RData"
#load(save.data.name)
data(ex1)
# # # # # load("~/Dropbox/R/unifrac_time_trials_2.RData")

# Threshold, in seconds. The maximum allowed time for one trial.
# Trials will not be conducted for that argument combination after
# this amount of time...
threshold <- 200 # seconds
# threshold <- 5 # seconds

# Number of cores to use for parallel trials
Ncores <- 8

################################################################################
# Create a vector of the numbers of species
# And then create a list of phyloseq objects for each one. Will loop on the list
################################################################################
# Define vector of the number of taxa to include
Nspecies <- c(seq(100, 500, 100), seq(600, nspecies(ex1), 300), nspecies(ex1) )
# Nspecies <- c(seq(100, 3000, 300) )

physeq_list <- vector("list", length(Nspecies))
names(physeq_list) <- as.character(Nspecies)

for( N in Nspecies ){
	########################################
	# Create random subset with N species
	########################################
	# Randomly subset N species from the original
	ex   <- prune_species(sample(species.names(ex1), N), ex1)
	# Pick a random root from the sub-sampled dataset
	tree <- root(tre(ex), sample(species.names(ex), 1), resolve.root=TRUE)
	OTU  <- otuTable(ex)
	# Rebuild to make sure you didn't screw up something.
	ex   <- phyloseq(OTU, tree)
	# Add to the list of physeq objects
	physeq_list[[as(N, "character")]] <- ex			
}

################################################################################
# Create table of combinations of function calls
################################################################################

UFcalls <- expand.grid(
	parallel=c(TRUE, FALSE),
	weighted=c(TRUE, FALSE),
	fast=c(TRUE, FALSE),
KEEP.OUT.ATTRS=FALSE)
# The following creates a list appropriate for do.call():
# UFcalls[1, , drop=TRUE]

################################################################################
# Loop over each incrementally more-complex data subset
# Stop tabulating a combination if run-time gets above a threshold
################################################################################
########################################	
# For each element of physeq_list, 
# track the time for each function call
# indicated by UFcalls
########################################
# Initialize data.frame of results
time.frame <- data.frame(matrix(NA, 0, 5, FALSE))

# index for keeping track of which call combinations in UFcalls are still below time threshold 
stillOK <- rep(TRUE, nrow(UFcalls))
# Loop over each N in Nspecies.
for( N in Nspecies ){ # N = Nspecies[1]
	if( sum(stillOK) > 0 ){
		for( i in which(stillOK) ){ # i<-1
			# Create list of arguments to UniFrac
			call.list <- c(physeq=physeq_list[[as.character(N)]], UFcalls[i, , drop=TRUE] )
			
			# Determine if should register parallel backend or sequential
			if( call.list$parallel){
				registerDoParallel(makeCluster(Ncores))
			}
			
			# Build, and call UniFrac with list of arguments
			runtime <- system.time(do.call("UniFrac", call.list))["elapsed"]
			
			# Store the values in growing data.frame
			time.frame <- rbind(time.frame, data.frame(species=N, runtime, UFcalls[i, , drop=TRUE]) )
			
			# Check if this runtime was below threshold, if not, mark it
			if( runtime > threshold ){
				# mark the combination as too long
				stillOK[i] <- FALSE
			}
		}		
	}
}
################################################################################
# Now run one more loop for picante::unifrac() time (unweighted parallel)
################################################################################
 # Initialize data.frame of results for picante
picante.frame <- data.frame(matrix(NA, 0, 2, FALSE))

# index for keeping track of which call combinations in UFcalls are still below time threshold 
stillOK <- TRUE
# Loop over each N in Nspecies.
for( N in Nspecies ){ # N = Nspecies[1]
	if( stillOK > 0 ){
		# Create list of arguments to UniFrac
		call.list <- c(physeq=physeq_list[[as.character(N)]], UFcalls[i, , drop=TRUE] )
		
		# Determine if should register parallel backend or sequential
		if( call.list$parallel){
			registerDoParallel(makeCluster(Ncores))
		}
		
		# calculate, score run-time
		exOTU <- as( t(otuTable(physeq_list[[as.character(N)]])), "matrix")
		extre <- tre(physeq_list[[as.character(N)]])
		runtime <- system.time( picante::unifrac(exOTU, extre) )["elapsed"]
		
		# Store the values in growing data.frame
		picante.frame <- rbind(picante.frame, data.frame(species=N, runtime) )
		
		# Check if this runtime was below threshold, if not, mark it
		if( runtime > threshold ){
			# stop loop
			stillOK <- FALSE
		}
	}
}
save(picante.frame, file=ps(save.data.name, "_picante.frame.RData"))
################################################################################
# That part takes a long time. Save the data before plotting. Plotting is quick
################################################################################
save.image(save.data.name)
# load(save.data.name)
# load(ps(save.data.name, "_picante.frame.RData"))
################################################################################
# Now build final data.frame of results
################################################################################
# Modify logical parameter values to informative text factor
time.frame[,"weighted"] <- factor(ifelse(time.frame[,"weighted"], "weighted", "unweighted"))
time.frame[,"fast"] <- factor(ifelse(time.frame[,"fast"], "fast", "orig"))
time.frame[,"parallel"] <- factor(ifelse(time.frame[,"parallel"], "parallel", "sequential"))

# Add an "Algorithm" column that summarizes the algorithm approach
time.frame.plot <- data.frame(time.frame, algorithm=paste(time.frame[, "fast"], time.frame[, "weighted"]) )

# add picante.frame to time.frame
picante.frame.plot <- data.frame(picante.frame, parallel="sequential", weighted="unweighted", fast="orig", algorithm="picante::unifrac")

# Finally, create a new separate, everything-combined data.frame for plotting
time.df <- rbind(time.frame.plot, picante.frame.plot)

################################################################################
# Now build ggplot of results.
################################################################################
# # # # # my.round <- function(x){
	# # # # # #x <- 0.00006476846534684645
    # # # # # format(round(as(x, "numeric")), scientific=FALSE)
# # # # # }

# # # # # TransInvNegLog10 <- Trans$new("InvNegLog10", f = function(x) 10^(-x), 
	# # # # # inverse = function(x) -log10(x), labels = function(x) x)

# # # # # TransInvNegLog10b <- Trans$new("InvNegLog10b", f = function(x) 10^(-x), 
	# # # # # inverse = function(x) -log10(x), labels = function(x) bquote(10^.(-log10(x)))
	# # # # # ) 

# Define your own log-scale transformation for pretty plotting.
TransSimpleLog10 <- Trans$new("SimpleLog10", f = function(x) log10(x), 
	inverse = function(x) 10^x, labels = function(x) round(10^x)
	)
	
ggplot(time.df, 
	aes(species, runtime, color = algorithm, shape=parallel )) + 
	# aes(species, log(runtime, 10), color = algorithm )) + 
	geom_path() +
	geom_point(size=3) +
	# facet_grid(weighted ~ parallel) +
	# facet_wrap(~parallel, 2) +
	# stat_smooth(method = "lm", formula = y ~ x) +
	# stat_smooth(method="loess") +	
	opts(title = ps("UniFrac performance in phyloseq: ", nsamples(ex), " samples, ", Ncores, " cores")) + 
	scale_colour_discrete(name = "Algorithm") + 
	# xlim(0, 600) + 
	# ylim(0, 100) +
	# scale_y_continuous("log10-runtime [seconds]")
	# scale_y_log10("runtime [seconds]", formatter="precision")
	# scale_y_log10("runtime [seconds]")
	scale_y_continuous("runtime [seconds]", trans= TransSimpleLog10)
	# scale_y_log10()
################################################################################
# save.image( ps(save.data.name, "_final_ggplot2_120208_12pm.RData") )