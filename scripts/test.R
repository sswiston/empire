############################################
#### LOADING ELLLIPSE EVOLUTION PACKAGE ####
############################################

library(devtools)
devtools::load_all("~/projects/ellipses/package")

library(ape)

######################
#### LOADING DATA ####
######################

outpath <- "~/projects/ellipses/simulations/output/test"

set.seed(2)
example_tree <- phytools::pbtree(b=.1,d=0,t=40,type="continuous")
tree <- dataTree$new(tree=example_tree)
tree$sim_data()
tree$save(filepath=outpath,prefix="true")
write.csv(tree$tip_data,file=paste0(outpath,"/tip_data.csv"))

##############
#### MCMC ####
##############

testing_conditions <- c("None") #testing_conditions <- c("under_prior","test_probs","record_proposals")
proposal_weights=list(sigma_x=5,sigma_y=5,sigma_r=5,sigma_s=5,sigma_a=5,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=5,kappa=5,W_d=1,W_m=1,W_c=1,W_h=1,V_r=3,V_s=3,V_a=3,tip=1)

mcmc <- make_MCMC(tree,proposal_weights=proposal_weights)
run_MCMC(mcmc,iterations=1000,moves_per_iteration=1,burnin=100,thinning=1,testing_conditions=testing_conditions,filepath=outpath)

######################
#### LOADING DATA ####
######################

outpath <- "~/projects/ellipses/skinks/output/test"

example_tree <- ape::read.tree("~/projects/ellipses/skinks/data/spheno.tre")
example_tree$tip.label <- gsub("_"," ",example_tree$tip.label)
data <- read.csv("~/projects/ellipses/skinks/data/ellipse_data.csv",header=TRUE)

##############
#### MCMC ####
##############

tree <- dataTree$new(tree=example_tree,tip_data=data)
mcmc <- make_MCMC(tree,
                  prior_root_x=distributions3::Normal(mean(data$x),sd(data$x)),
                  prior_root_y=distributions3::Normal(mean(data$y),sd(data$y)),
                  prior_root_r=distributions3::Normal(mean(data$r),sd(data$r)),
                  prior_root_s=distributions3::Normal(mean(data$s),sd(data$s)),
                  prior_root_a=distributions3::Normal(mean(data$a),sd(data$a)),
                  prior_sigma_x=distributions3::Uniform(0,10),
                  prior_sigma_y=distributions3::Uniform(0,10),
                  prior_sigma_r=distributions3::Uniform(0,10),
                  prior_sigma_s=distributions3::Uniform(0,10),
                  prior_sigma_a=distributions3::Uniform(0,10),
                  prior_mu=distributions3::Normal(mean(data$a),sd(data$a)),
                  prior_kappa=distributions3::Uniform(0,10),
                  proposal_weights=list(sigma_x=5,
                                        sigma_y=5,
                                        sigma_r=5,
                                        sigma_s=5,
                                        sigma_a=5,
                                        root_x=3,
                                        root_y=3,
                                        root_r=3,
                                        root_s=3,
                                        root_a=3,
                                        mu=5,
                                        kappa=5,
                                        W_d=1,
                                        W_m=1,
                                        W_c=1,
                                        W_h=1,
                                        V_r=3,
                                        V_s=3,
                                        V_a=3,
                                        tip=1))

testing_conditions <- c("None")
#testing_conditions <- c("under_prior","test_probs","record_proposals")

run_MCMC(mcmc,iterations=10000,moves_per_iteration=1,burnin=1000,thinning=10,testing_conditions=testing_conditions,filepath=outpath)

######################
#### LOADING DATA ####
######################

# args <- commandArgs(trailingOnly = TRUE)
# NUMBER <- args[1]
# CONDITION <- args[2]
#
# set.seed(NUMBER)
# print(paste0("Number: ",NUMBER))
# print(paste0("Condition: ",CONDITION))
# name <- paste0(CONDITION,".sim.",NUMBER)
# sim_directory <- paste0("./simulations/output/",name,"/")

sim_directory <- "~/projects/ellipses/simulations/output/test/"

sim_tree <- ape::read.tree(paste0(sim_directory,"true.tree.txt"))
sim_V <- read.csv(paste0(sim_directory,"true.V.tsv"),sep="\t")
sim_XY <- read.csv(paste0(sim_directory,"true.XY.tsv"),sep="\t")
sim_data <- cbind(sim_XY[1:length(sim_tree$tip.label),2:3],sim_V[1:length(sim_tree$tip.label),2:4])
sim_data$taxon <- sim_tree$tip.label

# Sampling half (rounded up) of taxa to remove due to extinction
num_total <- length(sim_tree$tip.label)
num_removed <- floor(num_total/2)
removed_taxa <- sample(1:num_total,num_removed)
trimmed_tree <- drop.tip(sim_tree,removed_taxa)
trimmed_data <- sim_data[-c(removed_taxa),]

tree <- dataTree$new(tree=trimmed_tree)
tree$tip_data <- trimmed_data

##############
#### MCMC ####
##############

proposal_weights=list(sigma_x=5,sigma_y=5,sigma_r=5,sigma_s=5,sigma_a=5,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=5,kappa=5,W_d=1,W_m=1,W_c=1,W_h=1,V_r=3,V_s=3,V_a=3,tip=1)

#print(paste0("Performing MCMC with Seed ",NUMBER))
mcmc <- make_MCMC(tree,proposal_weights=proposal_weights)
run_MCMC(mcmc,iterations=10000,moves_per_iteration=1,burnin=1000,thinning=10,filepath=paste0(sim_directory,"sensitivity"))
