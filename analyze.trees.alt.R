############################
###Analyze trees

N <- 10000
herit <- 1

time <- c(seq(0, 4, by = 0.001))

sel.intenses <- .005
n.locis <- c(100)
n_chromss <- c(200)

ts=commandArgs(trailingOnly=TRUE)[1]
t.offs =commandArgs(trailingOnly=TRUE)[2]
#####
input= commandArgs(trailingOnly=TRUE)[3]



phen_nums <- c(1:100) #1 number for each rep we want to do at each combination of parameters



pars <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)
# print(pars[1,])



err.array <- array(dim = c(length(time),5,dim(pars)[1]))   
print(dim(err.array))
err.std.array <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

err.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat_rent <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

err.array.rent.adj <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array.rent.adj <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat_rent_adj <- matrix(-1, nrow = dim(pars)[1], ncol = 12)


err.array.relate <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array.relate <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat_relate <- matrix(-1, nrow = dim(pars)[1], ncol = 12)


mat.true.phentrajs <- matrix(nrow = length(time), ncol = dim(pars)[1])




for(iter in 1:dim(pars)[1]){


	#####

	sel.intens <- pars[iter,1]
	n.loci <- pars[iter,2]
	n_chroms <- pars[iter,3]
	t <- pars[iter,4]
	t.off <- pars[iter,5]
	phen_num <- iter

	tree_dat <- sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_%s/out_%s/loci100_sintens0.005_N10000_nchr200_ton%s_toff%s_herit1_%s_alt.RData", ts, iter, ts, t.offs, iter)


    load(tree_dat)



	#print(paste("iter: ", iter, sep=""))
	# source("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/analyze_sim_true.R")
	# source("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/analyze_sim_rent.R")
	# source("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/analyze_sim_relate.R")
	# source("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/analyze_sim_rent_adj.R")
	#print (n.loci)
	source(sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/analyze_sim_%s.R", input))
	#print(paste("trial", as.character(iter), "complete."))
}



out_dir<- sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_%s/out_all_alt/", ts)
dir.create(out_dir)
save.image(paste(out_dir,"analyzed_trees_",input, ".RData", sep = ""))