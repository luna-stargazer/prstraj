
#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.
ts=commandArgs(trailingOnly=TRUE)[1]
t.offs=commandArgs(trailingOnly=TRUE)[2]
phen_nums=commandArgs(trailingOnly=TRUE)[3]
N <- 10000
herit <- 1
tout_dir<-sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_%s/temp_out_%s/", ts, phen_nums)
out_dir <- sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_%s/out_%s/", ts, phen_nums)

#nsamps <- 20
sel.intenses <- .005


n.locis <- c(100)
n_chromss <- c(200)


#n.locis <- c(10)
#n_chromss <- c(20)
#phen_nums <-start #1 number for each rep we want to do at each combination of parameters

traj.fn <- sprintf("%stemp.txt", tout_dir)
msout.fn <-sprintf("%sms_out.txt",tout_dir)
rent_in_fn <- sprintf("%srent_in.txt", tout_dir)


helper_fn <- "/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package
print("done sourcing help fn") 
ms_dir <- "/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/msseldir/"
rentplus_fn <- "/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/RentPlus.jar"
len_hap <- 200000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
sel_site <- 100000 #the position of the selected site in the haplotype
u <- 2e-8 #the neutral mutation rate per base pair/generation
r <- 2.5e-8 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers
sd.trait <- 1

time <- c(seq(0, 4, by = 0.001))
pars <- expand.grid(sel.intenses, n.locis, n_chromss, as.numeric(ts), as.numeric(t.offs), phen_nums)

print("start sim")
for(k in 1:dim(pars)[1]){
	sel.intens <- pars[k,1]
	n.loci <- pars[k,2]
	n_chroms <- pars[k,3]
	t <- pars[k,4]
	t.off <- pars[k,5]
	phen_num <- pars[k,6]
	source("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/pheno_sim_1iter.R")
	#source("../pheno_sim_1iter_norent.R")
	print(paste("trial", as.character(phen_num), "complete."))
}


save.image(paste(out_dir,"sim_trees_",phen_nums, ".RData", sep = ""))





