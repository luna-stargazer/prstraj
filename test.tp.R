#########
helper_fn <- "/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/helper_functions_coal_sel.R"
source(helper_fn)

on_list=c(0,0.02,0.04)
off_list=c(0,0,0.02)

for (j in 1:1){
ts=on_list[j]
t.offs=off_list[j]

departure=data.frame(matrix(nrow=100, ncol=4))
colnames(departure)<-c("RENT_fixed", "RENT_seg", "Relate_fixed", "Relate_seg")

for (iter in 1:2){

tree_dat <- sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_%s/out_%s/loci100_sintens0.005_N10000_nchr200_ton%s_toff%s_herit1_%s_alt.RData", ts, iter, ts, t.offs, iter)

load(tree_dat)

time=c(seq(0, 4, by = 0.001), 20)

times.c <- list()
lins.list <- list()
times.c.rent <- list()
lins.list.rent <- list()
times.c.rent.adj <- list()
lins.list.rent.adj <- list()
times.c.relate <- list()
lins.list.relate <- list()

for(i in 1:10){
	times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = TRUE, units_in = 4)
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)

	times.c.rent[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = TRUE, units_in = 2)
	lins.list.rent[[i]] <- times_to_lins(times.c.rent[[i]], time)

	times.c.rent.adj[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = TRUE, units_in = 4)
	lins.list.rent.adj[[i]] <- times_to_lins(times.c.rent.adj[[i]], time)

	times.c.relate[[i]] <- trees_to_times(relate.all.rescale[[i]], relate.anc.rescale[[i]], relate.der.rescale[[i]], time, sure.alt.is.derived = TRUE, units_in = 2)
	lins.list.relate[[i]] <- times_to_lins(times.c.relate[[i]], time)

}


trajs_neut_rent <- matrix(nrow = length(time), ncol = 10)
vars_neut_bin_rent <- matrix(nrow = length(time), ncol = 10)
for(i in 1:10){	
	trajs_neut_rent[,i] <- est_af_traj_neut(lins.list.rent[[i]])
	vars_neut_bin_rent[,i] <- est_af_var_neut_bin(lins.list.rent[[i]])
#	vars_neut_post_rent[,i] <- est_af_var_neut_post(lins.list.rent[[i]])
}


#neutral - relate
trajs_neut_relate <- matrix(nrow = length(time), ncol = 10)
vars_neut_bin_relate <- matrix(nrow = length(time), ncol = 10)
#vars_neut_post_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
for(i in 1:1:10){	
	trajs_neut_relate[,i] <- est_af_traj_neut(lins.list.relate[[i]])
	vars_neut_bin_relate[,i] <- est_af_var_neut_bin(lins.list.relate[[i]])
#	vars_neut_post_relate[,i] <- est_af_var_neut_post(lins.list.relate[[i]])
}

departure[iter,]=c(length(which(trajs_neut_rent[4002,]==1))/10, length(which(trajs_neut_rent[4002,]!=1 &trajs_neut_rent[4002,]!=0))/10, length(which(trajs_neut_relate[4002,]==1))/10, length(which(trajs_neut_relate[4002,]!=1 & trajs_neut_relate[4002,]!=0))/10)
print("iter successful")

}
assign(sprintf("departure%s", ts), departure)

}

save.image(paste("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/departureV2addt.alt", ".RData", sep = ""))



