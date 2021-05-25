on_list=c(0,0.02,0.04)
off_list=c(0,0,0.02)


for (j in 1:3){
ts=on_list[j]
t.offs=off_list[j]

departure=data.frame(matrix(nrow=100, ncol=4))
colnames(departure)<-c("RENT_fixed", "RENT_seg", "Relate_fixed", "Relate_seg")

for(iter in 1:100){


tree_dat <- sprintf("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_%s/out_%s/loci100_sintens0.005_N10000_nchr200_ton%s_toff%s_herit1_%s_alt.RData", ts, iter, ts, t.offs, iter)

load(tree_dat)

time=c(seq(0, 4, by = 0.001), 20)

rent_loci=c()
relate_loci=c()

trajs_neut_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
vars_neut_bin_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
#vars_neut_post_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
for(i in 1:length(rent_trees_list)){	
	trajs_neut_rent[,i] <- est_af_traj_neut(lins.list.rent[[i]])
	vars_neut_bin_rent[,i] <- est_af_var_neut_bin(lins.list.rent[[i]])
#	vars_neut_post_rent[,i] <- est_af_var_neut_post(lins.list.rent[[i]])
}


#neutral - relate
trajs_neut_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
vars_neut_bin_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
#vars_neut_post_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
for(i in 1:length(relate.all.rescale)){	
	trajs_neut_relate[,i] <- est_af_traj_neut(lins.list.relate[[i]])
	vars_neut_bin_relate[,i] <- est_af_var_neut_bin(lins.list.relate[[i]])
#	vars_neut_post_relate[,i] <- est_af_var_neut_post(lins.list.relate[[i]])
}


departure[iter,]=c(length(which(trajs_neut_rent[4002,]==1))/100, length(which(trajs_neut_rent[4002,]!=1 &trajs_neut_rent[4002,]!=0))/100, length(which(trajs_neut_relate[4002,]==1))/100, length(which(trajs_neut_relate[4002,]!=1 & trajs_neut_relate[4002,]!=0))/100)







}
assign(sprintf("departure%s", ts), departure)

}



save.image(paste("/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/departureV2addtp.alt", ".RData", sep = ""))


