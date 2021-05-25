#neutral - relate
trajs_neut_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
vars_neut_bin_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
#vars_neut_post_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
for(i in 1:length(relate.all.rescale)){	
	trajs_neut_relate[,i] <- est_af_traj_neut(lins.list.relate[[i]])
	vars_neut_bin_relate[,i] <- est_af_var_neut_bin(lins.list.relate[[i]])
#	vars_neut_post_relate[,i] <- est_af_var_neut_post(lins.list.relate[[i]])
}

print("here1")

trajs_neut_relate[time == 0,] <- (n_ders/n_chroms) #in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
vars_neut_bin_relate[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms
traj.phen.neut.relate <- 2 * trajs_neut_relate %*%  eff_sizes 
var.phen.neut.bin.relate <- 4 * vars_neut_bin_relate %*% eff_sizes^2
print("here2")
#var.phen.neut.post.relate <- 4 * vars_neut_post_relate %*% eff_sizes^2




#Method of moments from smoothed coalescent time estimates---relate.
trajs_mom_smoothtime_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
trajs_var_mom_smoothtime_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
for(i in 1:length(relate.all.rescale)){		
	trajs_mom_smoothtime_relate[,i] <- est_af_traj_mom.smoothtime(lins.list.relate[[i]], time)
	trajs_var_mom_smoothtime_relate[,i] <- est_af_var_mom.smoothtime(lins.list.relate[[i]], time*2*N)
}

print("here3")
traj.phen.mom_smoothtime_relate <- 2 * trajs_mom_smoothtime_relate %*%  eff_sizes 
var.phen.mom_smoothtime_relate <- 4 * trajs_var_mom_smoothtime_relate %*%  eff_sizes^2 
traj.phen.mom_smoothtime_relate[time == 0] <- traj.phen.neut.relate[time == 0]
var.phen.mom_smoothtime_relate[time == 0] <- var.phen.neut.bin.relate[time == 0]

print("here4")



#waiting time-based estimates and variance---relate
trajs_est_wt_l1_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
trajs_var_wt_l1_relate <- matrix(nrow = length(time), ncol = length(relate.all.rescale))
for(i in 1:length(relate.all.rescale)){
	wt.estvar.relate <- p_ests_wait(times.c.relate[[i]], time, ell.ref = 1, ell.alt = 1)
	trajs_est_wt_l1_relate[,i] <- wt.estvar.relate[,1]	
	trajs_var_wt_l1_relate[,i] <- wt.estvar.relate[,2]	
}
print("here5")
traj.phen.wt_l1_relate <- 2 * trajs_est_wt_l1_relate %*%  eff_sizes 
var.phen.wt_l1_relate <- 4 * trajs_var_wt_l1_relate %*%  eff_sizes^2 
traj.phen.wt_l1_relate[time == 0] <- traj.phen.neut.relate[time == 0]
var.phen.wt_l1_relate[time == 0] <- var.phen.neut.bin.relate[time == 0]
print("here6")


true.per.time <- numeric(0)
for(j in 1:length(time)){
	true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}
print("here7")

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
	true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}
print("here8")





#save errors, unscaled and scaled by estimated se, for each method at all times---relate.

err.neut.relate <- traj.phen.neut.relate - true.per.time
print("here9")
err.smoothmom.relate <- traj.phen.mom_smoothtime_relate - true.per.time
print("here10")
err.wt_l1.relate <- traj.phen.wt_l1_relate - true.per.time
print("here11")


err.neut.std.bin.relate <- err.neut.relate / sqrt(var.phen.neut.bin.relate)
print("here12")
err.smoothmom.std.relate <- err.smoothmom.relate / sqrt(var.phen.mom_smoothtime_relate)
print("here13")
err.wt_l1.std.relate <- err.wt_l1.relate / sqrt(var.phen.wt_l1_relate)
err.smoothmom.std.relate <- err.smoothmom.relate / sqrt(var.phen.mom_smoothtime_relate)
print("here14")


err.mat.relate <- cbind(err.neut.relate, err.smoothmom.relate, err.wt_l1.relate)
err.mat.std.relate <- cbind(err.neut.std.bin.relate, err.smoothmom.std.relate, err.wt_l1.std.relate)
print("here15")

err.array.relate[,,iter] <- err.mat.relate
err.std.array.relate[,,iter] <- err.mat.std.relate
mat.true.phentrajs[,iter] <- true.per.time


#Tests
qxtest_mat_relate[iter,1:3] <- Qx_test(trajs_neut_relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
print("here16")
qxtest_mat_relate[iter,4] <- Qx_test(trajs_neut_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
print("here17")
qxtest_mat_relate[iter,5:7] <- Qx_test(trajs_mom_smoothtime_relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
print("here18")
qxtest_mat_relate[iter,8] <- Qx_test(trajs_mom_smoothtime_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
print("here18")
qxtest_mat_relate[iter,9:11] <- Qx_test(trajs_est_wt_l1_relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
print("here19")
qxtest_mat_relate[iter,12] <- Qx_test(trajs_est_wt_l1_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
print("here20")

print("relate tree T_X statistic, number of timepoints, and permutation p")
print(Qx_test(trajs_neut_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000))
