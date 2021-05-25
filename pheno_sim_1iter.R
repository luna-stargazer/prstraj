#Doc Edge, 6/13/18
#Goal: Simulate phenotypic scores (possibly under selection) that are additive functions
# of genotypes at many loci. 


#Plan: simulate n_loci allele frequency trajectories.
#to do this, select an effect size from normal dist with mean 0 and sd
#equal to heritability over n_loci (from Martin et al). locus undergoes selection 
#according to effect size and selection gradient of trait (which changes through time). 
#Generate an allele frequency trajectory obeying the selection coefficient.
#Do this for all loci and save "phenotype" trajectory.
#Run through ms and save trees to estimate phenotype trajectory using coalescent times.
#Run haplotypes through RENT+ to see whether estimated trees recover phenotype trajectory.


trajs <- list()
sel.parts <- list()
eff_sizes <- numeric(0)
ss <- numeric(0)
curr.freqs <- numeric(0)


#need a constant that reflects a harmonic series from 1 to 2N-1, but only the
#terms where i/(2N) \in [.01, .99]. for N large this will be very close to the 
#analogous integral, which is ln(.99) - ln(.01) ~= 4.595
ser <- 1:(2*N-1)
harm.const <- sum(1/ser[ser >= .01*2*N & ser <= .99*2*N])



#Simulate and write n.loci allele-frequency trajectories. None of these should have fixed,
#so reject and do it again if fixed. Reject if minor allele frequency is < 0.01.
shift.achieved <- 0
while(shift.achieved == 0){
	pre.sel.freqs <- numeric(0)
	post.sel.freqs	<- numeric(0)
	for(i in 1:n.loci){
		fix <- 1
		while(fix == 1){
			eff_size <- rnorm(1,0,sqrt(herit * sd.trait^2 * harm.const /n.loci))
			s_locus <- eff_size * sel.intens / sd.trait #Charlesworth B3.7.7
			ss[i] <- s_locus
			p0 <- gen.neut.sfs.condit(1, 2*N, .01) #select from neutral sfs conditional on common.
			pre.sel.freqs[i] <- p0 
			sel.part <- sel.traj(p0, s = s_locus, N = N, t = t-t.off)
			post.sel.freqs[i] <- sel.part[nrow(sel.part),2]
			if(t.off > 0){
				post.drift <- neut.traj.time(p0 = sel.part[nrow(sel.part),2], N, t.off)
				sel.part <- rbind(sel.part, cbind(post.drift[,1] + t - t.off, post.drift[,2] ))			
			}
			if(sel.part[nrow(sel.part),2] >= 0.01 & sel.part[nrow(sel.part),2] <= 0.99){fix = 0}
		}
		eff_sizes[i] <- eff_size 
		curr.freqs[i] <- sel.part[nrow(sel.part),2]
		sel.parts[[i]] <- sel.part
	}	
	#Check whether the achieved shift is close to the target.
	trait.sd <- sqrt(sum(2 * curr.freqs * (1 - curr.freqs) *eff_sizes^2))
	sel.shift <- sum(2 * eff_sizes * (post.sel.freqs - pre.sel.freqs)) / trait.sd
	target.shift <- (sel.intens * sd.trait) * (t - t.off) * 2 * N
	if(target.shift*.95 <= sel.shift & sel.shift <= target.shift*1.05){shift.achieved <- 1}
#	print("trait.attempted")
}


for(i in 1:n.loci){
	sel.part <- sel.parts[[i]]
	driftup <- neut.traj(pre.sel.freqs[i], N, loss = TRUE)
	traj.fwd <- rbind(cbind(driftup[,1], rev(driftup[,2])), cbind(sel.part[-1,1] + max(driftup[,1]), sel.part[-1,2] ))
	traj <- traj.fwd
	traj[,2] <- rev(traj.fwd[,2])
	trajs[[i]] <- traj
}

eff_sizes_unscaled <- eff_sizes
eff_sizes <- eff_sizes / trait.sd #The SD of the polygenic score is set
#to be 1 in the present


rm(ser)
rm(harm.const)
rm(eff_size)
rm(s_locus)
rm(p0)
rm(sel.part)
#rm(post.drift)
rm(post.sel.freqs)
rm(pre.sel.freqs)


mat.trajs <- matrixify.list.of.trajs(trajs)

phen.traj <- as.numeric(2 * eff_sizes %*% t(mat.trajs))  
pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) )




n_ders <- numeric(n.loci)
ms_trees_list <- list()
rent_trees_list <- list()
ms_haplotypes_list <- list()
relate_trees_list <- list()


#for each locus, write trajectory, run ms to simulate a sample and tree,
#and run rent+ to infer tree. Save both the ms and the rent+ trees, as well
#as the number of derived alleles for each locus.
for(i in 1:n.loci){
	write.traj(traj.fn, trajs[[i]])
	curr.freq.der <- curr.freqs[i]
	n_der <- rbinom(1, n_chroms, curr.freq.der) #number of chroms with derived allele.
	if(n_der == 0){n_der <- n_der + 1}
	if(n_der == n_chroms){n_der <- n_der - 1}
	n_ders[i] <- n_der

	counter=0
	success=1
	while (success==1 & counter<=100){
	
	#run mssel
	ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(4*N*(1 - dbinom(0,len_hap, r))), " ", as.character(len_hap), " -t ", as.character(4*N* (1 - dbinom(0, len_hap, u))), " -T -L > ", msout.fn, sep = "")
	system(ms.string)
	
	
	#cut down mssel output to just the part rent+ needs
	process.string <- paste("sed '/positions: /,$!d' ", msout.fn, " | sed 's/.*positions: //' > ", rent_in_fn, sep = "")
	system(process.string)

	#run rent+
	rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_in_fn, " ", msout.fn, sep = "")
	success=system(rent.string)
	counter=counter+1
	#system(rent.string)
	}

	print(paste("loci",i, "attempt:", counter, sep=" "))

	ms_haplotypes_list[[i]] <- readLines(rent_in_fn) #save the haplotypes produced by ms (and fed to rent+) in a list entry.
	#the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
	#site at 0.05. The remaining entries are string haplotypes of 0s and 1s.
	

	#Read in trees at selected site in rent+
	rent_trees_fn <- paste(rent_in_fn, ".trees", sep = "")
	rent_trees <- readLines(rent_trees_fn)
	#which tree is for the selected site?
	rent_sel_ind <- grep(paste(as.character(sel_site), "\t", sep = ""), rent_trees)

	#get tree at selected site in newick format, both true and rent-estimated.
	rent_sel_tree_newick <- sub(paste(as.character(sel_site), "\t", sep = ""), "", rent_trees[rent_sel_ind])
	rent_trees_list[[i]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))
	
	#pull trees from ms
	msoutlines <- readLines(msout.fn)
	dists <- numeric(0)
	for(m in 1:(length(msoutlines) - 4)){
		suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
	}
	dists <- dists[!is.na(dists) & dists != Inf]
	posits.rs <- cumsum(dists)
	ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
	if(length(ms_to_extract)  == 0){ms_to_extract <- 1}
	ms_trees_list[[i]] <- read.tree(text = msoutlines[4 + ms_to_extract])		


	#generate relate output
	
	ms4relate.string <- sprintf("grep -v '\\[*\\]' %s/ms_out.txt| sed -e '3d;5,7d' > %s/ms.4relate.txt", tout_dir, tout_dir)
	
	system(ms4relate.string)

	con<-file(sprintf("%s/ms.4relate.txt", tout_dir))
	open(con)


    #extract positions of segregating loci
	exp=read.table(con,skip=4,nrow=1) 

	pos=as.numeric(exp[2:length(exp)])

	r <- 2.5e-8
	hap_len=200000
	P_r <- (1 - dbinom(0,hap_len, r))
	cM.dist <- 50 * log(1/(1 - 2*P_r))

    #actual position in bp/basepairs
	scaled_pos=pos*hap_len

	#calculate genetic map info
	gen_map=cM.dist*pos

	n=length(pos)
	recom_rate2=gen_map*10^6/scaled_pos

	map=cbind(scaled_pos, recom_rate2, gen_map)

	colnames(map)=c("pos","COMBINED_rate","Genetic_Map")

	write.table(map, sprintf("%s/recomb.map.txt",tout_dir), append = FALSE, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

        relate_dir="/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/relate_static"
	ms2relate.string <- sprintf("Rscript %s/ms2haps.R %s/ms.4relate.txt %s/relate 200000", relate_dir,tout_dir, tout_dir)

	system(ms2relate.string)

	#get rid of scientific notation in haplotype files by relate, othetrwise would throw errors in next step
	nosci.string<- paste("awk '{$3=sprintf(\"%.f\",$3)}1'", sprintf("%s/relate.haps > %s/relate.nosci.haps", tout_dir, tout_dir))

	system(nosci.string)


	relate.string <- sprintf("%s/bin/Relate --mode All -m 2e-8 -N 20000 --haps %s/relate.nosci.haps --sample %s/relate.sample --map %s/recomb.map.txt -o relate", relate_dir, tout_dir, tout_dir, tout_dir)     

	system(relate.string)


	# Reformat tree at site of interest(site under selection) to Newick format
	relate.tree <- sprintf("%s/bin/RelateExtract --mode TreeAtSNPAsNewick --anc %s/relate.anc --mut %s/relate.mut --bp_of_interest 100000 -o relate", relate_dir, tout_dir, tout_dir)

	system(relate.tree)

	relate_trees_list[[i]] <- read.tree(sprintf("%s/relate_at_100000.newick", tout_dir))




	print(paste("simulation locus", as.character(i), "of trait", as.character(k), "complete"))
}

#update tree tips to 1-based index 
rename.tips <- function(tree) {
    old<-as.numeric(tree$tip.label)
    new<-old+1
    tree$tip.label <- as.character(new)
    return(tree)
}


relate_trees_list=lapply(relate_trees_list, rename.tips)



#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.
for(i in 1:n.loci){
	if(is.null(names(ms_trees_list[[i]]))){
		cand.trees <- ms_trees_list[[i]]
		is.mono <- numeric(0)		
		for(j in 1:length(cand.trees)){
			der_tips <- which(cand.trees[[j]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
			is.mono[j] <- is.monophyletic(cand.trees[[j]], der_tips)
		}
		if(mean(is.mono) > 0){
			ind.to.take <- 	sample(which(is.mono == 1), 1)
		}		
		if(mean(is.mono) == 0){ #if none are monophyletic, choose at random
			ind.to.take <- sample(1:length(cand.trees), 1)
		}
		ms_trees_list[[i]] <- cand.trees[[ind.to.take]]
		#rent_trees_list[[i]] <- rent_trees_list[[i]][[ind.to.take]]
	}
}

for(i in 1:n.loci){
	if(is.null(names(rent_trees_list[[i]]))){
		cand.trees <- rent_trees_list[[i]]
		is.mono <- numeric(0)		
		for(j in 1:length(cand.trees)){
			der_tips <- which(cand.trees[[j]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
			is.mono[j] <- is.monophyletic(cand.trees[[j]], der_tips)
		}
		if(mean(is.mono) > 0){
			ind.to.take <- 	sample(which(is.mono == 1), 1)
		}		
		if(mean(is.mono) == 0){ #if none are monophyletic, choose at random
			ind.to.take <- sample(1:length(cand.trees), 1)
		}
		#ms_trees_list[[i]] <- cand.trees[[ind.to.take]]
		rent_trees_list[[i]] <- rent_trees_list[[i]][[ind.to.take]]
	}
}


#####
for(i in 1:n.loci){
	if(is.null(names(relate_trees_list[[i]]))){
		cand.trees <- relate_trees_list[[i]]
		is.mono <- numeric(0)		
		for(j in 1:length(cand.trees)){
			der_tips <- which(cand.trees[[j]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
			is.mono[j] <- is.monophyletic(cand.trees[[j]], der_tips)
		}
		if(mean(is.mono) > 0){
			ind.to.take <- 	sample(which(is.mono == 1), 1)
		}		
		if(mean(is.mono) == 0){ #if none are monophyletic, choose at random
			ind.to.take <- sample(1:length(cand.trees), 1)
		}
		#ms_trees_list[[i]] <- cand.trees[[ind.to.take]]
		relate_trees_list[[i]] <- relate_trees_list[[i]][[ind.to.take]]
	}
}



#split into ancestral and derived trees and retrieve coalescence times.
anc_trees_ms <- list()
der_trees_ms <- list()
ct_ms_all <- list()
ct_ms_anc <- list()
ct_ms_der <- list()

anc_trees_rent <- list()
der_trees_rent <- list()
ct_rent_all <- list()
ct_rent_anc <- list()
ct_rent_der <- list()


relate.anc.rescale <- list()
relate.der.rescale <- list()
ct_relate_all <- list()
ct_relate_anc <- list()
ct_relate_der <- list()



##rescaling branch length of relate trees
tmrcas.relate<-numeric(0)

for(i in 1:n.loci){
	tmrcas.relate[i]<- max(branching.times(relate_trees_list[[i]]))
	}

##assign copy of original relate tree to new rescaled tree
relate.all.rescale=relate_trees_list

for(i in 1:n.loci){
	relate.all.rescale[[i]]$edge.length <- relate.all.rescale[[i]]$edge.length*2/mean(tmrcas.relate)
}

for(i in 1:n.loci){
	anc_tips_rent <- which(rent_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_rent <- which(rent_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	anc_trees_rent[[i]] <- drop.tip(rent_trees_list[[i]], der_tips_rent)
	der_trees_rent[[i]] <- drop.tip(rent_trees_list[[i]], anc_tips_rent)

	anc_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	anc_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], der_tips_ms)
	der_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], anc_tips_ms)

	#########
	anc_tips_relate <- which(relate.all.rescale[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_relate <- which(relate.all.rescale[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	relate.anc.rescale[[i]] <- drop.tip(relate.all.rescale[[i]], der_tips_relate)
	relate.der.rescale[[i]] <- drop.tip(relate.all.rescale[[i]], anc_tips_relate)

	print(paste("round", as.character(i), "of tree processing complete"))
}



#########

times.c <- list()
lins.list <- list()
times.c.rent <- list()
lins.list.rent <- list()
times.c.rent.adj <- list()
lins.list.rent.adj <- list()
times.c.relate <- list()
lins.list.relate <- list()

for(i in 1:n.loci){
	times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = FALSE, units_in = 4)
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)

	times.c.rent[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = FALSE, units_in = 2)
	lins.list.rent[[i]] <- times_to_lins(times.c.rent[[i]], time)

	times.c.rent.adj[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = FALSE, units_in = 4)
	lins.list.rent.adj[[i]] <- times_to_lins(times.c.rent.adj[[i]], time)

	times.c.relate[[i]] <- trees_to_times(relate.all.rescale[[i]], relate.anc.rescale[[i]], relate.der.rescale[[i]], time, sure.alt.is.derived = FALSE, units_in = 2)
	lins.list.relate[[i]] <- times_to_lins(times.c.relate[[i]], time)

}



fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")
save.image(paste(out_dir, fn_str, ".RData", sep = ""))




times.c <- list()
lins.list <- list()
times.c.rent <- list()
lins.list.rent <- list()
times.c.rent.adj <- list()
lins.list.rent.adj <- list()
times.c.relate <- list()
lins.list.relate <- list()

#Infer coalescence times from tree data

for(i in 1:n.loci){
	times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = TRUE, units_in = 4)
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)

	times.c.rent[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = TRUE, units_in = 2)
	lins.list.rent[[i]] <- times_to_lins(times.c.rent[[i]], time)

	times.c.rent.adj[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = TRUE, units_in = 4)
	lins.list.rent.adj[[i]] <- times_to_lins(times.c.rent.adj[[i]], time)

	times.c.relate[[i]] <- trees_to_times(relate.all.rescale[[i]], relate.anc.rescale[[i]], relate.der.rescale[[i]], time, sure.alt.is.derived = TRUE, units_in = 2)
	lins.list.relate[[i]] <- times_to_lins(times.c.relate[[i]], time)

}


fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num),  "_alt", sep = "")
save.image(paste(out_dir, fn_str, ".RData", sep = ""))

