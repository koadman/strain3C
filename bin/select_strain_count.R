#!/usr/bin/env Rscript
#
# reads in a collection of strain count model outputs and estimates the number of strains
#
library(data.table)
library(loo)

get_cols<-function(post_sample_file, variable_prefix){
	post_df<-fread(file=post_sample_file, sep=",", skip=28, header=T)
	post_cols<-sapply(names(post_df),startsWith, variable_prefix)
	post_col_min <- min(which(post_cols==TRUE))
	post_col_max <- max(which(post_cols==TRUE))
	post_df[,post_col_min:post_col_max]
}

get_loo<-function(log_lik_csv_file){
	ll_cols <- get_cols(log_lik_csv_file, 'log_lik')
	log_liks <- as.matrix(ll_cols)
	lll<-array(log_liks,dim=c(nrow(log_liks),1,ncol(log_liks)))
	r_effs <- relative_eff(exp(lll), chain_id=1)
	loos <- loo(lll, r_eff = r_effs, cores = 2)
	loos
}

args <- commandArgs(TRUE)
post_fnames <- args[1:length(args)-1]
min_abundance <- args[length(args)]

loos<-lapply(post_fnames, get_loo)
names(loos)<-post_fnames
model_comp <- loo_compare(loos)
model_comp
best_model = which(row.names(model_comp)[1]==post_fnames)

abund_cols<-get_cols(args[best_model], 'abund')
abunds <- apply(abund_cols,2,sum) / nrow(abund_cols)
write.table(sum(abunds>min_abundance), file = "strain_count.txt", sep = "", row.names=F, col.names=F)

