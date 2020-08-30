#!/usr/bin/env Rscript
#
# reads in a collection of pointwise log likelihood matrices and uses PSIS-LOO to select the best fitting model
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

loos<-lapply(commandArgs(TRUE), get_loo)
names(loos)<-commandArgs(TRUE)
model_comp <- loo_compare(loos)
model_comp
best_model <- row.names(model_comp)[1]
write.table(best_model, file = "best_model.txt", sep = "", row.names=F, col.names=F, quote=F)
