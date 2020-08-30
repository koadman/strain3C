data {
    int V; // number of sites (nodes in graph)

    int L_hic; // number of links among sites (edges in graph)
    int hic_linkcounts[L_hic,4]; // counts of SNV co-observation. Genotypes are binary variables: mutant or not. Four possible combinations among two sites: 00,01,10,11
    int hic_linksites[L_hic,2]; // the indices of the two sites involved in each link

    int L_local; // number of links among sites (edges in graph)
    int local_linkcounts[L_local,4]; // counts of SNV co-observation. Genotypes are binary variables: mutant or not. Four possible combinations among two sites: 00,01,10,11
    int local_linksites[L_local,2]; // the indices of the two sites involved in each link

    int K; // number of components to factorize
    real<lower=0> abundance_prior;
    real subsample;
}

parameters {
    real<lower=0,upper=1> genotype[K,V]; // the genotype of the K strains. Ideally this would be a binary (or base 4) variable
    simplex[K] abundance;  // the abundance of each strain
}

model {
    for(k in 1:K){
        genotype[k] ~ beta(0.1,0.1);
    }
    abundance ~ exponential(abundance_prior);

    for( i in 1:L_hic ){
        vector[4] theta;
        theta = rep_vector(0,4);
        for( k in 1:K ){
            vector[4] t_k;
            // probability that strain k has genotype 0,0
            t_k[1] = (1-genotype[k,hic_linksites[i,1]]) * (1-genotype[k,hic_linksites[i,2]]);
            // same for genotype 0,1; then 1,0; then 1,1
            t_k[2] = (1-genotype[k,hic_linksites[i,1]]) * genotype[k,hic_linksites[i,2]];
            t_k[3] = genotype[k,hic_linksites[i,1]] * (1-genotype[k,hic_linksites[i,2]]);
            t_k[4] = genotype[k,hic_linksites[i,1]] * genotype[k,hic_linksites[i,2]];
            theta += abundance[k] * t_k; // weight the probabilities by strain abundance
        }
        target += multinomial_lpmf(hic_linkcounts[i] | theta / sum(theta));
    }

    for( i in 1:L_local ){
        vector[4] theta;
        theta = rep_vector(0,4);
        for( k in 1:K ){
            vector[4] t_k;
            // probability that strain k has genotype 0,0
            t_k[1] = (1-genotype[k,local_linksites[i,1]]) * (1-genotype[k,local_linksites[i,2]]);
            // same for genotype 0,1; then 1,0; then 1,1
            t_k[2] = (1-genotype[k,local_linksites[i,1]]) * genotype[k,local_linksites[i,2]];
            t_k[3] = genotype[k,local_linksites[i,1]] * (1-genotype[k,local_linksites[i,2]]);
            t_k[4] = genotype[k,local_linksites[i,1]] * genotype[k,local_linksites[i,2]];
            theta += abundance[k] * t_k; // weight the probabilities by strain abundance
        }
        target += multinomial_lpmf(local_linkcounts[i] | theta / sum(theta));
    }
}

generated quantities {
	vector[L_hic+L_local] log_lik;

    for( i in 1:L_hic ){
        vector[4] theta;
        theta = rep_vector(0,4);
        for( k in 1:K ){
            vector[4] t_k;
            // probability that strain k has genotype 0,0
            t_k[1] = (1-genotype[k,hic_linksites[i,1]]) * (1-genotype[k,hic_linksites[i,2]]);
            // same for genotype 0,1; then 1,0; then 1,1
            t_k[2] = (1-genotype[k,hic_linksites[i,1]]) * genotype[k,hic_linksites[i,2]];
            t_k[3] = genotype[k,hic_linksites[i,1]] * (1-genotype[k,hic_linksites[i,2]]);
            t_k[4] = genotype[k,hic_linksites[i,1]] * genotype[k,hic_linksites[i,2]];
            theta += abundance[k] * t_k; // weight the probabilities by strain abundance
        }
        log_lik[i] = multinomial_lpmf(hic_linkcounts[i] | theta / sum(theta));
    }

    for( i in 1:L_local ){
        vector[4] theta;
        theta = rep_vector(0,4);
        for( k in 1:K ){
            vector[4] t_k;
            // probability that strain k has genotype 0,0
            t_k[1] = (1-genotype[k,local_linksites[i,1]]) * (1-genotype[k,local_linksites[i,2]]);
            // same for genotype 0,1; then 1,0; then 1,1
            t_k[2] = (1-genotype[k,local_linksites[i,1]]) * genotype[k,local_linksites[i,2]];
            t_k[3] = genotype[k,local_linksites[i,1]] * (1-genotype[k,local_linksites[i,2]]);
            t_k[4] = genotype[k,local_linksites[i,1]] * genotype[k,local_linksites[i,2]];
            theta += abundance[k] * t_k; // weight the probabilities by strain abundance
        }
        log_lik[L_hic+i] = multinomial_lpmf(local_linkcounts[i] | theta / sum(theta));
    }
}
