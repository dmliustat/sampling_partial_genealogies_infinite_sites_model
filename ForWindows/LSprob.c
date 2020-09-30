#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tables.h"
#include "header.h"
#include "gamma.h"
#include "common.h"


extern long M; /*number of genes*/
extern long L; /*number of loci*/
extern double Mu;
extern double SEQLEN;
extern double P[2][2];
extern double THETA;
//extern long *mut_no, *rec_no; /*number of mutations or recombinations at each locus*/
//extern double *stopTime; /*Store the stopping time for each locus, should be <= tau*/
extern double NumLineage; /*The number of lineages to stop*/
extern double *positions;
//extern double empirical[2];
extern long Ne;
extern double adj_rate;
//extern double gamma[NODE_MAX];
extern double gamma2[2][NODE_MAX];


int Build(long* TYPE, long* type2node, long* ntype, const double RHO, const double THETA, double* logw, tsk_table_collection_t* tables, long* num_nodes, long* num_rec, long* num_mut, long* num_coal, long* num_nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node); /*success: return 0; otherwise, return 1*/

double CSDapprx(long *TYPE, long k, long *newseq, long nbranch, long ntype, long *nanc, const double RHO); /*calculate pi(alpha | H-alpha)*/

void FDupdate(long* TYPE, long* type2node, const double RHO, const double THETA, long k, long* ntype, long* nbranch, long* nanc, double* m, double* logw, tsk_table_collection_t* tables, int *rec_num, int *mut_num, int *coal_num, int* nonrecc, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node);

double rexp(double rate);

double CSDapprx_scaling(long* TYPE, long k, long* newseq, long nbranch, long ntype, long* nanc, const double RHO);





int Build(long* TYPE, long* type2node, long* ntype, const double RHO, const double THETA, double* logw, tsk_table_collection_t* tables, long* num_nodes, long* num_rec, long* num_mut, long* num_coal, long* num_nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node) {
    long *nanc, nbranch = M, *newseq; /*nanc: number of individuals ancestral at each locus */
    long *nanc_copy, *TYPE_tmp, ntype_copy, nlineage_copy, nbranch_copy; 
    long i, j, t, l, count, min, max, s1, s2;
    short int flag = 0; /*flag == 0, continue SIS; flag != 0, the stopping criteria is reached*/
    short int flag_MRCA = 0; /*flag_MRA == 0, the grand MRCA not found, continue loops*/
    //short int flag_H0 = 0; /*flag_H0 == 0, calculate the exponential rate for H_0*/
    double branchProb[TYPE_MAX], logH_tau=0.0, apprx;
    double u, tot, temp, m[2];
    unsigned long seed_copy[2]; /* store initial seed; print out if error */
	//int rec_num = 0, mut_num = 0, coal_num = 0;
	long *shuffle_index, * shuffle_seq, *TYPE_flat;
	double logPAC, *logPAC_1, logPAC_avg, PAC_temp;
	long nshuffle=1000;

    nanc = long_vec_init(L);
    nanc_copy = long_vec_init(L);
    newseq = long_vec_init(L);
	TYPE_tmp = long_vec_init((L + 1) * TYPE_MAX);
	TYPE_flat = long_vec_init((L + 1) * NumLineage);
	shuffle_seq = long_vec_init((L + 1) * NumLineage);
	shuffle_index = long_vec_init(NumLineage);
	logPAC_1 = dou_vec_init(nshuffle);


    for(i=0; i<L; i++) {
        *(nanc+i) = M;
    }
    /*for(i=0;i<L-1;i++) rec_no[i] = 0;*/

    *logw = 0.0; /*log SIS weight */
	*num_nodes = 0;
	*num_rec = 0;
	*num_mut = 0;
	*num_coal = 0;
	*num_nonrec = 0;

    *seed_copy = *seed;
    *(seed_copy+1) = *(seed+1);

    while((!flag) && (!flag_MRCA)){
        //for(i=0; i<L; i++) *(nanc_copy+i) = *(nanc+i);
        tot = 0;
        for(i=0; i<*ntype; i++){
            branchProb[i] = *(TYPE+(L+1)*i+L);
            //temp = (double)nbranch-1.0;
			temp = *(TYPE + (L + 1) * i + L) - 1.0;
            min = -1;
            max = -1; /* extreme ancestral sites*/
            for(j=0; j<L; j++){

				if ((*(TYPE + (L + 1) * i + j) > 0) && (*(mut_by_site + j) == 1)) {
					max = j;
					if (min == -1) min = j;
					temp += THETA;
				}

                /*
				if(*(TYPE+(L+1)*i+j) >= 0){
                    max = j;
                    if(min == -1) min = j;
                }
				if (*(mut_by_site + j) == 1) {
					temp += THETA;
				}
				*/
				
				/*
				if (*(TYPE + (L + 1) * i + j) >= 0) {
					max = j;
					if (min == -1) min = j;
					temp += THETA;
				}
				*/
            }
            /*turn site to the corresponding position */
            temp += RHO*(positions[max]-(min==0? 0 : positions[min-1]));   
            /*choose a branch with probability proportional to temp*/
            *(branchProb+i) *= temp;
			if (temp < 0) fprintf(stderr, "Probability of choosing branch is negative.\n");
            tot += *(branchProb+i);
        }
        if(tot <= 0){
            fprintf(stdout, "Error of non-positive probabilities for choosing the chromosome.\n");
            fprintf(stdout, "Initial seed %lu \t %lu \n\n", *(seed_copy), *(seed_copy+1));
            return(1);
        }


        /*start sampling a branch (type) out of the current configuration*/
        u = tot*runif();
        i = 0;
        while(u > 0){
            u -= *(branchProb+i);
            i++;
        }
        i--; /*i is the chosen type (chromosome)*/

      
       if(*(branchProb+i)<=0){
            fprintf(stdout, "Error. Non-possitive probability of choosing chromosome \n");
            fprintf(stdout, "Initial seed %lu \t %lu \n\n", *(seed_copy), *(seed_copy+1));
            return(1);
       }

		//fprintf(stderr, "Inter-event time %lf\n", omega);
        if(nbranch == NumLineage){
			//fprintf(stderr, "Stopping time reached. logw = %lf\n", *logw);
			fprintf(stderr, "\nStopping criterion is met.\n");
			flag = 1;
        }
        else{
			/*Given the chosen chromosome i, select the genetic event. Update configuration and logw*/
			FDupdate(TYPE, type2node, RHO, THETA, i, ntype, &nbranch, nanc, m, logw, tables, num_rec, num_mut, num_coal, num_nonrec, mut_by_site, mut_site, mut_state, mut_node);
		}

        /*check nbranch and ntype*/
        tot = 0;
        for(i=0; i<*ntype; i++) tot += *(TYPE+(L+1)*i+L);
		if ((long)tot != nbranch) {
			fprintf(stdout, "Error in the number of branches. nbranch = %ld, sum of # types = %ld\n", nbranch, (long)tot);
			return(1);
		}

        if(*ntype>TYPE_MAX){
            fprintf(stdout, "Number of distinct haplotypes exceeds the maximum\n");
            return(1);
        }

        /*flag_MRCA = 1 if the GMRCA is reached*/
        flag_MRCA = 1;
        for(j=0; j<L; j++){
            flag_MRCA *= (*(nanc+j)==1);
        }
    }
	*num_nodes = nbranch;
	fprintf(stderr, "num of nodes = %ld, recombinations = %ld; mutations = %ld; non-recurrent = %ld; coalescence = %ld\n", nbranch, *num_rec, *num_mut, *num_nonrec, *num_coal);
	//fprintf(stderr, "Before calculate H_tau, the log(weight) = %lf\n", *logw);

    /*calculate p(N_tau) and p(H_tau | N_tau); N_tau = nbranch*/
	//apprx = Normal_pdf(empirical, nbranch);
	//logH_tau += log(apprx);

	/*The product of approximate conditionals (PAC) depends on the order of the conditionals.
	Shuffle the order of sequences and take the average of PACs*/
	for (i = 0; i < NumLineage; i++) shuffle_index[i] = i;
	count = 0;
	for (i = 0; i < *ntype; i++) {
		for (j = 0; j < *(TYPE + (L + 1) * i + L); j++) {
			for (l = 0; l < L; l++) {
				*(TYPE_flat + (L + 1) * count + l) = *(TYPE + (L + 1) * i + l);
			}	
			*(TYPE_flat + (L + 1) * count + L) = 1;
			count += 1;
		}		
	}
	if (count != NumLineage) fprintf(stderr, "\nError shuffling sequences! count=%ld, nbranch=%ld\n", count, nbranch);

	logPAC = 0.0; /*log(PAC) for the first order*/
	
	for (i = 0; i < nshuffle; i++) {
		nbranch_copy = nbranch;
		nlineage_copy = NumLineage;
		PAC_temp = 0.0;

		/*shuffle the sequences*/
		shuffle(shuffle_index, NumLineage);
		for (j = 0; j < NumLineage; j++) {
			for (l = 0; l < L; l++) {
				shuffle_seq[(L+1)*j+l] = *(TYPE_flat + (L+1)*shuffle_index[j]+l);
			}	
			shuffle_seq[(L + 1) * j + L] = 1;
		}

		/*calculate PAC for the shuffled sequences*/
		for (j = NumLineage - 1; j > 0; j--) {
			for (l = 0; l < L; l++) *(nanc_copy + l) = *(nanc + l);
			for (l = 0; l < L; l++) newseq[l] = *(shuffle_seq + (L + 1) * j + l);
			apprx = CSDapprx_scaling(shuffle_seq, j, newseq, nbranch_copy, nlineage_copy, nanc_copy, RHO);
			PAC_temp += apprx;
			nbranch_copy--;
			nlineage_copy--;			
			for (l = 0; l < L; l++) {
				if (*(shuffle_seq + (L + 1) * j) + l >= 0) nanc_copy[l] -= 1;
				if (nanc_copy[l] < 0) fprintf(stderr, "\nError updating nanc_copy!\n");
			}
		}
		if (nlineage_copy != 1) fprintf(stderr, "\nError! Calculating PAC is wrong.\n");
		s1 = 0;
		s2 = 0;
		for (j = 0; j < L; j++) {
			if (*(shuffle_seq + (L + 1) * 0 + j) == 0) s1 += 1;
			if (*(shuffle_seq + (L + 1) * 0 + j) == 1) s2 += 1;
		}
		PAC_temp += (s1 * log(Mu) + s2 * log(1 - Mu));

		if (i == 0) {
			logPAC = PAC_temp;
			logPAC_1[i] = 0.0;
		}
		else {
			logPAC_1[i] = PAC_temp - logPAC;
			if ((logPAC_1[i] < -312) || (logPAC_1[i] > 313)) fprintf(stderr, "\nWarning! Difference between PACs is too large, %lf.\n", logPAC_1[i]);
		}
	}
	logPAC_avg = 0.0;
	for (i = 0; i < nshuffle; i++) logPAC_avg += exp(logPAC_1[i]);
	//fprintf(stderr, "\nThe averaged part = %lf, the log(PAC) for the first order = %lf\n", logPAC_avg/(double)nshuffle, logPAC);
	logPAC_avg = log(logPAC_avg / (double)nshuffle) + logPAC;
	//fprintf(stderr, "\nThe log(averaged PAC) = %lf\n", logPAC_avg);
	*logw += logPAC_avg;

	

	/*without shuffling*/
 //   ntype_copy = *ntype;
 //   for(i=0; i<(L+1)*ntype_copy; i++) TYPE_tmp[i] = *(TYPE+i);

	//for (i = *ntype - 1; i > 0; i--) {
	//	for (j = 0; j < *(TYPE + (L + 1) * i + L); j++) {
	//		for(t = 0; t < L; t++) newseq[t] = *(TYPE + (L + 1) * i + t);
	//		//apprx = CSDapprx(TYPE_tmp, i, newseq, nbranch, ntype_copy, nanc, RHO);
	//		//logH_tau += log(apprx);
	//		apprx = CSDapprx_scaling(TYPE_tmp, i, newseq, nbranch, ntype_copy, nanc, RHO);
	//		logH_tau += apprx;
	//		reduce_seq(TYPE_tmp, &ntype_copy, i);
	//		nbranch--;
	//		for (t = 0; t < L; t++) {
	//			if (*(TYPE + (L + 1) * i + t) >= 0) nanc[t] -= 1;
	//			if (nanc[t] < 0) fprintf(stderr, "\nError updating nanc!!\n");
	//		}

	//		if (nbranch == 1) break;
	//	}
	//}
	//if (nbranch > 1) {
	//	for (t = 0; t < L; t++) newseq[t] = *(TYPE + (L + 1) * 0 + t);
	//	while (nbranch > 1) {
	//		//apprx = CSDapprx(TYPE_tmp, i, newseq, nbranch, ntype_copy, nanc, RHO);
	//		//logH_tau += log(apprx);
	//		apprx = CSDapprx_scaling(TYPE_tmp, i, newseq, nbranch, ntype_copy, nanc, RHO);
	//		logH_tau += apprx;
	//		reduce_seq(TYPE_tmp, &ntype_copy, i);
	//		nbranch--;
	//		for (t = 0; t < L; t++) {
	//			if (*(TYPE + (L + 1) * i + t) >= 0) nanc[t] -= 1;
	//			if (nanc[t] < 0) fprintf(stderr, "Error updating nanc!!\n");
	//		}
	//	}
	//	s1 = 0;
	//	s2 = 0;
	//	for (j = 0; j < L; j++) {
	//		//fprintf(stderr, "%ld", *(TYPE_tmp + (L + 1) * 0 + j));
	//		if (*(TYPE_tmp + (L + 1) * 0 + j) == 0) s1 += 1;
	//		if (*(TYPE_tmp + (L + 1) * 0 + j) == 1) s2 += 1;
	//	}
	//}
	//else if (nbranch == 1) {
	//	s1 = 0;
	//	s2 = 0;
	//	for (j = 0; j < L; j++) {
	//		if(nanc[j] > 1) fprintf(stderr, "Error updating nanc!!!\n");
	//		//fprintf(stderr, "%ld", *(TYPE_tmp + (L + 1) * 0 + j));
	//		if (*(TYPE_tmp + (L + 1) * 0 + j) == 0) s1 += 1;
	//		if (*(TYPE_tmp + (L + 1) * 0 + j) == 1) s2 += 1;
	//	}
	//}
	//
	///*s1*log(Mu) + s2*log(1-Mu) is for the log of p(h_1), where Mu denotes the stationary distribution of allele 0; 
	//s1 denotes the number of loci with allele 0, and s2 allele 1 on haplotype h_1. */
 //   *logw += (logH_tau + s1*log(Mu) + s2*log(1-Mu));
	////* logw += (logr_tau - logr_0 + logH_tau);
	////fprintf(stderr, "logw=%lf\n", *logw);

    free(nanc);
    free(nanc_copy);
	free(TYPE_tmp);
	free(newseq);
	free(shuffle_index);
	free(shuffle_seq);
    
    return(0);
}





/*returns the Li & Stephens conditional approximation p(h_{t+1} | h_1, ..., h_t)*/
double CSDapprx(long *TYPE, long k, long *newseq, long nbranch, long ntype, long* nanc, const double RHO)
{
    long i, j, *index, index_len=0, d, * nanc_copy, a;
	int e;
    double pj, alpha_sum, *prop, *alpha, *newalpha; 

    index = long_vec_init(L);
	alpha = dou_vec_init(ntype);
	nanc_copy = long_vec_init(L);

	/*copy nanc and remove individual k from the sample*/
	for (i = 0; i < L; i++) {
		nanc_copy[i] = nanc[i];
		if (*(TYPE + (L + 1) * k + i) >= 0) nanc_copy[i] -= 1; 
	}
	//if (*(TYPE + k * (L + 1) + L) == 1) ntype--;

    /*remove one of type k*/
    nbranch--;
    *(TYPE+k*(L+1)+L) -= 1;

    /*set up index-loci ancestral in new*/
    for(i=0; i<L; i++){
        if(*(newseq +i) == 0 || *(newseq + i) == 1){
            *(index + index_len) = i;
            index_len++;
        }
    }

    if(index_len == 0){
        /*if new is non-ancestral on the whole chromosome*/
        *(TYPE + k*(L+1) + L) += 1;
        free(index);
		free(alpha);
		free(nanc_copy);
        return(1.0);
    }

    /*calculate prop and gamma */
    /*prop is the proportion of a chromosome in the current configuration
    that is of type 0 at each site. It's a vector of length index_len */
    prop = dou_vec_init(index_len);
    calcProp(prop, TYPE, ntype, index, index_len);

    if(index_len == 1){
        *(TYPE + k*(L+1) + L) += 1; 
        d = newseq[index[0]];
		for (i = 0; i < ntype; i++) {
			e = *(TYPE + i * (L + 1) + index[0]);
			if (e >= 0) {
				if (d == e) {
					//alpha[i] = gamma[nbranch] / (double)nbranch;
					alpha[i] = gamma2[d][nbranch] / (double)nbranch;
				}
				else {
					//alpha[i] = (1 - gamma[nbranch]) / (double)nbranch;
					alpha[i] = (1 - gamma2[d][nbranch]) / (double)nbranch;
				}
			}
			else {
				if (d == 0) {
					//alpha[i] = (prop[0] * gamma[*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					alpha[i] = (prop[0] * gamma2[d][*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = (prop[0] * gamma2[d][nbranch] + (1.0 - prop[0]) * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
				else {
					//alpha[i] = ((1.0 - prop[0]) * gamma[*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					alpha[i] = ((1.0 - prop[0]) * gamma2[d][*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = ((1.0 - prop[0]) * gamma2[d][nbranch] + prop[0] * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
			}
		}

        alpha_sum = 0.0;
		for (i = 0; i < ntype; i++) {
			alpha_sum += *(TYPE + (L + 1) * i + L) * alpha[i];
		}

        free(index);
        free(prop);
		free(alpha);
		free(nanc_copy);

        return(alpha_sum);
    }

    /*calculate alpha_1(x), for x = 1,...,k. alpha is of length ntype */
	newalpha = dou_vec_init(ntype);
    d = newseq[index[0]];
	if (*(nanc_copy + index[0]) >= NODE_MAX) a = NODE_MAX - 1;
	else a = *(nanc_copy + index[0]);
    for(i=0; i<ntype; i++){
		if (*(TYPE + i * (L + 1) + L) > 0) {
			e = *(TYPE + i * (L + 1) + index[0]);
			if (e >= 0) {
				if (d == e) {
					//alpha[i] = gamma[nbranch] / (double)nbranch;
					//alpha[i] = gamma2[d][nbranch] / (double)nbranch;
					alpha[i] = gamma2[d][a] / (double) * (nanc_copy + index[0]);
				}
				else {
					//alpha[i] = (1 - gamma[nbranch]) / (double)nbranch;
					//alpha[i] = (1 - gamma2[d][nbranch]) / (double)nbranch;
					alpha[i] = (1 - gamma2[d][a]) / (double) * (nanc_copy + index[0]);
				}
			}
			else {
				if (d == 0) {
					//alpha[i] = (prop[0] * gamma[*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = (prop[0] * gamma[*(nanc_copy + index[0])] + (1 - prop[0]) * (1 - gamma[*(nanc_copy + index[0])])) / *(nanc_copy + index[0]);
					//alpha[i] = (prop[0] * gamma2[d][*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					alpha[i] = (prop[0] * gamma2[d][a] + (1.0 - prop[0]) * (1.0 - gamma2[d][a])) / (double) * (nanc_copy + index[0]);
					//alpha[i] = (prop[0] * gamma2[d][nbranch] + (1.0 - prop[0]) * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
				else {
					//alpha[i] = ((1.0 - prop[0]) * gamma[*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = ((1 - prop[0]) * gamma[*(nanc_copy + index[0])] + prop[0] * (1 - gamma[*(nanc_copy + index[0])])) / *(nanc_copy + index[0]);
					//alpha[i] = ((1.0 - prop[0]) * gamma2[d][*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					alpha[i] = ((1.0 - prop[0]) * gamma2[d][a] + prop[0] * (1.0 - gamma2[d][a])) / (double) * (nanc_copy + index[0]);
					//alpha[i] = ((1.0 - prop[0]) * gamma2[d][nbranch] + prop[0] * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
			}
			/*if (e >= 0) {
				if (d == e) {
					alpha[i] = gamma2[d][a]/(double) * (nanc_copy + index[0]);
				}
				else {
					if (e == 1) alpha[i] = (1.0 - gamma2[d][a])/ (double)*(nanc_copy + index[0]);
					else alpha[i] = 0;
				}
			}*/

		}
    }
    alpha_sum = 0.0;
	for (i = 0; i < ntype; i++) {
		if (*(TYPE + i * (L + 1) + L) >= 0) {
			alpha_sum += *(TYPE + (L + 1) * i + L) * alpha[i];
		}
	}

    /*calculate alpha_j, for j=2,...,S */
    for(j = 1; j < index_len; j++){
        /*calculate recombination rate, p_j, in (A5) */
        pj = exp(-1.0 * (positions[*(index+j)] - positions[index[j-1]]) * adj_rate * RHO / (double)nbranch);
		if (*(nanc_copy + index[j]) >= NODE_MAX) a = NODE_MAX - 1;
		else a = *(nanc_copy + index[j]);

        d = newseq[index[j]];
		if (d != 0 && d != 1) fprintf(stderr, "newseq error! \n");
        for(i = 0; i < ntype; i++){
			if (*(TYPE + i * (L + 1) + L) > 0) {
				e = *(TYPE + i * (L + 1) + index[j]);

				if (e >= 0) {
					if (d == e) {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * gamma[nbranch];
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * gamma2[d][nbranch];
						newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * gamma2[d][a];

					}
					else {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (1.0 - gamma[nbranch]);
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (1.0 - gamma2[d][nbranch]);
						newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (1.0 - gamma2[d][a]);
					}
				}
				else {
					if (d == 0) {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (prop[j] * gamma[*(nanc_copy + index[j])] + (1.0 - prop[j]) * (1.0 - gamma[*(nanc_copy + index[j])]));
						newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (prop[j] * gamma2[d][a] + (1.0 - prop[j]) * (1.0 - gamma2[d][a]));
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (prop[j] * gamma2[d][nbranch] + (1.0 - prop[j]) * (1.0 - gamma2[d][nbranch]));

					}
					else {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * ((1.0 - prop[j]) * gamma[*(nanc_copy + index[j])] + prop[j] * (1.0 - gamma[*(nanc_copy + index[j])]));
						newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * ((1.0 - prop[j]) * gamma2[d][a] + prop[j] * (1.0 - gamma2[d][a]));
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * ((1.0 - prop[j]) * gamma2[d][nbranch] + prop[j] * (1.0 - gamma2[d][nbranch]));
					}
				}
				/*if (e >= 0) {
					if (d == e) {
						newalpha[i] = (pj * alpha[i] + (1 - pj) / (double)nbranch * alpha_sum) * gamma2[d][a];
					}
					else {
						if (e == 1) newalpha[i] = (pj * alpha[i] + (1 - pj) / (double)nbranch * alpha_sum) * (1.0 - gamma2[d][a]);
						else newalpha[i] = 0;
					}
				}*/

			}

        }
        alpha_sum = 0.0;
		for (i = 0; i < ntype; i++) {
			if (*(TYPE + i * (L + 1) + L) >= 0) {
				alpha_sum += *(TYPE + (L + 1) * i + L) * newalpha[i];
				alpha[i] = newalpha[i];
			}
		}

    }

    *(TYPE + k*(L+1) + L) += 1;

    free(index);
    free(prop);
	free(alpha);
	free(newalpha);
	free(nanc_copy);

    return(alpha_sum);
}




void FDupdate(long* TYPE, long* type2node, const double RHO, const double THETA, long k, long* ntype, long* nbranch, long* nanc, double* m, double* logw, tsk_table_collection_t* tables, int *rec_num, int *mut_num, int* coal_num, int* nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node)
{
    double tot, probtemp, probtemp2, *prob, *lweights, u, p_common; /*p_common=p(alpha | H-alpha) */
    long i, j, a, e, *newseq, *old, *rec1, *rec2, min, max, *TYPE_st;
    long onode, onode2, nnode, nnode2, edge, *edge_arr, nodeTime;
    double left, right; /* edge_sum = 0.0, rec_sum = 0.0 for edges*/
	int flag, flag_first = 0, flag_last = 0;
	char state[2];

    newseq = long_vec_init(L);
    old = long_vec_init(L);
    rec1 = long_vec_init(L);
    rec2 = long_vec_init(L);
	edge_arr = long_vec_init(L);
	TYPE_st = long_vec_init((L + 1) * TYPE_MAX);

    prob = dou_vec_init(K*L+L-1+TYPE_MAX);
    lweights = dou_vec_init(K*L+L-1+TYPE_MAX);

    /*onode is the index of the chosen/old node sampled from type k nodes*/
    onode = sample_node(TYPE, type2node, k);

    for(i=0; i<L; i++) {
        newseq[i] = *(TYPE+(L+1)*k+i);
        rec1[i] = -1;
        rec2[i] = newseq[i];
    }

	for (e = 0; e < (L + 1) * *ntype; e++) TYPE_st[e] = *(TYPE + e);

     /*calculate H(alpha | H-alpha)*/
    p_common = CSDapprx(TYPE, k, newseq, *nbranch, *ntype, nanc, RHO);


    /*mutation events, K*L in total (0 -- K*L-1)*/
    for(i=0; i<L; i++){
        /*on ancestral loci */
        if((*(TYPE+(L+1)*k+i) > 0) && (*(mut_by_site+i) == 1)){
            a = *(TYPE+(L+1)*k+i); /*a is the current allele of the chosen chromosome*/
            for(j=0; j<K; j++){
                if(P[j][a]>0){
					if (j == a) {
						prob[K * i + j] = THETA * p_common * P[j][a];
 						lweights[K * i + j] = 0;
					}
					else {
						newseq[i] = j; /*the new allele*/
						probtemp = CSDapprx(TYPE, k, newseq, *nbranch, *ntype, nanc, RHO);
						newseq[i] = a;
						
						lweights[K * i + j] = log(p_common) - log(probtemp);
						prob[K * i + j] = THETA * probtemp * P[j][a];
					}
                }
                else{
                    prob[K*i+j] = 0;
                    lweights[K*i+j] = -INFINITY;
                }
            }
        } 
        /*on non-ancestral loci */
        else{
            for(j=0; j<K; j++) {
                prob[K*i+j] = 0;
                lweights[K*i+j] = -INFINITY;
            }
        }      
    }

    /*recombination events, L-1 in total (K*L -- K*L+L-2)*/
    min = -1;
    max = -1; /*min and max are the extreme ancestral loci*/
    for(i=0; i<L; i++){
        prob[K*L+i] = 0;
        lweights[K*L+i] = -INFINITY;
        if(*(TYPE+(L+1)*k+i)>=0){
            max = i;
            if(min<0) min = i;
        }
    }
    for(i=min; i<max; i++){
        rec1[i] = rec2[i];
        rec2[i] = -1;
        probtemp = CSDapprx(TYPE, k, rec1, *nbranch, *ntype, nanc, RHO); /*p(beta | H-alpha) */
		/*add rec1 to the current configuration*/
		for (j = 0; j < L; j++) {
			if (rec1[j] >= 0) nanc[j] += 1;
		}
		j = add_s(TYPE, ntype, rec1, *nbranch);
		*nbranch += 1;
        probtemp2 = CSDapprx(TYPE, k, rec2, *nbranch, *ntype, nanc, RHO); /*p(gamma | H+beta-alpha) */
		/*remove rec1*/
		remove_s(TYPE, ntype, j, *nbranch);
		*nbranch -= 1;
		for (j = 0; j < L; j++) {
			if (rec1[j] >= 0) nanc[j] -= 1;
		}

		prob[K * L + i] = RHO * (positions[i + 1] - positions[i]) * probtemp * probtemp2;		
        lweights[K*L+i] = log(p_common) - log(probtemp) - log(probtemp2);
    }


    /*coalescence events, ntype in total (K*L+L-1 -- K*L+L-2+*ntype)*/
    for(i=0; i<*ntype; i++){
        /*coalescence of two different haplotypes*/
        if(i!=k){
            prob[K*L+L-1+i] = 0;
            lweights[K*L+L-1+i] = -INFINITY;
            for(j=0; j<L; j++){
                if(*(TYPE+(L+1)*i+j) == *(TYPE+(L+1)*k+j)){
                    newseq[j] = *(TYPE+(L+1)*i+j);
                }
                else{
                    /*if the two haplotypes are both ancestral but differen at site j, 
                    this event cannot occur*/
                    if(*(TYPE+(L+1)*i+j)>=0 && *(TYPE+(L+1)*k+j)>=0) j = L+1;
                    else{
                        if(*(TYPE+(L+1)*i+j) >= 0) newseq[j] = *(TYPE+(L+1)*i+j);
                        else newseq[j] = *(TYPE+(L+1)*k+j);
                    }                    
                }
                if(j<L) old[j] = *(TYPE+(L+1)*i+j);
            }
            if(j==L){
				//for (e = 0; e < (L + 1) * *ntype; e++) TYPE_st[e] = *(TYPE + e);
				/*remove beta=old from the current configuration*/
				for (e = 0; e < L; e++) {
					if (old[e] >= 0) nanc[e] -= 1;
				}
				e = remove_s(TYPE_st, ntype, i, *nbranch);
				*nbranch -= 1;
				if (e == 1 && i < k) {
					probtemp = CSDapprx(TYPE_st, k-1, newseq, *nbranch, *ntype, nanc, RHO); /*p(gamma | H-alpha-beta)*/
					probtemp2 = CSDapprx(TYPE_st, k-1, old, *nbranch, *ntype, nanc, RHO); /*p(beta | H-alpha-beta)*/
				}
				else {
					probtemp = CSDapprx(TYPE_st, k, newseq, *nbranch, *ntype, nanc, RHO); /*p(gamma | H-alpha-beta)*/
					probtemp2 = CSDapprx(TYPE_st, k, old, *nbranch, *ntype, nanc, RHO); /*p(beta | H-alpha-beta)*/
				}
                //probtemp = CSDapprx(TYPE, k, newseq, *nbranch, *ntype, nanc, RHO); /*p(gamma | H-alpha-beta)*/
                //probtemp2 = CSDapprx(TYPE, k, old, *nbranch, *ntype, nanc, RHO); /*p(beta | H-alpha-beta)*/
				/*add beta back; no need to change the configuration since TYPE is not modified*/
				*nbranch += 1;
				if(e == 1) *ntype += 1;
				for (e = 0; e < L; e++) {
					if (old[e] >= 0) nanc[e] += 1;
				}
				for (e = (L + 1) * i; e < (L + 1) * *ntype; e++) TYPE_st[e] = *(TYPE + e);
               
                lweights[K*L+L-1+i] = log(p_common) + log(probtemp2) - log(probtemp);
				prob[K * L + L - 1 + i] = *(TYPE + (L + 1) * i + L) * probtemp / probtemp2;
            }
        }
        /*a coalescence of two identical haplotypes*/
		/*have to decide if there is more than one type k*/
        else {            
            lweights[K*L+L-1+i] = log(p_common) + log((*(TYPE+(L+1)*i+L)==1) ? 0: 1);
			prob[K * L + L - 1 + i] = *(TYPE + (L + 1) * i + L) - 1.0;
        }
    }


    /*update the samples*/
    tot = 0.0;
	for (e = 0; e < (L * K + L + *ntype - 1); e++) {
		//if(*(prob + e) < 0) fprintf(stderr, "e %ld probability %f\n", e, *(prob + e));
		tot += *(prob + e);
	}

    if(tot <= 0){
        /**logq = 1; */
		//fprintf(stderr, "\ntot prob of events is <= 0!\n");
        *logw = 1;
        return;
    }

    u = runif()*tot;
    e = 0;
    while(u >= 0){
        u -= *(prob + e);
        e++;
    }
    /*e-1 is the chosen event*/
	if((*(prob + e - 1) > 0) && (*(lweights + e - 1) != -INFINITY)){
         *logw += *(lweights + e - 1);
		//fprintf(stderr, "\n log weight = %f\n", *logw);
    }
	
    else{
        fprintf(stdout, "\n e %ld, probability %f, log weight %f, tot %f\n", e, *(prob+e-1), *(lweights+e-1), tot);
        return;
    }
	

    if(e <= K*L){/*mutation*/
        a = (e - 1) % K; /*a is the new allele after mutation*/
        i = (e - a - 1) / K; /*i is the locus*/
		//fprintf(stderr, "\t mut");
		*mut_num += 1;

        m[0] = *(TYPE + (L + 1) * k + L)* P[a][*(TYPE + (L + 1) * k + i)]; /*m[0] is prior rate*/
        m[1] = (double)(L+i); /*m[1] is the type of transition; m[1]>=L is mutation*/

        /*calculate new type*/
        for(j=0; j<L; j++) newseq[j]=*(TYPE + (L + 1) * k + j);
		if (*(TYPE + (L + 1) * k + i) != a) {
			*(mut_by_site + i) -= 1;
			//fprintf(stderr, "_nonrec");
			*nonrec += 1;
			newseq[i] = a;
			remove_old(TYPE, type2node, ntype, k, onode, *nbranch);
			*nbranch -= 1;
			add_new(TYPE, type2node, ntype, newseq, onode, *nbranch);
			*nbranch += 1;
		}       

		mut_site[*mut_num - 1] = i;
		mut_state[*mut_num - 1] = a;
		mut_node[*mut_num - 1] = onode;

		//sprintf(state, "%d", a);
		//mut = tsk_mutation_table_add_row(&tables->mutations, i, onode, -1, state, 0, NULL, 0); 
		//if (mut < 0) fprintf(stderr, "Adding mutation failed.\n");
       
    }

    if(e>K*L && e<=(K*L+L-1)){/*recombination*/
        i = e-K*L;/*i is the breakpoint in [1, L-1]*/
		*rec_num += 1;
		nodeTime = *rec_num + *mut_num + *coal_num;

        m[0] = *(TYPE+(L+1)*k+L);
        m[1] = i;

        /*rec1 stores the alleles before the brakpoint; the alleles after the breakpoint are -1 in rec1*/
        /*rec2 stores the alleles after the brakpoint; the alleles before the breakpoint are -1 in rec2*/
        /*edge_arr: 0 non-ancestral; 1 rec1; 2 rec2*/

        for(j=0; j<i; j++){
            rec1[j] = *(TYPE+(L+1)*k+j);
            rec2[j] = -1;
			flag_first += rec1[j];
			flag_last += rec2[j];
            if(*(TYPE+(L+1)*k+j)<0) edge_arr[j] = 0;
			else {
				edge_arr[j] = 1;
			}
        }
        for(j=i; j<L; j++){
            rec1[j] = -1;
            rec2[j] = *(TYPE+(L+1)*k+j);
			flag_first += rec1[j];
			flag_last += rec2[j];
			if (*(TYPE + (L + 1) * k + j) < 0) edge_arr[j] = 0;
			else {
				edge_arr[j] = 2;
			}
        }
		
		if (flag_first == -L || flag_last == -L) {
			/*one of the new sequences is non-ancestral on all loci*/
			nnode = tsk_node_table_add_row(&tables->nodes, 0, nodeTime, TSK_NULL, TSK_NULL, NULL, 0);
			nnode2 = nnode;
		}
		else {
			nnode = tsk_node_table_add_row(&tables->nodes, 0, nodeTime, TSK_NULL, TSK_NULL, NULL, 0);
			nnode2 = tsk_node_table_add_row(&tables->nodes, 0, nodeTime, TSK_NULL, TSK_NULL, NULL, 0);
		}
		if (nnode < 0 || nnode2 < 0) fprintf(stderr, "Adding node failed.\n");

		/*store rec1 edges*/
		left = -1;
		flag = -1; /*flag=-1 haven't found the interval; flag=0 found left*/
		for (j = 0; j <= i; j++) {
			if (flag == -1 && j < i) {
				if (*(TYPE + (L + 1) * k + j) >= 0) {
					left = positions[j];
					flag = 0;
				}
			}			
			else if (flag == 0) {
				if (j == i) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
				}
				else if (*(TYPE + (L + 1) * k + j) < 0) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
					flag = -1;
				}
			}
		}

		/*store rec2 edges*/
		left = -1;
		flag = -1;
		for (j = i; j <= L; j++) {
			if (flag == -1 && j < L) {
				if (*(TYPE + (L + 1) * k + j) >= 0) {
					left = positions[j];
					flag = 0;
				}
			}
			else if (flag == 0) {
				if (j == L) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode2, onode);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
				}
				else if (*(TYPE + (L + 1) * k + j) < 0) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode2, onode);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
					flag = -1;
				}
			}
		}
		
		if (nnode != nnode2) {
			add_new(TYPE, type2node, ntype, rec1, nnode, *nbranch); /*nbranch is the current number of suquences (before adding or removing the node)*/
			*nbranch += 1;
			add_new(TYPE, type2node, ntype, rec2, nnode2, *nbranch);
			*nbranch += 1;
		}
		else if (flag_last == -1 * L) {
			/*rec2 is non-ancestral on the whole sequence*/
			add_new(TYPE, type2node, ntype, rec1, nnode, *nbranch);
			*nbranch += 1;
		}
		else{
			/*rec1 is non-ancestral on the whole sequence*/
			add_new(TYPE, type2node, ntype, rec2, nnode2, *nbranch);
			*nbranch += 1;
		}
		
        remove_old(TYPE, type2node, ntype, k, onode, *nbranch);
        *nbranch -= 1;
		/*fprintf(stderr, "Current nodes:\t");
		for (i = 0; i < *nbranch; i++) fprintf(stderr, "%ld\t", type2node[i]);
		fprintf(stderr, "\n");*/
    } 

    if(e>K*L+L-1){/*coalescence*/
		//fprintf(stderr, "\t coal");
		*coal_num += 1;
		nodeTime = *rec_num + *mut_num + *coal_num;
        i = e-(K*L+L); /*i is the type that coalesces with k*/
		for (j = 0; j < L; j++) {
			if((*(TYPE + (L + 1) * i + j) > 0) && (*(TYPE + (L + 1) * k + j) > 0)) *(mut_by_site + j) -= 1;
		}

        onode2 = sample_node(TYPE, type2node, i);
		while (onode == onode2) {
			onode2 = sample_node(TYPE, type2node, i);
		}
		if (onode > onode2) {
			onode = onode + onode2;
			onode2 = onode - onode2;
			onode = onode - onode2;
			k = k + i;
			i = k - i;
			k = k - i;
		}
        nnode = tsk_node_table_add_row(&tables->nodes, 0, nodeTime, TSK_NULL, TSK_NULL, NULL, 0);
		//fprintf(stderr, "Coalescence event is chosen. children %ld and %ld, parent %ld\n", onode, onode2, nnode);
		if (nnode < 0) fprintf(stderr, "Adding node failed.\n");
 
        /*edge_arr: 0 non-ancestral; 1 type k; 2 type i; 3: both*/
        for(j=0; j<L; j++){ 
            if(*(TYPE+(L+1)*k+j)<0 && *(TYPE+(L+1)*i+j)>=0){
                newseq[j] = *(TYPE+(L+1)*i+j);
                edge_arr[j] = 2;
            } 
            else if(*(TYPE+(L+1)*k+j)>=0 && *(TYPE+(L+1)*i+j)<0){
                newseq[j] = *(TYPE+(L+1)*k+j);
                edge_arr[j] = 1;
            }
            else if(*(TYPE+(L+1)*k+j)<0 && *(TYPE+(L+1)*i+j)<0){
                newseq[j] = -1;
                edge_arr[j] = 0;
            }
            else{ 
                /* *(TYPE+(L+1)*k+j)>=0 && *(TYPE+(L+1)*i+j)>=0 && they're identical on locus j */
				newseq[j] = *(TYPE+(L+1)*k+j);
                edge_arr[j] = 3;
                nanc[j] -= 1;
            }
        }

        /*add to the edge table*/
		/*within each given parent, the edges must be sorted by child index and then left*/
		/*store nnode1 edges*/
		left = -1;
		flag = -1; /*flag=-1 haven't found the interval; flag=0 found left*/
		for (j = 0; j <= L; j++) {
			if (flag == -1 && j < L) {
				if (edge_arr[j] == 1 || edge_arr[j] == 3) {
					left = positions[j];
					flag = 0;
				}
			}
			else if (flag == 0) {
				if (j == L) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
				}
				else if (edge_arr[j] == 0 || edge_arr[j] == 2) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
					flag = -1;
				}
			}
		}

		/*store nnode2 edges*/
		left = -1;
		flag = -1; /*flag=-1 haven't found the interval; flag=0 found left*/
		for (j = 0; j <= L; j++) {
			if (flag == -1 && j < L) {
				if (edge_arr[j] == 2 || edge_arr[j] == 3) {
					left = positions[j];
					flag = 0;
				}
			}
			else if (flag == 0) {
				if (j == L) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode2);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
				}
				else if (edge_arr[j] == 0 || edge_arr[j] == 1) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode2);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
					flag = -1;
				}
			}
		}

        m[0] = (*(TYPE+(L+1)*k+L) * (*(TYPE+(L+1)*i+L) - (double)(i==k)));
        m[1] = 0.0;
        add_new(TYPE, type2node, ntype, newseq, nnode, *nbranch);
        *nbranch += 1;
        if(*(TYPE+(L+1)*k+L)==1 && k<i) i--;
        remove_old(TYPE, type2node, ntype, k, onode, *nbranch);
        *nbranch -= 1;
        remove_old(TYPE, type2node, ntype, i, onode2, *nbranch);
        *nbranch -= 1;

    }

    free(newseq);
    free(old);
    free(rec1);
    free(rec2);
    free(prob);
    free(lweights);
	free(edge_arr);
	free(TYPE_st);
}



double rexp(double rate){
    double r;
    r=runif();
    return (-log(r)/rate);
}



double CSDapprx_scaling(long* TYPE, long k, long* newseq, long nbranch, long ntype, long* nanc, const double RHO)
{
	long i, j, * index, index_len = 0, d, * nanc_copy, a;
	int e;
	double pj, alpha_sum, * prop, * alpha_prime, * alpha_hat, logwt_scale; 

	index = long_vec_init(L);
	
	/*remove one of type k*/
	nbranch--;
	*(TYPE + k * (L + 1) + L) -= 1;

	/*set up index-loci ancestral in new*/
	for (i = 0; i < L; i++) {
		if (*(newseq + i) >= 0) {
			*(index + index_len) = i;
			index_len++;
		}
	}

	if (index_len == 0) {
		/*if new is non-ancestral on the whole chromosome*/
		*(TYPE + k * (L + 1) + L) += 1;
		free(index);
		return(0.0);
	}

	alpha_prime = dou_vec_init(ntype);
	nanc_copy = long_vec_init(L);

	/*copy nanc and remove individual k from the sample*/
	for (i = 0; i < L; i++) {
		nanc_copy[i] = nanc[i];
		if (*(TYPE + (L + 1) * k + i) >= 0) nanc_copy[i] -= 1;
	}
	//if (*(TYPE + k * (L + 1) + L) == 1) ntype--;


	/*calculate prop and gamma */
	/*prop is the proportion of a chromosome in the current configuration
	that is of type 0 at each site. It's a vector of length index_len */
	prop = dou_vec_init(index_len);
	calcProp(prop, TYPE, ntype, index, index_len);

	if (index_len == 1) {
		*(TYPE + k * (L + 1) + L) += 1;
		d = newseq[index[0]];
		for (i = 0; i < ntype; i++) {
			e = *(TYPE + i * (L + 1) + index[0]);
			if (e >= 0) {
				if (d == e) {
					//alpha[i] = gamma[nbranch] / (double)nbranch;
					alpha_prime[i] = gamma2[d][nbranch] / (double)nbranch;
				}
				else {
					//alpha[i] = (1 - gamma[nbranch]) / (double)nbranch;
					alpha_prime[i] = (1 - gamma2[d][nbranch]) / (double)nbranch;
				}
			}
			else {
				if (d == 0) {
					//alpha[i] = (prop[0] * gamma[*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					alpha_prime[i] = (prop[0] * gamma2[d][*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = (prop[0] * gamma2[d][nbranch] + (1.0 - prop[0]) * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
				else {
					//alpha[i] = ((1.0 - prop[0]) * gamma[*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					alpha_prime[i] = ((1.0 - prop[0]) * gamma2[d][*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = ((1.0 - prop[0]) * gamma2[d][nbranch] + prop[0] * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
			}
		}

		alpha_sum = 0.0;
		for (i = 0; i < ntype; i++) {
			alpha_sum += *(TYPE + (L + 1) * i + L) * alpha_prime[i];
		}
		logwt_scale = log(alpha_sum);

		free(index);
		free(prop);
		free(alpha_prime);
		free(nanc_copy);

		return(logwt_scale);
	}

	alpha_hat = dou_vec_init(ntype);

	/*calculate alpha_1(x), for x = 1,...,k. alpha is of length ntype */
	d = newseq[index[0]];
	if (*(nanc_copy + index[0]) >= NODE_MAX) a = NODE_MAX - 1;
	else a = *(nanc_copy + index[0]);
	for (i = 0; i < ntype; i++) {
		if (*(TYPE + i * (L + 1) + L) > 0) {
			e = *(TYPE + i * (L + 1) + index[0]);
			if (e >= 0) {
				if (d == e) {
					//alpha[i] = gamma[nbranch] / (double)nbranch;
					//alpha[i] = gamma2[d][nbranch] / (double)nbranch;
					alpha_prime[i] = gamma2[d][a] / (double) * (nanc_copy + index[0]);
				}
				else {
					//alpha[i] = (1 - gamma[nbranch]) / (double)nbranch;
					//alpha[i] = (1 - gamma2[d][nbranch]) / (double)nbranch;
					alpha_prime[i] = (1 - gamma2[d][a]) / (double) * (nanc_copy + index[0]);
				}
			}
			else {
				if (d == 0) {
					//alpha[i] = (prop[0] * gamma[*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = (prop[0] * gamma[*(nanc_copy + index[0])] + (1 - prop[0]) * (1 - gamma[*(nanc_copy + index[0])])) / *(nanc_copy + index[0]);
					//alpha[i] = (prop[0] * gamma2[d][*(nanc_copy + index[0])] + (1.0 - prop[0]) * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					alpha_prime[i] = (prop[0] * gamma2[d][a] + (1.0 - prop[0]) * (1.0 - gamma2[d][a])) / (double) * (nanc_copy + index[0]);
					//alpha[i] = (prop[0] * gamma2[d][nbranch] + (1.0 - prop[0]) * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
				else {
					//alpha[i] = ((1.0 - prop[0]) * gamma[*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma[*(nanc_copy + index[0])])) / (double)nbranch;
					//alpha[i] = ((1 - prop[0]) * gamma[*(nanc_copy + index[0])] + prop[0] * (1 - gamma[*(nanc_copy + index[0])])) / *(nanc_copy + index[0]);
					//alpha[i] = ((1.0 - prop[0]) * gamma2[d][*(nanc_copy + index[0])] + prop[0] * (1.0 - gamma2[d][*(nanc_copy + index[0])])) / (double)nbranch;
					alpha_prime[i] = ((1.0 - prop[0]) * gamma2[d][a] + prop[0] * (1.0 - gamma2[d][a])) / (double) * (nanc_copy + index[0]);
					//alpha[i] = ((1.0 - prop[0]) * gamma2[d][nbranch] + prop[0] * (1.0 - gamma2[d][nbranch])) / (double)nbranch;
				}
			}
			
		}
	}
	alpha_sum = 0.0;
	for (i = 0; i < ntype; i++) {
		if (*(TYPE + i * (L + 1) + L) > 0) {
			alpha_sum += *(TYPE + (L + 1) * i + L) * alpha_prime[i];
		}
	}
	for (i = 0; i < ntype; i++) {
		alpha_hat[i] = alpha_prime[i] / alpha_sum;
	}
	logwt_scale = log(alpha_sum);

	/*calculate alpha_j, for j=2,...,S */
	for (j = 1; j < index_len; j++) {
		/*calculate recombination rate, p_j, in (A5) */
		pj = exp(-1.0 * (positions[*(index + j)] - positions[index[j - 1]]) * adj_rate * RHO / (double)nbranch);
		if (*(nanc_copy + index[j]) >= NODE_MAX) a = NODE_MAX - 1;
		else a = *(nanc_copy + index[j]);

		d = newseq[index[j]];
		if (d != 0 && d != 1) fprintf(stderr, "newseq error!\n");
		for (i = 0; i < ntype; i++) {
			if (*(TYPE + i * (L + 1) + L) > 0) {
				e = *(TYPE + i * (L + 1) + index[j]);

				if (e >= 0) {
					if (d == e) {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * gamma[nbranch];
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * gamma2[d][nbranch];
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * gamma2[d][a];
					
					}
					else {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (1.0 - gamma[nbranch]);
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (1.0 - gamma2[d][nbranch]);
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * (1.0 - gamma2[d][a]);
						
					}
				}
				else {
					if (d == 0) {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (prop[j] * gamma[*(nanc_copy + index[j])] + (1.0 - prop[j]) * (1.0 - gamma[*(nanc_copy + index[j])]));
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * (prop[j] * gamma2[d][a] + (1.0 - prop[j]) * (1.0 - gamma2[d][a]));
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * (prop[j] * gamma2[d][nbranch] + (1.0 - prop[j]) * (1.0 - gamma2[d][nbranch]));
						

					}
					else {
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * ((1.0 - prop[j]) * gamma[*(nanc_copy + index[j])] + prop[j] * (1.0 - gamma[*(nanc_copy + index[j])]));
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * ((1.0 - prop[j]) * gamma2[d][a] + prop[j] * (1.0 - gamma2[d][a]));
						//newalpha[i] = (pj * alpha[i] + (1.0 - pj) / (double)nbranch * alpha_sum) * ((1.0 - prop[j]) * gamma2[d][nbranch] + prop[j] * (1.0 - gamma2[d][nbranch]));
						
					}
				}

			}
			//if (newalpha[i] <= 0)
				//fprintf(stderr, "Error! alpha[%ld] <= 0!\n", i);
		}
		alpha_sum = 0.0;
		for (i = 0; i < ntype; i++) {
			if (*(TYPE + i * (L + 1) + L) > 0) {
				alpha_sum += *(TYPE + (L + 1) * i + L) * alpha_prime[i];
			}
		}
		for (i = 0; i < ntype; i++) {
			alpha_hat[i] = alpha_prime[i] / alpha_sum;
		}
		logwt_scale += log(alpha_sum);
	}

	*(TYPE + k * (L + 1) + L) += 1;

	free(index);
	free(prop);
	free(alpha_prime);
	free(alpha_hat);
	free(nanc_copy);

	return(logwt_scale);
}












