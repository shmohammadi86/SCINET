#include <scinet.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_gene_activities(mat A, uvec samples, bool consider_baseline_expression = true, int thread_no = 4) {
	int K = samples.n_elem;
	for(int k = 0; k < K; k++)
		samples[k]--;
		
	mat Z = SCINET::compute_gene_activities(A, samples, consider_baseline_expression, thread_no);
	
	return Z;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_gene_activities_decoupled(mat archetypes, mat H, uvec samples, bool consider_baseline_expression = true, int thread_no = 4) {
	int K = samples.n_elem;
	for(int k = 0; k < K; k++)
		samples[k]--;
		
	mat Z = SCINET::compute_gene_activities_decoupled(archetypes, H, samples, consider_baseline_expression, thread_no);
	
	return Z;		
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_gene_activities_full(mat A, int thread_no = 4) {
	mat Z = SCINET::compute_gene_activities_full(A, thread_no);
	
	return Z;		
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat RIN_transform(mat A, int thread_no = 4) {
	mat Z = SCINET::RIN_transform(A, thread_no);
	
	return Z;		
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat subsample_gene_activities(mat A, uvec rows, uvec samples, int total_subsamples = 30, int cells_per_subsample = 10, int thread_no = 4, int seed = 0) {		
	int K = samples.n_elem;
	int M = rows.n_elem;
	
	for(int k = 0; k < K; k++)
		samples[k]--;
		
	for(int k = 0; k < M; k++)
		rows[k]--;
	
	
	mat subsampled_gene_activities = SCINET::subsample_gene_activities(A, rows, samples, total_subsamples, cells_per_subsample, thread_no, seed);

    	
	return(subsampled_gene_activities);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat subsample_gene_activities_decoupled(mat archetypes, mat H, uvec rows, uvec samples, int total_subsamples, int cells_per_subsample, int thread_no, int seed) {
	int K = samples.n_elem;
	int M = rows.n_elem;
	
	for(int k = 0; k < K; k++)
		samples[k]--;
		
	for(int k = 0; k < M; k++)
		rows[k]--;
	
	
	mat subsampled_gene_activities = SCINET::subsample_gene_activities_decoupled(archetypes, H, rows, samples, total_subsamples, cells_per_subsample, thread_no, seed);

    	
	return(subsampled_gene_activities);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
field<sp_mat> construct_cell_networks(sp_mat net, mat gene_activities, int thread_no) {

	field<sp_mat> res = SCINET::construct_cell_networks(net, gene_activities, thread_no);
    	
	return(res);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
field<mat> construct_cell_networks_noPrior(mat gene_activities, int thread_no) {

	field<mat> res = SCINET::construct_cell_networks_noPrior(gene_activities, thread_no);
    	
	return(res);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
field<sp_mat> construct_cell_networks_summary(sp_mat net, mat gene_activities, int total_subsamples = 30, int cells_per_subsample = 10, int seed = 0, int thread_no = 4) {

	field<sp_mat> res = SCINET::construct_cell_networks_summary(net, gene_activities, total_subsamples, cells_per_subsample, seed, thread_no);
    	
	return(res);
}
