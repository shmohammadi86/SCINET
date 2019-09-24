#ifndef DECODE_H
#define DECODE_H

//#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG

#include "armadillo"
using namespace arma;
using namespace std;

#include <omp.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <random>

void cumchi ( double *x, double *df, double *cum, double *ccum );
double r8_normal_01_cdf_inverse ( double p );
double nchoosek(int n, int k);
double r8_choose(int n, int k);

sp_mat read_from_mm(char *path);
sp_mat read_from_table(char *path);
sp_mat read_from_csv(char *path);
uvec read_uvec(char *path);

namespace SCINET {

	// Computes transformed gene activity scores for each sample
	mat compute_gene_activities(mat A, uvec samples, bool consider_baseline_expression, int thread_no);
	mat compute_gene_activities_decoupled(mat archetypes, mat H, uvec samples, bool consider_baseline_expression, int thread_no);
	mat compute_gene_activities_full(mat A, int thread_no);
		
	// Computes transformed gene activity scores using subsampling
	mat subsample_gene_activities(mat A, uvec rows, uvec samples, int total_subsamples, int cells_per_subsample, int thread_no, int seed);
	mat subsample_gene_activities_decoupled(mat archetypes, mat H, uvec rows, uvec samples, int total_subsamples, int cells_per_subsample, int thread_no, int seed);
	

	// Constructs one network per cell, using a reference network
	field<sp_mat> construct_cell_networks(sp_mat net, mat gene_activities, int thread_no);


	// Constructs one network per cell, without a reference network
	field<mat> construct_cell_networks_noPrior(mat gene_activities, int thread_no);


	// Construct cell networks and aggregates edge p-values
	field<sp_mat> construct_cell_networks_summary(sp_mat net, mat gene_activities, int total_subsamples, int cells_per_subsample, int seed, int thread_no);
}

#endif
