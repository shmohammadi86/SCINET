#include <scinet.h>
#include <cdflib.hpp>

/* 
 * From:Milton Abramowiz and Irene A Stegun. Handbook of Mathematical Functions.  National Bureau of Standards, 1964. 
 * Imp: http://finance.bi.no/~bernt/gcc_prog/
*/
double normalCDF(const double &z) {
	if (z >  6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };

	double b1 =  0.31938153; 
	double b2 = -0.356563782; 
	double b3 =  1.781477937;
	double b4 = -1.821255978;
	double b5 =  1.330274429; 
	double p  =  0.2316419; 
	double c2 =  0.3989423; 

	double a=fabs(z); 
	double t = 1.0/(1.0+a*p); 
	double b = c2*exp((-z)*(z/2.0)); 
	double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t; 
	n = 1.0-b*n; 
	if ( z < 0.0 ) n = 1.0 - n; 
	return n; 
}


vec FDR(vec pvals) {	
	pvals.replace(datum::nan, 1.0);

	unsigned int V = pvals.n_elem;
	double cV = 0, prev;
	for(int k = 1; k <= V; k++)
		cV += (1.0/k);

	uvec oidx;
	vec pval_adj;
	
	oidx = sort_index(pvals);
	pvals = sort(pvals); 

	pval_adj = zeros( V, 1 );		
	prev = 1.0;
	for (unsigned int kk = V-1; kk >= 1; kk--) { 
		pval_adj(oidx(kk)) = std::min(prev, pvals(kk)*V*cV/(kk+1));
		prev = pval_adj(oidx(kk));
	}
	pval_adj(oidx(0)) = std::min(prev, pvals(0)*V*cV);
			
	return pval_adj;
}


mat sampleUnif(int l, int m, double a, double b, int seed) {
	std::default_random_engine gen (seed);	
	std::uniform_real_distribution<double> unif(a, b);
	
	
	mat R(l, m);
	for (register int j = 0; j < m; j++) {
		for(register int i = 0; i < l; i++) {
			R(i, j) = unif(gen);
		}
	}
	return R;
}



namespace SCINET {
	mat compute_gene_activities_full(mat A, int thread_no = 4) {
		int M = A.n_rows;
		int N = A.n_cols;
		
		printf("M = %d, N = %d\n", M, N);
		
		printf("Normalizing columns (samples) ... ");
		mat Zr = zeros(M, N);
		for(int i = 0; i < N; i++) { // One sample at a time			
			vec v = A.col(i);
								
			uvec row_perm_forward = stable_sort_index(v);
			uvec row_perm = stable_sort_index(row_perm_forward);	
			vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
			
			vec v_RINT = zeros(size(p));
			for (int j = 0; j < p.n_elem; j++) {
				double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
				v_RINT(j) = norm_inv;
			}

			Zr.col(i) = v_RINT;						
		}
		printf("done\n");
		
		printf("Normalizing rows (genes) ... ");
		mat Zc = zeros(M, N);
		for(int i = 0; i < M; i++) { // One gene at a time			
			vec v = trans(A.row(i));
			
			uvec row_perm_forward = stable_sort_index(v);
			uvec row_perm = stable_sort_index(row_perm_forward);	
			vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
			
			vec v_RINT = zeros(size(p));
			for (int j = 0; j < p.n_elem; j++) {
				double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
				v_RINT(j) = norm_inv;
			}
			
			Zc.row(i) = trans(v_RINT);
		}	
		printf("done\n");
		

		mat Z = (Zr + Zc) / sqrt(2);

		return(Z);
	}
	


	mat compute_gene_activities(mat A, uvec samples, bool consider_baseline_expression = true, int thread_no = 4) {
		int M = A.n_rows;
		int N = A.n_cols;
		int K = samples.n_elem;
		printf("M = %d, N = %d, K = %d\n", M, N, K);

		mat Z = zeros(M, K);	
		
		vec row_RINT;
		if(consider_baseline_expression) {
			printf("Computing gene-specificity factor\n");	
			mat with_sample_A = A.cols(samples);
			
			vec row_means = mean(with_sample_A, 1);
			uvec row_perm_forward = stable_sort_index(row_means);
			uvec row_perm = stable_sort_index(row_perm_forward);	
			vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
			
			row_RINT = zeros(size(p));
			//#pragma omp parallel for num_threads(thread_no) 
			for (int j = 0; j < p.n_elem; j++) {
				double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
				row_RINT(j) = norm_inv;
			}
		}
	
		
		printf("Computing Rank-Based Inverse Normal Transformation\n");	

		#pragma omp parallel for num_threads(thread_no) 
		for(int i = 0; i < M; i++) { // One gene at a time
			
			rowvec v = A.row(i);
			
			uvec v_perm_forward = stable_sort_index(v);
			uvec v_perm = stable_sort_index(v_perm_forward);	
			vec p = (v_perm + ones(size(v_perm))) / (v_perm.n_elem + 1);		
			
			
			rowvec v_RINT(K);			
			for (int j = 0; j < K; j++) {				
				double norm_inv = r8_normal_01_cdf_inverse ( p(samples(j)) );
				v_RINT(j) = norm_inv;
			}
			
			if(consider_baseline_expression) {			
				Z.row(i) = (v_RINT + row_RINT(i)) / sqrt(2);
			}
			else {
				Z.row(i) = v_RINT;
			}			
		}
		
		printf("done\n");
		return(Z);
	}
	
	mat compute_gene_activities_decoupled(mat archetypes, mat H, uvec samples, bool consider_baseline_expression = true, int thread_no = 4) {
		int M = archetypes.n_rows;
		int N = H.n_cols;
		int K = samples.n_elem;
		mat Z = zeros(M, K);	
		printf("M = %d, N = %d, K = %d\n", M, N, K);
		
		vec row_RINT;
		if(consider_baseline_expression) {
			printf("Computing gene-specificity factor\n");	
			mat with_sample_A = archetypes*H.cols(samples);
			
			vec row_means = mean(with_sample_A, 1);
			uvec row_perm_forward = stable_sort_index(row_means);
			uvec row_perm = stable_sort_index(row_perm_forward);	
			vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
			
			row_RINT = zeros(size(p));
			
			//#pragma omp parallel for num_threads(thread_no) 
			for (int j = 0; j < p.n_elem; j++) {
				double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
				row_RINT(j) = norm_inv;
			}
		}
	
		
		printf("Computing Rank-Based Inverse Normal Transformation\n");	
		
		#pragma omp parallel for num_threads(thread_no) 
		for(int i = 0; i < M; i++) { // One gene at a time
			
			rowvec v = archetypes.row(i)*H;
			
			uvec v_perm_forward = stable_sort_index(v);
			uvec v_perm = stable_sort_index(v_perm_forward);	
			vec p = (v_perm + ones(size(v_perm))) / (v_perm.n_elem + 1);		
			
			
			rowvec v_RINT(K);			
			for (int j = 0; j < K; j++) {				
				double norm_inv = r8_normal_01_cdf_inverse ( p(samples(j)) );
				v_RINT(j) = norm_inv;
			}
			
			if(consider_baseline_expression) {			
				Z.row(i) = (v_RINT + row_RINT(i)) / sqrt(2);
			}
			else {
				Z.row(i) = v_RINT;
			}			
		}
		
		printf("done\n");
		return(Z);
	}
	
	mat subsample_gene_activities(mat A, uvec rows, uvec samples, int total_subsamples = 30, int cells_per_subsample = 10, int thread_no = 4, int seed = 0) {		
		umat subsamples = conv_to<umat>::from(round(sampleUnif(cells_per_subsample, total_subsamples, 0, samples.n_elem-1, seed)));
		
		
		mat Z = zeros(rows.n_elem, total_subsamples);	
		#pragma omp parallel for num_threads(thread_no) 		
		for(int i = 0; i < rows.n_elem; i++) { // One gene at a time			
			int r = rows(i);
			
			rowvec v = A.row(r);
			uvec v_perm_forward = stable_sort_index(v);
			uvec v_perm = stable_sort_index(v_perm_forward);	
			vec p = (v_perm + ones(size(v_perm))) / (v_perm.n_elem + 1);		
			
			
			rowvec v_RINT(samples.n_elem);			
			for (int j = 0; j < samples.n_elem; j++) {				
				double norm_inv = r8_normal_01_cdf_inverse ( p(samples(j)) );
				v_RINT(j) = norm_inv;
			}
			
			for(int j = 0; j < total_subsamples; j++) {
				Z(i, j) = sum(v_RINT(subsamples.col(j)));
			}
		}
		// Re-adjust the variance to become a standard normal distribution
		Z /= sqrt(cells_per_subsample);


		// Transform within-sample row means
		vec within_sample_row_means = zeros(A.n_rows);		
		for(int i = 0; i < A.n_rows; i++) {						
			rowvec v = A.row(i);
			within_sample_row_means(i) = mean(v(samples));
		}

		uvec row_perm_forward = stable_sort_index(within_sample_row_means);
		uvec row_perm = stable_sort_index(row_perm_forward);	
		vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);		
		
		p = p(rows);
		vec row_RINT = zeros(rows.n_elem);
		for (int j = 0; j < p.n_elem; j++) {
			double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
			row_RINT(j) = norm_inv;
		}

		// Incorporate row-factor (proportions?)
		for(int i = 0; i < rows.n_elem; i++) { // One gene at a time			
			for(int j = 0; j < total_subsamples; j++) {
				Z(i, j) = (Z(i, j) + row_RINT(i));
			}			
		}
		Z /= sqrt(2);
		
		printf("done\n");
		return(Z);
	}
	
	mat subsample_gene_activities_decoupled(mat archetypes, mat H, uvec rows, uvec samples, int total_subsamples = 30, int cells_per_subsample = 10, int thread_no = 4, int seed = 0) {		
		umat subsamples = conv_to<umat>::from(round(sampleUnif(cells_per_subsample, total_subsamples, 0, samples.n_elem-1, seed)));
		
		
		mat Z = zeros(rows.n_elem, total_subsamples);	
		#pragma omp parallel for num_threads(thread_no) 		
		for(int i = 0; i < rows.n_elem; i++) { // One gene at a time			
			int r = rows(i);
			
			rowvec v = archetypes.row(r)*H;
			uvec v_perm_forward = stable_sort_index(v);
			uvec v_perm = stable_sort_index(v_perm_forward);	
			vec p = (v_perm + ones(size(v_perm))) / (v_perm.n_elem + 1);		
			
			
			rowvec v_RINT(samples.n_elem);			
			for (int j = 0; j < samples.n_elem; j++) {				
				double norm_inv = r8_normal_01_cdf_inverse ( p(samples(j)) );
				v_RINT(j) = norm_inv;
			}
			
			for(int j = 0; j < total_subsamples; j++) {
				Z(i, j) = sum(v_RINT(subsamples.col(j)));
			}
		}
		// Re-adjust the variance to become a standard normal distribution
		Z /= sqrt(cells_per_subsample);


		// Transform within-sample row means
		vec within_sample_row_means = zeros(archetypes.n_rows);		
		for(int i = 0; i < archetypes.n_rows; i++) {						
			rowvec v = archetypes.row(i)*H;
			within_sample_row_means(i) = mean(v(samples));
		}

		uvec row_perm_forward = stable_sort_index(within_sample_row_means);
		uvec row_perm = stable_sort_index(row_perm_forward);	
		vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);		
		
		p = p(rows);
		vec row_RINT = zeros(rows.n_elem);
		for (int j = 0; j < p.n_elem; j++) {
			double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
			row_RINT(j) = norm_inv;
		}

		// Incorporate row-factor (proportions?)
		for(int i = 0; i < rows.n_elem; i++) { // One gene at a time			
			for(int j = 0; j < total_subsamples; j++) {
				Z(i, j) = (Z(i, j) + row_RINT(i));
			}			
		}
		Z /= sqrt(2);
		
		printf("done\n");
		return(Z);
	}
	
	
	field<mat> construct_cell_networks_noPrior(mat gene_activities, int thread_no = 4) {
		int M = gene_activities.n_rows;
		int K = gene_activities.n_cols;

		
		int perc = 0;
		int total_cells = 0;
		
		field<mat> networks(K);
		#pragma omp parallel for num_threads(thread_no)
		for(int k = 0; k < K; k++) {
			vec z = gene_activities.col(k);
			
			mat network = zeros(M, M);			
			for(int j = 0; j < M; j++) {
				double x = z(j);
				
				double z_j = z(j);
				for(int i = j+1; i < M; i++) {
					//double mean_stat = (z(i) + z(j)) / sqrt(2);
					//double pval = (1.0 - normalCDF(mean_stat));				
					
					double min_stat = std::min(z(i), z(j));								
					double tail_prob = (1.0 - normalCDF(min_stat));				
					double pval = tail_prob*tail_prob;
					
					network(i, j) = -log10(pval);
				}
				network = max(network, trans(network));
			}
			networks(k) = network;
						
			total_cells++;
			if(round(100*total_cells / K) > perc) {
				printf("%d %%\n", perc);
				perc++;
			}
			
		}
				
		return(networks);
	}



	// Returns value of edge weights per sample (each col corresponds to nnz encoded in subs)
	mat construct_cell_networks_driver(umat subs, mat gene_activities, int thread_no = 4) {
		int K = gene_activities.n_cols;
		int nnz = subs.n_cols;
		mat pvals_mat = zeros(nnz, K);
		
		int perc = 0;
		int total_cells = 0;
		
		#pragma omp parallel for num_threads(thread_no)
		for(int k = 0; k < K; k++) {
			vec z = gene_activities.col(k);
						
			vec pvals = ones(nnz);
			for(int l = 0; l < nnz; l++) {				
				int i = subs(0, l);
				int j = subs(1, l);
				//double mean_stat = (z(i) + z(j)) / sqrt(2);
				//double pval = (1.0 - normalCDF(mean_stat));				
				
				double min_stat = std::min(z(i), z(j));								
				double tail_prob = (1.0 - normalCDF(min_stat));				
				double pval = tail_prob*tail_prob;
				
				pvals(l) = pval;
			}
					
			pvals.transform( [](double val) { return (val < 1e-300?1e-300:val); } );
			
			pvals_mat.col(k) = pvals;

			total_cells++;
			if(round(100*total_cells / K) > perc) {
				perc = round(100*total_cells / K);
				printf("%d %%\n", perc);
			}
			
		}
		
		return pvals_mat;
	}

	// net: scaffold, A: gene activity scores
	field<sp_mat> construct_cell_networks(sp_mat net, mat gene_activities, int thread_no = 4) {
		int M = gene_activities.n_rows;
		int K = gene_activities.n_cols;
		field<sp_mat> networks(K);

		if(net.n_rows != M) {
			fprintf(stderr, "Error:: Number of genes in the gene activity score matrix (%d) doesn't match the number of nodes in the input network (%d)\n", gene_activities.n_rows, net.n_rows);
			return(networks);			
		}
		
		// Extract indices of nonzeros		
		net = trimatl(net);

		int nnz = net.n_nonzero;
		
		umat subs(2, nnz);
		int idx = 0;
		for(sp_mat::iterator it = net.begin(); it != net.end(); ++it) {
			subs(0, idx) = it.row();
			subs(1, idx) = it.col();
			idx++;
		}
  		
		mat pvals_mat = construct_cell_networks_driver(subs, gene_activities, thread_no);
		pvals_mat = -log10(pvals_mat);
		
		for(int k = 0; k < K; k++) {

			vec weights = pvals_mat.col(k);

			sp_mat network = sp_mat(subs, weights, M, M);
			network = (network + trans(network));
			network.diag().zeros();
			
			networks(k) = network;			
		}
				
		return(networks);
	}

	// net: scaffold, A: gene activity scores
	field<sp_mat> construct_cell_networks_summary(sp_mat net, mat gene_activities, int total_subsamples = 30, int cells_per_subsample = 10, int seed = 0, int thread_no = 4) {
		int M = gene_activities.n_rows;
		int K = gene_activities.n_cols;
		field<sp_mat> networks(K);

		if(net.n_rows != M) {
			fprintf(stderr, "Error:: Number of genes in the gene activity score matrix (%d) doesn't match the number of nodes in the input network (%d)\n", gene_activities.n_rows, net.n_rows);
			return(networks);			
		}
		
		// Extract indices of nonzeros		
		net = trimatl(net);

		int nnz = net.n_nonzero;
		
		umat subs(2, nnz);
		int idx = 0;
		for(sp_mat::iterator it = net.begin(); it != net.end(); ++it) {
			subs(0, idx) = it.row();
			subs(1, idx) = it.col();
			idx++;
		}
  		
		mat pvals_mat = construct_cell_networks_driver(subs, gene_activities, thread_no);
		
		mat aggregate_pvals = ones(pvals_mat.n_rows, total_subsamples);

		umat subsamples = conv_to<umat>::from(round(sampleUnif(cells_per_subsample, total_subsamples, 0, pvals_mat.n_cols-1, seed)));
		
		#pragma omp parallel for num_threads(thread_no)
		for(int i = 0; i < total_subsamples; i++) {

			mat subPvals = pvals_mat.cols(subsamples.col(i));
			
			// From: The harmonic mean p-value for combining dependent tests (Daniel J. Wilson, 2019)
			aggregate_pvals.col(i) = 1 / ((sum(1 / subPvals, 1) / cells_per_subsample));			
		}
		
		mat logPvals = -log10(aggregate_pvals);		
		for(int k = 0; k < K; k++) {

			vec weights = logPvals.col(k);

			sp_mat network = sp_mat(subs, weights, M, M);
			network = (network + trans(network));
			network.diag().zeros();
			
			networks(k) = network;			
		}
				
		return(networks);
	}
}
