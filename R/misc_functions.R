getLCC <- function(adj, genes) {
  G = igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  comps = igraph::components(G)
  comp_membership = comps$membership
  comp_sizes = comps$csize
  LCC_idx = which.max(comp_sizes)
  idx = which(comp_membership == LCC_idx)
  LCC = igraph::induced_subgraph(G, idx)
  igraph::V(LCC)$genes = genes[idx]
  return(LCC)
}

get.symPR.matrix <- function(adj, alpha_val) {
  d = Matrix::colSums(adj)
  isolated_vertices_mask = d == 0;
  d[isolated_vertices_mask] = 1;
  Dsqrt_inv = Matrix::Diagonal(n = length(d), x = 1 / sqrt(d))

  P = Dsqrt_inv %*% adj %*% Dsqrt_inv;
  
  I = Matrix::.symDiagonal(length(d))
  A = I - alpha_val*P
  
  return(A)
}

compute.genesets.compactness <- function(G, genesets, sample_no = 100, min.genes = 3) {
  nonisolated.vertices = which(igraph::degree(G) > 0)
  G = igraph::induced_subgraph(G, igraph::V(G)[nonisolated.vertices])
  
  if(is.null(names(genesets))) {
	  names(genesets) = 1:length(genesets)
  }

  set.seed(0)
  
  # From: "Systematic Evaluation of Molecular Networks for Discovery of Disease Genes"
  m=-0.02935302
  b=0.74842057
  log_edge_count = log10(igraph::ecount(G))
  alpha_val = m*log_edge_count+b

  # Symmetric transition matrix
  A = get.symPR.matrix(igraph::as_adjacency_matrix(G), alpha_val)
  Q = solve (as.matrix(A))

  genes = igraph::V(G)$name
  N = length(genes)
  
  specificity = igraph::V(G)$specificity # Z-score
  names(specificity) = genes


  compactness = array(0, length(genesets))  
  for(i in 1:length(genesets)) {
    message(sprintf('\t(%d/%d) %s ... ', i, length(genesets), names(genesets)[[i]]))
    
    geneset = genesets[[i]]
    shared_genes = intersect(geneset, genes)
    idx = match(shared_genes, genes)
    K = length(idx)
    
  	if(K < min.genes)
  		next
	
    e_gs = sparseVector(1 / K, idx, N);
    e_spec = sparseVector(specificity[idx], idx, N); 
    e_spec = e_spec / sum(e_spec)
    
    stat = Matrix::t(e_gs) %*% Q %*% e_spec
    
    rand_samples = replicate(sample_no, sample(N, K))

    
    E_gs_rand = sapply(1:sample_no, function(i) {
      rand_idx = rand_samples[,i]
      e_gs_rand = as.numeric(sparseVector(1 / K, rand_idx, N));
      return(e_gs_rand)} )
    
    E_spec_rand = sapply(1:sample_no, function(i) {
      rand_idx = rand_samples[,i]
      e_spec_rand = as.numeric(sparseVector(specificity[rand_idx], rand_idx, N));
      return(e_spec_rand)} )

    ss = colSums(E_spec_rand)
    ss[ss == 0] = 1
    E_spec_rand = scale(E_spec_rand, center = FALSE, scale = ss)
        
    rand_stats = diag(Matrix::t(E_gs_rand) %*% Q %*% E_spec_rand);
    
    compactness[i] = as.numeric((stat - mean(rand_stats)) / sd(rand_stats))

        
    message(sprintf('%.2f', compactness[i]))
    
  }
  
  return(compactness)
}


prioritize.geneset.LOO <- function(G, geneset, sample_no = 100) {
  # From: "Systematic Evaluation of Molecular Networks for Discovery of Disease Genes"
  nonisolated.vertices = which(igraph::degree(G) > 0)
  G = igraph::induced_subgraph(G, igraph::V(G)[nonisolated.vertices])
  genes = igraph::V(G)$name

  
  set.seed(0)
  
  m=-0.02935302
  b=0.74842057
  log_edge_count = log10(igraph::ecount(G))
  alpha_val = m*log_edge_count+b

  # Symmetric transition matrix
  A = get.symPR.matrix(igraph::as_adjacency_matrix(G), alpha_val)
  Q = solve (as.matrix(A))

  N = length(genes)
  
  specificity = igraph::V(G)$specificity # pre-normalized in [0, 1]
  names(specificity) = genes

  seeds = intersect(geneset, genes)
  seed_ids = match(seeds, genes) 

  gene.score = array(-Inf, length(seeds))
  rownames(gene.score) = seeds
  
  rest = setdiff(1:N, seed_ids)  
  for (i in 1:length(seeds)) {
    seed_id = seed_ids[i]
    print(i)
    idx = setdiff(seed_ids, seed_id)

    e_gs = as.numeric(sparseVector(specificity[idx], idx, N));
    e_gs = e_gs / sum(e_gs)

    scores = (1-alpha_val)*Q%*%e_gs
    mu = mean(scores[rest])
    sigma = sd(scores[rest])

    current_seed_name = seeds[i]
    gene.score[seeds[i]] = ( scores[seed_id] - mu ) / sigma
  }

  
  return(gene.score)
}



read.edgelist <- function(fname, header = TRUE, src.col = 1, dst.col = 2) {
  EdgeList = read.table(fname, sep = '\t', as.is = TRUE, header = header)
  src = EdgeList[, src.col]
  dst = EdgeList[, dst.col]
  
  nodes = union(src, dst)
  
  nV = length(nodes)
  ii = match(src, nodes)
  jj = match(dst, nodes)
  A = Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(nV, nV))
  G = igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = NULL)
  igraph::V(G)$name = nodes
  
  blacklist = c('UBC', 'SUMO1', 'SUMO2', 'SUMO3',
                igraph::V(G)$name[grep('^RPL|^RPS|^MRP', igraph::V(G)$name)])
  
  G = igraph::induced_subgraph(G, igraph::V(G)[!(igraph::V(G)$name %in% blacklist)])

  return(igraph::simplify(G))
}

pair.datasets <- function(G, activity.scores, selected.genes = NULL) {
  
  g1 = igraph::V(G)$name
  g2 = rownames(activity.scores)
  common.genes = intersect(g1, g2)
  if(! is.null(selected.genes) ) {
	  common.genes = intersect(selected.genes, common.genes)
  }
  
  idx = match(common.genes, g1)
  sub.net = igraph::induced_subgraph(G, igraph::V(G)[idx])
  
  idx = match(common.genes, g2)
  sub.activity.scores = activity.scores[idx, ]
  
  nV = length(common.genes)
  EL = igraph::as_edgelist(sub.net)
  idx = (match(EL[, 1], igraph::V(sub.net)$name) - 1) * nV +
        match(EL[, 2], igraph::V(sub.net)$name)
  igraph::E(sub.net)$idx = idx
  
  out = list(net = sub.net, activity.scores = sub.activity.scores, genes = common.genes)
  
  return(out)
}

#' Run the SCINET pipeline end to end
#'
#' Subsamples gene activity scores over the requested cells, restricts both the
#' reference network and the activity scores to their shared genes, and builds
#' one cell-type-specific network per subsample.
#'
#' @param G reference interactome as an igraph object whose vertices are named
#'   with gene symbols.
#' @param expression imputed expression profile, genes in rows.
#' @param samples column indices of the cells to profile.
#' @param genes optional subset of gene symbols to restrict to.
#' @param total_subsamples number of subsamples to draw.
#' @param cells_per_subsample cells drawn per subsample.
#' @param thread_no worker threads.
#' @param seed random seed for subsampling.
#'
#' @return a list of sparse adjacency matrices, one per subsample.
#' @export
run.SCINET <- function(G, expression, samples = NULL, genes = NULL,
                       total_subsamples = 100, cells_per_subsample = 10,
                       thread_no = 8, seed = 0) {
  if (is.null(samples) || length(samples) == 0) {
    stop("run.SCINET: provide `samples`, the column indices of the cells to profile.")
  }

  A <- as.matrix(expression)
  if (is.null(genes)) {
    rows <- seq_len(nrow(A))
  } else {
    rows <- match(toupper(genes), toupper(rownames(A)))
    rows <- rows[!is.na(rows)]
    if (length(rows) == 0) {
      stop("run.SCINET: none of `genes` are present in rownames(expression).")
    }
  }

  activity.scores <- subsample_gene_activities(
    A = A, rows = rows, samples = samples,
    total_subsamples = total_subsamples,
    cells_per_subsample = cells_per_subsample,
    thread_no = thread_no, seed = seed)
  rownames(activity.scores) <- rownames(A)[rows]

  paired.ds <- pair.datasets(G, activity.scores)
  net <- methods::as(igraph::as_adjacency_matrix(paired.ds$net), "TsparseMatrix")

  construct_cell_networks(net = net,
                          gene_activities = paired.ds$activity.scores,
                          thread_no = thread_no)
}

topo.spec <- function(G, sample_no = 100) {
  d = igraph::degree(G)
  w = igraph::strength(G)

  edge_weights = igraph::E(G)$weight
  N = length(edge_weights)
  
  UD = setdiff(sort(unique(d)), 0)
  mean_rand_stat = array(0, max(UD))  
  sd_rand_stat = array(0, max(UD))  
  for(deg in UD) {
    rand_samples = sample(edge_weights, sample_no*deg, replace=TRUE);
    start_idx = seq(from=1, to = length(rand_samples), by = deg);
    rand_sums = sapply(start_idx, function(i) {sum(rand_samples[i:(i+deg-1)])} )
    
    mean_rand_stat[deg] = mean(rand_sums);
    sd_rand_stat[deg] = sd(rand_sums);    
  }
  
  specificity = sapply( 1:length(w), function(i) {if(w[i] == 0) {x = -Inf;} else {x = (w[i] - mean_rand_stat[d[i]]) / sd_rand_stat[d[i]];}; return(x) } )
  
  specificity[is.na(specificity)] = -Inf
    
  names(specificity) = V(G)$name
  return(specificity)
}


trasnc.spec.ttest <- function(Arch.imputed.profile, Labels) {
	gene_names = rownames(Arch.imputed.profile)
	
	if(!is.factor(Labels))
		Labels = factor(Labels, levels = sort(unique(Labels)))
	
	UL = levels(Labels)
	
  transcriptional.gene.specificity.Z = array(0, dim = c(length(gene_names), length(UL)))
  rownames(transcriptional.gene.specificity.Z) = gene_names
  colnames(transcriptional.gene.specificity.Z) = UL
  
  for (celltype in UL) {
    message(celltype)
    
    Group = which(Labels %in% celltype)
    Null = which(!(Labels %in% celltype))

    A_group = Arch.imputed.profile[, Group]
    A_null = Arch.imputed.profile[, Null]
    
    delta_mean = rowMeans(A_group) - rowMeans(A_null);
    sigma1_sq = apply(A_group, 1, var)
    sigma2_sq = apply(A_null, 1, var)
    sigma_pooled = sqrt( (sigma1_sq / length(Group)) + (sigma2_sq / length(Null)) )
    z.stat = delta_mean / sigma_pooled
    
    transcriptional.gene.specificity.Z[, celltype] = z.stat
  }

	return(transcriptional.gene.specificity.Z)
}
