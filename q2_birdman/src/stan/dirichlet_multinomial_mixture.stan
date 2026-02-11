data {
  int<lower=1> N;              // number of samples
  int<lower=1> D;              // number of features (taxa)
  int<lower=1> K;              // number of clusters
  array[N, D] int y;           // count matrix [samples, features]
  vector[N] depth;             // sequencing depths
  
  // hyperparameters
  vector<lower=0>[K] alpha;    // prior for mixture weights
  vector<lower=0>[D] beta;     // prior for Dirichlet components
}

parameters {
  simplex[K] pi;               // mixing proportions
  array[K] simplex[D] phi;     // cluster-specific taxa distributions
  array[N] simplex[K] theta;   // sample-cluster assignments
}

model {
  // priors
  pi ~ dirichlet(alpha);
  for (k in 1:K) {
    phi[k] ~ dirichlet(beta);
  }
  
  // likelihood
  for (n in 1:N) {
    theta[n] ~ dirichlet(alpha);
    
    // multinomial likelihood with mixing
    vector[D] log_prob = log(theta[n]' * to_matrix(phi)) + depth[n];
    y[n] ~ multinomial(softmax(log_prob));
  }
}

generated quantities {
  // Compute cluster assignments and probabilities
  array[N] int<lower=1, upper=K> cluster_assignments;
  matrix[N, K] cluster_probs;
  
  for (n in 1:N) {
    // Get most likely cluster
    cluster_assignments[n] = categorical_rng(theta[n]);
    
    // Store probabilities
    cluster_probs[n] = to_row_vector(theta[n]);
  }
} 