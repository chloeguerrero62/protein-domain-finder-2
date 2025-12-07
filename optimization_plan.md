# Hyperparameter Optimization Strategy

## Method: Two-Stage Spectral Clustering

### Why Two-Stage Spectral?
- **Only viable unsupervised method** (Louvain gets 0% exact matches)
- **Promising preliminary results**: ~35% exact match, MAE ~1.2
- **Clear improvement path** through parameter tuning

### Parameters to Optimize

1. **Domain count estimation method**:
   - Silhouette score maximization
   - Eigengap heuristic
   - Consensus (median of multiple methods)

2. **Similarity kernel sigma**:
   - σ = median(distances) × factor
   - Test factors: [0.5, 1.0, 1.5, 2.0, 2.5]

3. **Search range for k estimation**:
   - max_domains: [6, 8, 10, 12]
   - Trade-off: larger range = better coverage, slower runtime

4. **k-NN graph construction**:
   - k: [5, 10, 15, 20]
   - Affects Louvain too, but main use is in Two-Stage

### Optimization Strategy

**Training set**: 30 proteins (stratified by n_domains)
**Validation metric**: Mean Absolute Error (MAE) in domain count
**Grid search**: ~80 parameter combinations
**Runtime**: ~20 minutes on 30 proteins

### Expected Improvement
- Baseline (current defaults): MAE ~1.2, 35% exact
- Optimized (tuned params): Target MAE <1.0, >40% exact
