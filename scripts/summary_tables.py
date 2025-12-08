"""
Generate formatted tables for presentation and paper
"""

import pandas as pd
import numpy as np

# Load results
method_comp = pd.read_csv('data/results/method_comparison.csv')
random_ctrl = pd.read_csv('data/results/random_controls.csv')
hyperparam = pd.read_csv('data/results/hyperparameter_grid_search.csv')

successful = method_comp[method_comp['success'] == True]

# =============================================================================
# TABLE 1: Method Comparison Summary
# =============================================================================

table1 = successful.groupby(['method', 'supervised']).agg({
    'pdb_chain': 'count',
    'exact_match': lambda x: (x == 1).sum(),
    'absolute_error': ['mean', 'std'],
    'n_predicted': 'mean'
}).round(2)

table1.columns = ['N', 'Exact Matches', 'MAE', 'MAE SD', 'Avg Predicted']
table1['Exact Match %'] = (table1['Exact Matches'] / table1['N'] * 100).round(1)

print("="*80)
print("TABLE 1: Method Comparison Summary")
print("="*80)
print(table1.to_string())
print()

table1.to_csv('data/results/table1_method_summary.csv')

# =============================================================================
# TABLE 2: Random Controls
# =============================================================================

table2 = random_ctrl.groupby('method').agg({
    'pdb_chain': 'count',
    'exact_match': lambda x: (x == 1).sum(),
    'absolute_error': ['mean', 'std']
}).round(2)

table2.columns = ['N', 'Exact Matches', 'MAE', 'MAE SD']
table2['Exact Match %'] = (table2['Exact Matches'] / table2['N'] * 100).round(1)

print("="*80)
print("TABLE 2: Random Controls")
print("="*80)
print(table2.to_string())
print()

table2.to_csv('data/results/table2_random_controls.csv')

# =============================================================================
# TABLE 3: Best Hyperparameters
# =============================================================================

best_params = hyperparam.nsmallest(1, 'mean_error')

print("="*80)
print("TABLE 3: Optimized Hyperparameters")
print("="*80)
print(f"Estimation Method: {best_params['estimation_method'].values[0]}")
print(f"Sigma Factor: {best_params['sigma_factor'].values[0]}")
print(f"Max Domains: {best_params['max_domains'].values[0]}")
print(f"k-NN: {best_params['k_graph'].values[0]}")
print(f"\nTraining Performance:")
print(f"  MAE: {best_params['mean_error'].values[0]:.2f}")
print(f"  Exact Match Rate: {best_params['exact_match_rate'].values[0]*100:.1f}%")
print()

# =============================================================================
# TABLE 4: Performance by Domain Count
# =============================================================================

two_stage = successful[successful['method'] == 'Two-Stage-Spectral']

table4 = two_stage.groupby('n_true').agg({
    'pdb_chain': 'count',
    'exact_match': lambda x: (x == 1).sum(),
    'absolute_error': ['mean', 'std']
}).round(2)

table4.columns = ['N Proteins', 'Exact Matches', 'MAE', 'MAE SD']
table4['Exact Match %'] = (table4['Exact Matches'] / table4['N Proteins'] * 100).round(1)

print("="*80)
print("TABLE 4: Performance by True Domain Count (Two-Stage Spectral)")
print("="*80)
print(table4.to_string())
print()

table4.to_csv('data/results/table4_performance_by_domain_count.csv')

print("All tables saved to data/results/")