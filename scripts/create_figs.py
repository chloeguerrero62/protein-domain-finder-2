import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

figures_dir = Path('data/results/figures')
figures_dir.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11

method_comp = pd.read_csv('data/results/method_comparison.csv')
random_ctrl = pd.read_csv('data/results/random_controls.csv')
hyperparam = pd.read_csv('data/results/hyperparameter_grid_search.csv')

successful = method_comp[method_comp['success'] == True]



####
# FIGURE 1: Method Comparison Summary

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Exact Match Rates
exact_rates = successful.groupby('method')['exact_match'].apply(
    lambda x: (x == 1).sum() / len(x) * 100
).reset_index()
exact_rates.columns = ['Method', 'Exact Match (%)']

# Color mapping for the three methods
def get_color(method):
    if 'Louvain' in method:
        return '#d62728'  # red - poor performance
    elif 'Two-Stage' in method:
        return '#2ca02c'  # green - good unsupervised
    else:  # Spectral-Graph
        return '#1f77b4'  # blue - supervised baseline

colors = [get_color(m) for m in exact_rates['Method']]
ax1.barh(exact_rates['Method'], exact_rates['Exact Match (%)'], color=colors, alpha=0.7)
ax1.set_xlabel('Exact Match Rate (%)')
ax1.set_title('A) Domain Count Accuracy', fontweight='bold')
ax1.set_xlim(0, 105)

# Panel B: Mean Absolute Error
mae = successful.groupby('method')['absolute_error'].mean().reset_index()
mae.columns = ['Method', 'MAE']

colors = [get_color(m) for m in mae['Method']]
ax2.barh(mae['Method'], mae['MAE'], color=colors, alpha=0.7)
ax2.set_xlabel('Mean Absolute Error (domains)')
ax2.set_title('B) Average Prediction Error', fontweight='bold')
ax2.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5)

# Panel C: Error Distribution (all three methods)
methods_subset = ['Louvain', 'Two-Stage-Spectral', 'Spectral-Graph']
data_subset = successful[successful['method'].isin(methods_subset)]
sns.violinplot(data=data_subset, y='method', x='absolute_error', ax=ax3, orient='h')
ax3.set_xlabel('Absolute Error (domains)')
ax3.set_ylabel('')
ax3.set_title('C) Error Distribution', fontweight='bold')

# Panel D: Method Type Comparison (Supervised vs Unsupervised)
supervised = successful[successful['supervised'] == True]
unsupervised = successful[successful['supervised'] == False]

method_type_mae = pd.DataFrame({
    'Type': ['Supervised\n(Spectral-Graph)', 'Unsupervised\n(Two-Stage)', 'Unsupervised\n(Louvain)'],
    'MAE': [
        supervised['absolute_error'].mean(),
        unsupervised[unsupervised['method'] == 'Two-Stage-Spectral']['absolute_error'].mean(),
        unsupervised[unsupervised['method'] == 'Louvain']['absolute_error'].mean()
    ]
})

colors = ['#1f77b4', '#2ca02c', '#d62728']
ax4.bar(method_type_mae['Type'], method_type_mae['MAE'], color=colors, alpha=0.7)
ax4.set_ylabel('Mean Absolute Error')
ax4.set_title('D) Supervised vs Unsupervised', fontweight='bold')
ax4.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig(figures_dir / 'fig1_comprehensive_comparison.png', bbox_inches='tight')
plt.close()

print("Figure 1: Comprehensive Method Comparison")


####
# FIGURE 2: Oracle vs No-Oracle Comparison


fig, ax = plt.subplots(figsize=(10, 6))

# Combine method comparison with random controls
oracle_methods = successful[successful['supervised'] == True].copy()
oracle_methods['category'] = 'Supervised\n(Oracle)'

two_stage = successful[successful['method'] == 'Two-Stage-Spectral'].copy()
two_stage['category'] = 'Two-Stage\n(No Oracle)'

louvain = successful[successful['method'] == 'Louvain'].copy()
louvain['category'] = 'Louvain\n(No Oracle)'

# Combine
combined = pd.concat([
    oracle_methods[['category', 'exact_match', 'absolute_error']],
    two_stage[['category', 'exact_match', 'absolute_error']],
    louvain[['category', 'exact_match', 'absolute_error']]
])

# Group and plot
summary = combined.groupby('category').agg({
    'exact_match': lambda x: (x == 1).sum() / len(x) * 100,
    'absolute_error': 'mean'
}).reset_index()

x = np.arange(len(summary))
width = 0.35

bars1 = ax.bar(x - width/2, summary['exact_match'], width, 
               label='Exact Match Rate (%)', alpha=0.7, color='steelblue')
ax2 = ax.twinx()
bars2 = ax2.bar(x + width/2, summary['absolute_error'], width,
                label='Mean Absolute Error', alpha=0.7, color='coral')

ax.set_ylabel('Exact Match Rate (%)', fontsize=12)
ax2.set_ylabel('Mean Absolute Error (domains)', fontsize=12)
ax.set_xlabel('Method Category', fontsize=12)
ax.set_title('Oracle Knowledge vs Clustering Quality', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(summary['category'])
ax.set_ylim(0, 110)

# Add legend
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.tight_layout()
plt.savefig(figures_dir / 'fig2_oracle_comparison.png', bbox_inches='tight')
plt.close()

print("Figure 2: Oracle Knowledge Analysis")

####
# FIGURE 3: Hyperparameter Optimization Heatmap

# Create heatmap showing sigma_factor vs max_domains for silhouette method
silhouette_results = hyperparam[hyperparam['estimation_method'] == 'silhouette']

# Pivot for heatmap
heatmap_data = silhouette_results.pivot_table(
    values='mean_error',
    index='sigma_factor',
    columns='max_domains',
    aggfunc='mean'
)

fig, ax = plt.subplots(figsize=(10, 6))
sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='RdYlGn_r', 
            cbar_kws={'label': 'Mean Absolute Error'}, ax=ax)
ax.set_title('Hyperparameter Optimization: Sigma Factor vs Max Domains',
             fontsize=14, fontweight='bold')
ax.set_xlabel('Max Domains', fontsize=12)
ax.set_ylabel('Sigma Factor', fontsize=12)

plt.tight_layout()
plt.savefig(figures_dir / 'fig3_hyperparameter_heatmap.png', bbox_inches='tight')
plt.close()

print("Figure 3: Hyperparameter Optimization Heatmap")

####
# FIGURE 4: Performance by Domain Count

two_stage_only = successful[successful['method'] == 'Two-Stage-Spectral']

# Group by true domain count
performance_by_n = two_stage_only.groupby('n_true').agg({
    'absolute_error': ['mean', 'std'],
    'exact_match': lambda x: (x == 1).sum() / len(x) * 100,
    'pdb_chain': 'count'
}).reset_index()

performance_by_n.columns = ['n_domains', 'mae_mean', 'mae_std', 'exact_match_pct', 'count']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: MAE by true domain count
ax1.errorbar(performance_by_n['n_domains'], performance_by_n['mae_mean'],
             yerr=performance_by_n['mae_std'], marker='o', linewidth=2,
             markersize=8, capsize=5, capthick=2)
ax1.set_xlabel('True Number of Domains', fontsize=12)
ax1.set_ylabel('Mean Absolute Error', fontsize=12)
ax1.set_title('A) Accuracy vs Protein Complexity', fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xticks(performance_by_n['n_domains'])

# Panel B: Exact match rate by domain count
bars = ax2.bar(performance_by_n['n_domains'], performance_by_n['exact_match_pct'],
               alpha=0.7, color='steelblue')
ax2.set_xlabel('True Number of Domains', fontsize=12)
ax2.set_ylabel('Exact Match Rate (%)', fontsize=12)
ax2.set_title('B) Success Rate by Complexity', fontweight='bold')
ax2.set_xticks(performance_by_n['n_domains'])
ax2.set_ylim(0, 100)

# Add count labels
for i, (bar, count) in enumerate(zip(bars, performance_by_n['count'])):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 2,
             f'n={int(count)}', ha='center', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / 'fig4_performance_by_complexity.png', bbox_inches='tight')
plt.close()

print("Figure 4: Performance by Domain Count")

print(f"\nAll figures saved to: {figures_dir}/")
