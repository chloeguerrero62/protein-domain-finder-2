import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_two_stage_performance():
    # Load data
    df = pd.read_csv('data/results/method_comparison.csv')
    
    # Filter for Two-Stage Spectral Method
    two_stage = df[df['method'] == 'Two-Stage-Spectral'].copy()
    
    # 1. Calculate Segmentation Categories
    perfect = two_stage[two_stage['n_predicted'] == two_stage['n_true']]
    under = two_stage[two_stage['n_predicted'] < two_stage['n_true']]
    over = two_stage[two_stage['n_predicted'] > two_stage['n_true']]
    
    n_total = len(two_stage)
    pct_perfect = (len(perfect) / n_total) * 100
    pct_under = (len(under) / n_total) * 100
    pct_over = (len(over) / n_total) * 100
    
    print(f"Total Proteins: {n_total}")
    print(f"Perfect Matches: {len(perfect)} ({pct_perfect:.1f}%)")
    print(f"Under-Segmented (Conservative): {len(under)} ({pct_under:.1f}%)")
    print(f"Over-Segmented (Aggressive): {len(over)} ({pct_over:.1f}%)")
    
    # 2. Extract Specific Case Studies mentioned in slide
    print(f"\n=== CASE STUDIES ===")
    targets = ['1BJR_E', '2J28_D', '1VK3_A']
    cols = ['pdb_chain', 'n_residues', 'n_true', 'n_predicted', 'exact_match']
    print(two_stage[two_stage['pdb_chain'].isin(targets)][cols].to_string(index=False))

    # 3. Generate the "Segmentation Breakdown" Chart
    plt.figure(figsize=(8, 5))
    categories = ['Perfect Match', 'Under-Segmented\n', 'Over-Segmented\n']
    counts = [len(perfect), len(under), len(over)]
    colors = ['#99ff99', '#ffcc99', '#ff9999'] # Green, Orange, Red
    
    bars = plt.bar(categories, counts, color=colors, edgecolor='black', alpha=0.8)
    
    # Add counts on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                 f'{int(height)}\n({height/n_total*100:.1f}%)',
                 ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.title('Two-Stage Spectral: Segmentation Behavior', fontsize=14, pad=15)
    plt.ylabel('Number of Proteins')
    plt.ylim(0, max(counts) * 1.2) # Add headroom for text
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('data/results/slide2_segmentation_breakdown.png', dpi=300)
    print("\n[Visual Generated]: data/results/slide2_segmentation_breakdown.png")

if __name__ == "__main__":
    analyze_two_stage_performance()