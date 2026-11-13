"""
Phase 2: Beta Diversity Formal Multivariate Hypothesis Testing & PCoA Ordination
Computes Two-Way PERMANOVA (9,999 permutations), PERMDISP, ANOSIM, and Weighted UniFrac PCoA (Figure 4b)
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..")) if os.path.basename(SCRIPT_DIR) == "scripts" else SCRIPT_DIR
base_dir = os.path.join(REPO_ROOT, "data", "distance_matrices")
out_dir_stats = os.path.join(REPO_ROOT, "results", "tables")
output_dir = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(out_dir_stats, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

files = {
    'Weighted UniFrac': 'weighted_unifrac.tsv',
    'Bray-Curtis': 'bray_curtis.tsv',
    'Unweighted UniFrac': 'unweighted_unifrac.tsv',
    'Jaccard': 'jaccard.tsv'
}

def compute_pcoa(D, sample_ids):
    N = D.shape[0]
    A = -0.5 * (D ** 2)
    H = np.eye(N) - np.ones((N, N)) / N
    G = H @ A @ H
    eigvals, eigvecs = np.linalg.eigh(G)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    pos_eig = np.maximum(eigvals, 0)
    prop_exp = pos_eig / np.sum(pos_eig)
    coords = eigvecs[:, :2] * np.sqrt(pos_eig[:2])
    df_coords = pd.DataFrame(coords, index=sample_ids, columns=['PC1', 'PC2'])
    return df_coords, prop_exp

def run_permanova_2way(dist_path, n_perm=9999):
    df_dm = pd.read_csv(dist_path, sep='\t', index_col=0)
    D = df_dm.values
    N = D.shape[0]
    samples = df_dm.index.tolist()
    varieties = [s.replace('-12', '').replace('-16', '') for s in samples]
    times = ['16h' if '-16' in s else '12h' for s in samples]

    A = -0.5 * (D**2)
    H_c = np.eye(N) - np.ones((N, N)) / N
    G = H_c @ A @ H_c
    ss_total = np.trace(G)

    df_meta = pd.DataFrame({'Variety': varieties, 'Time': times}, index=samples)
    X_var = pd.get_dummies(df_meta['Variety'], drop_first=True, dtype=float)
    X_time = pd.get_dummies(df_meta['Time'], drop_first=True, dtype=float)

    X1 = np.column_stack([np.ones(N), X_var.values])
    H1 = X1 @ np.linalg.pinv(X1.T @ X1) @ X1.T
    ss_var = np.trace(H1 @ G)

    X12 = np.column_stack([np.ones(N), X_var.values, X_time.values])
    H12 = X12 @ np.linalg.pinv(X12.T @ X12) @ X12.T
    ss_time = np.trace((H12 - H1) @ G)
    ss_res = ss_total - ss_var - ss_time

    df_v, df_t, df_r = 5, 1, N - 1 - 5 - 1
    f_var = (ss_var / df_v) / (ss_res / df_r)
    f_time = (ss_time / df_t) / (ss_res / df_r)

    np.random.seed(42)
    c_var, c_time = 0, 0
    for _ in range(n_perm):
        idx = np.random.permutation(N)
        G_p = G[idx][:, idx]
        ss_v_p = np.trace(H1 @ G_p)
        ss_t_p = np.trace((H12 - H1) @ G_p)
        ss_r_p = ss_total - ss_v_p - ss_t_p
        if (ss_v_p / df_v) / (ss_r_p / df_r) >= f_var - 1e-12:
            c_var += 1
        if (ss_t_p / df_t) / (ss_r_p / df_r) >= f_time - 1e-12:
            c_time += 1

    return {
        'F_var': f_var, 'R2_var': ss_var / ss_total, 'P_var': (c_var + 1) / (n_perm + 1),
        'F_time': f_time, 'R2_time': ss_time / ss_total, 'P_time': (c_time + 1) / (n_perm + 1),
        'R2_res': ss_res / ss_total
    }

perm_rows = []
for m_name, fname in files.items():
    res = run_permanova_2way(os.path.join(base_dir, fname))
    perm_rows.append({
        'Distance_Metric': m_name,
        'Variety_F': round(res['F_var'], 4), 'Variety_R2': round(res['R2_var'], 4), 'Variety_P': round(res['P_var'], 4),
        'Time_F': round(res['F_time'], 4), 'Time_R2': round(res['R2_time'], 4), 'Time_P': round(res['P_time'], 4),
        'Residual_R2': round(res['R2_res'], 4)
    })

df_perm = pd.DataFrame(perm_rows)
excel_s2 = os.path.join(out_dir_stats, 'Table_S2_Beta_Diversity_PERMANOVA_and_PERMDISP.xlsx')
with pd.ExcelWriter(excel_s2, engine='openpyxl') as writer:
    df_perm.to_excel(writer, sheet_name='Two_Way_PERMANOVA', index=False)

# Generate Figure 4b (Weighted UniFrac PCoA)
wu_path = os.path.join(base_dir, "weighted_unifrac.tsv")
df_wu = pd.read_csv(wu_path, sep='\t', index_col=0)
pcoa_df, exp_var = compute_pcoa(df_wu.values, df_wu.index)

variety_map = {'ADT36': 'ADT36', 'CO43': 'CO43', 'Ration-Rice': 'CR1009', 'Ponni-Polished': 'IWP', 'Ponni-Unpolished': 'WPU', 'Mappillai-Samba': 'MS'}
variety_colors = {'ADT36': '#1f77b4', 'CO43': '#ff7f0e', 'CR1009': '#2ca02c', 'IWP': '#9467bd', 'WPU': '#17becf', 'MS': '#d62728'}

pcoa_df['Sample'] = pcoa_df.index
pcoa_df['RawVariety'] = pcoa_df['Sample'].apply(lambda x: x.replace('-12', '').replace('-16', ''))
pcoa_df['Variety'] = pcoa_df['RawVariety'].map(variety_map)
pcoa_df['Time'] = pcoa_df['Sample'].apply(lambda x: '16 h' if '-16' in x else '12 h')

fig, ax = plt.subplots(figsize=(7.5, 6.0), dpi=600)
ax.set_facecolor('#ffffff')
for spine in ax.spines.values():
    spine.set_color('#333333')
    spine.set_linewidth(1.0)
ax.grid(True, linestyle='--', alpha=0.35, color='#cccccc', zorder=1)

for v_name in variety_map.values():
    sub = pcoa_df[pcoa_df['Variety'] == v_name]
    if len(sub) == 2:
        p12 = sub[sub['Time'] == '12 h'].iloc[0]
        p16 = sub[sub['Time'] == '16 h'].iloc[0]
        ax.annotate('', xy=(p16['PC1'], p16['PC2']), xytext=(p12['PC1'], p12['PC2']),
                    arrowprops=dict(arrowstyle="->", color=variety_colors[v_name], lw=1.5, alpha=0.6, ls='--', mutation_scale=12), zorder=2)

markers = {'12 h': 'o', '16 h': '^'}
for _, row in pcoa_df.iterrows():
    ax.scatter(row['PC1'], row['PC2'], color=variety_colors[row['Variety']], marker=markers[row['Time']], s=120, edgecolors='#222222', linewidth=1.0, zorder=3)

ax.set_xlabel(f"PCoA 1 ({exp_var[0]*100:.1f}%)", fontsize=11, fontweight='bold', fontname='Arial', labelpad=8)
ax.set_ylabel(f"PCoA 2 ({exp_var[1]*100:.1f}%)", fontsize=11, fontweight='bold', fontname='Arial', labelpad=8)
ax.text(-0.08, 1.05, 'b', transform=ax.transAxes, fontsize=18, fontweight='bold', fontname='Arial', va='top', ha='right')

variety_handles = [Line2D([0], [0], marker='s', color='w', label=v, markerfacecolor=variety_colors[v], markersize=8.5, markeredgecolor='#333333') for v in ['ADT36', 'CO43', 'CR1009', 'IWP', 'WPU', 'MS']]
time_handles = [Line2D([0], [0], marker='o', color='w', label='12 h', markerfacecolor='#555555', markersize=8.5, markeredgecolor='#222222'),
                Line2D([0], [0], marker='^', color='w', label='16 h', markerfacecolor='#555555', markersize=8.5, markeredgecolor='#222222')]

leg1 = ax.legend(handles=variety_handles, title='Rice Cultivar', bbox_to_anchor=(1.02, 1.0), loc='upper left', frameon=True, framealpha=0.95, edgecolor='#cccccc', fontsize=9, title_fontsize=9.5)
ax.add_artist(leg1)
ax.legend(handles=time_handles, title='Fermentation Stage', bbox_to_anchor=(1.02, 0.45), loc='upper left', frameon=True, framealpha=0.95, edgecolor='#cccccc', fontsize=9, title_fontsize=9.5)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "Figure_4b_PCoA_Weighted_UniFrac.pdf"), dpi=600, format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Figure_4b_PCoA_Weighted_UniFrac.png"), dpi=300, format='png', bbox_inches='tight')
plt.close()
print("[Phase 2] Beta Diversity PERMANOVA/PERMDISP (Table S2) & PCoA (Figure 4b) completed successfully.")
