import os
import pandas as pd
import numpy as np
from scipy import stats
import scikit_posthocs as sp
from statsmodels.stats.multitest import multipletests

# -------------------------------------------------------------
# 1. RAW DATA DEFINITION (Verified against 16S Amplicon Dataset)
# -------------------------------------------------------------
data = [
    {'Sample_ID': 'Ponni-Unpolished-12', 'Cultivar': 'White Ponni Unpolished (WPU)', 'Cultivar_Code': 'WPU', 'Time': '12 h', 'Observed': 45, 'Chao1': 45.000, 'Shannon': 0.727, 'Simpson': 0.302, 'Evenness': 0.191, 'Fisher': 4.360, 'Faith_PD': 1.006},
    {'Sample_ID': 'ADT36-12',            'Cultivar': 'ADT36',                       'Cultivar_Code': 'ADT36', 'Time': '12 h', 'Observed': 56, 'Chao1': 57.500, 'Shannon': 0.755, 'Simpson': 0.273, 'Evenness': 0.187, 'Fisher': 5.361, 'Faith_PD': 1.615},
    {'Sample_ID': 'Ration-Rice-12',      'Cultivar': 'CR1009 (Ration Rice)',         'Cultivar_Code': 'CR1009','Time': '12 h', 'Observed': 56, 'Chao1': 56.600, 'Shannon': 0.282, 'Simpson': 0.106, 'Evenness': 0.070, 'Fisher': 5.787, 'Faith_PD': 0.284},
    {'Sample_ID': 'Ponni-Polished-12',    'Cultivar': 'Improved White Ponni (IWP)',   'Cultivar_Code': 'IWP',   'Time': '12 h', 'Observed': 73, 'Chao1': 74.250, 'Shannon': 0.955, 'Simpson': 0.352, 'Evenness': 0.223, 'Fisher': 7.704, 'Faith_PD': 1.318},
    {'Sample_ID': 'Mappillai-Samba-12',   'Cultivar': 'Mapillai Samba (MS)',          'Cultivar_Code': 'MS',    'Time': '12 h', 'Observed': 75, 'Chao1': 75.000, 'Shannon': 1.571, 'Simpson': 0.638, 'Evenness': 0.364, 'Fisher': 7.958, 'Faith_PD': 1.597},
    {'Sample_ID': 'CO43-12',             'Cultivar': 'CO43',                        'Cultivar_Code': 'CO43',  'Time': '12 h', 'Observed': 54, 'Chao1': 54.000, 'Shannon': 2.045, 'Simpson': 0.825, 'Evenness': 0.513, 'Fisher': 5.379, 'Faith_PD': 1.337},
    
    {'Sample_ID': 'Ponni-Polished-16',    'Cultivar': 'Improved White Ponni (IWP)',   'Cultivar_Code': 'IWP',   'Time': '16 h', 'Observed': 42, 'Chao1': 42.000, 'Shannon': 0.311, 'Simpson': 0.100, 'Evenness': 0.083, 'Fisher': 4.322, 'Faith_PD': 0.515},
    {'Sample_ID': 'Ponni-Unpolished-16',  'Cultivar': 'White Ponni Unpolished (WPU)', 'Cultivar_Code': 'WPU', 'Time': '16 h', 'Observed': 79, 'Chao1': 79.000, 'Shannon': 1.474, 'Simpson': 0.591, 'Evenness': 0.337, 'Fisher': 8.764, 'Faith_PD': 1.633},
    {'Sample_ID': 'ADT36-16',            'Cultivar': 'ADT36',                       'Cultivar_Code': 'ADT36', 'Time': '16 h', 'Observed': 75, 'Chao1': 75.429, 'Shannon': 1.311, 'Simpson': 0.589, 'Evenness': 0.304, 'Fisher': 7.815, 'Faith_PD': 1.567},
    {'Sample_ID': 'CO43-16',             'Cultivar': 'CO43',                        'Cultivar_Code': 'CO43',  'Time': '16 h', 'Observed': 21, 'Chao1': 21.000, 'Shannon': 2.150, 'Simpson': 0.799, 'Evenness': 0.706, 'Fisher': 2.050, 'Faith_PD': 1.932},
    {'Sample_ID': 'Mappillai-Samba-16',   'Cultivar': 'Mapillai Samba (MS)',          'Cultivar_Code': 'MS',    'Time': '16 h', 'Observed': 51, 'Chao1': 51.000, 'Shannon': 2.031, 'Simpson': 0.789, 'Evenness': 0.516, 'Fisher': 5.023, 'Faith_PD': 1.006},
    {'Sample_ID': 'Ration-Rice-16',      'Cultivar': 'CR1009 (Ration Rice)',         'Cultivar_Code': 'CR1009','Time': '16 h', 'Observed': 75, 'Chao1': 75.000, 'Shannon': 1.104, 'Simpson': 0.402, 'Evenness': 0.256, 'Fisher': 7.795, 'Faith_PD': 1.354}
]

df = pd.DataFrame(data)
metrics = ['Observed', 'Chao1', 'Shannon', 'Simpson', 'Evenness', 'Fisher', 'Faith_PD']

# -------------------------------------------------------------
# 2. STATISTICAL COMPUTATION ENGINE
# -------------------------------------------------------------
# A. Temporal Analysis (12 h vs 16 h)
temporal_stats = []
for m in metrics:
    v12 = df[df['Time'] == '12 h'][m].values
    v16 = df[df['Time'] == '16 h'][m].values
    
    # Paired Wilcoxon signed-rank test
    wide = df.pivot(index='Cultivar_Code', columns='Time', values=m)
    w_stat, w_p = stats.wilcoxon(wide['12 h'], wide['16 h'])
    
    # Mann-Whitney U test (unpaired)
    mwu_stat, mwu_p = stats.mannwhitneyu(v12, v16, alternative='two-sided')
    
    temporal_stats.append({
        'Metric': m,
        'Mean_12h': f"{np.mean(v12):.3f} ± {np.std(v12, ddof=1):.3f}",
        'Median_12h': f"{np.median(v12):.3f}",
        'Mean_16h': f"{np.mean(v16):.3f} ± {np.std(v16, ddof=1):.3f}",
        'Median_16h': f"{np.median(v16):.3f}",
        'Wilcoxon_W': f"{w_stat:.1f}",
        'Wilcoxon_P': f"{w_p:.4f}",
        'MWU_U': f"{mwu_stat:.1f}",
        'MWU_P': f"{mwu_p:.4f}"
    })

df_temporal = pd.DataFrame(temporal_stats)

# B. Cultivar Analysis (Across 6 Varieties)
cultivar_stats = []
dunn_results = {}

for m in metrics:
    groups = [group[m].values for name, group in df.groupby('Cultivar_Code')]
    kw_stat, kw_p = stats.kruskal(*groups)
    
    # Dunn's post-hoc with Benjamini-Hochberg FDR
    dunn = sp.posthoc_dunn(df, val_col=m, group_col='Cultivar_Code', p_adjust='fdr_bh')
    dunn_results[m] = dunn
    
    cultivar_stats.append({
        'Metric': m,
        'Kruskal_Wallis_H': f"{kw_stat:.4f}",
        'Degrees_of_Freedom': len(groups) - 1,
        'Kruskal_Wallis_P': f"{kw_p:.4f}"
    })

df_cultivar = pd.DataFrame(cultivar_stats)

# C. Variety Summary Table (Mean ± SD across both time points)
variety_summary = []
for code, group in df.groupby('Cultivar_Code'):
    c_name = group['Cultivar'].iloc[0]
    row = {'Cultivar': c_name, 'Code': code}
    for m in metrics:
        vals = group[m].values
        row[m] = f"{np.mean(vals):.2f} ± {np.std(vals, ddof=1):.2f}"
    variety_summary.append(row)

df_variety_summary = pd.DataFrame(variety_summary)

# -------------------------------------------------------------
# 3. SAVE EXCEL SPREADSHEET
# -------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..')) if os.path.basename(SCRIPT_DIR) == 'scripts' else SCRIPT_DIR
out_dir_stats = os.path.join(REPO_ROOT, 'results', 'tables')
os.makedirs(out_dir_stats, exist_ok=True)
excel_path = os.path.join(out_dir_stats, 'Table_S1_Alpha_Diversity_Statistical_Significance.xlsx')

with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    df.to_excel(writer, sheet_name='Individual_Samples', index=False)
    df_temporal.to_excel(writer, sheet_name='Temporal_12h_vs_16h', index=False)
    df_cultivar.to_excel(writer, sheet_name='Cultivar_Kruskal_Wallis', index=False)
    df_variety_summary.to_excel(writer, sheet_name='Cultivar_Summaries', index=False)
    for m, d_df in dunn_results.items():
        d_df.to_excel(writer, sheet_name=f'Dunn_FDR_{m[:10]}')

print(f"Excel workbook generated: {excel_path}")

# -------------------------------------------------------------

print('[Phase 1] Alpha Diversity analysis and Table S1 completed successfully.')
