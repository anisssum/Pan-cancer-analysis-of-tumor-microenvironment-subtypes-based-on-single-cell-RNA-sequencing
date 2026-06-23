import os
import glob
import re
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

OUTDIR = 'results_liana_cluster'
FIGDIR = os.path.join(OUTDIR, 'figures')
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

def find_liana_files(path):
    files = list(Path(path).rglob("*.csv"))
    files = [str(p) for p in files if re.search(r"liana_results_.*\.csv$", os.path.basename(str(p)))]
    return sorted(list(set(files)))

cluster_df = pd.read_csv('cluster_df_464.csv', index_col=0)
cluster_map = cluster_df.set_index('Sample')['clusters'].to_dict()

liana_files = []

for folder in ["Korean", "TISCH_new", "3CA"]:
    path = os.path.join("/home/anya/flamingo/", folder, "liana_results")
    liana_files.extend(find_liana_files(path))

def read_liana_csv(path):
    # read using pandas, ensure expected columns exist
    df = pd.read_csv(path)
    # add metadata: dataset from filename if present
    fname = os.path.basename(path)
    m = re.match(r'liana_results_(.+)\.csv', fname)
    dataset = m.group(1) if m else fname
    # ensure sample column name exists: user showed a sample column named 'sample' or 'sample' value in file
    if 'sample' in df.columns:
        df.rename(columns={'sample': 'sample_name'}, inplace=True)
    elif 'Sample' in df.columns:
        df.rename(columns={'Sample': 'sample_name'}, inplace=True)
    return df

frames = []
for f in liana_files:
    try:
        df = read_liana_csv(f)
        frames.append(df)
    except Exception as e:
        print('Failed to read', f, e)

liana_all = pd.concat(frames, ignore_index=True, sort=False)
print('Combined LIANA rows:', len(liana_all))
del frames

sample_list = pd.read_csv("/home/anya/flamingo/samples_464.txt", header=None)[0].tolist()

len(set(sample_list).intersection(set(liana_all['sample_name']))) #429
missing = set(sample_list) - set(liana_all['sample_name'])
with open('missing_samples.txt', 'w') as f:
    f.write('\n'.join(missing))

liana_all['cluster_label'] = liana_all['sample_name'].map(cluster_map)

liana_all['lr_id'] = (
    liana_all['source'].astype(str) + '||' + liana_all['target'].astype(str) + '||' +
    liana_all['ligand_complex'].astype(str) + '||' + liana_all['receptor_complex'].astype(str)
)
liana_all.to_csv(os.path.join(OUTDIR, 'liana_all.csv'))

sample_groups = liana_all.groupby('sample_name')
per_sample = []
for sample, g in sample_groups:
    total_sum = g['lrscore'].sum(skipna=True)
    mean = g['lrscore'].mean(skipna=True)
    med = g['lrscore'].median(skipna=True)
    n_pairs = g['lr_id'].nunique()
    n_sig = (g['cellphone_pvals'] <= 0.05).sum()
    n_non_na = g['lrscore'].notna().sum()
    dataset = g['__dataset'].iloc[0]
    cluster_label = g['cluster_label'].iloc[0] if 'cluster_label' in g.columns else None
    per_sample.append({
        'sample_name': sample,
        'dataset': dataset,
        'cluster_label': cluster_label,
        'total_lrscore_sum': total_sum,
        'mean_lrscore': mean,
        'median_lrscore': med,
        'n_pairs': n_pairs,
        'n_sig_pairs': n_sig,
        'n_observed_pairs': n_non_na
    })
per_sample_df = pd.DataFrame(per_sample).set_index('sample_name')
del per_sample
per_sample_df.to_csv(os.path.join(OUTDIR, 'per_sample_summary.csv'))

# Create a sender->receiver matrix per sample: sum lrscore for each source|target
# We'll create a wide feature matrix where columns are 'src__tgt' pairs (limited to top N most variable to avoid exploding dimensionality)
pair_scores = liana_all.copy()
pair_scores['src_tgt'] = pair_scores['source'].astype(str) + '___' + pair_scores['target'].astype(str)
pair_agg = pair_scores.groupby(['sample_name', 'src_tgt'])['lrscore'].sum().reset_index()
pair_pivot = pair_agg.pivot(index='sample_name', columns='src_tgt', values='lrscore').fillna(0)

# Join with per_sample_df
features_df = per_sample_df.join(pair_pivot, how='left').fillna(0)

# ---------------------------
# Unsupervised visualization
# ---------------------------
feat_mat = features_df.drop(columns=['dataset', 'cluster_label'])
scaler = StandardScaler()
X = scaler.fit_transform(feat_mat)

pca = PCA(n_components=10)
pcs = pca.fit_transform(X)

umapper = umap.UMAP(n_components=2, random_state=0)
umap_coords = umapper.fit_transform(X)

viz_df = pd.DataFrame({
    'sample_name': feat_mat.index,
    'PC1': pcs[:, 0],
    'PC2': pcs[:, 1],
    'UMAP1': umap_coords[:, 0],
    'UMAP2': umap_coords[:, 1],
    'cluster_label': features_df['cluster_label'].astype(str),
    'dataset': features_df['dataset'].astype(str)
}).set_index('sample_name')

sns.set(style='whitegrid', context='talk')
plt.figure(figsize=(12, 10))
sns.scatterplot(data=viz_df, x='UMAP1', y='UMAP2', hue='cluster_label')
plt.title('UMAP of per-sample LR features')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'umap_per_sample.png'), dpi=300)
plt.close()

# ---------------------------
# Differential testing per LR across cluster labels
# ---------------------------
# We'll test lrscore aggregated to sample-level for each lr_id. For each sample, define lrscore for that lr_id as summed lrscore across source/target in that sample
lr_sample = liana_all.groupby(['lr_id', 'sample_name'])['lrscore'].sum().reset_index()
lr_wide = lr_sample.pivot(index='sample_name', columns='lr_id', values='lrscore').fillna(0)
sample_cluster = per_sample_df['cluster_label'].astype(str)
clusters = sample_cluster.loc[lr_wide.index]

# For each lr_id perform Kruskal-Wallis test across groups (non-parametric) if more than one group, else skip
from collections import defaultdict
pvals = {}
statdict = {}
for lr in lr_wide.columns:
    groups = [lr_wide.loc[clusters[clusters == cl].index, lr].values
              for cl in clusters.unique()]
    # need at least two groups with >0 variance
    try:
        if len(groups) < 2:
            continue
        stat, p = stats.kruskal(*groups)
        pvals[lr] = p
        statdict[lr] = stat
    except Exception as e:
        pvals[lr] = np.nan
        statdict[lr] = np.nan

pvals_series = pd.Series(pvals).dropna()
reject, pvals_corrected, _, _ = multipletests(pvals_series.values, method='fdr_bh')
df_de = pd.DataFrame({
    'lr_id': pvals_series.index,
    'pvalue': pvals_series.values,
    'p_adj': pvals_corrected,
    'reject_null': reject,
    'stat': [statdict.get(lr, np.nan) for lr in pvals_series.index]
})

# add a human readable column by splitting lr_id
split_cols = df_de['lr_id'].str.split('\|\|', expand=True)
split_cols.columns = ['source', 'target', 'ligand_complex', 'receptor_complex']
df_de = pd.concat([df_de, split_cols], axis=1)
df_de = df_de.sort_values('p_adj')
df_de.to_csv(os.path.join(OUTDIR, 'lr_differential_tests_by_cluster.csv'), index=False)


TOP_HITS = df_de[df_de['reject_null'] & df_de['p_adj'] < 0.01]
TOP_HITS.to_csv(os.path.join(OUTDIR, 'lr_p_adj_0_01.csv'), index=False)

# ---------------------------
# Enrichment of senders/receivers
# ---------------------------
hit_lrs = set(TOP_HITS['lr_id'])
liana_all['is_top_hit'] = liana_all['lr_id'].isin(hit_lrs)
# contingency table per source cell type
source_counts = liana_all.groupby('source')['is_top_hit'].agg(['sum', 'count'])
source_counts['not_hit'] = source_counts['count'] - source_counts['sum']
# fisher exact test per source
enrichment = []
total_hits = liana_all['is_top_hit'].sum()
total_not = (~liana_all['is_top_hit']).sum()
for src, row in source_counts.iterrows():
    a = int(row['sum'])
    b = int(row['not_hit'])
    c = int(total_hits - a)
    d = int(total_not - b)
    # fisher test
    try:
        odds, p = stats.fisher_exact([[a, b], [c, d]])
    except Exception:
        p = np.nan
        odds = np.nan
    enrichment.append({'source': src, 'n_hit': a, 'n_total': int(row['count']), 'oddsratio': odds, 'p': p})
enr_df = pd.DataFrame(enrichment).sort_values('p')
enr_df['p_adj'] = multipletests(enr_df['p'].fillna(1), method='fdr_bh')[1]
enr_df.to_csv(os.path.join(OUTDIR, 'sender_enrichment_topLRs.csv'), index=False)

# ---------------------------
# Predictive modeling: can we predict cluster label from LR features?
# ---------------------------
# encode cluster label
features_df['cluster_label_enc'] = features_df['cluster_label'].astype('category').cat.codes
X = features_df.drop(columns=['dataset', 'cluster_label', 'cluster_label_enc'])
y = features_df['cluster_label_enc']

clf = RandomForestClassifier(n_estimators=200, random_state=0)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
acc_scores = cross_val_score(clf, X, y, cv=cv, scoring='accuracy')
print('Cross-validated accuracy:', acc_scores.mean(), '+-', acc_scores.std())

clf.fit(X, y)
importances = pd.Series(clf.feature_importances_, index=X.columns).sort_values(ascending=False)
importances.to_csv(os.path.join(OUTDIR, 'rf_feature_importances.csv'))

# For binary problems you can compute AUC; for multiclass compute per-class AUC if needed
# Save model metrics
with open(os.path.join(OUTDIR, 'model_metrics.txt'), 'w') as fh:
    fh.write('cv_accuracy_mean\tcv_accuracy_std\n')
    fh.write(f"{acc_scores.mean()}\t{acc_scores.std()}\n")

plt.figure(figsize=(8, 20))
importances.head(50)[::-1].plot(kind='barh')
plt.title('Random Forest feature importances')
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'rf_feature_importances.png'), dpi=300)
plt.close()

# ---------------------------
# Mixed-effects modeling example
# ---------------------------
top_lr = df_de['lr_id'].iloc[0]
mdf = lr_wide[[top_lr]].reset_index().merge(per_sample_df[['cluster_label', 'dataset']], left_on='sample_name', right_index=True)
mdf = mdf.dropna()
# rename columns for formula
mdf = mdf.rename(columns={top_lr: 'lrscore', 'cluster_label': 'cluster', 'dataset': 'dataset'})
try:
    model = smf.mixedlm('lrscore ~ C(cluster)', mdf, groups=mdf['dataset'])
    res = model.fit(reml=False)
    with open(os.path.join(OUTDIR, 'mixedlm_topLR.txt'), 'w') as fh:
        fh.write(res.summary().as_text())
    print('MixedLM fit and saved summary')
except Exception as e:
    print('MixedLM failed:', e)


per_sample_df.to_csv(os.path.join(OUTDIR, 'per_sample_summary.csv'))
features_df.to_csv(os.path.join(OUTDIR, 'features_per_sample_full.csv'))
viz_df.to_csv(os.path.join(OUTDIR, 'viz_coords.csv'))

def plot_top_lrs_boxplots(topn=12):
    top = df_de['lr_id'].head(topn).tolist()
    long = lr_sample[lr_sample['lr_id'].isin(top)].copy()
    long = long.merge(per_sample_df[['cluster_label']], left_on='sample_name', right_index=True)
    for lr in top:
        f = os.path.join(FIGDIR, f'box_{re.sub(r"[^0-9A-Za-z]", "_", lr)}.png')
        plt.figure(figsize=(6, 4))
        sns.boxplot(data=long[long['lr_id'] == lr], x='cluster_label', y='lrscore')
        sns.stripplot(data=long[long['lr_id'] == lr], x='cluster_label', y='lrscore', color='k', size=3, jitter=True, alpha=0.6)
        plt.title(lr)
        plt.tight_layout()
        plt.savefig(f, dpi=200)
        plt.close()

plot_top_lrs_boxplots(12)
