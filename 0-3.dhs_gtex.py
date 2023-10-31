

# regulatory potential of cis-eQTL varaints from DHS

# cis-regulatory elements only
#
from pybedtools import BedTool
import pandas as pd

# Load the DHS data into a pandas dataframe
dhs_df = pd.read_csv("/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
                     sep='\t'  )

# Load the cCREs data into a pandas dataframe
ccres_df = pd.read_csv("/mnt/d/sandbox/0.gxe/dhs/GRCh38-cCREs.bed.gz",
                       sep='\t', header = None )
ccres_df.columns = ['seqname1', 'start1', 'end1', "id1", "id2", 'epigenetic_mark']


# Convert the dataframes into BedTool objects
dhs_bed = BedTool.from_dataframe(dhs_df)
ccres_bed = BedTool.from_dataframe(ccres_df)

# Intersect the two BedTool objects
overlap = dhs_bed.intersect(ccres_bed, wa=True, wb=True)

# Convert overlapping regions back into a pandas dataframe and write into a new file
overlap_df = overlap.to_dataframe( names= list( dhs_df.columns) +  list(ccres_df.columns ) ] )
overlap_df[ ["seqname", "start", "end", "mean_signal", "epigenetic_mark", "component"]].to_csv("/mnt/d/sandbox/0.gxe/dhs/overlap_DHS_cCREs.bed", sep='\t', index=False)



#
# longest transcript
# Filter gene_biotype for specified categories and shortest length type

# Load the gene information
gene_biotype = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype.csv' )
filtered_gene_biotype = gene_biotype[(gene_biotype['gene_category'].isin(["PCG", "lncRNA", "ncRNA", "pseudogene"])) &
                                     (gene_biotype['Length type'] == 'lonest')]

# Drop 'Length type' column as it's no longer needed
filtered_gene_biotype.drop(columns='Length type', inplace=True)

# Set 'Gene stable ID' as index for easy access
filtered_gene_biotype.set_index('Gene stable ID', inplace=True)

filtered_gene_biotype.gene_category.value_counts()




import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from intervaltree import Interval, IntervalTree

# Load DHS data and create an interval tree
dhs_df = pd.read_csv('/mnt/d/sandbox/0.gxe/dhs/overlap_DHS_cCREs.bed.gz',
                     sep="\t", usecols=[0,1,2],  compression='gzip'  )

dhs_df = pd.read_csv("/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
                     sep='\t'  )


dhs_id = {}
dhs_tree = {chrom: IntervalTree() for chrom in dhs_df['seqname'].unique()}
for idx, row in dhs_df.iterrows():
    dhs_tree[row['seqname']].addi(row['start'], row['end'])
    dhs_id[ ( row["seqname"], row['start'], row['end']) ] = row["identifier"]


# Define path and get all the gz files
path = '/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/'
all_files = glob.glob(os.path.join(path, "*.v8.signif_variant_gene_pairs.txt.gz"))


# For each tissue, process the eQTL data
for file in all_files :
    tissue = os.path.basename(file).split('.')[0]
    tissue_data = pd.read_csv(file, sep='\t')
    tissue_data['chrom'] = tissue_data['variant_id'].str.split('_').str[0]
    tissue_data['pos'] = tissue_data['variant_id'].str.split('_').str[1].astype(int)
    tissue_data['slope'] = tissue_data['slope'] # .abs()

    tissue_data_list = []
    print( tissue )
    for index, row in tissue_data.iterrows():
        region = 'DHS' if len(dhs_tree[row['chrom']].at(row['pos'])) > 0 else 'Non-DHS'
        if region == "DHS":
            tissue_data_list.append({ 'chrom':row["chrom"], 'pos':row["pos"], 'region': region, 'slope': row['slope'], 'gene_id':row["gene_id"].split(".")[0] })

    # Convert the list to a DataFrame and write it to a csv file
    tissue_data_df = pd.DataFrame(tissue_data_list)
    tissue_data_df.to_csv("/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/" + f'{tissue}_data.csv', index=False)





import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Define an empty DataFrame to store mean and CI data
mean_ci_data = pd.DataFrame(columns=['tissue', 'region', 'mean_abs_slope', 'ci_low', 'ci_high'])

# Iterate over all the CSV files
for file in glob.glob('/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/*_data.csv.gz'):
    tissue = os.path.basename(file).split('_data')[0]
    tissue_data = pd.read_csv(file, compression = "gzip")
    print( tissue )
    for region in ['DHS', 'Non-DHS']:
        region_data = tissue_data[tissue_data['region'] == region]['slope'].abs()
        mean_abs_slope = np.mean(region_data)
        ci_low, ci_high = stats.norm.interval(0.95, loc=mean_abs_slope, scale=stats.sem(region_data))

        #mean_ci_data = mean_ci_data.append({'tissue': tissue, 'region': region, 'mean_abs_slope': mean_abs_slope, 'ci_low': ci_low, 'ci_high': ci_high}, ignore_index=True)
        mean_ci_row = pd.DataFrame({'tissue': [tissue],
                                    'region': [region],
                                    'mean_abs_slope': [mean_abs_slope],
                                    'ci_low': [ci_low],
                                    'ci_high': [ci_high]})

        mean_ci_data = pd.concat([mean_ci_data, mean_ci_row], ignore_index=True)





import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Create two separate dataframes for DHS and non-DHS
dhs_data = mean_ci_data[mean_ci_data['region'] == 'DHS'].reset_index(drop=True)
non_dhs_data = mean_ci_data[mean_ci_data['region'] == 'Non-DHS'].reset_index(drop=True)

# Convert the confidence interval values into error values
dhs_data['ci'] = dhs_data['ci_high'] - dhs_data['mean_abs_slope']
non_dhs_data['ci'] = non_dhs_data['ci_high'] - non_dhs_data['mean_abs_slope']

# Sort the dataframes by mean_abs_slope in descending order
dhs_data = dhs_data.sort_values('mean_abs_slope', ascending=False)
non_dhs_data = non_dhs_data.reindex(dhs_data.index)  # Match the order of the DHS data

# Get the mean values and the errors for both DHS and non-DHS
dhs_means = dhs_data['mean_abs_slope'].values
dhs_errors = dhs_data['ci'].values
non_dhs_means = non_dhs_data['mean_abs_slope'].values
non_dhs_errors = non_dhs_data['ci'].values

# Combine the mean values and errors for DHS and non-DHS
means = np.array([dhs_means, non_dhs_means])
errors = np.array([dhs_errors, non_dhs_errors])

# Create a barplot with error bars
# Create a barplot with error bars
fig, ax = plt.subplots(figsize=(20,10))
bar_width = 0.3
opacity = 0.8
index = np.arange(len(dhs_means))

# Choose colors, e.g., darkgreen and darkorange
bar1 = ax.bar(index, dhs_means, bar_width, yerr=dhs_errors, capsize=2, alpha=opacity, color='darkgreen', label='DHS')
bar2 = ax.bar(index+bar_width, non_dhs_means, bar_width, yerr=non_dhs_errors, capsize=2, alpha=opacity, color='darkorange', label='non-DHS')

ax.set_xlabel('Tissue')
ax.set_ylabel('Mean Abs Slope')
ax.set_title('Mean Abs Slope by Tissue and Region')
ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(dhs_data['tissue'].values, rotation=90)  # Assumes tissues are in the same order in both datasets
ax.legend()
ax.set_xlim([min(index)-1, max(index)+1])

fig.subplots_adjust(left=0.05)  # Adjusts the left margin, making it closer to y-axis
fig.tight_layout()

plt.show()




##
import seaborn as sns

# Stacking the data
#dhs_data['region'] = 'DHS'
#non_dhs_data['region'] = 'non-DHS'
combined_data = pd.concat([dhs_data, non_dhs_data])

fig, ax = plt.subplots(figsize=(12,8))

ax = sns.barplot(data=combined_data, x='tissue', y='mean_abs_slope', hue='region', errorbar=None ,  capsize=0.1, palette=['darkgreen', 'darkorange'], ax=ax)

x_coords = [p.get_x() + 0.5*p.get_width() for p in ax.patches]

ax.errorbar(data= combined_data , x=x_coords, y='mean_abs_slope', yerr='ci', ls='', lw=8, color='black')

plt.title('Mean Abs Slope by Tissue and Region')
plt.xticks(rotation=90)
plt.tight_layout()

plt.show()


#
# DHS density per gene in obese-up and obese-down
#


# RNA-Seq
# lean Vs. obese patients
rnaseq_sig_genes = pd.read_csv( "/mnt/d/sandbox/0.gxe/rnaseq/rnaseq_sig_up_down.csv", index_col = 0)

# gene biotype and annotation
gene_biotype_o = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_obese.csv' , low_memory = False )
gene_biotype_df = gene_biotype_o[ gene_biotype_o["Length type"] == "average"]


# slope vs. tissue specificity
# from

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from sklearn.manifold import TSNE


data = []
for file in glob.glob('/mnt/d/sandbox/0.gxe/dhs/gtex_evariant_id/*.csv.gz'):
#for file in glob.glob('/mnt/d/sandbox/0.gxe/dhs//mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/*.csv.gz'):
    df = pd.read_csv(file, compression='gzip')
    df = df[df['region'] == 'DHS']  # keep only eQTLs in DHS regions
    tissue = os.path.basename(file).split('_data')[0]  # infer tissue from filename
    df['tissue'] = tissue
    data.append(df)

df = pd.concat(data, ignore_index=True)

# biotype
genecate = ["PCG", "lncRNA", "ncRNA", "pseudogene"][0]

#
f_df = df[( df['gene_id'].isin( gene_biotype_df[ gene_biotype_df["gene_category"] == genecate ][ "Gene stable ID"] ) ) ]

# dhs count gene
genes_up = rnaseq_sig_genes[ rnaseq_sig_genes[ "obese_up"] > 0 ].index.tolist()
genes_down = rnaseq_sig_genes[ rnaseq_sig_genes[ "obese_down"] < 0 ].index.tolist()


# up & down in obese
f_df_obese_up = f_df[  f_df[ "gene_id"].isin ( genes_up ) ]
f_df_obese_down = f_df[  f_df[ "gene_id"].isin ( genes_down ) ]


# obatin number of DHS sites per gene in each tissue from f_df dataset
# where gene_id is for gene, dhs_id is for DHS, tissue is for tissue name.
# there are two sets of genes, genes_up and genes_down, that have a set of gene_ids
# visualize distribution of DHS counts for genes_up and genes_down in each tissue.

# Count the number of DHS sites per gene in each tissue
dhs_counts = f_df.groupby(['gene_id', 'tissue']).size().reset_index(name='dhs_count')

# Merge DHS counts with genes_up and genes_down subsets
dhs_counts_up = dhs_counts[dhs_counts['gene_id'].isin(genes_up)]
dhs_counts_down = dhs_counts[dhs_counts['gene_id'].isin(genes_down)]

# side-by-side violin plot
dhs_counts_up['type'] = 'up'
dhs_counts_down['type'] = 'down'
dhs_counts_combined = pd.concat([dhs_counts_up, dhs_counts_down])


# dhs count for down-regulated in obese
# Calculate mean DHS count for each tissue in 'dhs_counts_up'
mean_dhs_counts_up = dhs_counts_up.groupby('tissue')['dhs_count'].mean().reset_index()
mean_dhs_counts_down = dhs_counts_down.groupby('tissue')['dhs_count'].mean().reset_index()


# Order the tissues by the mean DHS count
ordered_tissues = mean_dhs_counts_up.sort_values('dhs_count', ascending = False )['tissue']

# Create a categorical type for tissues ordered by the mean DHS count
dhs_counts_combined['tissue'] = pd.Categorical(dhs_counts_combined['tissue'], categories=ordered_tissues, ordered=True)


# Create the violin plot
fig, ax = plt.subplots(figsize=(14, 5))
sns.violinplot(x='tissue', y='dhs_count', hue='type', split=True,
data=dhs_counts_combined, inner=None, scale='width', palette=[sns.color_palette("muted")[1], sns.color_palette("muted")[9] ],
  ax=ax, bw = 0.02 )
ax.set_title('DHS counts for up-regulated and down-regulated genes across tissues')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.ylim( 0, 300)

plt.show()


# mean difference in dhs counts betwen up- and down-regulated gene sets
#

pval = []
for tissue in dhs_counts_up[ "tissue"].unique():
    dhs_count_up_tissue = dhs_counts_up[ dhs_counts_up[ "tissue"] == tissue ]["dhs_count"]
    dhs_count_down_tissue = dhs_counts_down[ dhs_counts_down[ "tissue"] == tissue ]["dhs_count"]
    t, p = scipy.stats.mannwhitneyu( dhs_count_up_tissue, dhs_count_down_tissue)
    print( tissue, np.mean( dhs_count_up_tissue), np.mean( dhs_count_down_tissue), p )
    pval.append( [ tissue, np.mean( dhs_count_up_tissue), np.mean( dhs_count_down_tissue), p ] )


pval_df = pd.DataFrame( pval)
pval_df.columns = [ "tissue", "m(up)", "m(down)", "pvalue"]

from statsmodels.stats.multitest import multipletests

# Correct for multiple testing
_, corrected_p_values, _, _ = multipletests( pval_df["pvalue"].tolist() , method='fdr_bh')

pval_df[ "adjpval"] = corrected_p_values

pval_df.to_csv( "dhs_count_per_pcg.csv")

#
# DHS density per gene.
# significantly higher in up-regulated genes than down-regulated genes
#

# Count the number of DHS sites per gene in each tissue
dhs_counts = f_df.groupby(['dhs_id', 'tissue']).size().reset_index(name='dhs_count')
dhstissue = dhs_counts[ "dhs_id"].value_counts()

# DHS shared by all 49 tissues

#
# obtain genome sequence for those coordinates of DHS that are shared by 49 tissues and up-regulated in obese.
#
dhs_df = pd.read_csv("/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
                     sep='\t'  )


dhsid_49 = dhstissue[ dhstissue == 49].index.tolist()
dhs_df_49 = dhs_df[ dhs_df["identifier"].astype( str ).isin( dhsid_49) ]


bedtools getfasta -fi /mnt/d/sandbox/hg38/hg38.fa  -bed dhs_shared_49tissues.csv -fo output.fasta

chr17   45719800        45720020        17.594248       1.60956270819672        305     45719910        45719890    45719950        Stromal A
chr15   66423240        66423460        15.68598        0.933201        1       66423350        66423350        66423350        Musculoskeletal
chr21   46203543        46203780        21.98976        0.1605475       2       46203680        46203652        46203708        Musculoskeletal
chr6    32625600        32625860        6.271863        0.656307        1       32625710        32625710        32625710        Lymphoid
chr17   12992400        12992760        17.240477       0.286159        1       12992630        12992630        12992630        Neural
chr22   23913080        23913320        22.523081       0.194876        1       23913190        23913190        23913190        Lymphoid
chr20   45812049        45812440        20.739353       1.05576666666667        3       45812230        45812097        45812344        Organ devel. / renal
chr20   45810245        45810480        20.739333       1.8928386875    32      45810370        45810312        45810394        Myeloid / erythroid
chr7    100212750       100212880       7.66587 0.270331        1       100212810       100212810       100212810       Digestive
chr10   45848557        45848720        10.408287       0.631972444827586       29      45848630        45848610        45848670        Myeloid / erythroid




# dhs density Vs aFC Vs expression FC for significntly higher in obese.
# we have a set of up-regulated genes in obese adipose tissue.
# some of those genes have strong effect size in other tissue.
# what would be implications of this sort of observation?


# tissue breadth of DHSs
#






#
f_df['variant'] = f_df['chrom'].astype(str) + ':' + f_df['pos'].astype(str)
f_df['variant_count'] = f_df.groupby('variant')['tissue'].transform('nunique')
f_df['gene_count'] = f_df.groupby('gene_id')['tissue'].transform('nunique')


def categorize(count):
    return count
    if count == 1:
        return 'tissue-specific'
    elif 1 < count <= 5:
        return 'high tissue-specific'
    elif 5 < count <= 10:
        return 'moderate tissue-specific'
    elif 30 < count <= 48:
        return 'low tissue-specific'
    else:
        return 'not tissue-specific'


df['variant_category'] = df['variant_count'].apply(categorize)
df['gene_category'] = df['gene_count'].apply(categorize)

mean_slope_by_variant_category = df.groupby('variant_category')['slope'].mean()
mean_slope_by_gene_category = df.groupby('gene_category')['slope'].mean()

sns.boxplot(x='variant_category', y='slope', data=df)
sns.boxplot(x='gene_category', y='slope', data=df)




#

# Pivot the DataFrame to have genes as rows and tissues as columns
f_df[ "abs_slope"] = f_df[ "slope"].abs()
# Binarize slope values: present (1) or absent (0)
#df['abs_slope'] = df['abs_slope'].apply(lambda x: 1 if x != 0 else 0)

# Now pivot to create a gene by variant matrix
#df_80 = df[ df["abs_slope"] >= 0.569118  ]  # top 20%
df_gene_variant = f_df.pivot_table(index='gene_id', columns='tissue', values='abs_slope', aggfunc='max')
# Drop NA values
df_gene_variant = df_gene_variant.fillna(0)

#df_gene_max = df_gene_max.fillna(0)
tsne = TSNE(n_components=2, random_state=0, perplexity=50, max_iter=2000, early_exag_coeff=12, stop_lying_iter=1000 )
tsne_result = tsne.fit_transform(df_gene_variant)



# color genes
tissues = df_gene_variant.columns

gene_tissue_mapping = df.groupby('gene_id').apply(lambda group: group.set_index('tissue')['slope'].to_dict()).to_dict()

# Convert it to a dictionary with a single tissue (with maximum slope value) for each gene
gene_tissue_mapping = {gene: max(tissues.items(), key=lambda item: item[1])[0] for gene, tissues in gene_tissue_mapping.items()}

# Create a color palette with enough distinct colors for each tissue
#color_palette = plt.get_cmap('tab20') # or any other colormap with sufficient distinct colors
color_palette = sns.color_palette("Spectral", 49)
#sns.color_palette("Spectral", as_cmap=True)

tissues = df_gene_variant.columns
# Map tissues to colors
tissue_to_color = dict(zip(tissues, color_palette))

# Map genes to tissue colors
gene_colors = [tissue_to_color[gene_tissue_mapping[gene]] for gene in df_gene_variant.index]

fig, ax = plt.subplots(figsize=(15,15))
plt.title('tSNE projection of the gene by tissue dataset', fontsize=24)
# Finally, you can use the 'colors' list to color your t-SNE plot
plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c=gene_colors)
plt.show()



fig, ax = plt.subplots(figsize=(15,15))
plt.title('tSNE projection of the gene by tissue dataset', fontsize=24)
# Finally, you can use the 'colors' list to color your t-SNE plot
#plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c=gene_colors)

# Create a scatter plot
scatter = plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c=gene_colors)

# Create legend handles manually
legend_handles = [mpatches.Patch(color=tissue_to_color[tissue], label=tissue) for tissue in tissues]

# Add the legend to the plot
plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')


plt.show()

# Map tissues to colors
#tissue_to_color = dict(zip(tissues, color_palette))


# color_mapping = {tissue: color_palette(i) for i, tissue in enumerate(tissues)}

# Create a dictionary where keys are gene names and values are tissues
# Create a dictionary where keys are gene names and values are tissues from original df
# Create a dictionary where keys are gene names and values are dictionaries with tissues and max slope values from the original df

# Create color mapping for each gene in df_gene_variant based on gene_tissue_mapping
#gene_color_mapping = {gene: gene_colors[tissue] for gene, tissue in gene_tissue_mapping.items()}

# Create a color list for the genes in df_gene_variant in the correct order
#colors = [gene_color_mapping[gene] for gene in df_gene_variant.index]



# gene x variants #
# genes into PCA  show if necessary
# Assuming the dataframe is named df
#
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from sklearn.manifold import TSNE

data = []
for file in glob.glob('/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/*.csv.gz'):
    df = pd.read_csv(file, compression='gzip')
    df = df[df['region'] == 'DHS']  # keep only eQTLs in DHS regions
    tissue = os.path.basename(file).split('_data')[0]  # infer tissue from filename
    df['tissue'] = tissue
    data.append(df)

df = pd.concat(data, ignore_index=True)

f_df = df[( df['gene_id'].isin( filtered_gene_biotype[ filtered_gene_biotype["gene_category"] == "PCG" ][ "Gene stable ID"] ) ) ]

f_df['variant'] = f_df['chrom'].astype(str) + ':' + f_df['pos'].astype(str)
#f_df['variant_count'] = f_df.groupby('variant')['tissue'].transform('nunique')
#f_df['gene_count'] = f_df.groupby('gene_id')['tissue'].transform('nunique')

# mark duplicates as True except for the first occurrence
duplicates = f_df.duplicated(subset=['variant'], keep=False)

# only keep rows where 'variant' is duplicated
f_df = f_df[duplicates]

# Define a lambda function to filter groups with less than or equal to 2 rows
filter_func = lambda group: len(group) >= 2

# Group the DataFrame by 'variant' and apply the filter function
f_df = f_df.groupby('variant').filter(filter_func)

df_gene_variant = f_df.pivot_table(index='gene_id', columns='variant', values='slope', aggfunc='max')







import umap
import umap.umap_ as umap
import matplotlib.pyplot as plt


# Perform UMAP
reducer = umap.UMAP()
umap_result = reducer.fit_transform(df_gene_variant.fillna(0))

# Plot UMAP result
fig, ax = plt.subplots(figsize=(15,15))
plot.scatter( umap_result[:, 0], umap_result[:, 1], alpha=0.1, c = colors )
plt.title('UMAP projection of the gene by tissue dataset', fontsize=24)
plt.show()



from sklearn.manifold import TSNE

tsne = TSNE(n_components=2, random_state=0)
tsne_result = tsne.fit_transform(df_gene_variant)

# Plot the t-SNE result
fig, ax = plt.subplots(figsize=(12,12))
plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c=colors)
#sns.scatterplot( x = tsne_result[:, 0], y = tsne_result[:, 1] , c = colors )
plt.title('tSNE projection of the gene by tissue dataset', fontsize=24)

plt.show()



# Then, use PCA to visualize tissue specificity of genes
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

pca = PCA(n_components=2)
pca_result = pca.fit_transform(pivot_df)

# Create a DataFrame from PCA results for plotting
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
pca_df['gene_id'] = pivot_df.index

# Plot the results
plt.figure(figsize=(10, 10))
sns.scatterplot(x='PC1', y='PC2', data=pca_df)
plt.title('PCA of Gene Slope Values Across Tissues')
plt.show()


import matplotlib.pyplot as plt








# Process each file
import glob
import numpy as np

# Collecting all file paths
file_paths = glob.glob('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz')

# Lists to store slope data for each group
dhs_slopes = []
non_dhs_slopes = []

# Process each file
for file_path in file_paths:
    # Read the data
    df = pd.read_csv(file_path, sep="\t", usecols=[0,7] )
    # Process each variant
    for index, row in df.iterrows():
        # Parse the chromosomal position of the variant
        chr, pos, _, _ , _ = row['variant_id'].split("_")
        pos = int(pos)

        # Categorize the variant and store the slope
        overlapping_regions = dhs_tree[chr].at(pos )
        if overlapping_regions:
            dhs_slopes.append(np.abs(row['slope']))
        else:
            non_dhs_slopes.append(np.abs(row['slope']))



# compare distribution of slopes between DHS and Non-DHS
#
# Prepare a figure to plot boxplots
path = "/mnt/d/sandbox/0.gxe/dhs/gtex_slopes/"
all_files = glob.glob( os.path.join( path, "*.csv.gz" ) )


fig, ax = plt.subplots(figsize=(15, 10))

# Loop through the files
for file in all_files:
    # Read tissue data
    tissue_data = pd.read_csv(file, sep=',', compression='gzip')
    tissue_data['abs_slope'] = tissue_data['slope'].abs()  # Calculate absolute value of slope

    # Initialize lists to store slope values for DHS and non-DHS regions
    dhs_slopes = tissue_data['abs_slope'][tissue_data['region' ] == "DHS"]
    non_dhs_slopes = tissue_data['abs_slope'][tissue_data['region' ] == "Non-DHS"]

    # Convert the slope lists to dataframes
    dhs_df = pd.DataFrame({'Tissue': [os.path.basename(file).split('.')[0]] * len(dhs_slopes),
                           'Region': ['DHS'] * len(dhs_slopes), 'Abs Slope': dhs_slopes})
    non_dhs_df = pd.DataFrame({'Tissue': [os.path.basename(file).split('.')[0]] * len(non_dhs_slopes),
                               'Region': ['Non-DHS'] * len(non_dhs_slopes), 'Abs Slope': non_dhs_slopes})

    # Combine the dataframes and plot the boxplot
    combined_df = pd.concat([dhs_df, non_dhs_df])
    sns.boxplot(x='Tissue', y='Abs Slope', hue='Region', data=combined_df, ax=ax)





import os
import glob
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns

# Load DHS data and create an interval tree for quick overlap checking
dhs_data = pd.read_csv('/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz', sep='\t',usecols=[0,1,2] )
dhs_data.columns = ['chrom', 'start', 'end']
dhs_data['chrom'] = dhs_data['chrom'].apply(str)
dhs_tree = {chrom: IntervalTree() for chrom in dhs_data['chrom'].unique()}

for index, row in dhs_data.iterrows():
    dhs_tree[row['chrom']].addi(row['start'], row['end'])

# For each tissue, process the eQTL data
for file in glob.glob('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz'):
    tissue = os.path.basename(file).split('.')[0]
    tissue_data = pd.read_csv(file, sep='\t')
    tissue_data['chrom'] = tissue_data['variant_id'].str.split('_').str[0]
    tissue_data['pos'] = tissue_data['variant_id'].str.split('_').str[1].astype(int)
    tissue_data['slope'] = tissue_data['slope'] # .abs()

    tissue_data_list = []
    print( tissue )
    for index, row in tissue_data.iterrows():
        region = 'DHS' if len(dhs_tree[row['chrom']].at(row['pos'])) > 0 else 'Non-DHS'
        tissue_data_list.append({'tissue': tissue, 'chrom':row["chrom"], 'pos':row["pos"], 'region': region, 'slope': row['slope'], 'gene_id':row["gene_id"].split(".")[0] })

    # Convert the list to a DataFrame and write it to a csv file
    tissue_data_df = pd.DataFrame(tissue_data_list)
    tissue_data_df.to_csv("/mnt/d/sandbox/0.gxe/dhs/gtex_slope/" + f'{tissue}_data.csv', index=False)
