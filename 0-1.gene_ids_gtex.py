
# check average transcript length
in_dir= "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/" + ".v8.egenes.txt.gz"
gl = "/mnt/d/sandbox/0.gxe/data/mart_human.txt.gz"

import pandas as pd
import gzip
import os
import re
import glob
import scipy.stats as stats


# Load the gene information

# genes tested in GTEx
gene_biotype = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_GTEx.csv' , low_memory = False )

# genes defined in ensembl
gene_biotype = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_obese.csv' , low_memory = False )

gtex_directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"


# Initialize a dictionary to store average transcript lengths for each tissue
average_transcript_lengths_tissue = {}

# Specify the directory containing the gene expression files
# directory = '/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/'

# Iterate over the files in the directory
# Set the directory path
# all genes test for cis-eQTL

# Get a list of all files in the directory
files = glob.glob(os.path.join(gtex_directory, "*.v8.egenes.txt.gz"))

# significant eVariants
#files = glob.glob(os.path.join(gtex_directory, "*.signif_variant_gene_pairs.txt.gz"))

# Create an empty list to store the results
results = []
# Iterate over the files
# average_transcript_lengths = pd.read_csv( "../average_transcript_lengths_GTEx_common.csv" )

glist = [ "PCG", "lncRNA", "ncRNA", "pseudogene"]

average_lengths_per_tissue = []
gtype = glist[ 0 ]


gene_biotype = gene_biotype_gtex[gene_biotype_gtex["Length type"] == "shortest"]

for file in os.listdir(gtex_directory):
    # if file.endswith(".v8.egenes.txt.gz"):
    # v8.signif_variant_gene_pairs.txt.gz
    if file.endswith( ".v8.signif_variant_gene_pairs.txt.gz"):
        tissue_name = file.split(".")[0]  # get the tissue name from the filename

        gene_expression = pd.read_csv(gtex_directory + file, sep="\t", compression='gzip')
        gene_expression['gene_id'] = gene_expression['gene_id'].str.split('.').str[0]  # remove the version number

        #merged = pd.merge(average_transcript_lengths, gene_expression, left_on='Gene stable ID', right_on='gene_id')
        merged = pd.merge( gene_biotype[ gene_biotype["gene_category"] == gtype ] , gene_expression, left_on='Gene stable ID', right_on='gene_id')
        merged_re = merged.drop_duplicates( subset = ["Gene stable ID"] )
        tlen = merged_re[ "Transcript length (including UTRs and CDS)"]
        mean_length = tlen.mean()

        # calculate the confidence interval
        ci_low, ci_high = stats.t.interval(0.95, len( tlen)-1,
                                            loc=mean_length,
                                            scale=stats.sem(tlen))

        average_lengths_per_tissue.append((tissue_name, mean_length, ci_low, ci_high))



# Convert the list to a DataFrame
average_lengths_per_tissue = pd.DataFrame(average_lengths_per_tissue, columns=["Tissue", "Transcript Length", "CI Low", "CI High"])
# Convert the results to a dataframe
# Print the results
#
average_lengths_per_tissue_sorted = average_lengths_per_tissue.sort_values(by='Transcript Length', ascending=False)

# overall transcrip legnths across 49 tissues
import seaborn as sns
import matplotlib.pyplot as plt

# Set the figure size
plt.figure(figsize=(12, 5))

# Create the barplot with error bars
plot = sns.barplot(x='Tissue', y='Transcript Length', data=average_lengths_per_tissue_sorted,
                   yerr=average_lengths_per_tissue_sorted['CI High'] - average_lengths_per_tissue_sorted['CI Low'])

plt.xlim( 0, 500 )
# Rotate the x-axis labels for better readability
plot.set_xticklabels(plot.get_xticklabels(), rotation=90)

# Display the plot
plt.show()




#
#  DHS and cis-eQTL varinats stat overview
#

from scipy.stats import chi2_contingency

egene_path = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"
egene_files = glob.glob(dir_path + 'v8.egenes.txt.gz')

dir_path = '/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/'
csv_files = glob.glob(dir_path + '*_data.csv.gz')

eqtlv_results = []

bty = ["PCG", "lncRNA", "ncRNA", "pseudogene" ]

egene_49tissue = []


for file in csv_files:
    tissue = file.split('/')[-1].split('_data')[0]  # Get the tissue name from the filename

    # egene
    #
    egene_df = pd.read_csv( egene_path + tissue + '.v8.egenes.txt.gz' , sep = "\t") # , usecols = [0, ] )

    # significant gene sets
    egene_df_sig = egene_df[  egene_df["qval"] <= 0.05 ]
    log2fc_bottom = egene_df_sig[ ( egene_df_sig[ "log2_aFC"] <= -np.log2( 1.5) ) ][ "gene_id"].tolist()
    log2fc_top = egene_df_sig[ ( egene_df_sig[ "log2_aFC"] >= np.log2( 1.5) )  ][ "gene_id"].tolist()
    egene_df_sig[ abs( egene_df_sig[ "log2_aFC"]) >= np.log2( 1.5 )  ]
    break
    gene_bottom = set( [x.split(".")[0] for x in log2fc_bottom ] )
    gene_top = set( [x.split(".")[0] for x in log2fc_top ] )

    numgene = max( [len( gene_bottom) , len( gene_top)]  )

        #
    # pick up the same number of top geneset
    egene_df_sig = egene_df_sig.copy()
    egene_df_sig[ "abs log2_aFC"] = abs( egene_df_sig[ "log2_aFC"] )
    egene_df_sig = egene_df_sig.sort_values( by = [  "abs log2_aFC" ] )

    afc_threshold_n = egene_df_sig[ egene_df_sig ["abs log2_aFC" ] <  np.log2( 1.5)  ][ "abs log2_aFC" ].quantile( [0.5, 0.6, 0.7, 0.8, 0.9, 0.95] )

    log2fc_mid_pos = set( egene_df_sig[ ( egene_df_sig[ "log2_aFC"] > 0 ) & ( egene_df_sig[ "log2_aFC"] <= afc_threshold_n[ 0.95]  ) ] [ "gene_id"].tolist() ) # [0:len(log2fc_top)]
    log2fc_mid_neg = set( egene_df_sig[ ( egene_df_sig[ "log2_aFC"] < 0 ) & ( egene_df_sig[ "log2_aFC"] >= -afc_threshold_n[ 0.95]  ) ] [ "gene_id"].tolist() ) # [-len(log2fc_bottom):]
    gene_mid =set( [x.split(".")[0] for x in log2fc_mid_neg ] + [x.split(".")[0] for x in log2fc_mid_pos ]  )

    # egene_sig = [x.split(".")[0] for x in  egene_df[ "gene_id"].tolist()]
    #print(  tissue, len(egene_df),  egene_df[ ( egene_df[ "log2_aFC"] > 0 )][0:len(log2fc_top)]["log2_aFC"].max(),  egene_df[ (egene_df[ "log2_aFC"] < 0 ) ][-len(log2fc_bottom):]["log2_aFC"].min()  )
    #print( tissue, len( log2fc_bottom), len( log2fc_top), len( log2fc_mid))
    df = pd.read_csv(file)
    df_dhs = df[df['region'] == 'DHS']
    df_non_dhs = df[df['region'] == 'Non-DHS']

    #eqtlv_results[ tissue ] = []
    for bt in bty :
        btgid = set( gene_biotype[ gene_biotype["gene_category"] == bt  ] ['Gene stable ID'].tolist()  )

        top_dhs = df_dhs['gene_id'].isin( gene_top & btgid  )
        top_non_dhs = df_non_dhs['gene_id'].isin( gene_top & btgid  )

        bottom_dhs = df_dhs['gene_id'].isin( gene_bottom & btgid )
        bottom_non_dhs = df_non_dhs['gene_id'].isin( gene_bottom & btgid  )

        mid_dhs = df_dhs['gene_id'].isin( gene_mid & btgid )
        mid_non_dhs = df_non_dhs['gene_id'].isin( gene_mid & btgid  )

        top_dhs_counts = top_dhs.sum()
        top_non_dhs_counts = top_non_dhs.sum()
        bottom_dhs_counts = bottom_dhs.sum()
        bottom_non_dhs_counts = bottom_non_dhs.sum()
        mid_dhs_counts = mid_dhs.sum()
        mid_non_dhs_counts = mid_non_dhs.sum()

        table_top = [ [ top_dhs_counts, top_non_dhs_counts] ,
                 [mid_dhs_counts, mid_non_dhs_counts] ]

        table_bottom = [ [ bottom_dhs_counts, bottom_non_dhs_counts] ,
                         [mid_dhs_counts, mid_non_dhs_counts] ]
        try:
            chi2, p_top, dof, expected = chi2_contingency(table_top)
            chi2, p_bottom, dof, expected = chi2_contingency(table_top)
        except ValueError:
            p_top, p_bottom = 1, 1
        #print( tissue, bt, len( gene_top), len( gene_bottom), len( gene_mid) , len( btgid ) )
        eqtlv_results.append( [ tissue, bt,  top_dhs_counts, top_non_dhs_counts, bottom_dhs_counts, bottom_non_dhs_counts , mid_dhs_counts, mid_non_dhs_counts, p_top, p_bottom ] )

        #print( tissue, bt, [ top_dhs_counts, top_non_dhs_counts], [ bottom_dhs_counts, bottom_non_dhs_counts] , [mid_dhs_counts, mid_non_dhs_counts], p_top, p_bottom )



eqtlv_df = pd.DataFrame( eqtlv_results )
eqtlv_df.columns = [ "tissue", "gene type", "top dhs", "top non-dhs", "bottom dhs", "bottom non-dhs", "mid dhs", "mid non-dhs", "p(top)", "p(bottom)"]





import matplotlib.pyplot as plt
import numpy as np


# Assuming you have dataframes top_df and bottom_df with columns 'within_DHS', 'outside_DHS' and 'pvalue'
# eqtlv_df_p = eqtlv_df[ eqtlv_df["gene type"] == "PCG" ]
eqtlv_df_p = eqtlv_df[ eqtlv_df["gene type"] == "PCG" ]

eqtlv_df_p = eqtlv_df_p.copy()
fdr_corrected_pvals = multi.multipletests( eqtlv_df_p["p(bottom)"], method='fdr_bh')[1]
eqtlv_df_p[ "-log10pbottom"] = -np.log10( fdr_corrected_pvals )

fdr_corrected_pvals = multi.multipletests( eqtlv_df_p["p(top)"], method='fdr_bh')[1]
eqtlv_df_p[ "-log10ptop"] = -np.log10( fdr_corrected_pvals )

eqtlv_df_p[ "pvaluecut"] = eqtlv_df_p[ "-log10ptop"] >= -np.log10(0.05)
top_df = eqtlv_df_p.sort_values( by=['pvaluecut', "top dhs" ] )  # Sort by p-value
#bottom_df = eqtlv_df_p.sort_values('-log10ptop')


fig, ax1 = plt.subplots(figsize=(6,10))
# Plotting the bars
ax1.barh(bar_l_up, top_df['top non-dhs'],  label='Non-DHS Top', alpha=0.9, color=col_up[1])
ax1.barh(bar_l_up, top_df['top dhs'], left=top_df['top non-dhs'],  label='DHS Top', alpha=0.9, color=col_up[3])


#ax1.set_xlabel("Number of significant cis-eQTL variants")
#ax1.set_ylabel("Tissues")

ax1.set_yticks([])
ax2.set_yticks(tick_pos)
ax2.set_yticklabels( top_df[ "tissue"] ,rotation=0,  horizontalalignment='center' )

# Move the y-axis to the right
ax2.yaxis.tick_right()


# Creating the secondary y-axis and plotting the p-values
ax2 = ax1.twiny()
ax2.plot(top_df['-log10ptop'], tick_pos, marker='o', color=col_up[3], label='-log10(adjusted p-values)')
ax2.set_xlabel('-log10(Adjusted P-value)')
ax2.axvline( x = -np.log10( 0.05 ) , color = col_up[3] , linestyle = "--")

# Set the ticks to be first names
ax1.yaxis.labelpad = 30

ax1.set_yticks(tick_pos)
ax1.set_yticklabels(  top_df[ "tissue"] )#  ,  horizontalalignment='center'  )
# Adding the legend
fig.legend(loc="lower left", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
plt.ylim( 0, 50 )
plt.tight_layout()
plt.show()





fig, ax1 = plt.subplots(figsize=(6,10))
# Plotting the bars


ax1.barh(bar_l_up, top_df['bottom non-dhs'],  label='Non-DHS Up-regulated', alpha=0.9, color=col_dn[1])
ax1.barh(bar_l_up, top_df['bottom dhs'], left=top_df['bottom non-dhs'],  label='DHS Up-regulated', alpha=0.9, color=col_dn[3])


#ax1.set_xlabel("Number of significant cis-eQTL variants")
#ax1.set_ylabel("Tissues")

ax1.set_yticks([])
ax2.set_yticks(tick_pos)
ax2.set_yticklabels( top_df[ "tissue"] )

# Move the y-axis to the right
ax2.yaxis.tick_right()


# Creating the secondary y-axis and plotting the p-values
ax2 = ax1.twiny()
ax2.plot(bottom_df['-log10pbottom'], tick_pos, marker='o', color='red', label='-log10(adjusted p-values)')
ax2.set_xlabel('-log10(Adj.P-value)')
ax2.axvline( x = -np.log10( 0.05 ) , color = "red", linestyle = "--")

# Set the ticks to be first names
ax1.set_yticks(tick_pos)
ax1.set_yticklabels(  top_df[ "tissue"] )
# Adding the legend
fig.legend(loc="lower left", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
plt.ylim( 0, 50 )
plt.tight_layout()
plt.show()



# eGene sharing network

[ "log2_aFC",  "gene_id" ]


import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from matplotlib.colors import ListedColormap
# Import necessary libraries
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Construct a dictionary that maps each tissue to a dictionary of gene_ids and their log2_aFC values
egene_path = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"
egene_files = glob.glob(egene_path + 'v8.egenes.txt.gz')

gene_dict = {}
for filename in glob.glob(gtex_directory + "*v8.egenes.txt.gz"):
    tissue = os.path.basename(filename).split(".")[0]
    df = pd.read_csv(filename, sep='\t')
    df = df[( ( df['log2_aFC'] >= np.log2(1.5) ) | (df['log2_aFC'] <= -np.log2(1.5)) ) & ( df["qval"] <= 0.05) ]
    gene_dict[tissue] = df.set_index('gene_id')['log2_aFC'].to_dict()

# Construct a matrix where rows are tissues, columns are gene_ids, and values are log2_aFC
aFC_matrix = pd.DataFrame(gene_dict).transpose()

# Replace NaNs with 0s
aFC_matrix.fillna(0, inplace=True)

# Calculate hierarchical clustering
linkage_matrix = linkage(squareform(pdist(aFC_matrix)), method='complete')

# Create a dendrogram ordering to use for the heatmap
dendro = dendrogram(linkage_matrix, no_plot=True)
dendro_order = [aFC_matrix.index[i] for i in dendro['leaves']]

# Reorder the aFC_matrix according to the dendrogram
aFC_matrix = aFC_matrix.loc[dendro_order, :]

# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 10))
cluster_map = sns.clustermap(aFC_matrix,  col_cluster=False, cmap='bwr', center=0,
    figsize=(7, 15), cbar = True, cbar_pos=(0, .1, .01, .1))


plt.show()



# gene shring network

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from matplotlib.colors import ListedColormap


gene_dict = {}
for filename in glob.glob(gtex_directory + "*v8.egenes.txt.gz"):
    tissue = os.path.basename(filename).split(".")[0]
    df = pd.read_csv(filename, sep='\t')
    df[ "Gene stable ID"] = [ x.split(".")[0] for x in df[ "gene_id" ] ]
    df_sig = df[( ( df['log2_aFC'].abs() >= np.log2(1.5) )  ) & ( df["qval"] <= 0.05) ]  # | (df['log2_aFC'] <= -np.log2(1.5))
    gene_dict[tissue] = df_sig.set_index('Gene stable ID')['log2_aFC'].to_dict()
    # aFC threshold



# Construct a matrix for the number of shared genes between tissues
shared_genes_matrix = pd.DataFrame(0, index=gene_dict.keys(), columns=gene_dict.keys())
for tissue1 in gene_dict.keys():
    for tissue2 in gene_dict.keys():
        shared_genes_matrix.loc[tissue1, tissue2] = len(set( gene_dict[tissue1].keys()).intersection( set( gene_dict[tissue2].keys() )) )

# for log2_aFC >= 1.5
cluster_map = sns.clustermap(shared_genes_matrix,  cmap='icefire', center=0,
    figsize=(14, 14), cbar = True, cbar_pos=(0, .1, .01, .1))

# for log2_aFC <= -1.5
cluster_map = sns.clustermap(shared_genes_matrix,  cmap="icefire_r", center=0,
    figsize=(14, 14), cbar = True, cbar_pos=(0, .1, .01, .1))


# Calculate hierarchical clustering
linkage_matrix = linkage(squareform(pdist(binary_matrix.T)), method='complete')

# Create a dendrogram ordering to use for the heatmap
dendro = dendrogram(linkage_matrix, no_plot=True)
dendro_order = [binary_matrix.columns[i] for i in dendro['leaves']]

# Reorder the binary_matrix and shared_genes_matrix according to the dendrogram
binary_matrix = binary_matrix[dendro_order]
shared_genes_matrix = shared_genes_matrix.loc[dendro_order, dendro_order]

# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 10))
cmap = ListedColormap(['white', 'blue'])
ax.matshow(binary_matrix.T, cmap=cmap)  # We transpose the binary_matrix so tissues are on the y-axis

# Set the ticks to be tissue names
ax.set_yticks(range(len(binary_matrix.columns)))
ax.set_yticklabels(binary_matrix.columns, rotation=0)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

plt.show()

# Print the shared_genes_matrix
# print(shared_genes_matrix)





# func.
# show sharing genes for a given set of genes
# options :
# select: PCG, lncRNA, ncRNA, pseudogene
# select: strong positive effect, negative effect, or both
#
# clustering with give gene set


# Merge DEG and gene_biotype dataframes
deg_files = glob.glob("/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv")
deg_dfs = [pd.read_csv(file, sep='\t') for file in deg_files]
deg_df = pd.concat(deg_dfs, ignore_index=True)


#
# p-value < 0.05
# gene set : FC > log2(1.5) , FC < -log2(1.5), -log2(1.5) < FC < log2(1.5)
#

gene_dict = {} #
gene_dict["obese_down"]={}
gene_dict["obese_up"]={}
gene_dict["obese_mid"] = {}


gene_category = merged_df[ "gene_category"] == "PCG"

alldegsets = []
for file in deg_files:
        deg_df = pd.read_csv(file, sep='\t')
        sid = file.split( "DEG_")[1].split(".")[0]
        # Find the column that starts with "log2(fold change)"
        fold_change_col = [col for col in deg_df.columns if col.startswith('log2(fold change)')][0]
        pval_col = [col for col in deg_df.columns if col.startswith('-log10(Pvalue)')][0]
        # Merge DEG and gene_biotype dataframes
        merged_df = pd.merge(deg_df, gene_biotype_df, left_on='Symbol', right_on='Gene name', how='inner')
        # fold change convertion
        # lean vs obese to obese vs lean
        merged_df[ fold_change_col ] = -merged_df[ fold_change_col ]

        #obtain up-regulated and down-regulated genes in obese
        merged_df_sig_neg = merged_df[ ( merged_df[ fold_change_col] <= -np.log2( 1.5 ) )  ]  # | ( merged_df[ fold_change_col] <= -log2( 1.5 ) )
        merged_df_sig_pos = merged_df[ ( merged_df[ fold_change_col] >= np.log2( 1.5 ) )  ]  # | ( merged_df[ fold_change_col] <= -log2( 1.5 ) )
        merged_df_sig_mid = merged_df[ ( merged_df[ fold_change_col].abs() < np.log2( 1.5 ) )   ]
        #
        # up-regulated genes
        gene_dict["obese_up"].update( merged_df_sig_pos.set_index('Gene stable ID')[fold_change_col ].to_dict() )
        gene_dict["obese_down"].update( merged_df_sig_neg.set_index('Gene stable ID')[fold_change_col ].to_dict() )
        gene_dict["obese_mid"].update( merged_df_sig_mid.set_index('Gene stable ID')[fold_change_col ].to_dict() )
        gene_dict[ sid ] =  merged_df_sig_pos.set_index('Gene stable ID')[fold_change_col ].to_dict()
        alldegsets += merged_df_sig_neg[ "Gene stable ID"].tolist()
        print( sid, len( merged_df_sig_neg), len( merged_df_sig_pos), len( merged_df_sig_mid) ) # sum( merged_df_sig_neg[fold_change_col].abs() >= np.log2( 1.5) ), sum( merged_df_sig_neg[pval_col] > -np.log10(0.05)  )  )


# del redundent gene
redun_genes = set( gene_dict["obese_down"].keys() ) & set( gene_dict["obese_up"] )
alldegsets = set( alldegsets )



for g in redun_genes:
    del gene_dict["obese_down"][ g ]
    del gene_dict["obese_up"][ g ]

rnaseq_df = pd.DataFrame( gene_dict )
rnaseq_df.fillna( 0 )
rnaseq_df.to_csv( "/mnt/d/sandbox/0.gxe/rnaseq/rnaseq_sig_up_down.csv")



# Construct a dictionary that maps each tissue to a dictionary of gene_ids and their log2_aFC values
egene_path = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"
egene_files = glob.glob(egene_path + 'v8.egenes.txt.gz')

numgenes = []

for filename in glob.glob(gtex_directory + "*v8.egenes.txt.gz"):
    tissue = os.path.basename(filename).split(".")[0]
    df = pd.read_csv(filename, sep='\t')
    df_up = df[( ( df['log2_aFC'] >= np.log2(1.5) )  ) & ( df["qval"] <= 0.05) ]  #  | ( df['log2_aFC'] <= -np.log2(1.5) )
    df_down = df[( ( df['log2_aFC'] <= -np.log2(1.5) )  ) & ( df["qval"] <= 0.05) ]  #  | ( df['log2_aFC'] <= -np.log2(1.5) )
    df[ "gid"] = [ x.split(".")[0] for x in df["gene_id"].tolist() ]
    numgenes.append( [tissue, len( df["gid"] ) ])
    df_gid_deg = df[ df["gid"].isin( alldegsets ) ]
    gene_dict[tissue] = df_gid_deg.set_index('gid')['log2_aFC'].to_dict()
    if "Adipose" in tissue:
        break

adi_gene_up = set( [x.split(".")[0] for x in df_up[ "gene_id"].tolist() ] )
adi_gene_down = set( [x.split(".")[0] for x in df_down[ "gene_id"].tolist() ] )
obese_gene_up = set ( list( gene_dict["obese_up"].keys() ) )
obese_gene_down = set ( list( gene_dict["obese_down"].keys() ) )


numgenes = pd.DataFrame( numgenes )



# clustering based on gene's FC from expression or effect size

# Construct a matrix where rows are tissues, columns are gene_ids, and values are log2_aFC
aFC_matrix = pd.DataFrame(gene_dict).transpose()

# Replace NaNs with 0s
aFC_matrix.fillna(0, inplace=True)

# Calculate hierarchical clustering
linkage_matrix = linkage(squareform(pdist(aFC_matrix)), method='complete')

# Create a dendrogram ordering to use for the heatmap
dendro = dendrogram(linkage_matrix, no_plot=True)
dendro_order = [aFC_matrix.index[i] for i in dendro['leaves']]

# Reorder the aFC_matrix according to the dendrogram
aFC_matrix = aFC_matrix.loc[dendro_order, :]


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 8))
sns.clustermap(shared_genes_matrix, cmap='viridis', standard_scale=1)
plt.title("Number of Shared Genes Across Tissues")
plt.show()


# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 10))
cluster_map = sns.clustermap(aFC_matrix,  col_cluster=False, cmap='bwr', center=0,
    figsize=(7, 15), cbar = True, cbar_pos=(0, .1, .01, .1))


plt.show()


tname = numgenes[0].tolist()

# Construct a matrix for the number of shared genes between tissues
shared_genes_matrix = pd.DataFrame(0, index=gene_dict.keys(), columns=gene_dict.keys())
for tissue1 in gene_dict.keys():
    for tissue2 in gene_dict.keys():
        if tissue1 in tname:
            inx = tname.index( tissue1 )
            tnum = numgenes[ 1 ].tolist()[ inx ]
            shared_genes_matrix.loc[tissue1, tissue2] = len(set( gene_dict[tissue1].keys()).intersection( set( gene_dict[tissue2].keys() )) )/tnum
        else:
            shared_genes_matrix.loc[tissue1, tissue2] = len(set( gene_dict[tissue1].keys()).intersection( set( gene_dict[tissue2].keys() )) )


selectindex = [ 'obese_up', "obese_down", 'GSE110729_genes', 'GSE141432_lean_NGT_genes',
       'GSE141432_lean_T2D_genes', 'GSE156906_MHO_genes',
       'GSE156906_MUO_genes', 'GSE162653_genes', 'GSE165932_genes',
       'GSE179455_obese_genes', 'GSE205668_genes', 'GSE55008_genes']


s_genes = shared_genes_matrix["obese"]
s_genes = pd.DataFrame( s_genes )

fig, ax = plt.subplots(figsize=(8,4))

sns.barplot( x = np.arange(len( s_genes[ 11:] )), y = "obese",  data = s_genes[ 11:].sort_values( by= ["obese"] , ascending = False) , color = "grey" )
ax.set_xticks(range(len( s_genes[ 11:] ) ) )
ax.set_xticklabels( s_genes[ 11:].sort_values( by= ["obese"] , ascending = False ).index, rotation=90)

plt.show()

# for log2_aFC >= 1.5
cluster_map = sns.clustermap( shared_genes_matrix[selectindex],  cmap='icefire', center=0,
    figsize=(6, 14), cbar = True, cbar_pos=(0, .1, .01, .1))

# for log2_aFC <= -1.5
cluster_map = sns.clustermap(shared_genes_matrix,  cmap="icefire_r", center=0,
    figsize=(14, 14), cbar = True, cbar_pos=(0, .1, .01, .1))


#
# chi2-square table
# up-regulated genes , sharing genes with strong positive effect size
# down-regulated genes, sharing genes with strong negative effect size
#

# contingency tabel for Fisher's exact test
#











#
#  DEGs from obese, metabolic diseases studies
#
# /mnt/d/sandbox/0.gxe/rnaseq/*.tsv

degfiles  = glob.glob( "/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv" )
# in this dir, "/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv"
# 1. should iterate each file
# 2. column name starts with "log2(fold change)" but end with different names in different study
# 3. you also test significance of mean differences between two subsets.
# so we print out file, means, and p-value

import pandas as pd
import glob


# Read all DEG files
deg_files = glob.glob("/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv")
deg_dfs = [pd.read_csv(file, sep='\t') for file in deg_files]
deg_df = pd.concat(deg_dfs, ignore_index=True)

# Load the gene_biotype DataFrame
gene_biotype_o = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_obese.csv' , low_memory = False )
#gene_biotype_gtex = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_GTEx.csv' )
# glist = set( gene_biotype_gtex["Gene stable ID"].tolist()) & set( gene_biotype_o[ "Gene stable ID"].tolist() )

gene_biotype_df = gene_biotype_o[ gene_biotype_o["Length type"] == "average"]


# Merge DEG and gene_biotype dataframes
deg_files = glob.glob("/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv")
pngenes = []
posgenes = []
neggenes = []
for file in deg_files:
        deg_df = pd.read_csv(file, sep='\t')
        # Find the column that starts with "log2(fold change)"
        fold_change_col = [col for col in deg_df.columns if col.startswith('log2(fold change)')][0]

        # Merge DEG and gene_biotype dataframes
        merged_df = pd.merge(deg_df, gene_biotype_df, left_on='Symbol', right_on='Gene name', how='inner')
        #merged_df_filtered = merged_df[ merged_df["Gene stable ID"].isin( pnglist )]
        # Split the DataFrame into two subsets
        df_positive_fc = merged_df[merged_df[fold_change_col] > np.log2(1.5 ) ]
        df_negative_fc = merged_df[merged_df[fold_change_col] < -np.log2(1.5) ]
        posgenes += list( df_positive_fc[  "Gene name"].unique() )
        neggenes += list( df_negative_fc[  "Gene name"].unique() )


        # Calculate mean transcript length
        mean_length_positive_fc = df_positive_fc['Transcript length (including UTRs and CDS)'].mean()
        mean_length_negative_fc = df_negative_fc['Transcript length (including UTRs and CDS)'].mean()

        # Perform t-test
        t_statistic, p_value = scipy.stats.mannwhitneyu(df_positive_fc['Transcript length (including UTRs and CDS)'].dropna(),
                                         df_negative_fc['Transcript length (including UTRs and CDS)'].dropna())

        # Print results
        print(f"For file {file}:")
        print("Mean transcript length for positive fold change: ", mean_length_positive_fc)
        print("Mean transcript length for negative fold change: ", mean_length_negative_fc)
        print("P-value for the difference in means: ", p_value)
        print()




#
# type 3 diabetes caused by a high fat diet
# diabetic dementia

# control up (down in obese)
png_pos = pd.DataFrame( posgenes ).value_counts()
png_pos2 = png_pos # [ png_pos >= 2 ]

# up in obese (down in control)
png_neg = pd.DataFrame( neggenes ).value_counts()
png_neg2 = png_neg # [ png_neg >= 2 ]

pnglist_pos = [x[0] for x in png_pos2.index.to_list() ]
pnglist_neg = [x[0] for x in png_neg2.index.to_list() ]

obese_dn = set( pnglist_pos ) - set( pnglist_neg)
obese_up = set( pnglist_neg ) - set( pnglist_pos)

#
gbiotype_up =  gene_biotype_df[ gene_biotype_df[ "Gene name"].isin( obese_up ) ]
gbiotype_dn =  gene_biotype_df[ gene_biotype_df[ "Gene name"].isin( obese_dn ) ]

# count total number of genes
#print ( len( obese_dn), len( obese_up) )
#(1378, 2598)

# up
PCG           1138
lncRNA         211
pseudogene      45
ncRNA           32

#down
PCG           1138
lncRNA         211
pseudogene      45
ncRNA           32



# check gene length
from scipy.stats import ttest_ind

# Read the gene_biotype DataFrame
deg_files = glob.glob("/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv")
#gene_biotype_df = pd.read_csv('gene_biotype_file_path', sep='\t')
gene_biotype_o = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_obese.csv' , low_memory = False )
gene_biotype_df = gene_biotype_o[ gene_biotype_o["Length type"] == "average"]



gbiotype_up =  pd.read_csv( "/mnt/d/sandbox/0.gxe/result/1.gbiotype_up.csv")  # gene_biotype_df[ gene_biotype_df[ "Gene name"].isin( obese_up ) ]
gbiotype_dn =  pd.read_csv( "/mnt/d/sandbox/0.gxe/result/1.gbiotype_dn.csv")  #gene_biotype_df[ gene_biotype_df[ "Gene name"].isin( obese_dn ) ]

#
sns.histplot( gbiotype_up[gbiotype_up[ "gene_category"] == "PCG"] [ "Transcript length (including UTRs and CDS)"] , color =  sns.color_palette()[1])
sns.histplot( gbiotype_dn[gbiotype_dn[ "gene_category"] == "PCG"] [ "Transcript length (including UTRs and CDS)"] , color = sns.color_palette()[0])

#
sns.histplot( gbiotype_up[gbiotype_up[ "gene_category"] == "lncRNA"] [ "Transcript length (including UTRs and CDS)"] , color =  sns.color_palette()[1])
sns.histplot( gbiotype_dn[gbiotype_dn[ "gene_category"] == "lncRNA"] [ "Transcript length (including UTRs and CDS)"] , color = sns.color_palette()[0])


gcate = "pseudogene"
scipy.stats.mannwhitneyu( gbiotype_up[gbiotype_up[ "gene_category"] == gcate ] [ "Transcript length (including UTRs and CDS)"] ,
gbiotype_dn[gbiotype_dn[ "gene_category"] == gcate ] [ "Transcript length (including UTRs and CDS)"] )



#
# cis-eQTL eVariants denstiy across tissues
# DHS density - acrpss phenotypes
#
"""datasets of cis-eQTL for 49 tissues *data.csv.gz in
/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/
header goes like: chrom,pos,region,maf,slope,gene_id
where there are cis-eQTL variants for gene_id
region column has information where varait is in DHS or Non-DHS

The, we have DEG sets in gbiotype_up and gbiotype_dn.
column name for "Gene stable ID" has gene id of DEGs.

Firstly, given DEG sets for up and dn, you compute tissue specificity by counting number fo those variants
in DHS or Non-DHS across tissues.

work fine. we need to obtain variants counts in "DHS" and "Non-DHS" for up and down regulated DEGs.
from this 2 by 2 table, we conduct chi-squre test.
we iterate over tissues in from  /mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/
There are 49 dataframe that starts with tissue name and end with data.csv.gz in this dir,"""





import pandas as pd
import glob
import numpy as np

# Load DEG sets
gbiotype_up =  pd.read_csv( "/mnt/d/sandbox/0.gxe/result/1.gbiotype_up.csv")  # gene_biotype_df[ gene_biotype_df[ "Gene name"].isin( obese_up ) ]
gbiotype_dn =  pd.read_csv( "/mnt/d/sandbox/0.gxe/result/1.gbiotype_dn.csv")  #gene_biotype_df[ gene_biotype_df[ "Gene name"].isin( obese_dn ) ]

# Combine upregulated and downregulated DEG sets
all_degs = pd.concat([gbiotype_up['Gene stable ID'], gbiotype_dn['Gene stable ID']])

# Load all the cis-eQTL datasets
cis_eqtl_files = glob.glob('/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/*data.csv.gz')
dfs = [ pd.read_csv(f, usecols=['region', 'gene_id']) for f in cis_eqtl_files]
cis_eqtl_df = pd.concat(dfs)

# Separate into DHS and non-DHS subsets
dhs_df = cis_eqtl_df[cis_eqtl_df['region'] == 'DHS']
non_dhs_df = cis_eqtl_df[cis_eqtl_df['region'] == 'Non-DHS']

# Count the occurrences of each DEG's variants in DHS and non-DHS regions
dhs_counts = dhs_df['gene_id'].value_counts()
non_dhs_counts = non_dhs_df['gene_id'].value_counts()

# Compute the tissue specificity for each gene
tissue_specificity_dhs = dhs_counts.reindex(all_degs, fill_value=0)
tissue_specificity_non_dhs = non_dhs_counts.reindex(all_degs, fill_value=0)



# distirbuiton of number of significant eQTL varinats
#

from scipy.stats import chi2_contingency

dir_path = '/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/'
csv_files = glob.glob(dir_path + '*_data.csv.gz')

chi2_results = {}

bty = ["PCG", "lncRNA", "ncRNA", "pseudogene" ][ 0 ]

for file in csv_files:
    tissue = file.split('/')[-1].split('_data')[0]  # Get the tissue name from the filename
    df = pd.read_csv(file)
    df_dhs = df[df['region'] == 'DHS']
    df_non_dhs = df[df['region'] == 'Non-DHS']

    dhs_counts_up = df_dhs['gene_id'].isin( gbiotype_up[ gbiotype_up["gene_category"] == bty ] ['Gene stable ID']).sum()
    dhs_counts_dn = df_dhs['gene_id'].isin(gbiotype_dn[ gbiotype_dn["gene_category"] == bty ] ['Gene stable ID']).sum()
    non_dhs_counts_up = df_non_dhs['gene_id'].isin(gbiotype_up[ gbiotype_up["gene_category"] == bty ]['Gene stable ID']).sum()
    non_dhs_counts_dn = df_non_dhs['gene_id'].isin(gbiotype_dn[ gbiotype_dn["gene_category"] == bty ]['Gene stable ID']).sum()

    # Construct the 2x2 table
    table = [[dhs_counts_up, dhs_counts_dn],
             [non_dhs_counts_up, non_dhs_counts_dn]]

    #table = [[dhs_counts_up/num_counts_up, dhs_counts_dn/num_counts_dn],
    #         [non_dhs_counts_up/num_counts_up, non_dhs_counts_dn/num_counts_dn]]

    # Perform the chi-square test
    try:
        chi2, p, dof, expected = chi2_contingency(table)
        chi2_results[tissue] = (chi2, p, dhs_counts_up, dhs_counts_dn, non_dhs_counts_up, non_dhs_counts_dn )
    except ValueError:
        chi2_results[tissue] = ( 1 , 1, dhs_counts_up, dhs_counts_dn, non_dhs_counts_up, non_dhs_counts_dn )


    # Perform the chi-square test
    # try:
    #    chi2, p, dof, expected = chi2_contingency(table)
    #    chi2_results[tissue] = (chi2, p, dhs_counts_up/num_counts_up, dhs_counts_dn/num_counts_dn, non_dhs_counts_up/num_counts_up, non_dhs_counts_dn/num_counts_dn )
    #except ValueError:
    #    chi2_results[tissue] = ( 1 , 1, dhs_counts_up/num_counts_up, dhs_counts_dn/num_counts_dn, non_dhs_counts_up/num_counts_up, non_dhs_counts_dn/num_counts_dn )


# Convert the results to a DataFrame for easy viewing
chi2_df = pd.DataFrame(chi2_results).T
chi2_df.columns=['Chi2', 'p-value', 'dhs_up', 'dhs_dn', 'non_dhs_up', 'non_dhs_dn' ]


# print(chi2_df)
# r_up = chi2_df[ "dhs_up"]/chi2_df[ "non_dhs_up"]
# r_down = chi2_df[ "dhs_dn"]/chi2_df[ "non_dhs_dn"]


# motified code and have additional columns for  dhs_counts_up, dhs_counts_dn, non_dhs_counts_up, non_dhs_counts_dn in chi2_df.
# let visualize those counts per tissue with stake barplot with x-axis to be tissue types.
# for each tissue, side-by-side barplots with one for [ dhs_counts_up, non_dhs_counts_up] and the other [dhs_counts_dn,  non_dhs_counts_dn ]
#

# visualize
import matplotlib.pyplot as plt

# Prepare the data
# Sort the DataFrame by 'DHS Up' column in ascending order
chi2_df["signum"]  = chi2_df[ "p-value"] < 0.05/49
#chi2_df[ "dhs_up_ratio"] = chi2_df['dhs_up']/chi2_df['non_dhs_up']
chi2_df[ "dhs_dn_ratio"] = chi2_df['dhs_dn'] + chi2_df['non_dhs_dn']
chi2_df = chi2_df.sort_values(by='dhs_up', ascending=False)
chi2_df = chi2_df.sort_values(by=[ "signum", "dhs_up", 'dhs_dn_ratio'], ascending=False)


#
bar1 = chi2_df['dhs_up']
bar2 = chi2_df['non_dhs_up']
bar3 = chi2_df['dhs_dn']
bar4 = chi2_df['non_dhs_dn']

# Positions of the left bar-boundaries
bar_l = [i+1 for i in range(len(chi2_df['dhs_up']))]

# Positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i for i in bar_l]

# Create the total score for each participant
totals = [i+j+k+l for i,j,k,l in zip(bar1, bar2, bar3, bar4)]
totals = [1] *49


# Create the percentage of the total score the pre_score value for each participant was
bar1_ratio = [i / j * 1 for  i,j in zip(bar1, totals)]
bar2_ratio = [i / j * 1 for  i,j in zip(bar2, totals)]
bar3_ratio = [i / j * 1 for  i,j in zip(bar3, totals)]
bar4_ratio = [i / j * 1 for  i,j in zip(bar4, totals)]


# let's sort x-axis tissue names by DHS up in ascending order
# Create a bar plot, in position bar_l


# let's sort x-axis tissue names by DHS up in ascending order
# Create a bar plot, in position bar_l
fig, ax = plt.subplots(figsize=(10,6))

col_up = sns.color_palette("rocket", 5)
col_dn = sns.color_palette("mako", 5 )



bar_l_up = np.asarray( bar_l ) - 0.12
ax.bar(bar_l_up, bar1_ratio, bottom=bar2_ratio, label='DHS Up-regulated', alpha=0.9, color=col_up[3] , width=0.4)
ax.bar(bar_l_up, bar2_ratio, label='Non-DHS Up-regulated', alpha=0.9, color=col_up[1], width=0.4)

bar_l_dn = np.asarray( bar_l ) + 0.12
ax.bar(bar_l_dn , bar3_ratio, bottom=bar4_ratio, label='DHS Down-regulated', alpha=0.9, color=col_dn[3], width=0.4, align='edge')
ax.bar(bar_l_dn , bar4_ratio, label='Non-DHS Down-regulated', alpha=0.9, color=col_dn[1], width=0.4, align='edge')

ax.set_ylabel("Number of significant cis-eQTL variants")
ax.set_xlabel("Tissues")

# add p-value
ax2 = ax.twinx()
ax2.plot(bar_l_up, -np.log10( chi2_df['p-value']), marker='o', color=col_up[4], label='-log10(p-values)')
ax2.set_ylabel('-log10(P-value)')
ax2.axhline( y =  -np.log10( 0.05/49), linestyle= "--", color = col_up[4])


# Set the ticks to be first names

# Set the ticks to be first names
ax.set_xticks(bar_l_up)
ax.set_xticklabels(chi2_df.index, rotation=90)


# Adding the legend
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

plt.xlim([min(tick_pos)-0.5, max(tick_pos)+0.5])

# Add a legend
plt.legend(loc='upper right')
#plt.grid(True)
plt.show()








# Sort the DataFrame by 'DHS Up' column in ascending order
chi2_df[ "dhs_up_ratio"] = chi2_df['dhs_up']/chi2_df['non_dhs_up']
chi2_df = chi2_df.sort_values(by='dhs_counts_up', ascending=True)

# Now the rest of the code remains the same
# Prepare the data
bar1 = chi2_df['dhs_counts_up']
bar2 = chi2_df['non_dhs_counts_up']
bar3 = chi2_df['dhs_counts_dn']
bar4 = chi2_df['non_dhs_counts_dn']

# Positions of the left bar-boundaries
bar_l = [i+1 for i in range(len(chi2_df['dhs_counts_up']))]

# Positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i for i in bar_l]

# Create the total score for each participant
totals = [i+j+k+l for i,j,k,l in zip(bar1, bar2, bar3, bar4)]

# Create the percentage of the total score the pre_score value for each participant was
bar1_ratio = [i / j * 100 for  i,j in zip(bar1, totals)]
bar2_ratio = [i / j * 100 for  i,j in zip(bar2, totals)]
bar3_ratio = [i / j * 100 for  i,j in zip(bar3, totals)]
bar4_ratio = [i / j * 100 for  i,j in zip(bar4, totals)]

# Create a bar plot, in position bar_l
plt.bar(bar_l, bar1_ratio, label='DHS Up', alpha=0.9, color='#019600', width=0.4)
plt.bar(bar_l, bar2_ratio, bottom=bar1_ratio, label='Non-DHS Up', alpha=0.9, color='#3C5F5A', width=0.4)

plt.bar(bar_l, bar3_ratio, label='DHS Down', alpha=0.9, color='#219AD8', width=0.4, align='edge')
plt.bar(bar_l, bar4_ratio, bottom=bar3_ratio, label='Non-DHS Down', alpha=0.9, color='#F1C40F', width=0.4, align='edge')

# Set the ticks to be first names
plt.xticks(tick_pos, chi2_df.index, rotation=90)
plt.ylabel("Percentage")
plt.xlabel("Tissues")

# Let the borders of the graphic
plt.xlim([min(tick_pos)-0.5, max(tick_pos)+0.5])
plt.ylim(-10, 110)

# Add a legend
plt.legend(loc='upper left')
plt.grid(True)
plt.show()





#
fig, ax1 = plt.subplots(figsize=(10,6))

# Plotting the bars
ax1.bar(bar_l_up, bar1_ratio, bottom=bar2_ratio, label='DHS Up-regulated', alpha=0.9, color=col_up[3] , width=0.85)
ax1.bar(bar_l_up, bar2_ratio, label='Non-DHS Up-regulated', alpha=0.9, color=col_up[1], width=0.85)

ax1.set_ylabel("Number of significant cis-eQTL variants")
ax1.set_xlabel("Tissues")

# Creating the secondary y-axis and plotting the p-values
ax2 = ax1.twinx()
ax2.plot(bar_l_up, -np.log10( chi2_df['p-value']), marker='o', color=col_up[2], label='-log10(p-values)')
ax2.set_ylabel('-log10(P-value)')
ax2.axhline( y =  -np.log10( 0.05/49), linestyle= "--", color = col_up[2])

# Set the ticks to be first names
ax1.set_xticks(bar_l_up)
ax1.set_xticklabels(chi2_df.index, rotation=90)


# Adding the legend
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

plt.xlim([min(tick_pos)-0.5, max(tick_pos)+0.5])
plt.tight_layout()
plt.show()

















# miRNA target

import os
import pandas as pd
import glob

gene_biotype = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype.csv' , low_memory = False )

gene_stable_ids = gene_biotype[gene_biotype["Gene type"] == "miRNA" ] ["Gene stable ID"].tolist()

# Directory path
gtex_directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"

# Create an empty list to hold all matching gene_ids from all dataframes
all_matching_ids = []

# Iterate over each file in the directory
for filename in glob.glob(os.path.join(gtex_directory, '*.v8.egenes.txt.gz')):
    # Load the dataframe
    df = pd.read_csv(filename, compression='gzip', header=0, sep='\t')

    # Filter the dataframe
    df['gene_id'] = df['gene_id'].str.split('.', expand=True)[0]
    matching_ids = df[df['gene_id'].isin(gene_stable_ids)]

    # Append the gene_ids to the all_matching_ids list
    all_matching_ids.extend(matching_ids['gene_id'].tolist())



all_mirna = gene_biotype[ gene_biotype["Gene stable ID"].isin( all_matching_ids ) ][ "Gene name"].tolist()
all_mirna = [x for x in all_mirna if pd.notna(x)]

# Filter out None or nan values


import pandas as pd

gb_df_dn = pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_dn.csv')  # adjust as necessary
gb_df_up =  pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_up.csv')

# DEGs in obese studies
obese_df = pd.read_csv( "/mnt/d/sandbox/0.gxe/rnaseq/rnaseq_sig_up_down.csv", index_col = 0)
obese_gene_up = obese_df[ obese_df["obese_up"] >= np.log2( 1.5)].index.tolist()
obese_gene_down = obese_df[ obese_df["obese_down"] <= -np.log2( 1.5)].index.tolist()


#
# mirna_df = pd.read_csv('/mnt/d/sandbox/0.gxe/data/miRDB_v6.0_prediction_result.txt.gz', sep='\t', header=None, names=['miRNA_ID', 'mRNA_ID', "score"])
mirna_hsa_df = pd.read_csv('/mnt/d/sandbox/0.gxe/data/miRDB_v6.0_hsa.txt', sep='\t', header=None, names=['miRNA_ID', 'mRNA_ID', "score"])
mirna_hsa_sig_df = mirna_hsa_df[ mirna_hsa_df["score"] >= 80]

gene_df = pd.read_csv('/mnt/d/sandbox/0.gxe/data/ens_nm_mart_export.txt.gz', sep='\t' ) #, header=None, names=['Gene_ID', 'mRNA_ID'])

mirna_all = ['hsa-' + x.lower().replace('mir', 'miR-') for x in all_mirna]
miRNAs_b35 =  [s + "-3p" for s in set( mirna_all) ] + [s + "-5p" for s in set( mirna_all) ]
mirna_all35 = set( list( mirna_all) + miRNAs_b35 )

# down miRNAs in obese
miRNAs_down_genename = gene_biotype_df[ gene_biotype_df['Gene stable ID'].isin(obese_gene_down) & ( gene_biotype_df["Gene type"]== 'miRNA')]['Gene name'].tolist()
miRNAs_down = ['hsa-' + x.lower().replace('mir', 'miR-') for x in miRNAs_down_genename]  # 16
miRNAs_b35_down =  [s + "-3p" for s in set( miRNAs_down) ] + [s + "-5p" for s in set( miRNAs_down) ]
mirna_dn35 = set( list( miRNAs_down) + miRNAs_b35_down )

#up miRNAs in obese
miRNAs_up_genename = gene_biotype_df[ gene_biotype_df['Gene stable ID'].isin(obese_gene_up) & ( gene_biotype_df["Gene type"]== 'miRNA')]['Gene name'].tolist()
miRNAs_up = ['hsa-' + x.lower().replace('mir', 'miR-') for x in miRNAs_up_genename]  #23
miRNAs_b35_up =  [s + "-3p" for s in set( miRNAs_up) ] + [s + "-5p" for s in set( miRNAs_up) ]
mirna_up35 = set ( list( miRNAs_up) + miRNAs_b35_up  )


# all targets
targets_all = mirna_hsa_sig_df[mirna_hsa_sig_df['miRNA_ID'].isin(mirna_all35)]['mRNA_ID'].tolist()
genes_all = mirna_hsa_sig_df[mirna_hsa_sig_df['RefSeq mRNA ID'].isin(targets_all)]['Gene stable ID'].tolist()

targets_all_dn_up =  mirna_hsa_sig_df[mirna_hsa_sig_df['miRNA_ID'].isin( set( mirna_all35 ) - set( miRNAs_dn35) -set( miRNAs_up35)  )]['mRNA_ID'].tolist()
genes_all_dn_up = mirna_hsa_sig_df[mirna_hsa_sig_df['RefSeq mRNA ID'].isin(targets_all_dn_up)]['Gene stable ID'].tolist()

# donw regualted miRNAs' target
targets_dn = mirna_hsa_sig_df[mirna_hsa_sig_df['miRNA_ID'].isin(miRNAs_dn35)]['mRNA_ID'].tolist()
genes_bypass = gene_df[gene_df['RefSeq mRNA ID'].isin(targets_dn)]['Gene stable ID'].tolist()

# up regualted miRNAs' target
targets_up = mirna_hsa_sig_df[mirna_hsa_sig_df['miRNA_ID'].isin(miRNAs_up35)]['mRNA_ID'].tolist()#
genes_silence = gene_df[gene_df['RefSeq mRNA ID'].isin(targets_up)]['Gene stable ID'].tolist()


pd.DataFrame( genes_bypass, columns =['Gene stable ID']).to_csv( "obese_dn_miRNA_target_gene.csv",  index=False)
pd.DataFrame( genes_silence, columns =['Gene stable ID']).to_csv( "obese_up_miRNA_target_gene.csv", index=False)

# genes bypass miRNA silencing
# remains PCG only all other types, lncRNA, pseudogenes filtered out
#
# filtered_up_ex_df = gbiotype_up[ ~gbiotype_up['Gene stable ID'].isin( set( genes_silence) )]

filtered_up_df = gbiotype_up[gbiotype_up['Gene stable ID'].isin( set( genes_bypass) - set( genes_silence) )]
filtered_up_df.to_csv( "gbiotype_up_non_miRNA_target.csv")

# number of genes sharing between targets





#
# GTEx with up-regualted genes passing miRNAs targeting
#
from scipy.stats import chi2_contingency

dir_path = '/mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/'
csv_files = glob.glob(dir_path + '*_data.csv.gz')

chi2_results = {}

bty = ["PCG", "lncRNA", "ncRNA", "pseudogene" ][ 0 ]

for file in csv_files:
    tissue = file.split('/')[-1].split('_data')[0]  # Get the tissue name from the filename
    df = pd.read_csv(file)
    df_dhs = df[df['region'] == 'DHS']
    df_non_dhs = df[df['region'] == 'Non-DHS']

    dhs_counts_up = df_dhs['gene_id'].isin( filtered_up_df[ filtered_up_df["gene_category"] == bty ] ['Gene stable ID']).sum()
    dhs_counts_dn = df_dhs['gene_id'].isin( gbiotype_dn[ gbiotype_dn["gene_category"] == bty ] ['Gene stable ID']).sum()
    non_dhs_counts_up = df_non_dhs['gene_id'].isin(filtered_up_df[ filtered_up_df["gene_category"] == bty ]['Gene stable ID']).sum()
    non_dhs_counts_dn = df_non_dhs['gene_id'].isin(gbiotype_dn[ gbiotype_dn["gene_category"] == bty ]['Gene stable ID']).sum()

    # Construct the 2x2 table
    table = [[dhs_counts_up, dhs_counts_dn],
             [non_dhs_counts_up, non_dhs_counts_dn]]

    # Perform the chi-square test
    try:
        chi2, p, dof, expected = chi2_contingency(table)
        chi2_results[tissue] = (chi2, p, dhs_counts_up, dhs_counts_dn, non_dhs_counts_up, non_dhs_counts_dn )
    except ValueError:
        chi2_results[tissue] = ( 1 , 1, dhs_counts_up, dhs_counts_dn, non_dhs_counts_up, non_dhs_counts_dn )


# Convert the results to a DataFrame for easy viewing
chi2_df = pd.DataFrame(chi2_results).T
chi2_df.columns=['Chi2', 'p-value', 'dhs_up', 'dhs_dn', 'non_dhs_up', 'non_dhs_dn' ]

chi2_df.to_csv( "eVariant_density_chi2_non_miRNA_target_pcg.csv")

# print(chi2_df)
# r_up = chi2_df[ "dhs_up"]/chi2_df[ "non_dhs_up"]
# r_down = chi2_df[ "dhs_dn"]/chi2_df[ "non_dhs_dn"]


# motified code and have additional columns for  dhs_counts_up, dhs_counts_dn, non_dhs_counts_up, non_dhs_counts_dn in chi2_df.
# let visualize those counts per tissue with stake barplot with x-axis to be tissue types.
# for each tissue, side-by-side barplots with one for [ dhs_counts_up, non_dhs_counts_up] and the other [dhs_counts_dn,  non_dhs_counts_dn ]
#

# visualize
import matplotlib.pyplot as plt

# Prepare the data
# Sort the DataFrame by 'DHS Up' column in ascending order

chi2_df["signum"]  = chi2_df[ "p-value"] < 0.05/49

chi2_df[ "dhs_up_ratio"] = chi2_df['dhs_up']/chi2_df['non_dhs_up']
chi2_df["signum"]  = chi2_df[ "p-value"] < 0.05/49

chi2_df = chi2_df.sort_values(by=[ 'signum', 'dhs_up'], ascending=False)

#
bar1 = chi2_df['dhs_up']
bar2 = chi2_df['non_dhs_up']
bar3 = chi2_df['dhs_dn']
bar4 = chi2_df['non_dhs_dn']

# Positions of the left bar-boundaries
bar_l = [i+1 for i in range(len(chi2_df['dhs_up']))]

# Positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i for i in bar_l]

# Create the total score for each participant
totals = [i+j+k+l for i,j,k,l in zip(bar1, bar2, bar3, bar4)]
totals = [1] *49


# Create the percentage of the total score the pre_score value for each participant was
bar1_ratio = [i / j * 1 for  i,j in zip(bar1, totals)]
bar2_ratio = [i / j * 1 for  i,j in zip(bar2, totals)]
bar3_ratio = [i / j * 1 for  i,j in zip(bar3, totals)]
bar4_ratio = [i / j * 1 for  i,j in zip(bar4, totals)]


# let's sort x-axis tissue names by DHS up in ascending order
# Create a bar plot, in position bar_l

# let's sort x-axis tissue names by DHS up in ascending order
# Create a bar plot, in position bar_l
fig, ax1 = plt.subplots(figsize=(12,8))

col_up = sns.color_palette("rocket", 5)
col_dn = sns.color_palette("mako", 5 )

#bar_l_up = np.asarray( bar_l ) # - 0.12
#plt.bar(bar_l_up, bar2_ratio, label='Within DHS', alpha=0.9, color=col_up[1], width=0.9 )
#plt.bar(bar_l_up, bar1_ratio, bottom=bar2_ratio, label='Outside of DHS', alpha=0.9, color=col_up[3] , width=0.9)

bar_l_up = np.asarray( bar_l ) #- 0.12

fig, ax1 = plt.subplots(figsize=(10,6))

# Plotting the bars
ax1.bar(bar_l_up, bar1_ratio, bottom=bar2_ratio, label='DHS Up-regulated', alpha=0.9, color=col_up[3] , width=0.85)
ax1.bar(bar_l_up, bar2_ratio, label='Non-DHS Up-regulated', alpha=0.9, color=col_up[1], width=0.85)

ax1.set_ylabel("Number of significant cis-eQTL variants")
ax1.set_xlabel("Tissues")

# Creating the secondary y-axis and plotting the p-values
ax2 = ax1.twinx()
ax2.plot(bar_l_up, -np.log10( chi2_df['p-value']), marker='o', color=col_up[2], label='-log10(p-values)')
ax2.set_ylabel('-log10(P-value)')
ax2.axhline( y =  -np.log10( 0.05/49), linestyle= "--", color = col_up[2])

# Set the ticks to be first names
ax1.set_xticks(bar_l_up)
ax1.set_xticklabels(chi2_df.index, rotation=90)


# Adding the legend
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

plt.xlim([min(tick_pos)-0.5, max(tick_pos)+0.5])
plt.tight_layout()
plt.show()








# log aFC
#
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from intervaltree import Interval, IntervalTree

fig, ax = plt.subplots(figsize=(12,8))

sns.histplot( afc_up.abs() , label = "Up-regulated PGC" , color = sns.color_palette()[3] , alpha = 0.9 ) # , stat = "density" )
sns.histplot( afc_dn.abs() , label = "Donw-regulated PGC" , color = sns.color_palette()[9] , alpha = 0.9 ) # , stat = "density"  )
sns.histplot( afc_mirna.abs(), bins = 50 , label = "Non-miRNA target PGC", color = sns.color_palette()[1] , alpha= 0.9) #, stat = "density" )

ax.legend()
plt.show()


#
# transcript length for total up-regulated genes Vs up-regulated genes passing miRNA silencing
#
gbiotype_dn = pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_dn.csv')  # adjust as necessary
gbiotype_up =  pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_up.csv')
filtered_up_df = pd.read_csv( "/mnt/d/sandbox/0.gxe/result/gbiotype_up_pass_miRNA.csv")
totalup_except = set(gbiotype_up[ "Gene name"] ) - set( filtered_up_df[ "Gene name"] )


nt = "PCG"

total_except = set( gbiotype_up["Gene name"] ) - set( filtered_up_df["Gene name"] )
gbiotype_up_except =  gbiotype_up[ gbiotype_up[ "Gene name" ].isin( list( totalup_except) ) ]

sns.histplot( gbiotype_up_except[ gbiotype_up_except[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"] , color = sns.color_palette()[3], label='Up-regulated PCGs', alpha = 0.9 ) # ,  stat="density" ) # , bins = 40 )
sns.histplot( gbiotype_dn[ gbiotype_dn[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"] , color = sns.color_palette()[9], label='Down-regulated PCGs', alpha = 0.9) # ,  stat="density" ) # , bins = 40 )
sns.histplot( filtered_up_df[ filtered_up_df[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"] , color =  sns.color_palette()[1], label="Non-miRNA target PCGs", alpha = 0.9) # ,  stat="density", bins = 25)
plt.legend()


# test differences in gene length
nt = "PCG"
up-regulated PCG Vs non-miRNA targ PCGs

print( gbiotype_up_except[gbiotype_up_except[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"].mean() ,
filtered_up_df[filtered_up_df[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"].mean() )
print( scipy.stats.mannwhitneyu( gbiotype_up_except[gbiotype_up_except[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"] ,
 filtered_up_df[ filtered_up_df[ "gene_category" ] ==  nt ] [ "Transcript length (including UTRs and CDS)"] ) )

print( gbiotype_dn[gbiotype_dn[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"].mean() ,
filtered_up_df[filtered_up_df[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"].mean() )
print( scipy.stats.mannwhitneyu( gbiotype_up_except[gbiotype_up_except[ "gene_category"] == nt ] [ "Transcript length (including UTRs and CDS)"] ,
 gbiotype_dn[ filtered_up_df[ "gene_category" ] ==  nt ] [ "Transcript length (including UTRs and CDS)"] ) )


# slope distribuitons
# in down-regulated PCGs
# in up-regulated genes except non-miRNA target PCG
# in non-miRNA target PCG
# eVariants density per gene in donw, up, non-miRNA target
import pandas as pd
from intervaltree import Interval, IntervalTree
from scipy.stats import chi2_contingency
import pandas as pd
import glob
import gzip


dhs_df = pd.read_csv("/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
                     sep='\t'  )


dhs_tree = {chrom: IntervalTree() for chrom in dhs_df['seqname'].unique()}
for idx, row in dhs_df.iterrows():
    dhs_tree[row['seqname']].addi(row['start'], row['end'])



gbiotype_dn = pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_dn.csv')  # adjust as necessary
gbiotype_up =  pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_up.csv')
filtered_up_df = pd.read_csv( "/mnt/d/sandbox/0.gxe/result/gbiotype_up_pass_miRNA.csv")


dir_path = '/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/'
csv_files = glob.glob(dir_path + '*.v8.egenes.txt.gz')

bty = ["PCG", "lncRNA", "ncRNA", "pseudogene" ][ 0 ]
abs_table = {}
for file in csv_files:
    tissue = file.split('/')[-1].split('.v8.')[0]  # Get the tissue name from the filename
    df = pd.read_csv(file, sep = "\t")
    df[ "gene_id_v"] = [ x.split(".")[0] for x in df[ "gene_id"] ]
    #df_dhs = df[df['region'] == 'DHS']
    #df_non_dhs = df[df['region'] == 'Non-DHS']
    afc_up = df[ df["gene_id_v"].isin( gbiotype_up[ gbiotype_up["gene_category"] == bty ] ['Gene stable ID'] ) ][ "log2_aFC" ]
    afc_dn = df[ df["gene_id_v"].isin( gbiotype_dn[ gbiotype_dn["gene_category"] == bty ] ['Gene stable ID'] ) ][ "log2_aFC" ]
    afc_mirna = df[ df["gene_id_v"].isin( filtered_up_df[ filtered_up_df["gene_category"] == bty ] ['Gene stable ID'] ) ][ "log2_aFC" ]

    # check DHS
    print( tissue, "DHS", np.mean( afc_up.abs() ), np.mean( afc_dn.abs() ), np.mean( afc_mirna.abs() ),
        scipy.stats.mannwhitneyu( afc_up.abs()  ,   afc_dn.abs()  ),
        scipy.stats.mannwhitneyu(  afc_up.abs() ,  afc_mirna.abs()  )



import os
import pandas as pd

# Load gene ID subsets
df_dn = pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_dn.csv')  # adjust as necessary
df_up =  pd.read_csv('/mnt/d/sandbox/0.gxe/result/gbiotype_up.csv')
df_miRNA = pd.read_csv( "/mnt/d/sandbox/0.gxe/result/gbiotype_up_pass_miRNA.csv")


df_mirna =  set( filtered_up_df[ filtered_up_df["gene_category"] == bty ] ['Gene stable ID'] )
df_up =  set( gbiotype_up[ gbiotype_up["gene_category"] == bty ] ['Gene stable ID'] ) -  df_mirna
df_dn = set( gbiotype_dn[ gbiotype_dn["gene_category"] == bty ] ['Gene stable ID'] )


# Lists to store results
results_up = []
results_dn = []
results_miRNA = []

# Directory with eGene files
dir_path = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"

# Iterate over eGene files
import os
import pandas as pd
from intervaltree import Interval, IntervalTree

# Load gene ID subsets
# Lists to store results
results_up = []
results_dn = []
results_miRNA = []

# Directory with eGene files
dir_path = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"

# Construct the interval tree for DHS regions
dhs_tree = IntervalTree()
for index, row in DHS_df.iterrows():  # Assuming DHS_df is the DataFrame containing DHS regions
    dhs_tree[row['start']:row['end']] = row['seqname']

dhs_tree = {chrom: IntervalTree() for chrom in dhs_df['seqname'].unique()}
for idx, row in dhs_df.iterrows():
    dhs_tree[row['seqname']].addi(row['start'], row['end'])


# Iterate over eGene files
for filename in os.listdir(dir_path):
    if filename.endswith(".v8.egenes.txt.gz"):
        tissue = filename.split('.')[0]  # Get tissue name from filename
        egenes_df = pd.read_csv(os.path.join(dir_path, filename), sep="\t")
        egenes_df['gene_id'] = egenes_df['gene_id'].apply(lambda x: x.split('.')[0])  # Remove version number from gene_id

        # For each subset, filter egenes_df and count DHS variants
        for subset_df, results in [(df_up, results_up), (df_dn, results_dn), (df_miRNA, results_miRNA)]:
            filtered_df = egenes_df[egenes_df['gene_id'].isin(subset_df['GeneID'])]
            dhs_count = 0
            for index, row in filtered_df.iterrows():
                if row['chrom'] in dhs_trees and dhs_trees[row['chrom']].search(row['variant_pos']):
                    dhs_count += 1
            non_dhs_count = len(filtered_df) - dhs_count
            results.append((tissue, dhs_count, non_dhs_count))

        break

# Convert results to DataFrames
results_up_df = pd.DataFrame(results_up, columns=['Tissue', 'DHS_Variants', 'Non_DHS_Variants'])
results_dn_df = pd.DataFrame(results_dn, columns=['Tissue', 'DHS_Variants', 'Non_DHS_Variants'])
results_miRNA_df = pd.DataFrame(results_miRNA, columns=['Tissue', 'DHS_Variants', 'Non_DHS_Variants'])




# ... (initial setup as before)

# Construct the interval tree for DHS regions
dhs_trees = {}  # One IntervalTree per chromosome
for index, row in dhs_df.iterrows():
    # If we haven't seen this chromosome before, create a new IntervalTree for it
    if row['seqname'] not in dhs_trees:
        dhs_trees[row['seqname']] = IntervalTree()
    # Add the DHS region to the appropriate IntervalTree
    dhs_trees[row['seqname']][row['start']:row['end']] = True

# Iterate over eGene files
for filename in os.listdir(dir_path):
    if filename.endswith(".v8.egenes.txt.gz"):
        # ... (same as before)

        # For each subset, filter egenes_df and count DHS variants
        for subset_df, results in [(df_up, results_up), (df_dn, results_dn), (df_miRNA, results_miRNA)]:
            filtered_df = egenes_df[egenes_df['gene_id'].isin(subset_df['GeneID'])]
            dhs_count = 0
            for index, row in filtered_df.iterrows():
                if row['chrom'] in dhs_trees and dhs_trees[row['chrom']].search(row['variant_pos']):
                    dhs_count += 1
            non_dhs_count = len(filtered_df) - dhs_count
            results.append((tissue, dhs_count, non_dhs_count))

# ... (convert results to DataFrames as before)


###  gene classification  ##
# 0: gene biotype DataFrame
# each group has subgroups of PCG, lncRNA, ncRNA, pseudogene (gene_category)
# Gene stable ID, Gene name, transcript length (including UTRs and CDS), Chrom, Gene start, Gene end
#
# group 1: aggregated obese genes set
# non-redundent DEGs in obese patients from 10 comparison
# aggregated, but excluded genes in HEG set and LEG set
# a) highly expressed genes in obese aggregated from 10 comparisons
# b) lowerly expressed genes in obese aggregated from 10 comparisons
#
# group 2: GTEx gene set
# strong eGene set in tissue
# tissue-specificity for each of 49 tissue
# genes having cis-eQTL variants associated with strong est aFC (effect size)
# a) genes having cis-eQTL variants with aFC > 1.5 for each tissue
# b) genes having cis-eQTL variants with aFC < -1.5 for each tissue
# c) 95% qunatile of cis-eQTL variants as reference
#
# group 3: DHS set bearing strong and weak cis-eQTL
# DHS having strongest cis-eQTL and its associated genes (group 2)
# a) DHS for group 2-a) genes
# b) DHS for group 2-b) genes
# c) DHS not having cis-eQTL showing tissue-spcificity : not a) and not b)
#
# group 4: DHS set for obese genes
# DHS having significant cis-eQTLs and its associated genes being obese gene set (group 1)
# a) highly expressed genes in obese
# b) lowely expressed genes in obese
#
# group 5: cis-eQTL variants for DHS sets
# cis-eQTL within group
# a: associated genes wth significant & large positive effect size (3a DHS set)
# b: associated genes with significant & large negative effect size (3b DHS set)
# c: associated genes with significant effect size genes, not 3a nor 3b (3c DHS set)
# d: associated genes highly expressed in obese (4a DHS set)
# e: associated genes lowly expressed in obese (4b DHS set)
#
# group 7: Polygenic Risk Scores
# a) variants with strong positive PRS (top 5%)
# b) variants with strong negative PRS (bottom 5%)
# c) variants with weak PRS ( 10~90% quantile)
#
# group 8: PRS associaed with eQTL-DHS
# 7a, 7b, 7c for 5a~5e DHS
#
# group 9: recent rare variants for ethnistic groups from gnomAD
# a) rare varaints within group 5a~5e
# b) common varaints within group 5a~5e
# c) LoF variants within group 5a~5e
#

# obese genes filtration
# aggregate first.
# then filter in genes that have one or more significant cis-eQTL variants in a least one tissue












    break
    # filtering MAF < 1%

    dhs_slopes_up = df_dhs[ df_dhs['gene_id'].isin( gbiotype_up_except[ gbiotype_up_except["gene_category"] == bty ] ['Gene stable ID']) ]
    dhs_slopes_dn = df_dhs[ df_dhs['gene_id'].isin( gbiotype_dn[ gbiotype_dn["gene_category"] == bty ] ['Gene stable ID']) ]
    non_dhs_slopes_up = df_non_dhs[ df_non_dhs['gene_id'].isin(gbiotype_up_except[ gbiotype_up_except["gene_category"] == bty ]['Gene stable ID']) ]
    non_dhs_slopes_dn = df_non_dhs[ df_non_dhs['gene_id'].isin(gbiotype_dn[ gbiotype_dn["gene_category"] == bty ]['Gene stable ID']) ]
    dhs_slopes_up_mirna = df_dhs[ df_dhs['gene_id'].isin( filtered_up_df[ "Gene stable ID"]  ) ]
    non_dhs_slopes_up_mirna = df_non_dhs[ df_non_dhs['gene_id'].isin( filtered_up_df[ "Gene stable ID"] ) ]

    s0, s1, s2 = dhs_slopes_up_mirna[ "slope"].abs() , dhs_slopes_up[ "slope"].abs() , dhs_slopes_dn[ "slope"].abs()
    print( tissue, "DHS", np.mean( s2[ s2 > 0] ), np.mean( s1[ s1 > 0] ), np.mean( s0[ s0 > 0 ] ),
        scipy.stats.mannwhitneyu( s2[ s2 > 0] , s1[ s1 > 0 ] )[1], scipy.stats.mannwhitneyu( s0[ s0 > 0] , s1[ s1 > 0 ] )[1] ,  )

    #print( tissue,  np.mean( s2[ s2 < 0] ), np.mean( s1[ s1 < 0] ), np.mean( s0[ s0 < 0 ] ),
    #    scipy.stats.mannwhitneyu( s2[ s2 < 0] , s1[ s1 < 0 ] )[1], scipy.stats.mannwhitneyu( s0[ s0 < 0] , s1[ s1 < 0 ] )[1]  )

    s0, s1, s2 = non_dhs_slopes_up_mirna[ "slope"].abs() , non_dhs_slopes_up[ "slope"].abs() , non_dhs_slopes_dn[ "slope"].abs()
    print( tissue, "Non-DHS",  np.mean( s2[ s2 > 0] ), np.mean( s1[ s1 > 0] ), np.mean( s0[ s0 > 0 ] ),
        scipy.stats.mannwhitneyu( s2[ s2 > 0] , s1[ s1 > 0 ] )[1], scipy.stats.mannwhitneyu( s0[ s0 > 0] , s1[ s1 > 0 ] )[1]  )


    print( dhs_slopes_up[ "gene_id"].value_counts().mean(),
      dhs_slopes_dn[ "gene_id"].value_counts().mean(),
      dhs_slopes_up_mirna[ "gene_id"].value_counts().mean(),
      non_dhs_slopes_up[ "gene_id"].value_counts().mean(),
      non_dhs_slopes_dn[ "gene_id"].value_counts().mean(),
      non_dhs_slopes_up_mirna[ "gene_id"].value_counts().mean() )

    #print( tissue,  np.mean( s2[ s2 < 0] ), np.mean( s1[ s1 < 0] ), np.mean( s0[ s0 < 0 ] ),
    #    scipy.stats.mannwhitneyu( s2[ s2 > 0] , s1[ s1 > 0 ] )[1], scipy.stats.mannwhitneyu( s0[ s0 > 0] , s1[ s1 > 0 ] )[1]  )

    break










# Iterate over all DEG files
degall = []
deg_files = glob.glob("/mnt/d/sandbox/0.gxe/rnaseq/*_genes.tsv")

for file in deg_files:
    # Load the DEG DataFrame
    deg_df = pd.read_csv(file, sep='\t')

    # Find the column that starts with "log2(fold change)"
    fold_change_col = [col for col in deg_df.columns if col.startswith('log2(fold change)')][0]

    # Merge DEG and gene_biotype dataframes
    merged_df = pd.merge(deg_df, gene_biotype_df, left_on='Symbol', right_on='Gene name', how='inner')
    merged_df = merged_df[ merged_df[ "gene_category"] == "lncRNA" ]
    degall = degall + list( merged_df[ "Gene name"])
    merged_df_filtered = merged_df[ merged_df[ "Gene stable ID"].isin( pnglist ) ]
    # Split the DataFrame into two subsets
    df_positive_fc = merged_df_filtered[merged_df_filtered[fold_change_col] > 0]
    df_negative_fc = merged_df_filtered[merged_df_filtered[fold_change_col] < 0]

    # Calculate mean transcript length
    mean_length_positive_fc = df_positive_fc['Transcript length (including UTRs and CDS)'].mean()
    mean_length_negative_fc = df_negative_fc['Transcript length (including UTRs and CDS)'].mean()

    # Perform t-test
    t_statistic, p_value = ttest_ind(df_positive_fc['Transcript length (including UTRs and CDS)'].dropna(),
                                     df_negative_fc['Transcript length (including UTRs and CDS)'].dropna())

    # Print results
    if p_value <= 0.05:
        print(f"For file {file}:", len( df_positive_fc), len( df_negative_fc) )
        print("Mean transcript length for positive fold change: ", mean_length_positive_fc)
        print("Mean transcript length for negative fold change: ", mean_length_negative_fc)
        print("P-value for the difference in means: ", p_value)
        print()


#
# gene epxression breadth
#

ncount = pd.DataFrame( degall)[0].value_counts()
genesmorethan2 = ncount[ ncount>1].index
pnglist = [x[0] for x in png.index.to_list() ]


lean vs ( obese, Diabetic Obese, obeseNGT(normal glucose tolerance), obeseT2D,
non-diabetic obese Vs. Diabetic obese
metabolic healthy obese, metabolic disease obese


DEG_GSE151760_genes obese vs lean
GSE179455 T2D vs lean
GSE92724 cont vs diab
GSE106289 pre-weight loss vs pose-weright lose


DEG_GSE106289.tsv                 DEG_GSE141432_lean_T2D_genes.tsv  DEG_GSE179455_obese_genes.tsv
DEG_GSE110729.tsv                 DEG_GSE145412_genes.tsv           DEG_GSE185957_obese_post_genes.tsv
DEG_GSE121344_obeseD_genes.tsv    DEG_GSE151760_genes.tsv           DEG_GSE185957_obese_pre_genes.tsv
DEG_GSE121344_obeseND_genes.tsv   DEG_GSE156906_MHO_genes.tsv       DEG_GSE205668_genes.tsv
DEG_GSE126169_genes.tsv           DEG_GSE156906_MUO_genes.tsv       DEG_GSE55008.tsv
DEG_GSE129398_genes.tsv           DEG_GSE161042_genes.tsv           DEG_GSE92724.tsv
DEG_GSE132831_genes.tsv           DEG_GSE162653_genes.tsv           DEG_GSE95640.tsv
DEG_GSE137631_prerygb_genes.tsv   DEG_GSE165932_genes.tsv           arc
DEG_GSE141432_lean_NGT_genes.tsv  DEG_GSE179455_T2D_genes.tsv



# shortest
# Filter gene_biotype for specified categories and shortest length type
filtered_gene_biotype = gene_biotype[(gene_biotype['gene_category'].isin(["PCG", "lncRNA", "ncRNA", "pseudogene"])) &
                                     (gene_biotype['Length type'] == 'shortest')]

# Drop 'Length type' column as it's no longer needed
filtered_gene_biotype.drop(columns='Length type', inplace=True)

# Set 'Gene stable ID' as index for easy access
filtered_gene_biotype.set_index('Gene stable ID',v inplace=True)

# test difference btw [ "Whole_Blood"] and others in PCG
# Load Whole_Blood egenes

filtered_gene_biotype_pcg = filtered_gene_biotype[ filtered_gene_biotype["gene_category"] == "lncRNA" ]


whole_blood_egenes = pd.read_csv('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz', sep='\t')
whole_blood_egenes['gene_id'] = whole_blood_egenes['gene_id'].str.split('.').str[0]

# Filter for genes present in filtered_gene_biotype

whole_blood_egenes = whole_blood_egenes[ whole_blood_egenes['gene_id'].isin( filtered_gene_biotype_pcg.index)]

# Add 'Transcript length' to whole_blood_egenes
whole_blood_egenes['Transcript length'] = whole_blood_egenes['gene_id'].map(filtered_gene_biotype_pcg['Transcript length (including UTRs and CDS)'])


import glob
import scipy.stats as stats

# Get list of eQTL files
eqtl_files = glob.glob('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/*.v8.egenes.txt.gz')

# Remove Whole_Blood from list
eqtl_files.remove('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz')

# Initialize results list
results = []


adipose_egenes = pd.read_csv('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt.gz', sep='\t')
adipose_egenes['gene_id'] = adipose_egenes['gene_id'].str.split('.').str[0]
adipose_egenes = adipose_egenes[adipose_egenes['gene_id'].isin(filtered_gene_biotype_pcg.index)]
adipose_egenes['Transcript length'] = adipose_egenes['gene_id'].map( filtered_gene_biotype_pcg['Transcript length (including UTRs and CDS)'])

# Iterate over files
for file in eqtl_files:
    # Load egenes
    tissue_egenes = pd.read_csv(file, sep='\t')
    tissue_egenes['gene_id'] = tissue_egenes['gene_id'].str.split('.').str[0]

    # Filter for genes present in filtered_gene_biotype
    tissue_egenes = tissue_egenes[tissue_egenes['gene_id'].isin(filtered_gene_biotype_pcg.index)]

    # Add 'Transcript length' to tissue_egenes
    tissue_egenes['Transcript length'] = tissue_egenes['gene_id'].map( filtered_gene_biotype_pcg['Transcript length (including UTRs and CDS)'])

    # Perform t-test between Whole_Blood and tissue
    t_stat, p_val = stats.ttest_ind(whole_blood_egenes['Transcript length'], tissue_egenes['Transcript length'])

    # Add results to list
    results.append({'Tissue': file.split('/')[-1].split('.v8.egenes.txt.gz')[0], 'T-Stat': t_stat, 'P-Value': p_val })


from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Prepare the DataFrame (egene_data)
# ...
eqtl_files = glob.glob('/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/*.v8.egenes.txt.gz')

# Calculate p-values for each pair of tissues

# Calculate p-values for each pair of tissues
p_values = {}

for i in range(len(eqtl_files)):
    for j in range(i+1, len(eqtl_files)):
        file1 = eqtl_files[ i ]
        tissue_name1 = file1.split("/")[-1].split(".")[0]
        tissue_egenes1 = pd.read_csv(file1, sep='\t')
        tissue_egenes1['gene_id'] = tissue_egenes1['gene_id'].str.split('.').str[0]
        # Filter for genes present in filtered_gene_biotype
        tissue_egenes1 = tissue_egenes1[tissue_egenes1['gene_id'].isin(filtered_gene_biotype_pcg.index)]
        # Add 'Transcript length' to tissue_egenes
        tissue_egenes1['Transcript length'] = tissue_egenes1['gene_id'].map( filtered_gene_biotype_pcg['Transcript length (including UTRs and CDS)'])

        file2 = eqtl_files[ j ]
        tissue_name2 = file2.split("/")[-1].split(".")[0]
        tissue_egenes2 = pd.read_csv(file2, sep='\t')
        tissue_egenes2['gene_id'] = tissue_egenes2['gene_id'].str.split('.').str[0]
        # Filter for genes present in filtered_gene_biotype
        tissue_egenes2 = tissue_egenes2[tissue_egenes2['gene_id'].isin(filtered_gene_biotype_pcg.index)]
        # Add 'Transcript length' to tissue_egenes
        tissue_egenes2['Transcript length'] = tissue_egenes2['gene_id'].map( filtered_gene_biotype_pcg['Transcript length (including UTRs and CDS)'])

        _, p_value = mannwhitneyu(tissue_egenes1["Transcript length"], tissue_egenes2["Transcript length"])
        p_values[(tissue_name1, tissue_name2)] = p_value



# Correct for multiple testing
_, corrected_p_values, _, _ = multipletests(list(p_values.values()), method='fdr_bh')

# Map corrected p-values back to tissue pairs
corrected_p_values = {pair: p for pair, p in zip(p_values.keys(), corrected_p_values)}

# Now corrected_p_values contains the corrected p-values for each pair of tissues
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Assuming `corrected_p_values` is a dictionary where keys are a pair of tissues
# and values are the corresponding corrected p-values

# Generate a DataFrame from corrected_p_values
tissues = []
for i in range(len(eqtl_files)):
    file1 = eqtl_files[ i ]
    tissue_name1 = file1.split("/")[-1].split(".")[0]
    tissues.append( tissue_name1 )

p_value_df = pd.DataFrame(index=tissues, columns=tissues)

for pair, p_value in corrected_p_values.items():
    tissue1, tissue2 = pair
    p_value_df.loc[tissue1, tissue2] = p_value
    p_value_df.loc[tissue2, tissue1] = p_value

# Convert to float type
p_value_df = p_value_df.astype(float)

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(p_value_df, dtype=bool))

# Plot the heatmap
plt.figure(figsize=(20, 20))
sns.heatmap(-log10( p_value_df), mask=mask, cmap='plasma', square=True, annot=True, fmt=".2g")

plt.title('Pairwise p-values of transcript length comparisons between tissues')
plt.xlabel('Tissue')
plt.ylabel('Tissue')
plt.show()


# test difference btw ["Muscle_Skeletal", "Testis"] in lncRNA
#

#





# tissues specified
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gzip

# Path to the directory
directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"

# Load gene info data
gene_info_path = "/mnt/d/sandbox/0.gxe/data/mart_human.txt.gz"
gene_info = pd.read_csv(gene_info_path, sep="\t", usecols=[0, -1], names=["gene_id", "transcript_length"], compression='gzip')
gene_info = pd.read_csv( "/mnt/d/sandbox/0.gxe/gtex/average_transcript_lengths_GTEx_shared.csv" )


# Select transcript length data for Testis and Muscle_Skeletal
t1_lengths = []
t2_lengths = []
t3_lengths = []

genestested = []
# Iterate over the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".v8.egenes.txt.gz"):
        tissue_name = filename.replace(".v8.egenes.txt.gz", "")
        if tissue_name in ['Cells_Cultured_fibroblasts', 'Uterus', 'Kidney_Cortex' ]:
            # Load gene expression data
            gene_expression_path = os.path.join(directory, filename)
            gene_expression = pd.read_csv(gene_expression_path, sep="\t", usecols=[0], names=["gene_id"], compression='gzip')
            gene_expression['gene_id'] = gene_expression['gene_id'].str.split('.').str[0]

            # Merge the dataframes
            merged = pd.merge(average_transcript_lengths, gene_expression, left_on='Gene stable ID', right_on='gene_id')

            if tissue_name == 'Cells_Cultured_fibroblasts':
                t1_lengths.extend(merged['Average transcript length'].tolist())
            elif tissue_name == 'Uterus':
                t2_lengths.extend(merged['Average transcript length'].tolist())
            elif tissue_name == 'Kidney_Cortex':
                t3_lengths.extend(merged['Average transcript length'].tolist())

# Convert lists to pandas Series
t1 = pd.Series(t1_lengths)
t2 = pd.Series(t2_lengths)
t3 = pd.Series(t3_lengths)

# Plot histograms
plt.figure(figsize=(12, 6))
sns.histplot(t1, color='blue', kde=True, label='Cells_Cultured_fibroblasts')
sns.histplot(t2, color='red', kde=True, label='Uterus')
sns.histplot(t3, color='green', kde=True, label='Kidney_Cortex')

plt.title('Distribution of Transcript Lengths')
plt.xlabel('Transcript Length')
plt.ylabel('Frequency')
plt.legend()
plt.show()



# compare distribution for skeletal muscle and testis
import matplotlib.pyplot as plt
import seaborn as sns

# Extract transcript lengths for Testis and Muscle_Skeletal
testis_lengths = gene_expression.loc[gene_expression['Tissue'] == 'Testis', 'Transcript Length']
muscle_skeletal_lengths = gene_expression.loc[gene_expression['Tissue'] == 'Muscle_Skeletal', 'Transcript Length']

# Plot histograms
plt.figure(figsize=(12, 6))
sns.histplot(testis_lengths, color='blue', kde=True, label='Testis')
sns.histplot(muscle_skeletal_lengths, color='red', kde=True, label='Muscle_Skeletal')

plt.title('Distribution of Transcript Lengths')
plt.xlabel('Transcript Length')
plt.ylabel('Frequency')
plt.legend()
plt.show()


#
# save table of genes tested in GTEx
#
import pandas as pd
import os
import glob

# Set the directory
directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"

# Load the gene info data
gene_info = pd.read_csv('/mnt/d/sandbox/0.gxe/data/mart_human.txt.gz', sep='\t')


# Initialize an empty list to store all gene IDs
all_gene_ids = []
genestested = []
# Iterate over the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".v8.egenes.txt.gz"):
        tissue_name = filename.replace(".v8.egenes.txt.gz", "")
        if 1: # tissue_name in ['Testis', 'Muscle_Skeletal', "Spleen"]:
            # Load gene expression data
            gene_expression_path = os.path.join(directory, filename)
            gene_expression = pd.read_csv(gene_expression_path, sep="\t", usecols=[0], names=["gene_id"], compression='gzip')
            gene_expression['gene_id'] = gene_expression['gene_id'].str.split('.').str[0]
            merged = pd.merge(average_transcript_lengths, gene_expression, left_on='Gene stable ID', right_on='gene_id')

            # Add the gene IDs to the list
            all_gene_ids.extend(merged['gene_id'].unique().tolist())

# Remove duplicates from the list
all_gene_ids = list(set(all_gene_ids))

# Select the rows with gene IDs in all_gene_ids
all_genes_average_lengths = average_transcript_lengths[average_transcript_lengths['Gene stable ID'].isin(all_gene_ids)]

# Save average_transcript_lengths to a CSV file
all_genes_average_lengths.to_csv("../average_transcript_lengths_GTEx_all.csv", index=False)
all_gene_ids = pd.read_csv( "../average_transcript_lengths_GTEx_all.csv" )



# shared by all tissues
import pandas as pd
import os
import matplotlib.pyplot as plt

directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"
gene_counts = {}

# Iterate over the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".v8.egenes.txt.gz"):
        df = pd.read_csv(directory + filename, sep='\t', usecols=['gene_id'])
        for gene in df['gene_id']:
            gene = gene.str.split('.').str[0]
            if gene in gene_counts:
                gene_counts[gene] += 1
            else:
                gene_counts[gene] = 1

# Convert the dictionary to a DataFrame
df_counts = pd.DataFrame.from_dict(gene_counts, orient='index', columns=['count'])

# Count the number of genes in each number of tissues
tissue_counts = df_counts['count'].value_counts().sort_index()

# Plot the results
tissue_counts.plot(kind='bar')
plt.xlabel('Number of Tissues')
plt.ylabel('Number of Genes')
plt.title('Number of Genes in Each Number of Tissues')
plt.show()

#
genestested_0 = set( all_gene_ids )

# Iterate over the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".v8.egenes.txt.gz"):
        tissue_name = filename.replace(".v8.egenes.txt.gz", "")
        if 1: # tissue_name in ['Testis', 'Muscle_Skeletal', "Spleen"]:
            # Load gene expression data
            gene_expression_path = os.path.join(directory, filename)
            gene_expression = pd.read_csv(gene_expression_path, sep="\t", usecols=[0], names=["gene_id"], compression='gzip')
            gene_expression['gene_id'] = gene_expression['gene_id'].str.split('.').str[0]
            merged = pd.merge(average_transcript_lengths, gene_expression, left_on='Gene stable ID', right_on='gene_id')

            # Add the gene IDs to the list
            genestested_0 = genestested_0.intersection(  set( merged['gene_id'].unique().tolist() ) )
            #all_gene_ids.extend(merged['gene_id'].unique().tolist())


# Filter the DataFrame
filtered_df = df_counts[df_counts['count'] >= 46]

# Get the row names (index)
shared_gene_ids = filtered_df.index.tolist()


# Select the rows with gene IDs in all_gene_ids
shared_genes_average_lengths = average_transcript_lengths[average_transcript_lengths['Gene stable ID'].isin(shared_gene_ids)]

shared_genes_average_lengths= shared_genes_average_lengths[ shared_genes_average_lengths['Average transcript length'] <= 10000]
# Save average_transcript_lengths to a CSV file
shared_genes_average_lengths.to_csv("../average_transcript_lengths_GTEx_common.csv", index=False)



# Extract transcript lengths for Testis and Muscle_Skeletal

# Plot histograms
plt.figure(figsize=(12, 6))
sns.histplot(shared_genes_average_lengths["Average transcript length"], color='red', kde=True, label='Common')

plt.title('Distribution of Transcript Lengths')
plt.xlabel('Transcript Length')
plt.ylabel('Frequency')
plt.legend()
plt.show()



# other than shared genes
shared_genes_average_lengths = pd.read_csv( "../average_transcript_lengths_GTEx_common.csv" )

average_lengths_per_tissue = []
for file in os.listdir(gtex_directory):
    if file.endswith(".v8.egenes.txt.gz"):
        tissue_name = file.split(".")[0]  # get the tissue name from the filename

        gene_expression = pd.read_csv(file, sep="\t", compression='gzip')
        gene_expression['gene_id'] = gene_expression['gene_id'].str.split('.').str[0]  # remove the version number

        #merged = pd.merge(average_transcript_lengths, gene_expression, left_on='Gene stable ID', right_on='gene_id')
        gene_expression_specific = gene_expression[~gene_expression['gene_id'].isin(gene_info['Gene stable ID'])]
        merged = pd.merge( gene_expression_specific, all_gene_ids,  left_on='gene_id', right_on= 'Gene stable ID')

        mean_length = merged["Average transcript length"].mean()

        # calculate the confidence interval
        ci_low, ci_high = stats.t.interval(0.95, len(merged["Average transcript length"])-1,
                                            loc=mean_length,
                                            scale=stats.sem(merged["Average transcript length"]))

        average_lengths_per_tissue.append((tissue_name, mean_length, ci_low, ci_high))


# Convert the list to a DataFrame
average_lengths_per_tissue = pd.DataFrame(average_lengths_per_tissue, columns=["Tissue", "Average Transcript Length", "CI Low", "CI High"])

# Convert the results to a dataframe
#average_lengths_per_tissue = pd.DataFrame(results)

# Sort the DataFrame by 'Average Transcript Length' in descending order
average_lengths_per_tissue_sorted = average_lengths_per_tissue.sort_values(by='Average Transcript Length', ascending=False)

# Print the first few rows of the data
print(average_lengths_per_tissue_sorted.head())


# Print the results
for tissue, avg_length in average_transcript_lengths.items():
    print(f'Average transcript length for {tissue}: {avg_length}')


# overall transcrip legnths across 49 tissues
import seaborn as sns
import matplotlib.pyplot as plt

# Set the figure size
plt.figure(figsize=(15, 6))

# Create the barplot with error bars
plot = sns.barplot(x='Tissue', y='Average Transcript Length', data=average_lengths_per_tissue_sorted,
                   yerr=average_lengths_per_tissue_sorted['CI High'] - average_lengths_per_tissue_sorted['CI Low'])

# Rotate the x-axis labels for better readability
plot.set_xticklabels(plot.get_xticklabels(), rotation=90)

# Display the plot
plt.show()



# genes associed with significant tissue specific eQTLs

slope , pval_beta,
qval: 5th from end



# Iterate over the files in the directory
import os
import glob
import pandas as pd
import scipy.stats as stats

# Set the directory path
gtex_directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"

# Get a list of all files in the directory
files = glob.glob(os.path.join(directory_path, "*.v8.egenes.txt.gz"))

# Create an empty list to store the results
results = []


# Iterate over the files
average_transcript_lengths = pd.read_csv( "../average_transcript_lengths_GTEx_common.csv" )

average_lengths_per_tissue = []
for file in os.listdir(gtex_directory):
    if file.endswith(".v8.egenes.txt.gz"):
        tissue_name = file.split(".")[0]  # get the tissue name from the filename

        gene_expression = pd.read_csv(file, sep="\t", compression='gzip')
        gene_expression['gene_id'] = gene_expression['gene_id'].str.split('.').str[0]  # remove the version number
        gene_expression = gene_expression[ gene_expression[ "pval_nominal"] <= gene_expression[ "pval_nominal_threshold"]]

        #merged = pd.merge(average_transcript_lengths, gene_expression, left_on='Gene stable ID', right_on='gene_id')
        merged = pd.merge(average_transcript_lengths , gene_expression, left_on='Gene stable ID', right_on='gene_id')

        mean_length = merged["Average transcript length"].mean()

        # calculate the confidence interval
        ci_low, ci_high = stats.t.interval(0.95, len(merged["Average transcript length"])-1,
                                            loc=mean_length,
                                            scale=stats.sem(merged["Average transcript length"]))

        average_lengths_per_tissue.append((tissue_name, mean_length, ci_low, ci_high))


# Convert the list to a DataFrame
average_lengths_per_tissue = pd.DataFrame(average_lengths_per_tissue, columns=["Tissue", "Average Transcript Length", "CI Low", "CI High"])

# Convert the results to a dataframe
#average_lengths_per_tissue = pd.DataFrame(results)

# Print the first few rows of the data
print(average_lengths_per_tissue.head())


# Print the results
for tissue, avg_length in average_transcript_lengths.items():
    print(f'Average transcript length for {tissue}: {avg_length}')

average_lengths_per_tissue_sorted = average_lengths_per_tissue.sort_values(by='Average Transcript Length', ascending=False)

# overall transcrip legnths across 49 tissues
import seaborn as sns
import matplotlib.pyplot as plt

# Set the figure size
plt.figure(figsize=(15, 6))

# Create the barplot with error bars
plot = sns.barplot(x='Tissue', y='Average Transcript Length', data=average_lengths_per_tissue_sorted,
                   yerr=average_lengths_per_tissue_sorted['CI High'] - average_lengths_per_tissue_sorted['CI Low'])

# Rotate the x-axis labels for better readability
plot.set_xticklabels(plot.get_xticklabels(), rotation=90)

# Display the plot
plt.show()


# mouse. etiologic brain cell types for obesity
# Timshel et al. 2020
set( brain[ "Symbol"] )  & set( pnglist_neg)
{'ACHE', 'CEP19', 'HIP1R', 'LEP', 'LEPR',
'PCSK1', 'RAPGEF3', 'SIM1', 'ZBTB7B', 'ZFHX3'}

set( brain[ "Symbol"] )  & set( pnglist_neg)
 {'GIPR', 'POMC'}
