# DEGs's transcript length, DHS density, miRNAs


# DHS density, cis-eQTL density, per gene
# Transcript length

g_cate = [ "PCG", "lncRNA" ]
results = []
for bt_inx in range( 2 ):
    bt = g_cate[ bt_inx ]
    grouping = [ "Very Low", "Low", "Moderate", "High", "Very High"]
    for j in grouping :
        gset_up, gset_down = deg_ai_in[ bt][ j ]
        evar_up = obtain_evar_geneset( evar_df, gset_up )
        evar_down = obtain_evar_geneset( evar_df, gset_down )

        num_var_nodhs = sum( evar_up["dhs_id"] == "N" )
        num_var_dhs = sum( evar_up[ "dhs_id"] != "N")
        print( bt, j, len( gset_up), len( gset_down), len( evar_up), len( evar_down ) , num_var_nodhs, num_var_dhs)

        for r in [ "Up", "Down"] :
            data = evar_down[~evar_down['dhs_id'].str.startswith('N')]
            if r == "Up":
                data = evar_up[~evar_up['dhs_id'].str.startswith('N')]
            # Calculate average number of eQTL per gene
            num_eQTL_per_gene = data['gene_id'].value_counts()

            # Calculate average number of DHS per gene
            num_unique_dhs_per_gene = data.groupby('gene_id')['dhs_id'].nunique()

            # Match the Gene stable IDs with their transcript lengths
            gene_transcript_lengths = gene_biotype_df[gene_biotype_df['Gene stable ID'].isin(data['gene_id'])]['Transcript length (including UTRs and CDS)']
            # Iterate over genes in this group and tissue
            for gene_id in data['gene_id'].unique():
                results.append({
                    'Gene type': bt,
                    'Regulation': r ,
                    'Group': j ,
                    'Gene ID': gene_id,
                    'Num eQTL per gene': num_eQTL_per_gene[gene_id],
                    'Num DHS per gene': num_unique_dhs_per_gene[gene_id],
                    'Transcript Length': gene_transcript_lengths[gene_biotype_df['Gene stable ID'] == gene_id].values[0]
                })



# transcript length Vs number of DHS per gene
#


cond = (results_df["Gene type"] == "PCG" ) & ( results_df["Regulation"] == "Up" )
grouped_means = results_df[ cond  ].groupby('Group').mean()
grouped_sems = results_df[ cond ].groupby('Group').sem()

# Extract the values into variables
x_means = grouped_means['Transcript Length']
y_means = grouped_means['Num DHS per gene']
x_sems = grouped_sems['Transcript Length']
y_sems = grouped_sems['Num DHS per gene']

# Create the plot
fig, ax = plt.subplots( figsize = (6,5) )

# Scatter plot for mean values
label_order=['Very Low', "Low" , 'Moderate', 'High', "Very High"]


palette_up =  ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"]  # Up-regulated
palette_down = ['#08519c','#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']  # down-regulated



for inx in range( 5) : # grouped_means.index:
    group = label_order[ 4 - inx ]
    ax.errorbar(x_means.loc[group], y_means.loc[group],
                xerr=x_sems.loc[group], yerr=y_sems.loc[group],
                marker='o', linestyle='', label=group , color = palette_up[ 4- inx ] )


# Add labels, title, and legend
ax.set_xlabel('Average Transcript Length')
ax.set_ylabel('Average Num DHS per gene')
ax.set_title('Average Transcript Length vs Average Num DHS per gene by Group')
ax.legend()

# Show the plot
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig( "Fig.6A.Transcript_length_up_PCG.pdf" )
plt.show()



# UTRs' length
# legnths of 5,3'UTR,CDS, transcript, genic region
pathtomart = "/mnt/d/sandbox/1.gxe/data/"

# Step 1: Read in the datasets
df_utr_raw = pd.read_csv(pathtomart + 'mart_human_utr.txt', sep='\t')
df_genic_raw = pd.read_csv( pathtomart + 'mart_human_genic.txt', sep='\t')

df_utr = df_utr_raw[ df_utr_raw[ "Transcript stable ID"].isin( transcript_ids )]
# Calculate 5' UTR length and 3' UTR length
df_utr['5 UTR length'] = abs( df_utr['5\' UTR end'] - df_utr['5\' UTR start']) + 1
df_utr['3 UTR length'] = abs( df_utr['3\' UTR end'] - df_utr['3\' UTR start'] ) + 1
df_genic['Gene length'] = abs( df_genic['Gene end (bp)'] - df_genic['Gene start (bp)'] ) + 1


subsets = { 'Very Low': VL, 'Low': L, 'Moderate': M,  'High': H,  'Very High': VH,}
subsets = {}
for j in grouping + ["No cis-eQTL"] :
    gset_up, gset_down = deg_ai_in["PCG"][ j ]
    subsets[ j ]= gset_up

# Step 3: Calculate average and SEM for each subset and each feature
features = ['5 UTR length', '3 UTR length', 'CDS Length', 'Transcript length (including UTRs and CDS)', 'Gene length']


results = {'AI': [], 'Feature': [], 'Average': [] }

for subset_name, subset_genes in subsets.items():
    # Filter for genes in the subset
    subset_df_utr = df_utr[df_utr['Gene stable ID'].isin(subset_genes)]
    subset_df_genic = df_genic[df_genic['Gene stable ID'].isin(subset_genes)]

    # Deduplicate entries in subset_df_utr
    features_utr = [
        '5 UTR length',
        '3 UTR length',
        'CDS Length',
        'Transcript length (including UTRs and CDS)'
    ]
    for feature in features_utr:
        subset_df_utr[feature] = subset_df_utr.groupby('Transcript stable ID')[feature].transform('first')
    subset_df_utr = subset_df_utr.drop_duplicates(subset=['Transcript stable ID', 'Gene stable ID'])

    # Deduplicate entries in subset_df_genic for gene length
    subset_df_genic['Gene length'] = subset_df_genic.groupby('Transcript stable ID')['Gene length'].transform('first')
    subset_df_genic = subset_df_genic.drop_duplicates(subset=['Transcript stable ID', 'Gene stable ID'])

    # Group by Gene ID and average the features for UTR dataframe
    grouped_subset_df_utr = subset_df_utr.groupby('Gene stable ID').agg({
        '5 UTR length': 'mean',
        '3 UTR length': 'mean',
        'CDS Length': 'mean',
        'Transcript length (including UTRs and CDS)': 'mean'
    }).reset_index()

    # Group by Gene ID and average the features for genic dataframe
    grouped_subset_df_genic = subset_df_genic.groupby('Gene stable ID').agg({
        'Gene length': 'mean'
    }).reset_index()

    # Calculate the overall averages and SEM for the subset
    for feature in features:
        if feature == 'Gene length':
            avg = grouped_subset_df_genic[feature]
            sem = grouped_subset_df_genic[feature]
        else:
            avg = grouped_subset_df_utr[feature]
            sem = grouped_subset_df_utr[feature]

        results['AI'].extend( [subset_name]*len( avg) )
        results['Feature'].extend([ feature]*len(avg))
        results['Average'].extend(avg)


results_df = pd.DataFrame(results)

results_down_df = pd.DataFrame(results)


results_ud_df = pd.concat(  results_up_df, results_down_df )




plt.figure(figsize=(7, 6))


palette =  ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"] +  ['#08519c','#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']  # down-regulated



bar = sns.barplot( data = results_ud_df[ results_ud_df[ "Feature"] == "5 UTR length"], x = "AI", y = "Average" ,  hue= "Regulation")
for i,thisbar in enumerate(bar.patches):
    print( i, thisbar )
    # Set a different hatch for each bar
    thisbar.set_color( palette[i])
    thisbar.set_edgecolor("k")

plt.ylabel( "Average length (bp)" )
plt.xlabel( "Allelic Imbalance Level" )
plt.title("5'-UTR length by AI Level and Regulation")
sns.despine( offset=10, trim = True)
plt.show()



#

mirt = pd.read_csv( "/mnt/d/sandbox/1.gxe/data/hsa_MTI.csv")

gtype = "miRNA"
obese_up_gene = set( alldegsets["up"]) - set( alldegsets[ "com"] )
obese_down_gene = set( alldegsets["down"]) - set( alldegsets[ "com"] )

mirna_df =  gene_biotype_df[ (gene_biotype_df["Gene type"] == gtype) & (gene_biotype_df[ "Gene stable ID"].isin( obese_up_gene | obese_down_gene)) ]


# down miRNAs in obese
miRNAs_down_genename =  mirna_df[ mirna_df['Gene stable ID'].isin(obese_down_gene) ]['Gene name'].tolist()
miRNAs_down = ['hsa-' + x.lower().replace('mir', 'miR-') for x in miRNAs_down_genename]  # 16
miRNAs_b35_down =  [s + "-3p" for s in set( miRNAs_down) ] + [s + "-5p" for s in set( miRNAs_down) ]
mirna_dn35 = set( list( miRNAs_down) + miRNAs_b35_down )

#up miRNAs in obese
miRNAs_up_genename = mirna_df[ mirna_df['Gene stable ID'].isin( obese_up_gene) ]['Gene name'].tolist()
miRNAs_up = ['hsa-' + x.lower().replace('mir', 'miR-') for x in miRNAs_up_genename]  #23
miRNAs_b35_up =  [s + "-3p" for s in set( miRNAs_up) ] + [s + "-5p" for s in set( miRNAs_up) ]
mirna_up35 = set ( list( miRNAs_up) + miRNAs_b35_up  )

all_target_up_genes = set()
all_target_down_genes = set()
mirt_up_mirs_targetgenes = {}
for m in mirna_up35:
    tmp = mirt[ mirt[ "miRNA"].isin( [m] ) ][ "Target Gene"]
    if len( set( tmp ) ) > 0 :
        mirt_up_mirs_targetgenes[ m ] = set( tmp )
        all_target_up_genes.update( set( tmp ))

mirt_down_mirs_targetgenes = {}
for m in mirna_dn35:
    tmp = mirt[ mirt[ "miRNA"].isin( [m] ) ][ "Target Gene"]
    if len( set( tmp ) ) > 0 :
        mirt_down_mirs_targetgenes[ m ] = set( tmp )
        all_target_down_genes.update( set( tmp ))


for k in deg_ai_in[ bt ]:
    up, down = deg_ai_in[ bt ][ k ]
    up_gid = set( gene_biotype_df[ (  gene_biotype_df["Gene stable ID"].isin( up ) ) ][ "Gene name"].to_list() )
    down_gid = set( gene_biotype_df[ ( gene_biotype_df["Gene stable ID"].isin( down ) ) ][ "Gene name"].to_list() )
    print( k, len( up_gid & all_target_up_genes ), len( down_gid & all_target_up_genes ),
            len( up_gid & (all_target_down_genes - all_target_up_genes )), len( down_gid & ( all_target_down_genes - all_target_up_genes ) )  )



# make heatmap


# Initialize an empty DataFrame to store the results

# Example dicts
# mirna_up35 = {'mirna_1': ['gene_1', 'gene_2'], 'mirna_2': ['gene_2', 'gene_3']}
# deg_ai_in = {'Low': ['gene_1', 'gene_3'], 'High': ['gene_2', 'gene_4']}


# for network
term = { 'Very Low' : 'VL',
'Low' : 'L',
'Moderate' : "M",
'High' : 'H',
'Very High' : 'VH',
'No cis-eQTL' : 'N' }

mirna2gene = []
geneanno = []

# up-regulated miRNA dnd down-regulated PCG
result_df = pd.DataFrame()

for ai_level, ai_genes in deg_ai_in[ "PCG"].items():
    for mirna_id, mirna_genes in mirt_up_mirs_targetgenes.items():
        down_gid = set( gene_biotype_df[  gene_biotype_df["Gene stable ID"].isin( ai_genes[1] ) ][ "Gene name"].to_list() )
        shared_genes = set( down_gid ).intersection(set(mirna_genes))
        result_df = result_df.append({
            'AI_Level': ai_level,
            'mirna_id': mirna_id,
            'Shared_Count': len(shared_genes)
        }, ignore_index=True)
        if ai_level in [ "Very Low", "Low", "Moderate"] :#  mirna_id in [ "hsa-miR-3612", "hsa-miR-30a-5p", 'hsa-miR-6778-3p', 'hsa-miR-4649-3p', 'hsa-miR-22-3p', 'hsa-miR-650']:
            for g in shared_genes:
                mirna2gene.append( (mirna_id, g ))
                geneanno.append( ( term[ ai_level] +"_Down" ))


# Pivot the DataFrame to make it suitable for a heatmap
heatmap_data = result_df.pivot("mirna_id", "AI_Level", "Shared_Count").fillna(0)

# Perform clustering and create the heatmap
sns.clustermap(heatmap_data,  cmap="mako", figsize=(4,10))
plt.savefig( "Fig.6G.heatmap_up_miRNA_targets.pdf" )
plt.show()




# down-regulated miRNA and up-target PCG
result_df = pd.DataFrame()

for ai_level, ai_genes in deg_ai_in[ "PCG"].items():
    for mirna_id, mirna_genes in mirt_down_mirs_targetgenes.items():
        up_gid = set( gene_biotype_df[  gene_biotype_df["Gene stable ID"].isin( ai_genes[0] ) ][ "Gene name"].to_list() )
        shared_genes = set( up_gid ).intersection(set(mirna_genes))
        result_df = result_df.append({
            'AI_Level': ai_level,
            'mirna_id': mirna_id,
            'Shared_Count': len(shared_genes)
        }, ignore_index=True)


# Pivot the DataFrame to make it suitable for a heatmap
heatmap_data = result_df.pivot("mirna_id", "AI_Level", "Shared_Count").fillna(0)

# Perform clustering and create the heatmap
sns.clustermap(heatmap_data, figsize=(4,9))
plt.savefig( "Fig.6G.heatmap_down_miRNA_targets.pdf" )

plt.show()




#



# Variables
# Variables
directory = "/mnt/d/sandbox/1.gxe/result/"
regulations = ["up", "down"]

# Helper function to read and process data

def read_and_process(filename):
    df = pd.read_csv(filename)
    regulation = filename.split("_")[ 2].split("-")[0]
    print( filename, regulation )
    df["Regulation"] = regulation
    return df[[ "Regulation", "term_name", "negative_log10_of_adjusted_p_value", "intersection_size", "intersections"]]

# Read in all CSVs
all_data = []
for file in os.listdir(directory):
    if file.endswith("mirna.csv"):
        all_data.append(read_and_process(directory + file))

combined_df = pd.concat(all_data)


combined_df["Regulation"] = pd.Categorical(combined_df["Regulation"], categories=regulations, ordered=True)
combined_df = combined_df.sort_values(by=["Regulation", "negative_log10_of_adjusted_p_value"], ascending=[False, False])
combined_df.columns = ['Regulation', 'term_name', '-log10(adjusted P-value)',
       'Number of genes', 'intersections']
# Pl

# Plot
combined_df = combined_df.sort_values(by=["Regulation", "-log10(adjusted P-value)"], ascending=[True, False])

# Plot
if 1:
    subset = combined_df[combined_df["Regulation"] != 1]
    plt.figure(figsize=(2, 4))
    sns.scatterplot(data=subset, x="Regulation", y="term_name", size="Number of genes",
                    hue="-log10(adjusted P-value)", sizes=(10, 300), palette="viridis")
    plt.title(f"GO Enrichment for targets of up- & {reg}-regulated miRNAs")
    plt.xlim( -0.5, 1.5)
    plt.xlabel("Regulation")
    plt.ylabel("GO Term")
    plt.legend(title="Significance", bbox_to_anchor=(1, 1))
    #plt.tight_layout()
    plt.show()




# annotation to miRNAs' targets

mirtar_df = pd.read_csv( "mirna2gene_df.csv default node.csv" )

tmp =[]
for ai_level, ai_genes in deg_ai_in[ "PCG"].items():
    for mirna_id, mirna_genes in mirt_up_mirs_targetgenes.items():
        down_gid = set( gene_biotype_df[  gene_biotype_df["Gene stable ID"].isin( ai_genes[1] ) ][ "Gene name"].to_list() )
        up_gid = set( gene_biotype_df[  gene_biotype_df["Gene stable ID"].isin( ai_genes[0] ) ][ "Gene name"].to_list() )
        for inx, row in mirtar_df.iterrows():
            if row[ "name"] in down_gid:
                tmp.append( (row["name"] , "down_" + term[ ai_level] ) )
            if row[ "name"] in up_gid:
                tmp.append( (row["name"] , "up_" + term[ ai_level] ) )

tmp_df = pd.DataFrame( tmp )
tmp_df.to_csv( "miR_target_gene_annot.csv")
