# cis-eQTL defines regulatory DHS associated with eGenes stratified by AI levels



# DHS density, cis-eQTL density, per gene
# Transcript length

results = []
for bt_inx in ["PCG", "lncRNA"]:
    grouping = [ "Very Low", "Low", "Moderate", "High", "Very High"]
    group = ['Very Low', 'Low', 'Moderate', 'High', 'Very High']
    for j in range( 5 ):
        if bt_inx == "PCG":
            subset_genes = PCGs[ group[ j ]  ]
        elif bt_inx == "lncRNA":
            subset_genes = lncRs[ group[ j ] ]
        evar = obtain_evar_geneset( evar_df, subset_genes )
        evar_nr = evar[ [ "chrom", "pos", "dhs_id", "gene_id"] ].drop_duplicates()

        num_var_nodhs = sum( evar["dhs_id"] == "N" )
        num_var_dhs = sum( evar[ "dhs_id"] != "N")
        data = evar_nr[~evar_nr['dhs_id'].str.startswith('N')]

        # Calculate average number of eQTL per gene
        num_eQTL_per_gene = data['gene_id'].value_counts()

        # Calculate average number of DHS per gene
        num_unique_dhs_per_gene = data.groupby('gene_id')['dhs_id'].nunique()

        # Match the Gene stable IDs with their transcript lengths
        gene_transcript_lengths = gene_biotype_df[gene_biotype_df['Gene stable ID'].isin( data['gene_id'] )]['Transcript length (including UTRs and CDS)']
            # Iterate over genes in this group and tissue
        for gene_id in data['gene_id'].unique():
            results.append({
                'Gene type': bt_inx,
                'Group': grouping[j ],
                'Gene ID': gene_id,
                'Num eQTL per gene': num_eQTL_per_gene[gene_id],
                'Num DHS per gene': num_unique_dhs_per_gene[gene_id],
                'Transcript Length': gene_transcript_lengths[gene_biotype_df['Gene stable ID'] == gene_id].values[0]
            })


# Convert the results into a DataFrame
results_df = pd.DataFrame(results)



# Figure 2C Number of cis-eQTL and
# Figure 2D Number of DHS

bt = "lncRNA"
# [ results_df["Gene type"] == bt ]
ax = sns.barplot( data = results_df, y = "Num eQTL per gene",  x = "Gene type", hue = "Group" , palette=[ col_set[4 - u] for u in range(5)] ,
                 hue_order=[ group[4-u] for u in range(5)]  )
sns.despine( offset=10, trim=True)
plt.ylabel( "Average number of eQTL per gene" )
#plt.legend( labels = [group[4-i] for i in range( 5)] , )


# Figure 2E transcript length

sns.barplot(  data = results_df, y= "Transcript Length", hue = "Group", x = "Gene type" , palette=col_set)
plt.legend([],[],  frameon=False)
sns.despine( offset=10, trim=True)

for i in group:
    print( i,  results_df[ (results_df["Gene type"] ==bt) & (results_df["Group"] ==i)][ "Transcript Length"].mean() )



# Figure 2F

# legnths of 5,3'UTR,CDS, transcript, genic region
pathtomart = "/mnt/d/sandbox/1.gxe/data/"

# Step 1: Read in the datasets
df_protein_coding = pd.read_csv( pathtomart + "mart_human_prot.txt", sep = "\t")
# fileter-in protein-coding transcrts
transcript_ids = df_protein_coding[ df_protein_coding["Transcript type"] == "protein_coding" ][ "Transcript stable ID"].to_list()

df_utr_raw = pd.read_csv(pathtomart + 'mart_human_utr.txt', sep='\t')
df_genic_raw = pd.read_csv( pathtomart + 'mart_human_genic.txt', sep='\t')

df_utr = df_utr_raw[ df_utr_raw[ "Transcript stable ID"].isin( transcript_ids )]
df_genic = df_genic_raw[ df_genic_raw[ "Transcript stable ID"].isin( transcript_ids )]

# Calculate 5' UTR length and 3' UTR length
df_utr['5 UTR length'] = abs( df_utr['5\' UTR end'] - df_utr['5\' UTR start']) + 1
df_utr['3 UTR length'] = abs( df_utr['3\' UTR end'] - df_utr['3\' UTR start'] ) + 1
df_genic['Gene length'] = abs( df_genic['Gene end (bp)'] - df_genic['Gene start (bp)'] ) + 1




subsets = { 'Very Low': VL, 'Low': L, 'Moderate': M,  'High': H,  'Very High': VH,}

# Step 3: Calculate average and SEM for each subset and each feature
features = ['5 UTR length', '3 UTR length', 'CDS Length']


results = {'AI': [], 'Feature': [], 'Average': [] }

for i in range( 5 ):
    subset_genes = PCGs[ group[ j ]  ]
    # Filter for genes in the subset
    subset_df_utr = df_utr[df_utr['Gene stable ID'].isin(subset_genes)]
    subset_df_genic = df_genic[df_genic['Gene stable ID'].isin(subset_genes)]

    # Deduplicate entries in subset_df_utr
    features_utr = [
        '5 UTR length',
        '3 UTR length',
        'CDS Length',
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
    }).reset_index()

    # Calculate the overall averages and SEM for the subset
    for feature in features:
        avg = grouped_subset_df_utr[feature]
        sem = grouped_subset_df_utr[feature]

        results['AI'].extend( [ group[ j] ]*len( avg) )
        results['Feature'].extend([ feature]*len(avg))
        results['Average'].extend(avg)


results_df = pd.DataFrame(results)


# Plotting

# Plotting
plt.figure(figsize=(8, 5))
g = sns.barplot(data=results_df, x='Feature', y='Average', hue='AI', errorbar="sd", capsize=.3,  palette= col_set)


g.set_xticklabels(['5 UTR', '3 UTR', 'CDS ' ])

plt.ylabel('Average length in each feature (bp) ')
plt.yscale( "log")
plt.title('Average Lengh of Feature across AI groups')
sns.despine( offset=10, trim = True)
#plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
