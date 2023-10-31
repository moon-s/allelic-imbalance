# DEG


deg_files = glob.glob("/mnt/d/sandbox/1.gxe/rnaseq/DEG_*_genes.tsv")

# DEG in obese
adis = ["SWA_OB_GSE110729", "WAT_NGT_GSE141432", "WAT_T2D_GSE141432",
        "SAT_MHO_GSE156906", "SAT_MUO_GSE156906", "WAT_OB_GSE162653",
        "OAT_IR_GSE165932", "EAT_OB_GSE179455", "SAT_OB_GSE205668",
        "OAT_OB_GSE55008"]
gene_dict = {} #
#gene_category = gene_biotype_df[ "gene_category"] == "PCG"
alldegsets = {"up":[], "down":[]}
for i in range( len( deg_files)) :
        file = deg_files[ i ]
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
        merged_df_sig_down = merged_df[ ( merged_df[ fold_change_col] <= -np.log2( 1.5 ) )  ]  # | ( merged_df[ fold_change_col] <= -log2( 1.5 ) )
        merged_df_sig_up = merged_df[ ( merged_df[ fold_change_col] >= np.log2( 1.5 ) )  ]  # | ( merged_df[ fold_change_col] <= -log2( 1.5 ) )
        #
        compid = adis[ i]
        gene_dict[ compid ] = { "down": merged_df_sig_down.set_index('Gene stable ID')[fold_change_col ].index.tolist(),
                            "up": merged_df_sig_up.set_index('Gene stable ID')[fold_change_col ].index.tolist() }
        alldegsets[ "down"].extend( merged_df_sig_down.set_index('Gene stable ID')[fold_change_col ].index.tolist() )
        alldegsets[ "up"].extend( merged_df_sig_up.set_index('Gene stable ID')[fold_change_col ].index.tolist() )
        # up-regulated genes
        #gene_dict["obese_up"].update( merged_df_sig_pos.set_index('Gene stable ID')[fold_change_col ].to_dict() )
        #gene_dict["obese_down"].update( merged_df_sig_neg.set_index('Gene stable ID')[fold_change_col ].to_dict() )
        #gene_dict[ sid ] =  merged_df_sig_pos.set_index('Gene stable ID')[fold_change_col ].to_dict()
        #alldegsets += merged_df_sig_neg[ "Gene stable ID"].tolist()
        print( sid, compid, len( gene_dict[ compid]["up" ] ), len( gene_dict[ compid] ["down"] ) )  # sum( merged_df_sig_neg[fold_change_col].abs() >= np.log2( 1.5) ), sum( merged_df_sig_neg[pval_col] > -np.log10(0.05)  )  )



# Get unique list of all genes
all_genes = list({gene for study in gene_dict.values() for direction in study.values() for gene in direction})
alldegsets[ "com"] = set( alldegsets["up"]) &  set( alldegsets["down"])

obese_up_gene = set( alldegsets["up"]) - set( alldegsets[ "com"] )
obese_down_gene = set( alldegsets["down"]) - set( alldegsets[ "com"] )



# for merged set of eGenes between two tissues

# for merged set of eGenes between two tissues
bt = "PCG"
g_cate = [ "PCG", "lncRNA" ]
deg_ai_in = {}
for bt in g_cate:
    deg_ai_in[ bt ] = {}
    deg_ai_all = set()
    gbio = set( gene_biotype_df[ gene_biotype_df["gene_category"] ==bt ][ "Gene stable ID"] )
    if bt == "PCG":
        for ai in PCGs:
            deg_ai_in[ bt ][ ai ] = ( obese_up_gene & set( PCGs[  ai ]),   obese_down_gene & set( PCGs[ ai  ]) )
            print( bt, ai, len( deg_ai_in[ bt ][ ai ][0] ), len( deg_ai_in[ bt ][ ai ][1] ) )
            deg_ai_all.update(  obese_up_gene & set( PCGs[ ai]) )
            deg_ai_all.update(  obese_down_gene & set( PCGs[ ai]) )

    deg_ai_in[ bt ][ "No cis-eQTL"] = (  ( obese_up_gene - deg_ai_all) & gbio ,  (obese_down_gene - deg_ai_all) & gbio  )
    ai = "No cis-eQTL"
    print( bt, ai, len( deg_ai_in[ bt ][ ai ][0] ), len( deg_ai_in[ bt ][ ai ][1] ) )
    if bt == "lncRNA":
        for ai in lncRs:
            deg_ai_in[ bt ][ ai ] = ( obese_up_gene & set( lncRs[  ai ]),   obese_down_gene & set( lncRs[ ai  ]) )
            print( bt, ai, len( deg_ai_in[ bt ][ ai ][0] ), len( deg_ai_in[ bt ][ ai ][1] ) )
            deg_ai_all.update(  obese_up_gene & set( lncRs[ ai]) )
            deg_ai_all.update(  obese_down_gene & set( lncRs[ ai]) )
    ai = "No cis-eQTL"
    deg_ai_in[ bt ][ "No cis-eQTL"] = (  ( obese_up_gene - deg_ai_all) & gbio ,  (obese_down_gene - deg_ai_all) & gbio  )
    print( bt, ai, len( deg_ai_in[ bt ][ ai ][0] ), len( deg_ai_in[ bt ][ ai ][1] ) )




deg_count = []
for bt in deg_ai_in:
    for d in deg_ai_in[ "PCG" ]:
        print( bt, d,  len( deg_ai_in[ bt][d][0]) , len( deg_ai_in[ bt][ d][1] ))
        deg_count.append( (bt, d,  len( deg_ai_in[ bt][d][0]) , len( deg_ai_in[ bt][ d][1] )) )

deg_count = pd.DataFrame( deg_count)
deg_count.columns = ( "Gene type", "AI level", "Up-regulated", "Down-regulated" )



bt = "PCG"
# Convert the dataframe from wide form to long form
deg_count_long = pd.melt(deg_count[ deg_count["Gene type"] == bt ], id_vars=["AI level"], value_vars=["Up-regulated", "Down-regulated"],
                         var_name="Regulation", value_name="Count")


# Define custom color palette
# Very low to very high, and no cis-eQTL
palette = {
    "Up-regulated": ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"],
    "Down-regulated": ['#08519c','#3182bd'  ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']
}



deg_count_long["color"] = deg_count_long.apply(lambda row: palette[row["Regulation"]][deg_count["AI level"].tolist().index(row["AI level"])], axis=1)

# Plot
plt.figure(figsize=(6, 4))
bar = sns.barplot(x="AI level", y="Count", hue="Regulation", data=deg_count_long, edgecolor="white", palette=['#fd8d3c', '#6baed6' ])#  palette=deg_count_long["color"].unique())

for i,thisbar in enumerate(bar.patches):
    # Set a different hatch for each bar
    thisbar.set_color( deg_count_long["color"][i])
    thisbar.set_edgecolor("k")
plt.legend( title= "DEG in obese ATs")
plt.xlabel( "Allelic Imbalance Level" )
plt.title("Number of DE PCGs in ATs from individuals with obesity by AI Level")
sns.despine( offset=10, trim=True)
plt.tight_layout()
plt.savefig( "Fig.3B.DEG_number_AI_PCG.pdf")
plt.show()


# DEG PPI
# closeness centrality


ppi_genes = pd.read_csv( "./result/top 10 percent hubgenes default node.csv" )
bt = "PCG"
gtypes = []
bt = "PCG"
gtypes = []
for g in ppi_genes[ "shared name"]:
    for k in deg_ai_in[ bt]:
        a = deg_ai_in[ bt][ k ]
        gn_up = gene_biotype_gtex[ gene_biotype_gtex["Gene stable ID"].isin(  a[0] )][ "Gene name"].to_list()
        gn_down = gene_biotype_gtex[ gene_biotype_gtex["Gene stable ID"].isin(  a[1] )][ "Gene name"].to_list()
        if g in gn_up :
            g_type = "Up & " + k
            gtypes.append( g_type )
        elif g in gn_down:
            g_type = "Down & " + k
            gtypes.append( g_type )

ppi_genes[ "gene_ai_level" ] = gtypes
# Split 'gene_ai_level' into two columns 'Regulation' and 'AI level'
ppi_genes[['Regulation', 'AI level']] = ppi_genes['gene_ai_level'].str.split(' & ', expand=True)


# Plot
plt.figure(figsize=(6, 4))
bar = sns.barplot( y='ClosenessCentrality', x='AI level', hue='Regulation', data=ppi_genes, order=["Very Low", 'Low', 'Moderate', 'High', 'Very High', 'No cis-eQTL'],
            palette=['#fd8d3c', '#6baed6' ])
for i,thisbar in enumerate(bar.patches):
    # Set a different hatch for each bar
    thisbar.set_color( deg_count_long["color"][i])
    thisbar.set_edgecolor("k")
plt.title("Closeness Centrality of DEGs by AI Level")
plt.legend( [],[], frameon=False)
sns.despine( offset=10, trim=True)
plt.tight_layout()
plt.show()



# hub genes from PPI
plot_data = ppi_genes[["name", "gene_ai_level", "NumberOfUndirectedEdges", "Regulation"]]

print( sum( plot_data["Regulation"] == "Up"),  sum( plot_data["Regulation"] == "Down"),
      np.quantile( plot_data[ plot_data["Regulation"] == "Up" ]["NumberOfUndirectedEdges"], 1 - 0.02882),
np.quantile( plot_data[ plot_data["Regulation"] == "Down" ]["NumberOfUndirectedEdges"], 1- 0.0681))
plot_data = ppi_genes[["name", "gene_ai_level", "NumberOfUndirectedEdges", "Regulation"]]
ndown = 62
nup = 152.1
plot_data = plot_data[ ( plot_data["NumberOfUndirectedEdges"] >= nup ) & (plot_data[ "Regulation"] == "Up" ) ]

plot_data = plot_data.sort_values("NumberOfUndirectedEdges", ascending=False )


palette = {
    "Up-regulated": ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"],
    "Down-regulated": ['#08519c', '#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']
}


# Define the color mapping
colors = plot_data['gene_ai_level'].map({
    "Up & No cis-eQTL": "#feedde",
    "Up & Very High": "#fdd0a2",
    "Up & High": "#fdae6b",
    "Up & Moderate": "#fd8d3c",
    "Up & Low": "#e6550d",
    "Up & Very Low": "#a63603",

    "Down & No cis-eQTL": "#eff3ff",
    "Down & Very High": "#c6dbef",
    "Down & High": "#9ecae1",
    "Down & Moderate": "#6baed6",
    "Down & Low": "#3182bd",
    "Down & Very Low": "#08519c",

})

# 3. Plot
plt.figure(figsize=(4, 5))
bars = sns.barplot(data=plot_data, x="NumberOfUndirectedEdges", y="name"  , palette=colors )
plt.title("Top 10 hub up-regulated PCGs")
plt.xlabel("Number of interacting genes")
plt.ylabel("Hub gene")
plt.xticks(rotation=90)  # rotate x-axis labels for better readability
plt.tight_layout()
sns.despine( offset=10, trim = True)
plt.savefig( "Fig.3G.top10_hub_up_PCGS.pdf")
