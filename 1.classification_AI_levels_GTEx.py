#
egene_path = "/mnt/d/sandbox/1.gxe/gtex/GTEx_Analysis_v8_eQTL/"

label_order=['Very Low', "Low" , 'Moderate', 'High', "Very High"]

palette_up =  ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"]  # Up-regulated
palette_down = ['#08519c','#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']  # down-regulated

col_set = [ sns.color_palette( "mako" , n_colors= 12 )[x*2+3] for x in range( 5) ]


# by AI group
gene_biotype_gtex = pd.read_csv('/mnt/d/sandbox/1.gxe/data/gene_biotype_nonredun.csv' , low_memory = False )
gene_biotype_df = gene_biotype_gtex[gene_biotype_gtex["Length type"] == "average"]




# counting tissues per gene
# adipose S.V. aFC averag
def obtain_qval_adipo_df( qval ):
    qval = 0.05 # default
    allgenes = {}
    adipo_afc = {}
    for file in  glob.glob(egene_path + '*egenes.txt.gz') :
        tissue = file.split('/')[-1].split(".")[0]  # Get the tissue name from the filename
        egene_df = pd.read_csv( file , sep = "\t") # , usecols = [0, ] )
        egene_df_sig = egene_df[  ( egene_df["qval"] <= qval )].copy() #   &  ( abs( egene_df[ "log2_aFC"] ) > 2.5  ) & ( abs( egene_df[ "log2_aFC"] ) < 4.5 )  ]
        egene_df_sig[ "Gene stable ID"] = egene_df_sig['gene_id'].str.split('.').str[0]
        #if "Adipose" not in tissue:
        for inx, row in egene_df_sig.iterrows():
            g = row[ "Gene stable ID"]
            a = row[ "log2_aFC"]
            if g not in allgenes:
                allgenes[ g] =  []
            allgenes[ g].append( tissue )
        if "Adipose" in tissue:
            print( tissue, len( egene_df_sig) )
            for inx, row in egene_df_sig.iterrows():
                g, a = row [ [ "Gene stable ID", "log2_aFC"] ]
                if g not in adipo_afc:
                    adipo_afc[ g ] = []
                adipo_afc[ g ].append( abs( a)  ) # abs value of aFC

    # averaging aFC if gene shared by two adipose
    adipo_afc_df = {}
    for g in adipo_afc:
        adipo_afc_df[ g ] = ( np.mean( adipo_afc[ g] ), len( allgenes[ g ] ) )

    adipo_afc_df = pd.DataFrame.from_dict( adipo_afc_df , orient="index")
    adipo_afc_df.columns = [ "log2_aFC", "nTissue" ]
    return adipo_afc_df


#
# group by aFC quantiles (20% interval)
def quantile_group( adipo_afc_df, bt ):
    gtex_gene_type = set( gene_biotype_df[ gene_biotype_df["gene_category"] == bt ][ "Gene stable ID" ].to_list() )
    adipo_afc_df_filter =  adipo_afc_df[ adipo_afc_df.index.isin( gtex_gene_type ) ]

    aFC_values = adipo_afc_df_filter[ "log2_aFC" ]
    quantiles = np.quantile( aFC_values, [0.2, 0.4, 0.6, 0.8])

    # Assign genes to groups based on quantiles
    grouping = [ "Very Low", "Low", "Moderate", "High", "Very High"]

    gene_groups = {'Very Low': [], 'Low': [], 'Moderate': [], 'High': [], 'Very High': []}
    aFC_groups = {'Very Low': [], 'Low': [], 'Moderate': [], 'High': [], 'Very High': []}

    for gene, row in adipo_afc_df_filter.iterrows():
        aFC = row[ "log2_aFC"]
        if aFC <= quantiles[0]:
            gene_groups['Very Low'].append(gene)
            aFC_groups['Very Low'].append(aFC)
        elif aFC <= quantiles[1]:
            gene_groups['Low'].append(gene)
            aFC_groups['Low'].append(aFC)
        elif aFC <= quantiles[2]:
            gene_groups['Moderate'].append(gene)
            aFC_groups['Moderate'].append(aFC )
        elif aFC <= quantiles[3]:
            gene_groups['High'].append(gene)
            aFC_groups['High'].append(aFC )
        else:
            gene_groups['Very High'].append(gene)
            aFC_groups['Very High'].append(aFC)
    return gene_groups, aFC_groups

adipo_df = obtain_qval_adipo_df(0.05)

PCGs, pcg_afcs = quantile_group(adipo_df, "PCG")
lncRs, lnc_afcs = quantile_group(adipo_df, "lncRNA")



# total number of non-redundant cis-eQTL within DHS
#

# Figure 1B
data_evariant = []
tissue_list = []
for file in glob.glob('/mnt/d/sandbox/1.gxe/dhs/gtex_evariant_id/*.csv.gz')[0:2]:
#for file in glob.glob('/mnt/d/sandbox/0.gxe/dhs//mnt/d/sandbox/0.gxe/dhs/gtex_slope_all/*.csv.gz'):
    df = pd.read_csv(file, compression='gzip')
    tissue = os.path.basename(file).split('_data')[0]  # infer tissue from filename
    df['tissue'] = tissue
    tissue_list.append( tissue )
    data.append(df)

evar_df = pd.concat(data_evarint, ignore_index=True)

# obtain evariants and dhs associaed genes for a given tissue
def obtain_evar_geneset( evar_df, geneset ):
    # tissue = "Adipose_Subcutaneous"
    evariants_df_tissue  = evar_df[( evar_df['gene_id'].isin( geneset ) ) ][ [ "chrom","pos", "dhs_id", "slope", "gene_id"] ]
    return evariants_df_tissue.drop_duplicates()


#
grouping = [ "Very Low", "Low", "Moderate", "High", "Very High"]

#
# Distribtuion of allelic Fold Change
# x mean aFC, y mean number of cis-QTL per gene
num_evar_lists = []
evar_lists = []
mean_sem = []
for bt_inx in ["PCG", "lncRNA" ]:

    genesets, afcsets = quantile_group(adipo_df, bt_inx )
    for j in range( 5):
        group = grouping[ j ]
        gset = genesets[ group ]
        evar = obtain_evar_geneset( evar_df, gset )
        num_var_per_gene = evar["gene_id"].value_counts().mean(), evar["gene_id"].value_counts().sem()
        mean_afc = np.mean( pcg_afcs[group ] ), scipy.stats.sem( pcg_afcs[ group ] )
        print( bt_inx, group, num_var_per_gene, mean_afc )
        mean_sem.append( bt_inx, group, num_var_per_gene[0], num_var_per_gene[1], mean_afc[0], mean_afc[1] )

mean_df = pd.DataFrame( mean_sem )
mean_df.columns = [ "Gene type", "AI level", "Average number of cis-eQTL per gene", "cis-eQTL sem", "Average aFC", "aFC sem" ]


# Extract the values into variables
genetype = "PCG"
genetype = "lncRNA"
grouped_means = mean_df[ mean_df["Gene type"] == genetype ].groupby('AI level').mean()
#grouped_sems = mean_df[  mean_df["Gene type"] == "PCG" ].groupby('AI level').sem()

x_means = grouped_means['Average aFC']
y_means = grouped_means['Average number of cis-eQTL per gene']
x_sems = grouped_means['aFC sem']
y_sems = grouped_means['cis-eQTL sem']

# Create the plot
fig, ax = plt.subplots( figsize = (5,4) )

# Scatter plot for mean values
label_order=['Very Low', "Low" , 'Moderate', 'High', "Very High"]

palette_up =  ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"]  # Up-regulated
palette_down = ['#08519c','#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']  # down-regulated

col_set = [ sns.color_palette( "mako" , n_colors= 12 )[x*2+3] for x in range( 5) ]

for inx in range( 5) : # grouped_means.index:
    group = label_order[ 4 - inx ]
    ax.errorbar(x_means.loc[group], y_means.loc[group],
                xerr=x_sems.loc[group], yerr=y_sems.loc[group],
                marker='o', linestyle='', label=group , color = col_set[ 4- inx ] )

# Add labels, title, and legend
ax.set_xlabel('Average allelic fold change (aFC)')
ax.set_ylabel('Average number of cis-eQTL per gene')
ax.set_title('Protein-coding genes (PCG)' ) # Average allic vs Average Num DHS per gene by Group')
ax.set_title('Long non-coding RNA (lncRNA)' ) # Average allic vs Average Num DHS per gene by Group')

ax.legend( title = "AI level")
plt.savefig( "Fig.1C.average_num_cis-eQTL_%s.pdf" %( genetype ) )
# Show the plot
sns.despine(offset=10, trim=False);

plt.show()


# Figure 1F
# Total number of cis-eQTL within DHS
#

num_evar_lists = []
evar_lists = []
for bt_inx in ["PCG", "lncRNA" ]:
    genesets, afcsets = quantile_group(adipo_df, bt_inx )
    for j in range( 5):
        group = grouping[ j ]
        gset = genesets[ group ]
        evar = obtain_evar_geneset( evar_df, gset )
        num_var_nodhs = sum( evar["dhs_id"] == "N" )
        num_var_dhs = sum( evar[ "dhs_id"] != "N")
        num_dhs = len( set( evar[ evar[ "dhs_id"] != "N"]["dhs_id"] ) )
        num_evar_lists.append( (bt_inx, group, num_var_nodhs, num_var_dhs, num_dhs, len( gset ) ) )

        evar[ "Effect size"] = group
        evar[ "Gene type"] = bt_inx
        evar_lists.append( evar.copy() )

combined = pd.concat( evar_lists )
num_evar_df = pd.DataFrame( num_evar_lists )
num_evar_df.columns = ["Gene type", "AI level", "Number of cis-eQTL outside of DHS", "Number of cis-eQTL within DHS", "Number of DHS", "Number of gene"]
combined[ "slope"] = combined[ "slope"].abs()

# Draw a violin plot on the current axis
axs = axs.flatten()
fig, axs = plt.subplots(1, 2, figsize=(9, 4))


col_set = [ sns.color_palette( "mako" , n_colors= 12 )[x*2+3] for x in range( 5) ]
#col_set = [ cols[ i] for i in range( 5)  ]

fig, ax = plt.subplots( figsize = ( 5, 4))
sns.violinplot(x="Effect size", y='slope',  data=combined[ combined['Gene type'] == "PCG"] , palette=col_set ) #, ax=axs[0], palette= col_set)
plt.xlabel( "Magnitude of Allelic Imbalance")
plt.ylabel( "Effect size of cis-eQTLs")
plt.title( "Protein-coding gene (PCG)")
sns.despine( offset=-10, trim=True )
fig.savefig( "fig.1D.effect_size.PCG.pdf")


fig, ax = plt.subplots( figsize = ( 5, 4))
sns.violinplot(x="Effect size", y='slope',  data=combined[ combined['Gene type'] == "lncRNA"] , palette=col_set ) #, ax=axs[0], palette= col_set)
plt.xlabel( "Magnitude of Allelic Imbalance")
plt.ylabel( "Effect size of cis-eQTLs")
plt.title( "Long non-coding RNA (lncRNA)")
sns.despine( offset=-10, trim=True )
fig.savefig( "fig.1E.effect_size.lncRNA.pdf")

# difference of absolute value of slopes between very low and low groups.
#
filt1 = (combined['Gene type'] == "PCG") & (combined['Gene type'] == "Very Low")
filt2 = (combined['Gene type'] == "PCG") & (combined['Gene type'] == "Low")
x1 = combined[ filt1]["slope"]
x2 = combined[ filt2]["slope"]
scipy.stats.manhut

# distribution of slopes of cis-eQTL associated with tissue-specific genes (aFC > 2 or aFC <0.5 compared to control gene set)


# pleiotropy by allelic imbalnce levels
#
# tissue breadth
# number of shared genes bewteen aidpose tisssue and other tissues
#


#

tb = "lncRNA"
tb = "PCG"
# group 1: very low, to 5: very high
ai_labels = ['Very Low', 'Low', 'Moderate', "High", "Very High" ]

t_breadth_lnc = {}
t_breadth_pcg = {}

egene_49tissue = [] # {}
for file in  glob.glob(egene_path + '*egenes.txt.gz') :
    tissue = file.split('/')[-1].split(".")[0]  # Get the tissue name from the filename
    egene_df = pd.read_csv( file , sep = "\t") # , usecols = [0, ] )
    egene_df_sig = egene_df[  ( egene_df["qval"] <= 0.05 )].copy() #   &  ( abs( egene_df[ "log2_aFC"] ) > 2.5  ) & ( abs( egene_df[ "log2_aFC"] ) < 4.5 )  ]
    egene_df_sig[ "Gene stable ID"] = egene_df_sig['gene_id'].str.split('.').str[0]
    if "Adipose" not in tissue:
        for j in range( 5 ):
            group = ['Very Low', 'Low', 'Moderate', 'High', 'Very High'][ j ]
            subset_genes = PCGs[ group ]
            if group not in t_breadth_pcg:
                t_breadth_pcg[ group ] = {}
                t_breadth_lnc[ group ] = {}
            tmp =  set( egene_df_sig[ "Gene stable ID"].to_list() ) &  set( subset_genes )
            for g in tmp:
                if g not in t_breadth_pcg[ group ]:
                    t_breadth_pcg[ group ][ g ] = 0
                t_breadth_pcg[ group ][ g ] += 1
            #egene_49tissue.append( [ "PCG", group,  adipo_afc_df[ adipo_afc_df.index.isin( tmp["Gene stable ID"] )][ "nTissue"].mean(), tmp[ "log2_aFC"].abs().mean(),  len( tmp ) ) ])  # egene_df_sig[ "log2_aFC"] # )

            subset_genes = lncRs[ group ]
            tmp =  set( egene_df_sig[ "Gene stable ID"].to_list() ) &  set( subset_genes )
            for g in tmp:
                if g not in t_breadth_lnc[ group ]:
                    t_breadth_lnc[ group ][ g ] = 0
                t_breadth_lnc[ group ][ g ] += 1





# tissue breadth
#
x: AI level
y: average number of tissues


df = []
for group in t_breadth_pcg:
    for g in t_breadth_pcg[ group ]:
        df.append( ( group, "PCG", g, t_breadth_pcg[ group][ g ] ))
    for g in t_breadth_lnc[ group ]:
        df.append( ( group, "lncRNA", g, t_breadth_lnc[ group][ g ] ))

df = pd.DataFrame( df )

df.columns = [ "AI Level", "Gene type", "Gene stable ID", "Number of Tissue"]


#df = pd.read_csv(url, index_col='ID')

# Theme


#df = pd.read_csv(url, index_col='ID')

#df = pd.read_csv(url, index_col='ID')

# Theme
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':1})
palette = col_set # [ col_set[4- i] for i in range( 5) ]  # sns.color_palette("Set2", 12)

plt.figure(figsize= (5, 4) )

g = sns.FacetGrid(df[df["Gene type" ] == "PCG"], palette=palette, row="AI Level", hue="AI Level", aspect=5, height=1.,
    row_order=  [grouping[4-i] for i in range(5)])

# map df - Kernel Density Plot of IMDB Score for each Language
g.map_dataframe(sns.kdeplot, x="Number of Tissue", fill=True, alpha=1)
g.map_dataframe(sns.kdeplot, x="Number of Tissue", color='black')

# function to draw labels
def label(x, color, label):
    ax = plt.gca() #get current axis
    ax.text(0, .2, label, color='black', fontsize=13,
            ha="left", va="center", transform=ax.transAxes)
# iterate grid to plot labels
#g.map(label, "AI Level")

# adjust subplots to create overlap
g.fig.subplots_adjust(hspace=-.4)

# remove subplot titles
g.set_titles("")
plt.xlim( 0, 48)
# remove yticks and set xlabel
g.set( xlabel="Tissue Breadth")
# remove left spine
g.despine(left=True)
# set title
#plt.suptitle('Long non-coding RNA (lncRNA)', y=0.98)
plt.suptitle('Protein-coding gene (PCG)', y=0.98)
#plt.tight_layout()
plt.savefig('Fig.1I.ridgeplot_PCG.pdf')




# Figure 2 AB
# total number of cis-eQTLs withi DHS
#
numvar = []
for bt in ["PCG", "lncRNA"]:
    for g in grouping:
        tmp = combined[ (combined["Gene type"]==bt ) & ( combined['Effect size'] == g ) ]
        n1 = sum( tmp[ "dhs_id"] != "N" )
        numvar.append( ( bt, g, n1 , n1/len( tmp ) ) )

numvar = pd.DataFrame( numvar )


# Select only the "dhs" columns
numvar.columns = [ "Gene type", "AI Level" , "n_var_dhs" , "ratio" ]

# Create a figure

# Create a figure
plt.figure(figsize= (5, 4) )

# Create a barplot
sns.barplot(data=numvar[ numvar["Gene type"] == "PCG"] , x='AI Level', y='n_var_dhs', errorbar=None, palette= col_set )

plt.title("Protein-coding genes")
plt.ylabel("Total Number of cis-eQTLs within DHS")
# plt.xticks(rotation=90)  # Rotate x-axis labels for better visibility
sns.despine( offset=10, trim = True)
plt.savefig( "Fig1.F.Total_cis-eQTLs_inDHS.pdf")



# Create a figure

# Create a figure
plt.figure(figsize= (5, 4) )

# Create a barplot
sns.barplot(data=numvar[ numvar["Gene type"] == "lncRNA"] , x='AI Level', y='n_var_dhs', errorbar=None, palette= col_set )

plt.title("Long non-coding RNAs")
plt.ylabel("Total Number of cis-eQTLs within DHS")
# plt.xticks(rotation=90)  # Rotate x-axis labels for better visibility
sns.despine( offset=10, trim = True)
plt.savefig( "Fig1.G.Total_cis-eQTLs_inDHS_lncRNA.pdf")
