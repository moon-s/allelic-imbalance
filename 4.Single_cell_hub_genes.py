
# DimPlot( human, reduction = "umap" )
library(dplyr)
library(Seurat)
library(patchwork)

human = readRDS( "human_adipocytes.rds")

human = readRDS( "human_all.rds")
mouse = readRDS( "mouse_all.rds")

human = readRDS( "human_all.rds")

ctx =  levels( mouse)

for (x in  ctx ){
    xtable = FindMarkers( mouse, ident.1 = "HFD", ident.2 = "Chow", group.by ="diet", subset.ident = x , p_val_adj = 0.05, logfc.threshold = 0.1  )
    write.csv( xtable [ xtable ["p_val_adj"] <= 0.05, ], paste( "./cell_deg_mouse0.1/",  gsub( "/", "_", x) , ".csv") )
}

write.csv( xtable, "hfd_deg_cells.csv")


human = readRDS( "human_all.rds")
ctx =  levels( human)

for (x in  ctx ){
    xtable = FindMarkers( human, ident.1 = "40-50", ident.2 = "20-30", group.by ="bmi.range", subset.ident = x , p_val_adj = 0.05, logfc.threshold = 0.1  )
    write.csv( xtable [ xtable ["p_val_adj"] <= 0.05, ], paste( "./cell_deg_human0.1/",  gsub( "/", "_", x) , ".csv") )
}


DimPlot( human, features =  c('TNF','IL6','PTPRC'),  reduction="umap", raster = FALSE, split.by = "bmi.range" )


plots <- FeaturePlot(human, features =  c('CD4','PTPRC', 'ITGAM'), order = T , raster = FALSE, split.by = "bmi.range")
plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5)))
CombinePlots(plots = plots)

# PPI
up-regulated PCGs
FeaturePlot(human, "CD4", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c(option = "inferno", alpha = 0.7)
FeaturePlot(human, "PTPRC", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c(option = "inferno", alpha = 0.7)
FeaturePlot(human, "ITGAM", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c(option = "inferno", alpha = 0.7)

down-regulated PCGs
FeaturePlot(human, "PPARG", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "APOE", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "IGF1", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "PTEN", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "ESR1", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "UBA52", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "ADIPOQ", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "CEBPB", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)

# Ligand-receptor
FeaturePlot(human, "FN1", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c(option = "inferno", alpha = 0.7)
FeaturePlot(human, "GRM7", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c(option = "inferno", alpha = 0.7)

FeaturePlot(human, "WNT5A", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)
FeaturePlot(human, "VEGFA", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)

FeaturePlot(human, "TSHR", pt.size = 0.5, split.by = "bmi.range")  & scale_color_viridis_c( alpha = 0.7)



# GO enrichment
# biological process


# Variables
directory = "/mnt/d/sandbox/1.gxe/result/GO/"
ai_levels = ["Very Low", "Low", "Moderate", "High", "Very High", "No cis-eQTL"]  # Assuming these are the AI levels
regulations = ["Up", "Down"]

# Helper function to read and process data
def read_and_process(filename):
    df = pd.read_csv(filename)
    df = df.head( 5)
    ai_level, regulation, _ = filename.split("_", 2)
    df["AI level"] = ai_level.split( "/")[-1]
    df["Regulation"] = regulation
    return df[["AI level", "Regulation", "term_name", "negative_log10_of_adjusted_p_value", "intersection_size", "intersections"]]

# Read in all CSVs
all_data = []
for file in os.listdir(directory):
    if file.endswith(".csv"):
        all_data.append(read_and_process(directory + file))

combined_df = pd.concat(all_data)

combined_df.columns = ['AI level', 'Regulation', 'term_name',
       '-log10(adj. p-value)', 'intersection_size',
       'intersections']



# Pl

# Plot
ai_order = ["VL", "L", "M", "H", "VH", "N"]

combined_df["AI level"] = pd.Categorical(combined_df["AI level"], categories=ai_order, ordered=True)
combined_df = combined_df.sort_values(by=["AI level", "-log10(adj. p-value)"], ascending=[True, False])


for reg in regulations:
    subset = combined_df[combined_df["Regulation"] == reg]
    fig, ax = plt.subplots()
    fig.set_size_inches(4, 7, forward=True)
    #plt.figure(figsize=(3, 6))

    sns.scatterplot(data=subset, x="AI level", y="term_name", size="intersection_size", hue="-log10(adj. p-value)",
                    sizes=(10, 400), palette="viridis")
    plt.title(f"{reg}-regulated genes")
    plt.xlabel("AI Level")
    plt.ylabel("GO Term")
    plt.legend(title="Significance", bbox_to_anchor=(1, 1))

    sns.despine( offset=15, trim=True)

    plt.savefig( "Fig.5AB.BP_enrichment_%s.pdf" % (reg) )

    #plt.tight_layout()

    plt.show()




# Ligand-receptor


# ligand-receptor interaction
lr_data = pd.read_csv( "./data/human_lr_pair.txt", sep="\t")


#ligand- receptor pair.

go_deg_pairs = []
go_degs_annot = []

# subset of DEGs in significant GO terms
deg_go_total = []
for k in deg_go:
    degset = deg_go[ k ]
    for inx, row in lr_data.iterrows():
        l_id = row[  "ligand_ensembl_gene_id" ]
        r_id = row[  "receptor_ensembl_gene_id" ]

        if len( set( [l_id, r_id] ) & degset )  >= 1:
                go_deg_pairs.append( ( row[ "ligand_gene_symbol"], row[ "receptor_gene_symbol"] ) )
                if l_id in degset:
                    go_degs_annot.append( ( row[ "ligand_gene_symbol"], k[0], k[1] , "L"   )  )
                    deg_go_total.add( row[ "ligand_gene_symbol"] )
                if r_id in degset:
                    go_degs_annot.append( ( row[ "receptor_gene_symbol"], k[0], k[1] , 'R'  )  )
                    deg_go_total.add( row[ "receptor_gene_symbol"] )



# all DEGs # all DEGs
term = { 'Very Low' : 'VL',
'Low' : 'L',
'Moderate' : "M",
'High' : 'H',
'Very High' : 'VH',
'No cis-eQTL' : 'N' }


all_deg_pairs = []
all_degs_annot = []
all_deg_total = set([])
tmp = []
bt = "PCG"
for k in deg_ai_in[  bt ]:
    up, down = deg_ai_in[ bt ][ k ]
    for inx, row in lr_data.iterrows():
        l_id = row[  "ligand_ensembl_gene_id" ]
        r_id = row[  "receptor_ensembl_gene_id" ]
        if len( set( [ l_id, r_id ])  & ( set( up) | set( down ) ))  >= 1:
            all_deg_pairs.append( ( row[ "ligand_gene_symbol"], row[ "receptor_gene_symbol"] ) )
            count = set([])
            if l_id in up:
                all_degs_annot.append( ( row[ "ligand_gene_symbol"], term[ k ], "Up" , 'L'))
                all_deg_total.add( row[ "ligand_gene_symbol"] )
            if r_id in up:
                all_degs_annot.append( ( row[ "receptor_gene_symbol"], term[ k ], "Up", 'R' ))
                all_deg_total.add( row[ "receptor_gene_symbol"] )
            if l_id in down:
                all_degs_annot.append( ( row[ "ligand_gene_symbol"], term[ k ], "Down" , "L"))
                all_deg_total.add( row[ "ligand_gene_symbol"] )
            if r_id in down:
                all_degs_annot.append( ( row[ "receptor_gene_symbol"], term[ k ], "Down" , "R"))
                all_deg_total.add( row[ "receptor_gene_symbol"] )



for inx, row in lr_data.iterrows():
    l_id = row[  "ligand_gene_symbol" ]
    r_id = row[  "receptor_gene_symbol" ]
    if len( set([ l_id , r_id ]) & all_deg_total ) == 1:
        if l_id not in all_deg_total:
            all_degs_annot.append( ( row[ "ligand_gene_symbol"], "No DEG", "No DEG" , "L"))
        elif r_id not in all_deg_total:
            all_degs_annot.append( ( row[ "receptor_gene_symbol"], "No DEG", "No DEG" , "R"))




apd_pd = pd.DataFrame( set( all_deg_pairs) )
apd_pd.to_csv( "all_deg_pairs.csv")
apd_annot = pd.DataFrame( set( all_degs_annot) )
apd_annot.to_csv( "all_deg_lig_rec_annot.csv")



# Define custom color palette
# number of ligands

lr_sets12 = {}
num_lr = []
for ai in term.values():
    for r in [ "Up", "Down"]:
        k = (ai , r )
        tmp = apd_annot[ ( apd_annot[ 1] == ai ) & ( apd_annot[ 2] == r ) ]
        l_genes = tmp[ tmp[ 3] == "L"][ 0 ]
        r_genes = tmp[ tmp[3 ] == "R"][ 0 ]
        lr_sets12[ k ] = ( l_genes, r_genes )
        print( k, len( l_genes), len( r_genes ) )
        num_lr.append( ( k[0], k[1] +"-regulated", len( l_genes), len( r_genes) ) )

ai6 = {"VL":1, "L":2, "M":3, "H":4, "VH":5, "N":6 }

for k0 in lr_sets12:
    k0_l = lr_sets12[ k0][0]
    for k1 in lr_sets12:
        k1_r = lr_sets12[ k1][1]
        n = sum( apd_pd[0].isin( k0_l) & apd_pd[1].isin( k1_r))
        if n > 0 :
            print ( k0[1] + "_" + str( ai6[ k0[0] ] ) , k1[1] + "_" + str( ai6[ k1[0] ] ) , n  )
            #if "Down_2" ==  k1[1] + "_" + str( ai6[ k1[0] ] ):
            #    print( k1_r)

num_lr_df = pd.DataFrame( num_lr )
num_lr_df.columns = [ "AI level", "Regulation", "Number of Ligand", "Number of Receptor"]

palette = {
    "Up-regulated": ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"],
    "Down-regulated": ['#08519c','#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']
}

ai6 = {"VL":1, "L":2, "M":3, "H":4, "VH":5, "N":6 }

num_lr_df["color"] = num_lr_df.apply(lambda row: palette[row["Regulation"]][ ai6[ row["AI level"] ] -1], axis=1)


plt.figure(figsize=(6, 4))
bar = sns.barplot( data = num_lr_df , y = "Number of Ligand", x = "AI level", hue = "Regulation" , palette=['#fd8d3c', '#6baed6' ] )

for i,thisbar in enumerate(bar.patches):
    # Set a different hatch for each bar
    thisbar.set_color( deg_count_long["color"][i])
    thisbar.set_edgecolor("k")

plt.xlabel( "Allelic Imbalance Level" )
plt.title("Number of Ligands by AI Level and Regulation")
plt.legend([],[], frameon=False)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig( "Fig.5D.number_of_ligand.pdf")


plt.show()



# number of receptors
lr_sets12 = {}
num_lr = []
for ai in term.values():
    for r in [ "Up", "Down"]:
        k = (ai , r )
        tmp = apd_annot[ ( apd_annot[ 1] == ai ) & ( apd_annot[ 2] == r ) ]
        l_genes = tmp[ tmp[ 3] == "L"][ 0 ]
        r_genes = tmp[ tmp[3 ] == "R"][ 0 ]
        lr_sets12[ k ] = ( l_genes, r_genes )
        print( k, len( l_genes), len( r_genes ) )
        num_lr.append( ( k[0], k[1] +"-regulated", len( l_genes), len( r_genes) ) )

ai6 = {"VL":1, "L":2, "M":3, "H":4, "VH":5, "N":6 }

for k0 in lr_sets12:
    k0_l = lr_sets12[ k0][0]
    for k1 in lr_sets12:
        k1_r = lr_sets12[ k1][1]
        n = sum( apd_pd[0].isin( k0_l) & apd_pd[1].isin( k1_r))
        if n > 0 :
            print ( k0[1] + "_" + str( ai6[ k0[0] ] ) , k1[1] + "_" + str( ai6[ k1[0] ] ) , n  )
            #if "Down_2" ==  k1[1] + "_" + str( ai6[ k1[0] ] ):
            #    print( k1_r)

num_lr_df = pd.DataFrame( num_lr )
num_lr_df.columns = [ "AI level", "Regulation", "Number of Ligand", "Number of Receptor"]

palette = {
    "Up-regulated": ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"],
    "Down-regulated": ['#08519c','#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']
}

ai6 = {"VL":1, "L":2, "M":3, "H":4, "VH":5, "N":6 }

num_lr_df["color"] = num_lr_df.apply(lambda row: palette[row["Regulation"]][ ai6[ row["AI level"] ] -1], axis=1)


plt.figure(figsize=(6, 4))
bar = sns.barplot( data = num_lr_df , y = "Number of Receptor", x = "AI level", hue = "Regulation" , palette=['#fd8d3c', '#6baed6' ] )

for i,thisbar in enumerate(bar.patches):
    # Set a different hatch for each bar
    thisbar.set_color( deg_count_long["color"][i])
    thisbar.set_edgecolor("k")

plt.xlabel( "Allelic Imbalance Level" )
plt.title("Number of Receptors by AI Level and Regulation")
plt.legend([],[], frameon=False)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig( "Fig.5E.number_of_receptor.pdf")


plt.show()






#
# top 5%
tmp = []
lr_hubs = pd.read_csv( "./result/all_deg_pairs_network_anlaysis.csv")

for i in deg_ai_in[ "PCG"] :
    g = deg_ai_in["PCG"][ i ][0 ] | deg_ai_in["PCG"][ i ][1 ]
    print( i, len( g))
    #if "ENSG00000114251" in g:
    #    print ( i , "IN" )
    gn = gene_biotype_df[ gene_biotype_df[ "Gene stable ID"].isin(  g  ) ][ "Gene name"].to_list()
    a = lr_hubs[  lr_hubs[ "name"].isin( gn )].copy()
    a[ "AI level new"] = i
    for inx, row in a.iterrows():
        tmp.append( row )
    #print( i, a[ [ "AI level", "name"] ] )

tmp_df = pd.DataFrame( tmp )


lr_hubs = tmp_df.sort_values(by=["Outdegree", "4"], ascending=[False, False])
l_hubs = lr_hubs[ (lr_hubs[ "2"] != "No DEG") &( lr_hubs[ "4"] == "L")]

lr_hubs = tmp_df.sort_values(by=["Indegree", "4"], ascending=[False, False])
r_hubs = lr_hubs[ (lr_hubs[ "2"] != "No DEG") &( lr_hubs[ "4"] == "R")]

l_hubs.head()

#const colors = d3.schemeCategory10;
#const colors = d3.schemeCategory10;
palette = {
    "Up-regulated": ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"],
    "Down-regulated": ['#08519c', '#3182bd'    ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']
}

term = {  'VL' : 'Very Low' ,
'L' : 'Low' ,
"M" : 'Moderate' ,
'H' : 'High' ,
'VH' : 'Very High' ,
'N'  :'No cis-eQTL' }


# ligand hub genes
plot_data = l_hubs.sort_values( "Outdegree", ascending = False).head( 12)
tmp = []
for inx, row in plot_data.iterrows():
    tmp.append( row["2"] + " & " +  row["AI level new"]  )

plot_data[ "gene_ai_level"] = tmp


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
# 3. Plot
plt.figure(figsize=(3.5, 5))
bars = sns.barplot(data=plot_data, x="Outdegree", y="name"  , palette=colors , edgecolor="k")
plt.title("Top 5% hub ligands")
plt.xlabel("Number of Directed Edges")
plt.ylabel("Hub Ligand")
plt.xticks(rotation=90)  # rotate x-axis labels for better readability
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.show()



plot_data = r_hubs.sort_values( "Indegree", ascending = False).head( 15 )
tmp = []
for inx, row in plot_data.iterrows():
    tmp.append( row["2"] + " & " +  row["AI level new"]  )

plot_data[ "gene_ai_level"] = tmp



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
# 3. Plot
plt.figure(figsize=(3.5, 5))
bars = sns.barplot(data=plot_data, x="Indegree", y="name"  , palette=colors , edgecolor="k")
plt.title("Top 5% hub receptors")
plt.xlabel("Number of Directed Edges")
plt.ylabel("Hub Receptor")
plt.xticks(rotation=90)  # rotate x-axis labels for better readability
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.show()



#DEG transcript length
