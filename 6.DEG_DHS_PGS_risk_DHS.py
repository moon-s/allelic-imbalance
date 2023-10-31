# PGS risk variants from PGS catalog
#

# GWAS + PRS relative risk
#
gwas_df = pd.read_csv('/mnt/d/sandbox/1.gxe/gwas/gwas_catalog_updated.tsv', sep='\t', index=False)
gwas_dhs_hit = gwas_df['DHS_Identifier'].value_counts().index.astype(str).to_list()[1:]


g_cate = [ "PCG", "lncRNA" ]
dhs_lists = {}
res = []
dhs_all = set()
for bt_inx in range( 2 ):
    bt = g_cate[ bt_inx ]
    grouping = [ "Very Low", "Low", "Moderate", "High", "Very High"]
    dhs_lists[ bt ] = {}
    for j in grouping :
        gset_up, gset_down = deg_ai_in[ bt][ j ]
        evar_up = obtain_evar_geneset( evar_df, gset_up )
        evar_down = obtain_evar_geneset( evar_df, gset_down )
        data_up = evar_up[~evar_up['dhs_id'].str.startswith('N')]
        data_down = evar_down[~evar_down['dhs_id'].str.startswith('N')]
        dhs_up = data_up[ "dhs_id"].unique()
        dhs_down = data_down[ "dhs_id"].unique()
        if j not in dhs_lists[ bt ]:
            dhs_lists[ bt ][ j ] = [ set(), set() ]
        dhs_lists[ bt ][ j ][0 ].update( set( dhs_up ) )
        dhs_lists[ bt ][ j ][1 ].update( set( dhs_down ) )
        dhs_all.update( set( dhs_up) & set( gwas_dhs_hit) )
        dhs_all.update( set( dhs_down) & set( gwas_dhs_hit) )
        print( bt, j, set( dhs_up) & set( gwas_dhs_hit), set( dhs_down) & set( gwas_dhs_hit) )

#
pgsid_togo = [ ['PGS003494', 	'Body measurement'	,'Birth weight'] ,
['PGS002850', 	'Body measurement'	 ,'BMI'] ,
['PGS002133', 	'Body measurement'	 ,'Body fat percentage'] ,
['PGS002967', 	'Body measurement'	 ,'Height'] ,
['PGS002781', 	'Lipid or lipoprotein measurement'	 ,'HDL cholesterol'] ,
['PGS003029', 	'Lipid or lipoprotein measurement'	 ,'LDL'] ,
['PGS003139', 	'Lipid or lipoprotein measurement'	 ,'Total Cholesterol'] ,
['PGS003149', 	'Lipid or lipoprotein measurement'	 ,'Triglycerides'] ,
['PGS002027', 	'Metabolic disorder'	 ,'Diabetic retinopathy'] ,
['PGS002033',   'Metabolic disorder'	 ,'Overweight' ],
['PGS002098', 	'Metabolic disorder'	 ,'Recent poor appetite'] ,
['PGS002025',	'Metabolic disorder'	 ,'T1D'] ,
['PGS003092', 	'Metabolic disorder'	 ,'T2D'] ,
['PGS002237', 	'Other disease'	 ,'Chronic Kidney disease'] ,
['PGS000748', 	'Other disease'	 ,'Coronary artery disease'] ,
['PGS003584', 	'Other disease'	 ,'Major depressive disorder'] ,
['PGS003209',   'Other disease'  , 'Sleep apnea'],
 ]

pgsid_togo = pd.DataFrame( pgsid_togo )

# DHSs associated with subset of genes
bt = "PCG"
subsets_nr = {} #  { 'Very Low': dhs_subsets[0], 'Low': dhs_subsets[1], 'Moderate': dhs_subsets[2],  'High': dhs_subsets[3],  'Very High': dhs_subsets[4] }
for k in dhs_lists[ "PCG"]:
    subsets_nr[ k ] = dhs_lists[ "PCG" ][ k ][ 1]  # up-regulated 0: , down-regulated : 1


#
# low/high impact varaitns
#
risk_fc = []
directory2pgs = "/mnt/d/sandbox/1.gxe/dhs/prs_dhsid/"

for inx, row in pgsid_togo.iterrows() :
#if file.endswith(".csv.gz"):
    pgs_id, trait_category, trait = row
    print( bt,pgs_id, trait, trait_category )
        # Extract trait_reported from the header
        # Load dataset with 'region' as the index
    file = pgs_id + "_data.csv.gz"
    df_tissue = pd.read_csv(directory2pgs + file, sep=",", compression='gzip', skiprows=3, index_col='region', dtype = {'region': str, "hm_chr": str})

    q100 = np.quantile( df_tissue["effect_weight"] , [ 0.05, 0.1, 0.25, 0.33, 0.66, 0.75, 0.9, 0.95]   )

        # Remove 'Non-DHS' variants
    df_nondhs = df_tissue[df_tissue.index == 'Non-DHS']
    df = df_tissue[df_tissue.index != 'Non-DHS']

    for ai_level in subsets_nr:
        effect_weights_nodhs = df_nondhs[ "effect_weight"]
        effect_weights_nodhs_low = effect_weights_nodhs[ effect_weights_nodhs <= q100[2] ]
        effect_weights_nodhs_high = effect_weights_nodhs[ effect_weights_nodhs >= q100[5] ]
            # Select 'effect_weight' for the regions in each DHS set
        effect_weights_dhs = df.loc[df.index.isin( subsets_nr[ ai_level ]), 'effect_weight']
        effect_weights_dhs_low = effect_weights_dhs[ effect_weights_dhs <= q100[2] ]
        effect_weights_dhs_high = effect_weights_dhs[ effect_weights_dhs >= q100[5] ]

        pvalue_test_low_nodhs_vs_dhs = scipy.stats.ttest_ind(effect_weights_nodhs_low, effect_weights_dhs_low, nan_policy='omit').pvalue
        pvalue_test_high_nodhs_vs_dhs = scipy.stats.ttest_ind(effect_weights_nodhs_high, effect_weights_dhs_high, nan_policy='omit').pvalue

        diff_low =  np.mean(  effect_weights_dhs_low) - np.mean( effect_weights_nodhs_low)
        diff_high =  np.mean(  effect_weights_dhs_high) - np.mean( effect_weights_nodhs_high)

        nlow_nodhs , nhigh_nodhs = len( effect_weights_nodhs_low), len( effect_weights_nodhs_high)
        nlow_dhs , nhigh_dhs = len( effect_weights_dhs_low), len( effect_weights_dhs_high)
        ctable = ( [nlow_dhs , nhigh_dhs ], [ nlow_nodhs , nhigh_nodhs] )
        if ( nlow_dhs + nhigh_dhs >= 10 ) & ( nlow_dhs* nhigh_dhs > 0 ):
            try:
                _, pvalue_low_vs_high_risk , _, _ = chi2_contingency( ctable )
            except ValueError:
                pvalue_low_vs_high_risk = 1

            risk_fc.append( ( trait_category, trait ,ai_level, diff_low ,  diff_high, pvalue_test_low_nodhs_vs_dhs, pvalue_test_high_nodhs_vs_dhs,
              nlow_nodhs,  nhigh_nodhs,  nlow_dhs, nhigh_dhs , pvalue_low_vs_high_risk ) )



risk_df = pd.DataFrame( risk_fc)
risk_df.columns = ["Trait category", 'Trait', 'AI Level', 'Mean Diff(1Q)', 'Mean Diff(4Q)', 'P(1Q)', 'P(4Q)',
       'n(low_ctl)', 'n(high_ctl)', 'n(low_test)', 'n(high_test)', 'p(odd)']


risk_df .head()
import statsmodels.stats.multitest as multi
pcombined = []
aialpha = []
alpha_dict = { 'Very Low': 1, 'Low': 0.8, 'Moderate': 0.6,  'High': 0.4,  'Very High': 0.2}


for i, row in risk_df.iterrows():
    p1, p2 = row[ "P(1Q)"] , row[ "P(4Q)"]
    p12 = scipy.stats.combine_pvalues( [ p1, p2 ] ).pvalue
    #print( p12)
    pcombined.append( p12 )
    alph = alpha_dict[ row[ "AI Level"]]
    aialpha.append( alph )

results_corrected = multi.multipletests(  np.asarray( pcombined), method='fdr_bh')[1]
risk_df[ "-log10(adj.P(1Q,4Q))"] = -np.log10( pcombined )
risk_df[ "AI_alpha"] = aialpha

podd_corrected = multi.multipletests(  np.asarray( risk_df[ "p(odd)" ]), method='fdr_bh')[1]
risk_df[ "-log10(adj.p(Odd))"] = -np.log10( podd_corrected )

risk_df = risk_df.fillna(0)


#
# at least 10 variants
#

risk_fc = pd.read_csv( "risk_ratio_lncRNA_up.csv")
plt.figure(figsize=( 4, 4) )

p = sns.scatterplot( data = risk_fc, x = "Mean Diff(1Q)", y = "Mean Diff(4Q)", hue = "Trait", alpha = risk_fc["AI_alpha"], size = risk_fc["-log10(adj.P(1Q,4Q))"], palette=sns.color_palette("husl", 15)[2:])
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
p.axhline( 0 , ls = "--", color = "grey" )
p.axvline( 0 , ls = "--", color = "grey" )



df_so = pd.DataFrame(np.random.randint(0,100,size=(20, 4)), columns=list('ABCD'))
scatter_so=sns.lmplot(x='C', y='D', data=df_so,
           fit_reg=False,y_jitter=0, scatter_kws={'alpha':0.2})

scatter_so = sns.lmplot( data = risk_fc, x = "Mean Diff(1Q)", y = "Mean Diff(4Q)", hue = "Trait",  size = risk_fc["-log10(adj.P(1Q,4Q))"]*20, palette=sns.color_palette("husl", 19)[2:])

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
p.axhline( 0 , ls = "--", color = "grey" )
p.axvline( 0 , ls = "--", color = "grey" )


def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], str(point['val']))

label_point(risk_fc["Mean Diff(1Q)"], risk_fc["Mean Diff(4Q)"],  risk_fc["trait"], plt.gca())



 # GWAS Catalog
gwas_df = pd.read_csv('/mnt/d/sandbox/1.gxe/gwas/gwas_catalog_updated.tsv', sep='\t')

gwas_df[ gwas_df["DHS_Identifier"].astype(str).isin( gwas_dhs_hit) ][ "MAPPED_TRAIT"].value_counts()

# gwas_df[ gwas_df["DHS_Identifier"].astype(str).isin( gwas_dhs_hit) ][ "MAPPED_TRAIT"].value_counts()
gwas_dhs_hit = gwas_df['DHS_Identifier'].value_counts().index.astype(str).to_list()
gwas_3df = gwas_df[ [ "DHS_Identifier", "MAPPED_TRAIT", "SNPS"] ]

g_cate = [ "PCG", "lncRNA" ]
dhs_lists = {}
dhs_no_lists = {}
res = []
dhs_all = set()
for bt_inx in range( 2 ):
    bt = g_cate[ bt_inx ]
    grouping = [ "Very Low", "Low", "Moderate", "High", "Very High"]
    dhs_lists[ bt ] = {}
    dhs_no_lists[ bt ] = {}
    for j in grouping :
        gset_up, gset_down = deg_ai_in[ bt][ j ]
        evar_up = obtain_evar_geneset( evar_df, gset_up )
        evar_down = obtain_evar_geneset( evar_df, gset_down )
        data_up = evar_up[~evar_up['dhs_id'].str.startswith('N')]
        data_down = evar_down[~evar_down['dhs_id'].str.startswith('N')]
        dhs_up = data_up[ "dhs_id"].unique()
        dhs_down = data_down[ "dhs_id"].unique()

        gwas_df_up = gwas_3df[ gwas_3df[ "DHS_Identifier"].astype(str).isin( dhs_up ) ]
        gwas_df_down = gwas_3df[ gwas_3df[ "DHS_Identifier"].astype(str).isin( dhs_down ) ]

        gwas_dhs_up = gwas_df_up[ "DHS_Identifier"]
        gwas_dhs_down =  gwas_df_down[ "DHS_Identifier"]

        gwas_snp_up = gwas_df_up[ "SNPS"].unique()
        gwas_snp_down = gwas_df_down[ "SNPS"].unique()


        if j not in dhs_lists[ bt ]:
            dhs_lists[ bt ][ j ] = [ set(), set() ]
            dhs_no_lists[ bt ][ j ] = [ set(), set() ]
        dhs_lists[ bt ][ j ][0 ].update( set( gwas_dhs_up) )
        dhs_lists[ bt ][ j ][1 ].update( set( gwas_dhs_down) )
        dhs_no_lists[ bt ][ j ][0 ].update( set( dhs_up) - set( gwas_dhs_up) )
        dhs_no_lists[ bt ][ j ][1 ].update( set( dhs_down) - set( gwas_dhs_down) )


        dhs_all.update( set( gwas_dhs_up ) )
        dhs_all.update( set( gwas_dhs_down)  )
        res.append( ( bt, j, "Up-regulated", len( gwas_snp_up ), len( set( gwas_dhs_up) ), len( set( dhs_up)) ))
        res.append(   (bt, j, "Down-regulated", len( gwas_snp_down ), len(  set( gwas_dhs_down) ), len( set( dhs_down))  ) )
        print(  ( bt, j, len( gwas_snp_up ), len( set( gwas_dhs_up) ), len( set( dhs_up)),
                     len( gwas_snp_down ), len(  set( gwas_dhs_down) ), len( set( dhs_down))  )  )




res_df = pd.DataFrame( res )
res_df.columns = [ "Gene type", "AI level", "Number of risk variants(up)", 'Number of risk DHS(up)', 'Total DHS(up)',
                  "Number of risk variants(down)", 'Number of risk DHS(down)', 'Total DHS(down)']

# bar chart 1 -> top bars (group of 'smoker=No')
bar1 = sns.barplot(x="AI level",  y="Number of risk DHS(up)", data=res_df, color='darkblue')

# bar chart 2 -> bottom bars (group of 'smoker=Yes')
bar2 = sns.barplot(x="AI level", y="Total DHS(up)", data=res_df, color='lightblue')
# add legend
top_bar = mpatches.Patch(color='darkblue', label='GWAS hit')
bottom_bar = mpatches.Patch(color='lightblue', label='No GWAS')
plt.legend(handles=[top_bar, bottom_bar])



# number non-redundant SNPs associaed with complex traits
# Very low to very high, and no cis-eQTL
palette = {
    "Up-regulated": ['#a63603','#e6550d' ,'#fd8d3c' ,'#fdae6b','#fdd0a2',  "#feedde"],
    "Down-regulated": ['#08519c','#3182bd'  ,'#6baed6'  ,'#9ecae1' ,'#c6dbef', '#eff3ff']
}



res_df["color"] = res_df.apply(lambda row: palette[row["Regulation"]][res_df["AI level"].tolist().index(row["AI level"])], axis=1)

# Plot
plt.figure(figsize=(6, 4))
bar = sns.barplot(x="AI level", y="Number of risk variants", hue="Regulation", data=res_df, edgecolor="white", palette=['#fd8d3c', '#6baed6' ])#  palette=deg_count_long["color"].unique())

for i,thisbar in enumerate(bar.patches):
    # Set a different hatch for each bar
    thisbar.set_color( deg_count_long["color"][i])
    thisbar.set_edgecolor("k")
plt.legend( title= "DEG in obese ATs")
plt.xlabel( "Allelic Imbalance Level" )
plt.title("Number of GWAS-associated variants within DHSs in DEGs")
sns.despine( offset=10, trim=True)
plt.tight_layout()
plt.savefig( "Fig.3B.DEG_number_AI_PCG.pdf")
plt.show()

#
# test proprotion - risk susceptablity
#
#PCG ( up-gene, down-gene), (up-SNPs, down-SNPs
scipy.stats.chi2_contingency( [ [ 62, 53  ], [172, 92  ]]  )
scipy.stats.chi2_contingency( [ [ 165, 104  ], [ 487, 328  ]]  )
scipy.stats.chi2_contingency( [ [ 356, 161  ], [1270, 609  ]]  )
scipy.stats.chi2_contingency( [ [ 470, 246  ], [2591, 976  ]]  )
scipy.stats.chi2_contingency( [ [ 530, 254  ], [2556, 1335  ]]  )

#lncRNA
scipy.stats.chi2_contingency( [ [ 33, 33  ], [  69, 108  ]]  )
scipy.stats.chi2_contingency( [ [ 27, 18  ], [ 118, 70  ]]  )
scipy.stats.chi2_contingency( [ [ 39, 25  ], [ 68,  102 ]]  )
scipy.stats.chi2_contingency( [ [ 27, 15  ], [ 57,  28 ]]  )
scipy.stats.chi2_contingency( [ [ 36, 28  ], [ 998, 85  ]]  )


# number of DHS with GWAS-associaed SNPs to total DHS

# Plot

res_df[ "Fraction"] = res_df["Number of risk DHS"]/res_df["Total DHS"]
plt.figure(figsize=(5, 4))
bar = sns.barplot(x="AI level", y="Number of DHS", hue="Regulation", data=res_df[ res_df["Gene type"] =="PCG"] , edgecolor="white", palette=['#fd8d3c', '#6baed6' ], errorbar=None)#  palette=deg_count_long["color"].unique())

pal =palette[ "Up-regulated"] + palette["Down-regulated"]
for i,thisbar in enumerate(bar.patches):
    # Set a different hatch for each bar
    print( thisbar)
    thisbar.set_color( pal[i])
    thisbar.set_edgecolor("k")

plt.legend( title= "DEG in obese ATs")
plt.ylabel( "Number of GWAS-associated DHSs" )
sns.despine( offset=10, trim=True)

plt.xlabel( "Allelic Imbalance Level" )
plt.title("Protein-coding genes")
sns.despine( offset=10, trim=True)
plt.tight_layout()
plt.savefig( "Fig.7C.Number_GWAS_DHS_PCG.pdf")

# DHSs in PCG
scipy.stats.chi2_contingency( [ [  157,  1016 - 157  ], [ 89 , 1160 - 89   ]]  )
scipy.stats.chi2_contingency( [ [ 416  , 3246 - 416   ], [ 306, 2962 - 306  ]]  )
scipy.stats.chi2_contingency( [ [ 1134 , 9389 - 1134  ], [ 535, 5492 -535  ]]  )
scipy.stats.chi2_contingency( [ [ 1999 , 14157 - 1999   ], [ 894, 8084  ]]  )
scipy.stats.chi2_contingency( [ [  2048, 17550   ], [ 1084, 10567  ]]  )

# DHSs in lncRNA
scipy.stats.chi2_contingency( [ [ 68 , 719 -68   ], [ 100 ,  1139 - 100  ]]  )
scipy.stats.chi2_contingency( [ [ 105  , 631 - 105   ], [ 59 , 582 - 59   ]]  )
scipy.stats.chi2_contingency( [ [ 67 ,  1202 - 67   ], [ 97, 911-97  ]]  )
scipy.stats.chi2_contingency( [ [ 54 ,  936 - 54   ], [ 27, 545 - 27  ]]  )
scipy.stats.chi2_contingency( [ [ 553  , 2453 - 553    ], [ 81  ,  1321 - 81  ]]  )




# difference in per-allele effect weight between GWAS-associed DHS and other DHS for Q1




 #
 pgsid_togo = [ ['PGS003494', 	'Body measurement'	,'Birth weight'] ,
 ['PGS002850', 	'Body measurement'	 ,'BMI'] ,
 ['PGS002133', 	'Body measurement'	 ,'Body fat percentage'] ,
 ['PGS002967', 	'Body measurement'	 ,'Height'] ,
 ['PGS002781', 	'Lipid or lipoprotein measurement'	 ,'HDL cholesterol'] ,
 ['PGS003029', 	'Lipid or lipoprotein measurement'	 ,'LDL'] ,
 ['PGS003139', 	'Lipid or lipoprotein measurement'	 ,'Total Cholesterol'] ,
 ['PGS003149', 	'Lipid or lipoprotein measurement'	 ,'Triglycerides'] ,
 ['PGS002027', 	'Metabolic disorder'	 ,'Diabetic retinopathy'] ,
 ['PGS002033',   'Metabolic disorder'	 ,'Overweight' ],
 ['PGS002098', 	'Metabolic disorder'	 ,'Recent poor appetite'] ,
 ['PGS002025',	'Metabolic disorder'	 ,'T1D'] ,
 ['PGS003092', 	'Metabolic disorder'	 ,'T2D'] ,
 ['PGS002237', 	'Other disease'	 ,'Chronic Kidney disease'] ,
 ['PGS000748', 	'Other disease'	 ,'Coronary artery disease'] ,
 ['PGS003584', 	'Other disease'	 ,'Major depressive disorder'] ,
 ['PGS003209',   'Other disease'  , 'Sleep apnea'],
  ]

 pgsid_togo = pd.DataFrame( pgsid_togo )

 # DHSs associated with subset of genes
 bt = "PCG"
 subsets_nr = {} #  { 'Very Low': dhs_subsets[0], 'Low': dhs_subsets[1], 'Moderate': dhs_subsets[2],  'High': dhs_subsets[3],  'Very High': dhs_subsets[4] }
 for k in dhs_lists[ "PCG"]:
     subsets_nr[ k ] = dhs_lists[ "PCG" ][ k ][ 1]  # up-regulated 0: , down-regulated : 1


 #
 # low/high impact varaitns
 #
risk_fc = []
directory2pgs = "/mnt/d/sandbox/1.gxe/dhs/prs_dhsid/"
results = []
for inx, row in pgsid_togo.iterrows() :
 #if file.endswith(".csv.gz"):
    pgs_id, trait_category, trait = row
    print( "***",pgs_id, trait, trait_category )
         # Extract trait_reported from the header
         # Load dataset with 'region' as the index
    file = pgs_id + "_data.csv.gz"
    df_tissue = pd.read_csv(directory2pgs + file, sep=",", compression='gzip', skiprows=3, index_col='region', dtype = {'region': str, "hm_chr": str})

    q100 = np.quantile( df_tissue["effect_weight"] , [ 0.05, 0.1, 0.25, 0.33, 0.66, 0.75, 0.9, 0.95]   )
         # Remove 'Non-DHS' variants
    df_nondhs = df_tissue[df_tissue.index == 'Non-DHS']
    df = df_tissue[df_tissue.index != 'Non-DHS']
    for bt in dhs_lists:
        for ai_level in dhs_lists[ bt ]:
            effect_weights_nodhs = df_nondhs[ "effect_weight"]
            ew_dhs_gwas_up = df.loc[df.index.isin( [ str(x) for x in dhs_lists[bt][ ai_level ][0]] ), 'effect_weight']
            ew_dhs_no_up = df.loc[df.index.isin( [ str(x) for x in dhs_no_lists[bt][ ai_level ][0]] ), 'effect_weight']

            ew_dhs_gwas_up = ew_dhs_gwas_up[ ew_dhs_gwas_up > q100[5]]
            ew_dhs_no_up = ew_dhs_no_up[ ew_dhs_no_up > q100[5]]

            pvalue_test_dhs_gwas_vs_no_up = scipy.stats.ttest_ind(ew_dhs_gwas_up, ew_dhs_no_up, nan_policy='omit').pvalue
            means =  np.mean(  ew_dhs_gwas_up ) , np.mean( ew_dhs_no_up)
            sems =   np.std( ew_dhs_gwas_up  ) , np.std( ew_dhs_no_up )
            nums =   len( ew_dhs_gwas_up  ) , len( ew_dhs_no_up )
            results.append( ( bt, "Q4", ai_level, 'Up-regulated', means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_up,pgs_id, trait_category, trait ))
            if pvalue_test_dhs_gwas_vs_no_up < 0.05:
                if means[0] > means[ 1]:
                    print( bt, "Q4", ai_level, 'Up-regulated',  means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_up )

            ew_dhs_gwas_down = df.loc[df.index.isin( [ str(x) for x in dhs_lists[bt][ ai_level ][1]] ), 'effect_weight']
            ew_dhs_no_down = df.loc[df.index.isin( [ str(x) for x in dhs_no_lists[bt][ ai_level ][1]] ), 'effect_weight']

            ew_dhs_gwas_down = ew_dhs_gwas_down[ ew_dhs_gwas_down > q100[5]]
            ew_dhs_no_down = ew_dhs_no_down[ ew_dhs_no_down > q100[5]]

            # Select 'effect_weight' for the regions in each DHS set
            means =  np.mean(  ew_dhs_gwas_down ) , np.mean( ew_dhs_no_down)
            sems =   np.std( ew_dhs_gwas_down  ) , np.std( ew_dhs_no_down )
            nums =   len( ew_dhs_gwas_down  ) , len( ew_dhs_no_down )

            pvalue_test_dhs_gwas_vs_no_down = scipy.stats.ttest_ind(ew_dhs_gwas_down, ew_dhs_no_down, nan_policy='omit').pvalue
            results.append( ( bt, "Q4", ai_level,"Down-regulated", means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_down,pgs_id, trait_category, trait ))
            if pvalue_test_dhs_gwas_vs_no_down < 0.05:
                if means[ 0] > means[ 1]:
                    print( bt, "Q4", ai_level,"Down-regulated",  means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_down )



for inx, row in pgsid_togo.iterrows() :
 #if file.endswith(".csv.gz"):
    pgs_id, trait_category, trait = row
    print( "***",pgs_id, trait, trait_category )
         # Extract trait_reported from the header
         # Load dataset with 'region' as the index
    file = pgs_id + "_data.csv.gz"
    df_tissue = pd.read_csv(directory2pgs + file, sep=",", compression='gzip', skiprows=3, index_col='region', dtype = {'region': str, "hm_chr": str})

    q100 = np.quantile( df_tissue["effect_weight"] , [ 0.05, 0.1, 0.25, 0.33, 0.66, 0.75, 0.9, 0.95]   )

         # Remove 'Non-DHS' variants
    df_nondhs = df_tissue[df_tissue.index == 'Non-DHS']
    df = df_tissue[df_tissue.index != 'Non-DHS']
    for bt in dhs_lists:
        for ai_level in dhs_lists[ bt ]:
            effect_weights_nodhs = df_nondhs[ "effect_weight"]
            ew_dhs_gwas_up = df.loc[df.index.isin( [ str(x) for x in dhs_lists[bt][ ai_level ][0]] ), 'effect_weight']
            ew_dhs_no_up = df.loc[df.index.isin( [ str(x) for x in dhs_no_lists[bt][ ai_level ][0]] ), 'effect_weight']

            ew_dhs_gwas_up = ew_dhs_gwas_up[ ew_dhs_gwas_up < q100[2]]
            ew_dhs_no_up = ew_dhs_no_up[ ew_dhs_no_up < q100[2]]

            pvalue_test_dhs_gwas_vs_no_up = scipy.stats.ttest_ind(ew_dhs_gwas_up, ew_dhs_no_up, nan_policy='omit').pvalue
            means =  np.mean(  ew_dhs_gwas_up ) , np.mean( ew_dhs_no_up)
            sems =   np.std( ew_dhs_gwas_up  ) , np.std( ew_dhs_no_up )
            nums =   len( ew_dhs_gwas_up  ) , len( ew_dhs_no_up )

            results.append( ( bt, "Q1", ai_level, 'Up-regulated',  means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_up,pgs_id, trait_category, trait ))
            if pvalue_test_dhs_gwas_vs_no_up < 0.05:
                if means[0] < means[ 1]:
                    print( bt, "Q1", ai_level, 'Up-regulated',  means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_up )

            ew_dhs_gwas_down = df.loc[df.index.isin( [ str(x) for x in dhs_lists[bt][ ai_level ][1]] ), 'effect_weight']
            ew_dhs_no_down = df.loc[df.index.isin( [ str(x) for x in dhs_no_lists[bt][ ai_level ][1]] ), 'effect_weight']

            ew_dhs_gwas_down = ew_dhs_gwas_down[ ew_dhs_gwas_down < q100[2]]
            ew_dhs_no_down = ew_dhs_no_down[ ew_dhs_no_down < q100[2]]

            # Select 'effect_weight' for the regions in each DHS set

            pvalue_test_dhs_gwas_vs_no_down = scipy.stats.ttest_ind(ew_dhs_gwas_down, ew_dhs_no_down, nan_policy='omit').pvalue
            means =  np.mean(  ew_dhs_gwas_down ) , np.mean( ew_dhs_no_down )
            sems =   np.std( ew_dhs_gwas_down  ) , np.std( ew_dhs_no_down )
            nums =   len( ew_dhs_gwas_down  ) , len( ew_dhs_no_down )

            results.append( ( bt, "Q1", ai_level,"Down-regulated",  means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_down,pgs_id, trait_category, trait ))
            if pvalue_test_dhs_gwas_vs_no_down < 0.05:
                if means[ 0] < means[ 1]:
                    print( bt, "Q1", ai_level,"Down-regulated", means[0], sems[0], nums[0], means[1], sems[1], nums[1], pvalue_test_dhs_gwas_vs_no_down )





results_df = pd.DataFrame( results )
results_df.columns = [ "Gene type", "Quantile", "AI level", "Regulation", "Average PRS_GWAS", "Average PRS no GWAS",
                      "SEM PRS_GWAS", "SEM PRS no GWAS", "p-value", "pgs_id", "Trait category", "Trait" ]


results_df_ = results_df.fillna( 1)
multi.multipletests(  results_df_["p-value"] , method='fdr_bh')[1]
sum( results_df_["p-value"] == None )
results_df[ "-log10(Adj.p-value"] = -np.log10( multi.multipletests(  results_df_["p-value"] , method='fdr_bh')[1] )




*** PGS003584 Major depressive disorder Other disease
PCG Low Up-regulated (0.006973740965909091, 0.005550714014760147) 1.3132649335321967e-20
PCG High Up-regulated (0.006828129577464788, 0.0059822097444633725) 3.3215466316726337e-25
PCG Very High Up-regulated (0.006258394294478527, 0.0052874271506849315) 4.527376485388252e-42
PCG Very High Down-regulated (0.006258394294478527, 0.0052874271506849315) 3.8954135083373465e-18
lncRNA Very Low Up-regulated (0.004567397777777778, 0.004017902414860681) 0.031683723466584574
lncRNA Very Low Down-regulated (0.004567397777777778, 0.004017902414860681) 0.014416649563395953
lncRNA Very High Up-regulated (0.0072922064961832055, 0.006810478971088435) 4.52124870867291e-06
lncRNA Very High Down-regulated (0.0072922064961832055, 0.006810478971088435) 0.020544990914052774

*** PGS003584 Major depressive disorder Other disease
PCG Very Low Down-regulated (-0.0048908084375, -0.004622692942942943) 0.001118989227059698
PCG Low Up-regulated (-0.005851962576419214, -0.00493997400388727) 3.213335952585568e-06
PCG Moderate Up-regulated (-0.005283311245059289, -0.004982601891402715) 0.0013454661288568664
PCG High Up-regulated (-0.007203271103658537, -0.00602270999821013) 5.803180295547001e-28
PCG High Down-regulated (-0.007203271103658537, -0.00602270999821013) 0.0005158864298859066
PCG Very High Up-regulated (-0.006176795314043754, -0.005313189035131123) 7.949687523889683e-26
PCG Very High Down-regulated (-0.006176795314043754, -0.005313189035131123) 1.6267038269872199e-38
lncRNA Very High Up-regulated (-0.007130915745310958, -0.006559921628555493) 0.000169335243348951




#

def DH_GWAS_SMD_plot( bt, quant, regul ):
    bt, quant, regul = "PCG", "Q1", "Up-regulated"

    plt.figure(figsize=(5, 10))
    condi =  ( results_df[ "Gene type"] == bt ) & ( results_df[ "Quantile"] ==quant ) & ( results_df[ "Regulation"] =="Up-regulated" )
    bar = sns.barplot( data = results_df[condi  ],
              y = "Trait", x = "SMD", hue = "AI level", palette=palette["Up-regulated"] )
              #plt.xscale( "symlog")

    pointsize = [x*10 if x > -np.log10( 0.05) else 0 for x in results_df[condi][ "-log10(Adj.p-value"]  ]

    x_coords = [p.get_x() +  p.get_width() for p in bar.patches]
    y_coords = [ p.get_y() + p.get_height() for p in bar.patches]

    pvalue_coords = []
    for inx, thisbar in  enumerate(bar.patches):
        #print( inx, thisbar)
        w = thisbar.get_width()
        pos = thisbar.get_y() + thisbar.get_height()*0.5
        row = results_df[condi][ results_df[ condi][ "SMD"].isin( [w] ) ]
        pvalue_coords.append( ( pos, row["AI level"].tolist()[0] , float( row[ "-log10(Adj.p-value)"]) ))

    pvalue_coords_df = pd.DataFrame( pvalue_coords )
    pvalue_coords_df.columns = [ "pos", "AI level", "-log10(Adj.p-value)" ]
    plt.tick_params(axis = 'y', left = True, right = False, labelleft = True, labelright = False)

    pvalue_coords_df[ "pointsize"] =  [x*10 if x > -np.log10( 0.05) else 0 for x in pvalue_coords_df[ "-log10(Adj.p-value)"]  ]

    print( pvalue_coords_df.head(), sum( pvalue_coords_df[ "-log10(Adj.p-value)"] > 1) )

    sns.scatterplot( data= pvalue_coords_df, y = "pos", x = 0.5, hue="AI level", size= "pointsize" , sizes = (1, 70), edgecolor= "k", linewidth=0.5,  legend=False, palette=palette["Up-regulated"])

    #print( x_coords[0:3],  results_df[ condi & ( results_df["-log10(Adj.p-value)"] > -np.log10( 0.05) ) ] )
    plt.legend( [], [], frameon=False)
    plt.tight_layout()
    sns.despine( offset=10, trim=True)
    plt.xlim( -0.7, 0.7)

    plt.savefig( "Fig.7E.DHS_GWAS_SMD_%s_%s_%s.pdf" % (bt, quant, regul ) )


# childhood BMI
eg = [ 'LEPR', 'GLP1R', 'PCSK1', 'KLF14' ]
KLF14 = ( "down-regulated", "High AI" )

#
# thermogenesis
# Very low AI
FTO = "Very Low"



imprinted_DEGs. p-value = 0.0107
