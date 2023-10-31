
# make csv table of DARC-eQTL stratified by AI level
#

# load GWAS table

gwas_df = pd.read_csv('/mnt/d/sandbox/1.gxe/gwas/gwas_catalog_updated.tsv', sep='\t')
gwas_dhs_hit = gwas_df['DHS_Identifier'].value_counts().index.astype(str).to_list()[1:]

nr_snps = []
for inx, row in gwas_dhs_number_gwas.items():
    n = len( gwas_df[ gwas_df['DHS_Identifier'] == inx ]["SNPS"].unique())
    nr_snps.append( ( inx, n ) )

nr_snps = pd.DataFrame( nr_snps )
nr_snps[ 0].astype( str )
nr_snps.columns = ["DHS ID" , "count"]

# DHSs in each AI group

num_evar_lists = []
evar_lists = []
gene_lists = []
for bt_inx in ["PCG", "lncRNA" ]:
    genesets, afcsets = quantile_group(adipo_df, bt_inx )
    for j in range( 5):
        group = grouping[ j ]
        gset = genesets[ group ]
        evar = obtain_evar_geneset( evar_df, gset )
        num_var_nodhs = sum( evar["dhs_id"] == "N" )
        num_var_dhs = sum( (evar[ "dhs_id"] != "N")  )
        num_var_gwas_dhs = sum( evar[ "dhs_id"].isin( set( nr_snps[ "DHS ID"] ) )  )
        num_dhs = len( set( evar[ evar[ "dhs_id"] != "N"]["dhs_id"] ) )
        num_dhs_gwas = len( set( evar[ evar[ "dhs_id"] != "N"]["dhs_id"] ) & set( nr_snps[ "DHS ID"]) )
        num_evar_lists.append( (bt_inx, group, num_var_nodhs, num_var_dhs,num_var_gwas_dhs, num_dhs, num_dhs_gwas, len( gset ) ) )
        for g in evar["gene_id"].unique():
            num_var_nodhs1 = sum( ( evar["dhs_id"] == "N") & ( evar[ "gene_id"] == g ) )
            num_var_dhs1 = sum( (evar[ "dhs_id"] != "N")  & ( evar[ "gene_id"] == g ) )
            num_var_gwas_dhs1 = sum( evar[ "dhs_id"].isin( set( nr_snps[ "DHS ID"] ) ) & ( evar[ "gene_id"] == g ) )
            num_dhs1 = len( set( evar[ ( evar[ "dhs_id"] != "N" ) & ( evar[ "gene_id"] == g ) ]["dhs_id"] ) )
            num_dhs_gwas1 = len( set( evar[ ( evar[ "dhs_id"] != "N") & ( evar[ "gene_id"] == g ) ]["dhs_id"] ) & set( nr_snps[ "DHS ID"]) )
            gene_lists.append( (bt_inx, group, num_var_nodhs1, num_var_dhs1, num_var_gwas_dhs1, num_dhs1, num_dhs_gwas1, len( gset ) ) )


num_evar_lists = pd.DataFrame( num_evar_lists )
gene_lists = pd.DataFrame( gene_lists)

#GWAS-DHS-cis-eQTL
#DARC-eQTL

sns.histplot( gene_lists[ 3]/gene_lists[ 5] )
sns.histplot( gene_lists[ 4]/gene_lists[ 6] )
for bt_inx in ["PCG", "lncRNA" ]:
    for j in range( 5):
        group = grouping[ j ]
        tmp = gene_lists[ ( gene_lists[0] == bt_inx ) & ( gene_lists[ 1] == group ) ]
        n0, n1 = tmp[ 3]/tmp[ 5], tmp[ 4]/tmp[ 6]
        print( bt_inx, group, np.mean( n0[ n0 > 0] ), np.mean( n1[ n1 > 0 ] ) ,
            scipy.stats.mannwhitneyu( n0[ n0 > 0] , n1[ n1 > 0]  , nan_policy="omit"   ).pvalue )


# test each gene - chi-X2 test
for inx, row in gene_lists.items():
    n3 = row[ 3:7 ]
    break


# test DARC-DHS ratio to total DHS for each genes
# chi-squre test
# a. number of total cis-eQTL
# b. number of cis-eQTL within total DHSs
# c. number of cis-eQTL within DHSs associated with GWAS
#
