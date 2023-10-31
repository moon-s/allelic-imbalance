#
# PRS loci within out side of DHS
#

import gzip
import os
import glob
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns


# Load DHS data and create an interval tree for quick overlap checking
dhs_data = pd.read_csv('/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz', sep='\t',usecols=[0,1,2, 3] )
dhs_data.columns = ['chrom', 'start', 'end', "identifier" ]
dhs_data['chrom'] = dhs_data['chrom'].apply(str).str.split(pat = "r" , expand = True)[ [1 ] ]
dhs_tree = {chrom: IntervalTree() for chrom in dhs_data['chrom'].unique()}


dhs_id = {}
for index, row in dhs_data.iterrows():
    dhs_tree[row['chrom']].addi(row['start'], row['end'])
    dhs_id[ ( row["chrom"], row['start'], row['end']) ] = row["identifier"]




ntrait = 0
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/*.txt.gz'):
    with gzip.open(filename, 'rt') as f:
        # Read the header lines
        #trait_name = f.readline().split('=')[1].strip()
        #num_variants = int(f.readline().split('=')[1].strip())
        num = 0
        for l in f:
            if "#variants_number=" in l:
                num = int( l.split("=")[1].strip() )
            elif "#trait_reported=" in l:
                trait = l.split( "=")[1].strip()
            elif "#pgs_id=" in l:
                pgs_id = l.split( "=")[1].strip()
            elif "HmPOS_match_pos" in l:
                break

        if ( num > 100) : #  & (num < 100000) :
            headlines = "#pgs_id=%s\n" % ( pgs_id )
            headlines += "#varinats_number=%d\n" % ( num )
            headlines += "#trait_reported=%s\n" % ( trait )
            ntrait += 1
            # Load the file into a dataframe
            pgs_scores = pd.read_csv(f, sep='\t', comment='#', low_memory=False, dtype = {"hm_chr":str} )
            #print( pgs_id )
            pgs_scores = pgs_scores[  ["hm_chr", "hm_pos", "effect_weight"] ] # .drop(columns=["chr_name", "effect_weight",	"chr_position",	"effect_allele",	"other_allele",	"hm_source",	"hm_rsID",	"hm_inferOtherAllele"])
            prs_data_list = []
            print( pgs_id, trait, num )
            pgs_scores[ "hm_chr"] = pgs_scores[ "hm_chr"].astype("str")
            for index, row in pgs_scores.iterrows():
                if row["hm_chr"] in dhs_tree.keys():
                    if np.isnan( row[ "hm_pos"] ):
                        pass
                    else:
                        region = "DHS" if len(dhs_tree[row['hm_chr']].at(row['hm_pos'])) > 0 else 'Non-DHS'
                        d_id = "N"
                        if region == "DHS":
                            se = list( dhs_tree[row['hm_chr']].at(row['hm_pos']) )
                            d_id = dhs_id[ (row["hm_chr"], se[0][0], se[0][1] )]
                        prs_data_list.append({ 'hm_chr':row["hm_chr"], 'pos':row["hm_pos"], 'effect_weight': row['effect_weight'] , "dhs_id":d_id })
                        #break
            prs_file = "/mnt/d/sandbox/0.gxe/dhs/prs_dhsid/" + f'{pgs_id}_data.csv'
            with gzip.open(prs_file , 'w') as f:
                f.write( headlines )
            prs_data_list = pd.DataFrame(prs_data_list)
            prs_data_list[ "pos"] = prs_data_list[ "pos"].astype( "int" )
            prs_data_list.to_csv( prs_file , index=False, mode = "a")
            break




# check
ntrait = 0
dhsfiles = [ x.split("/")[-1].split("_")[0] for x in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/dhsScores/*.csv.gz') ]
prsfiles = [ x.split("/")[-1].split("_")[0] for x in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/*.txt.gz') ]

for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/*.txt.gz'):
    with gzip.open(filename, 'rt') as f:
        # Read the header lines
        #trait_name = f.readline().split('=')[1].strip()
        #num_variants = int(f.readline().split('=')[1].strip())
        num = 0
        for l in f:
            if "#variants_number=" in l:
                num = int( l.split("=")[1].strip() )
            elif "#trait_reported=" in l:
                trait = l.split( "=")[1].strip()
            elif "#pgs_id=" in l:
                pgs_id = l.split( "=")[1].strip()
            elif "HmPOS_match_pos" in l:
                break

        if (num > 10000) : # & ( filename in dhsfiles):
            if prs_id not in dhsfiles:
                print( pgs_id, trait , num )


# longevity associated PRS

filenames = [ "PGS000906_hmPOS_GRCh38.txt.gz",   # longevity
"PGS000318_hmPOS_GRCh38.txt.gz",   # all-cause mortality female
"PGS000319_hmPOS_GRCh38.txt.gz"]    # all-cause mortality male

pgs_id = "PGS000906"
prs_file = "/mnt/d/sandbox/0.gxe/pgsall/dhsScores/" + f'{pgs_id}_data.csv.gz'
with gzip.open(prs_file , 'r') as f:
    pgs_scores = pd.read_csv(f, sep='\t', comment='#', low_memory=False, dtype = {"hm_chr":str} )



ntrait = 0
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/*.txt.gz'):
    with gzip.open(  filename, 'rt') as f:
        # Read the header lines
        #trait_name = f.readline().split('=')[1].strip()
        #num_variants = int(f.readline().split('=')[1].strip())
        num = 0
        for l in f:
            if "#variants_number=" in l:
                num = int( l.split("=")[1].strip() )
            elif "#trait_reported=" in l:
                trait = l.split( "=")[1].strip()
            elif "#pgs_id=" in l:
                pgs_id = l.split( "=")[1].strip()
            elif "HmPOS_match_pos" in l:
                break

        if ( num <= 10000 )  & (num > 1000) :
            headlines = "#pgs_id=%s\n" % ( pgs_id )
            headlines += "#varinats_number=%d\n" % ( num )
            headlines += "#trait_reported=%s\n" % ( trait )
            ntrait += 1
            # Load the file into a dataframe
            with gzip.open(  filename, 'rt') as f:
                pgs_scores = pd.read_csv(f, sep='\t', comment='#', low_memory=False, dtype = {"hm_chr":str} )
            #print( pgs_id )
            pgs_scores = pgs_scores[  ["hm_chr", "hm_pos", "effect_weight"] ] # .drop(columns=["chr_name", "effect_weight",	"chr_position",	"effect_allele",	"other_allele",	"hm_source",	"hm_rsID",	"hm_inferOtherAllele"])
            prs_data_list = []
            print( pgs_id, trait )
            pgs_scores[ "hm_chr"] = pgs_scores[ "hm_chr"].astype("str")
            for index, row in pgs_scores.iterrows():
                if row["hm_chr"] in dhs_tree.keys():
                    if np.isnan( row[ "hm_pos"] ):
                        pass
                    else:
                        region = 'DHS' if len(dhs_tree[row['hm_chr']].at(row['hm_pos'])) > 0 else 'Non-DHS'
                        prs_data_list.append({ 'hm_chr':row["hm_chr"], 'pos':row["hm_pos"], 'region': region, 'effect_weight': row['effect_weight']  })
                        #break
            prs_file = "/mnt/d/sandbox/0.gxe/pgsall/dhsScores/" + f'{pgs_id}_data.csv'
            with open(prs_file , 'w') as f2:
                f2.write( headlines )
            prs_data_list = pd.DataFrame(prs_data_list)
            prs_data_list[ "pos"] = prs_data_list[ "pos"].astype( "int" )
            prs_data_list.to_csv( prs_file , index=False, mode = "a")
            break



# obtain DHS ID

ntrait = 0
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/dhsScores/*.txt.gz'):
    with gzip.open(filename, 'rt') as f:
        # Read the header lines
        #trait_name = f.readline().split('=')[1].strip()
        #num_variants = int(f.readline().split('=')[1].strip())
        num = 0
        headlines = ""
        for l in f:
            if "#" == l[0]:
                headlines += l
                if headlines.count( "#") >= 3:
                    break
        if 1 : #  & (num < 100000) :
            ntrait += 1
            # Load the file into a dataframe
            pgs_scores = pd.read_csv(filename, sep=',', comment='#', low_memory=False, dtype = {"hm_chr":str} )
            #print( pgs_id )
            prs_data_list = []
            print( pgs_id, trait, num )
            pgs_scores_nondhs = pgs_scores[ pgs_scores[ "region"] == "Non-DHS"]
            pgs_scores_dhs = pgs_scores[ pgs_scores[ "region"] == "DHS"]
            for index, row in pgs_scores_dhs.iterrows():
                #if row["hm_chr"] in dhs_tree.keys():
                    #region = "DHS" if len(dhs_tree[row['hm_chr']].at(row['hm_pos'])) > 0 else 'Non-DHS'
                    #d_id = "N"
                    #if region == "DHS":
                    se = list( dhs_tree[row['hm_chr']].at(row['hm_pos']) )
                    d_id = dhs_id[ (row["hm_chr"], se[0][0], se[0][1] )]
                    prs_data_list.append( d_id )
            break
            prs_file = "/mnt/d/sandbox/0.gxe/dhs/prs_dhsid/" + f'{pgs_id}_data.csv'
            with gzip.open(prs_file , 'w') as f:
                f.write( headlines )
            prs_data_list = pd.DataFrame(prs_data_list)
            prs_data_list[ "pos"] = prs_data_list[ "pos"].astype( "int" )
            prs_data_list.to_csv( prs_file , index=False, mode = "a")
            break



#
# compare effect weights between DHS and non-DHS
#
# whole distribution
#

results = []
filenames = [ "PGS000906_data.csv.gz",   # longevity
"PGS000318_data.csv.gz",   # all-cause mortality female
"PGS000319_data.csv.gz",
"PGS000014_data.csv.gz"]    # all-cause mortality male


# ratio of number of variants in DHS and Non-DHS
numvar = []
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/dhsScores/*.csv.gz') :
#for filename in filenames :
    #filename = '/mnt/d/sandbox/0.gxe/pgsall/dhsScores/' + filename
    with gzip.open(filename, 'rt') as f:
        trait = ""
        #for l in f:
        #    if "#variants_number=" in l:
        #        num = int( l.split("=")[1].strip() )
        #    elif "#pgs_id=" in l:
        #        pgs_id = l.split( "=")[1].strip()
        #    elif "#trait_reported=" in l:
        #        trait = l.split( "=")[1].strip()
        #        break
        prs_id = os.path.basename(filename ).split('_data')[0]
        prs_scores = pd.read_csv(f, sep=',', comment='#', low_memory=False, dtype = {"hm_chr":str} )
        n0 = sum(  prs_scores[ "region"] == "Non-DHS" )
        n1 = sum(  prs_scores[ "region"] == "DHS" )
        numvar.append( {"PRS_ID": prs_id, "nratio":float(n0)/n1  })

numvar_df = df.DataFrame( numvar )

filename = "/mnt/d/sandbox/0.gxe/pgsall/dhsScores/PGS001356_data.csv.gz"
with gzip.open(filename, 'rt') as f:
    prs_scores = pd.read_csv(f, sep=',', comment='#', low_memory=False, dtype = {"hm_chr":str} )



#
# faction of DHS regions
#
dhs_data = pd.read_csv('/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz', sep='\t',usecols=[0,1,2] )
dhs_data.columns = ['chrom', 'start', 'end']
dhs_data['chrom'] = dhs_data['chrom'].apply(str).str.split(pat = "r" , expand = True)[ [1 ] ]



# gnomAD
# testing colcalizaiton
#
dhs_treechr = {}
for chr in dhs_tree:
    dhs_treechr[ "chr" + chr ] = dhs_tree[ chr ]


filen = "/mnt/d/sandbox/0.gxe/data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"
numvar = []

import pandas as pd

# Define the column names. You might need to adjust this depending on the exact structure of your data
column_names = ["#CHROM", "POS"]

# Read the gnomAD file, skipping comment lines
gnomad_df = pd.read_csv(filen, comment='#', sep="\t", usecols=[0,1] )


##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count for samples">
##INFO=<ID=AC_afr,Number=A,Type=Integer,Description="Alternate allele count for samples of African-American/African ancestry">
##INFO=<ID=AC_afr_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of African-American/African ancestry">
##INFO=<ID=AC_afr_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of African-American/African ancestry">
##INFO=<ID=AC_amr,Number=A,Type=Integer,Description="Alternate allele count for samples of Latino ancestry">
##INFO=<ID=AC_amr_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of Latino ancestry">
##INFO=<ID=AC_amr_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of Latino ancestry">
##INFO=<ID=AC_asj,Number=A,Type=Integer,Description="Alternate allele count for samples of Ashkenazi Jewish ancestry">
##INFO=<ID=AC_asj_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of Ashkenazi Jewish ancestry">
##INFO=<ID=AC_asj_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of Ashkenazi Jewish ancestry">
##INFO=<ID=AC_eas,Number=A,Type=Integer,Description="Alternate allele count for samples of East Asian ancestry">
##INFO=<ID=AC_eas_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of East Asian ancestry">
##INFO=<ID=AC_eas_jpn,Number=A,Type=Integer,Description="Alternate allele count for samples of Japanese ancestry">
##INFO=<ID=AC_eas_kor,Number=A,Type=Integer,Description="Alternate allele count for samples of Korean ancestry">
##INFO=<ID=AC_eas_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of East Asian ancestry">
##INFO=<ID=AC_eas_oea,Number=A,Type=Integer,Description="Alternate allele count for samples of Other East Asian ancestry">
##INFO=<ID=AC_female,Number=A,Type=Integer,Description="Alternate allele count for female samples">
##INFO=<ID=AC_fin,Number=A,Type=Integer,Description="Alternate allele count for samples of Finnish ancestry">
##INFO=<ID=AC_fin_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of Finnish ancestry">
##INFO=<ID=AC_fin_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of Finnish ancestry">
##INFO=<ID=AC_male,Number=A,Type=Integer,Description="Alternate allele count for male samples">
##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="Alternate allele count for samples of Non-Finnish European ancestry">
##INFO=<ID=AC_nfe_bgr,Number=A,Type=Integer,Description="Alternate allele count for samples of Bulgarian (Eastern European) ancestry">
##INFO=<ID=AC_nfe_est,Number=A,Type=Integer,Description="Alternate allele count for samples of Estonian ancestry">
##INFO=<ID=AC_nfe_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of Non-Finnish European ancestry">
##INFO=<ID=AC_nfe_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of Non-Finnish European ancestry">
##INFO=<ID=AC_nfe_nwe,Number=A,Type=Integer,Description="Alternate allele count for samples of North-Western European ancestry">
##INFO=<ID=AC_nfe_onf,Number=A,Type=Integer,Description="Alternate allele count for samples of Other Non-Finnish European ancestry">
##INFO=<ID=AC_nfe_seu,Number=A,Type=Integer,Description="Alternate allele count for samples of Southern European ancestry">
##INFO=<ID=AC_nfe_swe,Number=A,Type=Integer,Description="Alternate allele count for samples of Swedish ancestry">
##INFO=<ID=AC_oth,Number=A,Type=Integer,Description="Alternate allele count for samples of Other ancestry">
##INFO=<ID=AC_oth_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of Other ancestry">
##INFO=<ID=AC_oth_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of Other ancestry">
##INFO=<ID=AC_popmax,Number=A,Type=Integer,Description="Allele count in the population with the maximum AF">
##INFO=<ID=AC_raw,Number=A,Type=Integer,Description="Alternate allele count for samples, before removing low-confidence genotypes">
##INFO=<ID=AC_sas,Number=A,Type=Integer,Description="Alternate allele count for samples of South Asian ancestry">
##INFO=<ID=AC_sas_female,Number=A,Type=Integer,Description="Alternate allele count for female samples of South Asian ancestry">
##INFO=<ID=AC_sas_male,Number=A,Type=Integer,Description="Alternate allele count for male samples of South Asian ancestry">

##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in samples">
##INFO=<ID=vep,Number=.,Type=String,
variant_type=snv;vep=C
Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|
BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|
VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|
AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|
CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info">

bcftools view -f .,PASS gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz > filtered.vcf
bcftools view -i 'INFO/AN>113173' gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz > filtered_AN.vcf
bcftools query -f '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/vep\n]' gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz > filtered_AN_AC_vep.tsv


with gzip.open(filen, 'rt') as f:
    trait = ""
    numvar = [ 0, 0]
    for l in f:
        if l[0] != "#":
            chrom, pos = l.split()[0:2]
            if chrom in dhs_treechr.keys():
                region = 1 if len(dhs_treechr[chrom ].at( int( pos) ) ) > 0 else 0
                numvar[ region ] += 1





import os
import gzip
import re

def parse_info(info_str, ids_to_extract):
    info_dict = {}
    info_parts = info_str.split(';')
    for part in info_parts:
        key, _, value = part.partition('=')
        if key in ids_to_extract:
            if key == 'vep':
                value = process_vep(value)
            info_dict[key] = value
    return info_dict

Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position| 13
Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE| 25
HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF| 40
AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF| 56
CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info

def process_vep(vep_str):
    vep_entries = vep_str.split(',')
    processed_entries = []
    for entry in vep_entries:
        vep_values = entry.split('|')
        processed_entry = f"{vep_values[1]}|{vep_values[4]}|{vep_values[55]}|{vep_values[64]}|{vep_values[66]}"
        processed_entries.append(processed_entry)
        #if len( vep_values[55]) > 0:
        #    print( vep_values[55])
    return ','.join(processed_entries)

def main():
    input_file = '/mnt/d/sandbox/0.gxe/data/filtered_AN.vcf'  # 2.1.1 liftover_grch38
    output_file = '/mnt/d/sandbox/0.gxe/gnomad/filtered_AN_AC.tsv.gz'
    ids_to_extract = ['AC_afr', 'AC_amr', 'AC_asj', 'AC_eas','AC_nfe', 'AC_sas',  'AC_eas_jpn', 'AC_eas_kor', 'AC_eas_oea', 'AC_fin', 'AC_nfe_bgr', 'AC_nfe_est', 'AC_nfe_nwe', 'AC_nfe_onf', 'AC_nfe_seu', 'AC_nfe_swe', 'vep']

    descriptions = {}

    with open(input_file, 'rt') as file, gzip.open(output_file, 'wt') as output:
        for line in file:
            if line.startswith('##'):
                match = re.search(r'##INFO=<ID=(\w+),', line)
                if match and match.group(1) in ids_to_extract:
                    desc_match = re.search(r'Description="([^"]+)"', line)
                    if desc_match:
                        descriptions[match.group(1)] = desc_match.group(1)
            elif line.startswith('#CHROM'):
                break

        for id in ids_to_extract:
            print(f"{id}: {descriptions.get(id, 'No description available')}")
        print()

        output.write('\t'.join([ "chrom",  "pos",   "snp_id" ] + ids_to_extract ) + '\n')

        for line in file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                info_str = parts[7]
                info_dict = parse_info(info_str, ids_to_extract)
                row =  [parts[0], parts[1], parts[2] ]  + [info_dict.get(id, '.') for id in ids_to_extract]
                output.write('\t'.join(row) + '\n')
                break

if __name__ == '__main__':
    main()




#



import pandas as pd
from intervaltree import Interval, IntervalTree

# Load the DHS regions
dhs_df = pd.read_csv("/mnt/d/sandbox/0.gxe/dhs/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz", sep="\t")

# Create empty dictionaries to store interval trees for each chromosome
trees = {"chr"+str(i): IntervalTree() for i in range(1, 23)}  # for autosomes only

# Build interval trees
for row in dhs_df.itertuples():
    chrom = str(row.seqname)
    if chrom in trees:
        # Add the interval and the identifier to the interval tree of the corresponding chromosome
        trees[chrom].addi(row.start, row.end, str( row.identifier) )

# Load the locus data
locus_df = pd.read_csv("/mnt/d/sandbox/0.gxe/gnomad/filtered_AN_AC.tsv.gz", sep="\t")

# Add a new column for the DHS identifier
locus_df['DHS_id'] = 'N'

# Check each locus if it's located within DHS
for row in locus_df.itertuples():
    chrom = str(row.chrom)
    if chrom in trees:
        # Query the interval tree for overlapping DHS regions
        overlapping_dhs = trees[chrom].at(row.pos)
        if overlapping_dhs:
            # If there are overlapping DHS regions, add their identifiers to the 'DHS_id' column
            locus_df.at[row.Index, 'DHS_id'] = ','.join(interval.data for interval in overlapping_dhs)

# Save the updated locus data
locus_df.to_csv("/mnt/d/sandbox/0.gxe/gnomad/filtered_AN_AC_DHS.tsv.gz", sep="\t", index=False)

chr1    729395  rs1214697418    0       1       0       0       0       0       0       0       0       0  downstream_gene_variant|ENSG00000229376|||,upstream_gene_var
iant|ENSG00000230021|||,upstream_gene_variant|ENSG00000224956|||,regulatory_region_variant||||,regulatory_re
gion_variant||||        1.10263,1.10262










f = "/mnt/d/sandbox/0.gxe/data/human_lr_pair.txt"
lrpair = pd.read_csv( f )

        ndhscut = ( prs_scores[ "region"] == "Non-DHS" ) & ( prs_scores["effect_weight"] > 0 )
        dhscut =  ( prs_scores[ "region"] == "DHS" ) & ( prs_scores["effect_weight"] > 0  )
        p = scipy.stats.mannwhitneyu( prs_scores[ dhscut ][ "effect_weight"] ,
              prs_scores[ ndhscut ][ "effect_weight"] )
        results.append( {"p-val":p[1], "mdiff":  np.mean( prs_scores[ dhscut ][ "effect_weight"] ) -
             np.mean( prs_scores[ ndhscut ][ "effect_weight"] ), "prsid":prs_id + "_P",
             "trait": trait } )
        print( prs_id, trait, sum( ndhscut), sum( dhscut ) )

        ndhscut = ( prs_scores[ "region"] == "Non-DHS" ) & ( prs_scores["effect_weight"] < 0 )
        dhscut  = ( prs_scores[ "region"] == "DHS" ) & ( prs_scores["effect_weight"] < 0  )
        p = scipy.stats.mannwhitneyu( prs_scores[ dhscut ][ "effect_weight"] ,
             prs_scores[ ndhscut ][ "effect_weight"] )
        results.append( {"p-val":p[1], "mdiff":  np.mean( prs_scores[ dhscut ][ "effect_weight"] ) -
             np.mean( prs_scores[ ndhscut ][ "effect_weight"] ), "prsid":prs_id + "_N",
             "trait": trait } )
        print( prs_id, trait, sum( ndhscut), sum( dhscut ) )


results_df = pd.DataFrame( results)


diff_h = results_df[ ( ( results_df["mdiff"] > 0 ) & ( results_df[ "p-val"] < 0.05/1697 ) ) ]
diff_l = results_df[ ( ( results_df["mdiff"] < 0 ) & ( results_df[ "p-val"] < 0.05/1697 ) ) ]




results = []
for index, row in results_df.iterrows():
    p_val, d, prsid, trait  = row[ "p-val"], row[ "mdiff"], row[ "prsid"] , row[ "trait" ]
    if ( p_val <=  0.05/1697 ) :
        f = '/mnt/d/sandbox/0.gxe/pgsall/dhsScores/*_data.csv.gz' % ( prsid )
        with gzip.open(filename, 'rt') as f:
            trait = ""
            num = 0
            for l in f:
                if "#variants_number=" in l:
                    num = int( l.split("=")[1].strip() )
                elif "#pgs_id=" in l:
                    pgs_id = l.split( "=")[1].strip()
                elif "#trait_reported=" in l:
                    trait = l.split( "=")[1].strip()
                    break
        break



#
# compare effect weights between DHS and non-DHS
#
# test the whole distribution for weight > 0
# test half and half if median is about 0
#



results = []
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/dhsScores/*.csv.gz'):
    with gzip.open(filename, 'rt') as f:
        prs_id = os.path.basename(filename).split('_data')[0]
        pgs_id = ""
        trait = ""
        for l in f:
            if "#variants_number=" in l:
                num = int( l.split("=")[1].strip() )
            elif "#pgs_id=" in l:
                pgs_id = l.split( "=")[1].strip()
                if trait != "":
                    break
            elif "#trait_reported=" in l:
                trait = l.split( "=")[1].strip()
                if pgs_id != "":
                    break

        prs_scores = pd.read_csv(f, sep=',', comment='#', low_memory=False, dtype = {"hm_chr":str} )
        if sum( prs_scores[ "effect_weight"] < 0 ) > 1000:
            #
            # test for weight < 0
            threD =  ( prs_scores["region"] == "DHS" ) & ( prs_scores[ "effect_weight"] < 0 )
            threND = ( prs_scores[ "region"] == "Non-DHS") & ( prs_scores[ "effect_weight"] < 0 )
            p = scipy.stats.mannwhitneyu( prs_scores[ threD ] [ "effect_weight"] ,
                         prs_scores[ threND ][ "effect_weight" ] )
            results.append( {"p-val":p[1], "mdiff":  np.mean( prs_scores[ threD ][ "effect_weight"] ) -
                         np.mean( prs_scores[ threND][ "effect_weight"] ), "prsid":prs_id,
                         "trait": trait, "nVar" : sum( prs_scores[ "effect_weight"]<0), "dist": -1  } )

            # test for weight > 0
            threD =  ( prs_scores["region"] == "DHS" ) & ( prs_scores[ "effect_weight"] > 0 )
            threND = ( prs_scores[ "region"] == "Non-DHS") & ( prs_scores[ "effect_weight"] > 0 )
            p = scipy.stats.mannwhitneyu( prs_scores[  threD ][ "effect_weight"] ,
                         prs_scores[ threND ][ "effect_weight"] )
            results.append( {"p-val":p[1], "mdiff":  np.mean( prs_scores[ threD ][ "effect_weight"] ) -
                         np.mean( prs_scores[ threND][ "effect_weight"] ), "prsid":prs_id,
                         "trait": trait, "nVar" : sum( prs_scores[ "effect_weight"] >0 ), "dist": 1  } )
            # break


results_df = pd.DataFrame( results )



results_dfcut =  results_df[ results_df[ "p-val"] < 0.05/3366 ]

log10pval = -np.log10( results_dfcut[ "p-val" ] )
logdiff = -np.log10( abs( results_dfcut[ "mdiff"] )  )

sns.scatterplot( x = log10pval, y = logdiff )



#
f = '/mnt/d/sandbox/0.gxe/pgsall/dhsScores/PGS000039_data.csv.gz'
prs_scores = pd.read_csv(f, sep=',', comment='#', low_memory=False, dtype = {"hm_chr":str} )
sns.histplot(   )
sns.histplot(   )


sns.histplot( prs_scores[ prs_scores[ "region"] == "Non-DHS"][ "effect_weight"]  )
sns.histplot( prs_scores[ prs_scores[ "region"] == "DHS"][ "effect_weight" ]  )



from scipy import stats

df1 = prs_scores[ ( prs_scores[ "region"] == "Non-DHS") & (prs_scores[ "effect_weight"] < 0 )]
df2 = prs_scores[ ( prs_scores[ "region"] == "DHS") & (prs_scores[ "effect_weight"] < 0 ) ]

# The second returned argument is the lambda parameter that was used
df1['boxcox_effect_weight'], _ = stats.boxcox( abs( df1['effect_weight']) )
df2['boxcox_effect_weight'], _ = stats.boxcox( abs( df2['effect_weight']) )

sns.histplot(df1['boxcox_effect_weight'], bins=30, kde=True, stat='density', label='Non-DHS', color='blue')
sns.histplot(df2['boxcox_effect_weight'], bins=30, kde=True, stat='density', label='DHS', color='red')

plt.legend()
plt.show()



#
sns.ecdfplot(df1['effect_weight'], label='Region 1', color='blue')
sns.ecdfplot(df2['effect_weight'], label='Region 2', color='red')

plt.legend()
plt.show()

#
# DHS associed with complex traits
# DHS density in genes for 4 categoreis
# PCG, lncRNA, ncRNA, pseudogene (control sets)
#






mean_ci_data = pd.DataFrame(columns=['prs_id', "q" ,'region', 'effect_weight', 'ci_low', 'ci_high'])
pvalue_data = pd.DataFrame(columns=['prs_id', "q", "pval", "num_nondhs", "num_dhs"] )


        qrange = np.arange(0.95, 1.01, 0.01)
        pgs_scores_weights = pgs_scores[  "effect_weight" ]
        for q in range(5):
            q_list = ( pgs_scores_weights.quantile( qrange[ q] ) < pgs_scores_weights ) & ( pgs_scores_weights < pgs_scores_weights.quantile( qrange[ q + 1 ] ) )
            pgs_scores_q =  pgs_scores[ q_list  ]
            pgs_region = {}
            for region in ['DHS', 'Non-DHS']:
                #mean_ci_data = mean_ci_data.append({'tissue': tissue, 'region': region, 'mean_abs_slope': mean_abs_slope, 'ci_low': ci_low, 'ci_high': ci_high}, ignore_index=True)
                region_data = pgs_scores_q[pgs_scores_q['region'] == region]['effect_weight']  # .abs()
                #region_data80 = region_data[ region_data > region_data.quantile( 0.2 ) ]
                pgs_region[ region ] = region_data
                mean_effect_weight = np.mean( region_data )
                ci_low, ci_high = stats.norm.interval(0.95, loc=mean_effect_weight, scale=stats.sem( region_data))
                mean_ci_row = pd.DataFrame({'prs_id': [prs_id],
                                    'q' : [ qrange[ q] ] ,
                                    'region': [region],
                                    'effect_wei  ght': [mean_effect_weight],
                                    'ci_low': [ci_low],
                                    'ci_high': [ci_high]})
                mean_ci_data = pd.concat([mean_ci_data, mean_ci_row], ignore_index=True)
            pval = scipy.stats.ttest_ind( pgs_region[ "Non-DHS" ]  ,  pgs_region[ "DHS" ] )
            pval_row = pd.DataFrame( {'prs_id': [prs_id],
                                'q': [qrange[ q] ],
                                'pval': [pval],
                                'num_nondhs' : [ len( pgs_region[ "Non-DHS"] ) ] ,
                                 'num_dhs' : [ len( pgs_region[ "DHS"] ) ] }  )
            pvalue_data = pd.concat( [ pvalue_data, pval_row], ignore_index =True)
            #pgs_scores_quantile = pgs_scores[ "effect_weight"].quantile( q = c( 0.05, 0.1, 0.15, 0.2) )



#
#
# visualize p-values for each quantiles in each PRS trait
#
##
import seaborn as sns

# Stacking the data
#dhs_data['region'] = 'DHS'
#non_dhs_data['region'] = 'non-DHS'

pvals = pvalue_data["pval"]
pvalue_data[ "-log10(pval)"] = [ -np.log10( p[1] ) for p in pvals ]
r = mean_ci_data[ mean_ci_data["region"]== "DHS"]["effect_weight"].reset_index(drop=True)/mean_ci_data[ mean_ci_data["region"]== "Non-DHS"]["effect_weight"].reset_index(drop=True)
pvalue_data[ "ratio"]= np.log2( r )



fig, ax = plt.subplots(figsize=(20,10))
ax = sns.barplot(data=pvalue_data, x="q", y='pval', errorbar=None )

#x_coords = [p.get_x() + 0.5*p.get_width() for p in ax.patches]
#ax.errorbar(data= combined_data , x=x_coords, y='log10pval', yerr='ci', ls='', lw=8, color='black')

plt.title('-log10( P-value) for each quantile of PRS distribution')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()


pvals = pvalue_data["pval"]
pvalue_data_nonzero = pvalue_data[ x != 0 ]
x = abs( pvalue_data[ "num_dhs"] - pvalue_data[ "num_nondhs" ]  )
# pvalue_data[ "-log(abs(difference))"] =   np.log10( x ] )


fig, ax = plt.subplots(figsize = (9, 6))
# Add scatterplot
sns.jointplot( x = "-log10(pval)", y = "ratio", hue = "q", data = pvalue_data )



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
        break
    break
    # Convert the list to a DataFrame and write it to a csv file
    tissue_data_df = pd.DataFrame(tissue_data_list)
    tissue_data_df.to_csv("/mnt/d/sandbox/0.gxe/dhs/gtex_slope/" + f'{tissue}_data.csv', index=False)
