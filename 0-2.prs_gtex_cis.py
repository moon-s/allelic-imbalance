# PGS

# Iterate over the files in the directory
import os
import glob
import pandas as pd
import scipy.stats as stats

# Set the directory path
# gtex_directory = "/mnt/d/sandbox/0.gxe/gtex/GTEx_Analysis_v8_eQTL/"
# Get a list of all files in the directory
# files = glob.glob(os.path.join(directory_path, "*.v8.egenes.txt.gz"))



# SNP in genic/up&down-stream/UTR/exonic/

# check genic region :
# -> average transcript length

average_transcript_lengths = pd.read_csv( "/mnt/d/sandbox/0.gxe/gtex/average_transcript_lengths_GTEx_common.csv" )
genic = pd.read_csv( "/mnt/d/sandbox/0.gxe/data/mart_human_genic.txt" , sep = "\t" , low_memory=False)

atl_merged = pd.merge(average_transcript_lengths, genic,  on='Gene stable ID')
# atl - genic : 156289

atl_merged = atl_merged.drop(columns=["Transcript stable ID", "Transcript length (including UTRs and CDS)"])
atl_merged = atl_merged.drop_duplicates()

atl_merged.to_csv("average_transcript_lengths_common_genic.csv", index=False)

# Iterate over the files
pgc_dir = "/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/"
files = glob.glob(os.path.join(pgc_dir,"*.txt.gz") )


# check number SNP in each PRS set :
# more than 10000 for 5%
# 1000 genes assigned
#
# top and bottom 5%, 10%, 15%, 20%, 25% from distribution of effect weights
# assignment: postive weight, negative weight, mixed weights
#

import pandas as pd
import glob
import gzip


# Iterate over all files in the directory
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/dhsScores/*.csv.gz'):
    comments = []
    data = []

    with gzip.open(filename, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                comments.append(line)
            else:
                data.append(line)

    # Filter the data lines
    dhs_data = [line for line in data if ',DHS,' in line]
    # Assign the effect_weight to the matching genes in the alt_merged dataframe
    nondhs_data = [line for line in data if ',Non-DHS,' in line ]
    head = data[0]
    if head[0] == ",":
        data[0] = head[1:]
        dhs_data = [ ",".join( line.split(",")[1:] ) for line in dhs_data]
        nondhs_data = [",".join( line.split(",")[1:] ) for line in nondhs_data]
    dhs_data_id = []
    print ( comments )
    for row in dhs_data:
            row1 = row.split(",")
            chrom, pos = row1[0], int( row1[1])
            se = list( dhs_tree[ chrom ].at( pos ) )
            d_id = dhs_id[ ( chrom , se[0][0], se[0][1] )]
            row1[ 2] = str( d_id )
            dhs_data_id.append( ",".join( row1 ) )

            # Create a subset dataframe of the PGS score table that includes the updated gene IDs and transcript lengths
            # Save the subset dataframe
            #subset.to_csv(f"{trait_name}_subset.csv", index=False)
    prs_file = "/mnt/d/sandbox/0.gxe/dhs/prs_dhsid/" + filename.split("/")[-1][:-3]
    with open( prs_file , 'w') as f:
        for comment in comments:
            f.write(comment)
        f.write( data[0])
        for line in dhs_data_id:
            f.write(line)
        for line in nondhs_data:
            f.write(line)



# DHS density given gene sets
# up-regulated genes in obese
# down-regulated genes in obese
#

# # slopes of eVariants in the DHS associaed with up-regulated genes in obese
# expected tissue-specificity of up-regulated genes in obese



import pandas as pd
import glob
import gzip

# Load the alt_merged dataframe
# alt_merged = pd.read_csv('alt_merged.csv')
alt_merged = pd.read_csv( "/mnt/d/sandbox/0.gxe/gtex/average_transcript_lengths_common_genic.csv" )

#pgs_id=PGS000001
#trait_reported=Breast Cancer
#variants_number=77
# Create a new empty DataFrame to store the rows of pgs_scores that fall within a genic region

# Iterate over all files in the directory
pgs_scores_shared = pd.DataFrame()


ntrait = 0
for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/*.txt.gz'):
    with gzip.open(filename, 'rt') as f:
        # Read the header lines
        #trait_name = f.readline().split('=')[1].strip()
        #num_variants = int(f.readline().split('=')[1].strip())
        for l in f:
            if "#variants_number=" in l:
                num = int( l.split("=")[1].strip() )
            elif "#trait_reported=" in l:
                trait = l.split( "=")[1].strip()
            elif "#pgs_id=" in l:
                pgs_id = l.split( "=")[1].strip()
            elif "HmPOS_match_pos" in l:
                break
        if num > 100000:
            ntrait += 1
            # Load the file into a dataframe
            pgs_scores = pd.read_csv(f, sep='\t', comment='#', low_memory=False )
            #print( pgs_id )
            pgs_scores = pgs_scores[  ["hm_chr", "hm_pos"] ] # .drop(columns=["chr_name", "effect_weight",	"chr_position",	"effect_allele",	"other_allele",	"hm_source",	"hm_rsID",	"hm_inferOtherAllele"])
            if pgs_scores_shared.empty:
                pgs_scores_shared = pgs_scores
            else:
                pgs_scores_shared = pd.concat( [pgs_scores_shared,  pgs_scores] )
                pgs_scores_shared = pgs_scores_shared.drop_duplicates()
                print( pgs_id, len( pgs_scores_shared) )


pgs_scores_shared.to_csv( "/mnt/d/sandbox/0.gxe/pgsall/scoringFiles_loci.csv")

# iterate from pgs_scores_shared
genic_rows = []  # list to store rows that fall within a genic region

for chrom in alt_merged[ "hm_chr"]:
    pgs_scores_chr = pgs_scores[ pgs_scores[ "hm_chr"] == chrom ]
    alt_merged_chr = alt_merged[ alt_merged['Chromosome/scaffold name'] == chrom ]
    for idx, row in pgs_scores_chr.iterrows():
        pos = row['hm_pos']

        # Get the corresponding genic region in alt_merged
        tmp =  ( alt_merged_chr['Gene start (bp)'] <= pos ) &
                              ( alt_merged_chr['Gene end (bp)'] >= pos )
        # If the genic region is not empty, add the row from pgs_scores to the list
        if not tmp.empty:
            genic_rows.append(row)



# Concatenate all the rows into a new DataFrame
pgs_scores_genic = pd.concat(genic_rows, axis=1).transpose()

# Reset the index of pgs_scores_genic
pgs_scores_genic.reset_index(drop=True, inplace=True)


            # If the genic region is not empty, add the row from pgs_scores to pgs_scores_genic
            if not genic_region.empty:
                a = row.to_frame( ).T
                a.reset_index(drop=True, inplace=True)
                        genic_region.reset_index(drop=True, inplace=True)
                        ag = a.join( genic_region)
                        #pgs_scores_genic = pgs_scores_genic.append( ag )
                        if pgs_scores_genic.empty:
                            pgs_scores_genic = ag
                        else:
                            pgs_scores_genic.loc[ len( pgs_scores_genic) ] = ag.loc[0]

            # Reset the index of pgs_scores_genic
            pgs_scores_genic.reset_index(drop=True, inplace=True)





import pandas as pd
from intervaltree import Interval, IntervalTree

# Load the data
pgs_scores_shared = pd.read_csv('/path/to/pgs_scores_shared.csv')
alt_merged = pd.read_csv('/path/to/alt_merged.csv')

# Create a dictionary of interval trees

pgs_chrom = pgs_scores_shared[ "hm_chr"].unique()

trees = dict()
genelen = dict()
for chrom in alt_merged['Chromosome/scaffold name'].unique():
    # Subset the data for the current chromosome

    chrom_data = alt_merged[alt_merged['Chromosome/scaffold name'] == chrom]

    # Create an interval tree for the current chromosome
    tree = IntervalTree()

    for index, row in chrom_data.iterrows():
        # Add each gene as an interval to the tree
        tree.add(Interval(row['Gene start (bp)'], row['Gene end (bp)'], row[ "Gene stable ID" ]))
        genelen[ row["Gene stable ID"] ] = row[ "Average transcript length" ]

    # Add the tree to our dictionary of trees
    trees[chrom] = tree

# Create a new dataframe to store the genes that overlap with our positions
genes_overlapping = pd.DataFrame()


# Create a new dataframe to store the genes that overlap with our positions
genes_overlapping = pd.DataFrame()
pgs_scores_shared = pgs_scores_shared.dropna(subset=['hm_pos'])
pgs_scores_shared['hm_pos'] = pgs_scores_shared['hm_pos'].astype(int)



# genes_overlapping = pd.DataFrame(columns=['chr', 'pos' ])
overlapping_records = []
chroms =  alt_merged['Chromosome/scaffold name'].unique()
for index, row in pgs_scores_shared.iterrows():
    # Look up the overlapping genes for each position in our position data
    if row['hm_chr'] in chroms:
        overlapping = trees[row['hm_chr']].overlap(row['hm_pos'], row['hm_pos']+1 )

        # Add the overlapping genes to our new dataframe
        for interval in overlapping:
            overlapping_records.append( {'chr': row['hm_chr'], 'pos': row['hm_pos'] } )



# Write the new dataframe to a CSV file
genes_overlapping = pd.DataFrame(overlapping_records)
genes_overlapping.to_csv('prs_loci_genic.csv', index=False)



# add transcript length

trees = dict()
genelen = dict()
for chrom in alt_merged['Chromosome/scaffold name'].unique():
    # Subset the data for the current chromosome

    chrom_data = alt_merged[alt_merged['Chromosome/scaffold name'] == chrom]

    # Create an interval tree for the current chromosome
    tree = IntervalTree()

    for index, row in chrom_data.iterrows():
        # Add each gene as an interval to the tree
        tree.add( Interval(row['Gene start (bp)'], row['Gene end (bp)'], row[ "Gene stable ID" ] ) )
        genelen[ row["Gene stable ID"] ] = row[ "Average transcript length" ]

    # Add the tree to our dictionary of trees
    trees[chrom] = tree




# genes_overlapping = pd.DataFrame(columns=['chr', 'pos' ])
overlapping_records_len = []
chroms =  alt_merged['Chromosome/scaffold name'].unique()
for index, row in genes_overlapping.iterrows():
    # Look up the overlapping genes for each position in our position data
    if row['chr'] in chroms:
        overlapping = trees[row['chr']].overlap(row['pos'], row['pos']+1 )

        # Add the overlapping genes to our new dataframe
        if len( overlapping) == 1:
            for interval in overlapping:
                geneid = interval[2]
                overlapping_records_len.append( {'chr': row['chr'], 'pos': row['pos'], 'Gene ID':geneid, "Average transcript length":genelen[ geneid ] } )
        else:
            mlen = max([ genelen[x] for x[-1] in overlapping ] )
            for interval in overlapping:
                geneid = interval[2]
                if genelen[ geneid ] == mlen:
                    overlapping_records_len.append( {'chr': row['chr'], 'pos': row['pos'], 'Gene ID':geneid, "Average transcript length":genelen[ geneid ] } )


# Write the new dataframe to a CSV file
genes_overlapping_len = pd.DataFrame(overlapping_records_len)
genes_overlapping_len.to_csv('/mnt/d/sandbox/0.gxe/pgsall/prs_loci_genic_len.csv', index=False)




# transcript length distribuiton for quantiles of PRS variants with 5% intervals
genes_overlapping_len = pd.read_csv( "/mnt/d/sandbox/0.gxe/pgsall/prs_loci_genic_len.csv" )

for filename in glob.glob('/mnt/d/sandbox/0.gxe/pgsall/ScoringFiles/*.txt.gz'):
    with gzip.open(filename, 'rt') as f:
        # Read the header lines
        #trait_name = f.readline().split('=')[1].strip()
        #num_variants = int(f.readline().split('=')[1].strip())
        for l in f:
            if "#variants_number=" in l:
                num = int( l.split("=")[1].strip() )
            elif "#trait_reported=" in l:
                trait = l.split( "=")[1].strip()
            elif "#pgs_id=" in l:
                pgs_id = l.split( "=")[1].strip()
            elif "HmPOS_match_pos" in l:
                break
        if num > 100000:
            # Load the file into a dataframe
            pgs_scores = pd.read_csv(f, sep='\t', comment='#', low_memory=False )
            #print( pgs_id )
            pgs_scores = pgs_scores[  ["hm_chr", "hm_pos"] ] # .drop(columns=["chr_name", "effect_weight",	"chr_position",	"effect_allele",	"other_allele",	"hm_source",	"hm_rsID",	"hm_inferOtherAllele"])
            break
