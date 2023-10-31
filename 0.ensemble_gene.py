# Gene from ENSEMBL
# filtering
# genes on autosomes
# multiple transcripts for a gene
# Longest, Shortest, Average
# MIR genes from noncording

# Import necessary libraries
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Load the data
gene_data = pd.read_csv('/mnt/d/sandbox/0.gxe/data/mart_human_prot.txt', sep='\t', low_memory=False)

# Filter out genes that are not on autosomes
# includes autosomes only
gene_data = gene_data[gene_data['Chromosome/scaffold name'].notna() & gene_data['Chromosome/scaffold name'].str.isdigit()]
gene_data = gene_data.rename(columns = {'Chromosome/scaffold name':'Chrom', "Gene start (bp)":"Gene start", "Gene end (bp)":"Gene end", })


# transript types
pcg = [ 'protein_coding' ]
igg = [ 'IG_V_gene', 'IG_C_gene',  'IG_D_gene', 'IG_J_gene' ]
trg = [ 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', 'TR_C_gene' ]
lncRNA = [ 'lncRNA'  ]
ncRNA = ['snoRNA', 'misc_RNA', 'rRNA', 'miRNA', 'snRNA',  'sRNA',  'scaRNA', 'ribozyme', 'vault_RNA',  'scRNA']
pseudogene = [  'processed_pseudogene',   'unprocessed_pseudogene', 'transcribed_unprocessed_pseudogene',
 'IG_V_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene',
 'TR_V_pseudogene', 'translated_processed_pseudogene',  'TR_J_pseudogene', 'rRNA_pseudogene',
  'unitary_pseudogene', 'pseudogene', 'IG_C_pseudogene',  'translated_unprocessed_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene' ]
etc = [ 'protein_coding_CDS_not_defined', 'protein_coding_LoF', 'retained_intron', 'nonsense_mediated_decay',  'TEC', 'non_stop_decay' , 'artifact' ]


"""
ncRNA: A non-coding gene short ones
miRNA: A small RNA (~22bp) that silences the expression of target mRNA.
miscRNA: Miscellaneous RNA. A non-coding RNA that cannot be classified.
piRNA: An RNA that interacts with piwi proteins involved in genetic silencing.
rRNA: The RNA component of a ribosome.
siRNA: A small RNA (20-25bp) that silences the expression of target mRNA through the RNAi pathway.
snRNA: Small RNA molecules that are found in the cell nucleus and are involved in the processing of pre messenger RNAs
snoRNA: Small RNA molecules that are found in the cell nucleolus and are involved in the post-transcriptional modification of other RNAs.
tRNA: A transfer RNA, which acts as an adaptor molecule for translation of mRNA.
vaultRNA: Short non coding RNA genes that form part of the vault ribonucleoprotein complex.

TR gene: T cell receptor gene that undergoes somatic recombination, annotated in collaboration with IMGT http://www.imgt.org/.
TR C gene: Constant chain T cell receptor gene that undergoes somatic recombination before transcription
TR D gene: Diversity chain T cell receptor gene that undergoes somatic recombination before transcription
TR J gene: Joining chain T cell receptor gene that undergoes somatic recombination before transcription
TR V gene: Variable chain T cell receptor gene that undergoes somatic recombination before transcription
"""


def categorize_transcript_type(transcript_type):
    pcg = [ 'protein_coding' ]
    igg = [ 'IG_V_gene', 'IG_C_gene',  'IG_D_gene', 'IG_J_gene' ]
    trg = [ 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', 'TR_C_gene' ]
    lncRNA = [ 'lncRNA' ]
    ncRNA = ['snoRNA', 'misc_RNA', 'rRNA', 'miRNA', 'snRNA', 'sRNA', 'scaRNA', 'ribozyme', 'vault_RNA', 'scRNA']
    pseudogene = ['processed_pseudogene', 'unprocessed_pseudogene', 'transcribed_unprocessed_pseudogene', 'IG_V_pseudogene',
                  'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene', 'TR_V_pseudogene',
                  'translated_processed_pseudogene', 'TR_J_pseudogene', 'rRNA_pseudogene', 'unitary_pseudogene', 'pseudogene',
                  'IG_C_pseudogene', 'translated_unprocessed_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene']

    if transcript_type in pcg:
        return 'PCG'
    elif transcript_type in igg:
        return 'IGG'
    elif transcript_type in trg:
        return 'TRG'
    elif transcript_type in lncRNA:
        return 'lncRNA'
    elif transcript_type in ncRNA:
        return 'ncRNA'
    elif transcript_type in pseudogene:
        return 'pseudogene'
    else:
        return 'other'

# Apply the function to the 'Transcript type' column
gene_data['gene_category'] = gene_data['Transcript type'].map(categorize_transcript_type)

# considering
gene_categories = ["PCG", "lncRNA", "ncRNA", "pseudogene"]
gene_data_filtered = gene_data[ gene_data ['gene_category'].isin(gene_categories)]


# Create separate DataFrames for the longest, shortest, and average transcript length per gene
grouped = gene_data_filtered.groupby(['Gene stable ID', "gene_category" ] ) # "Gene name", "Gene type" , 'Transcript type', "Chrom", "Gene start", "Gene end" ])

average_transcripts = grouped['Transcript length (including UTRs and CDS)'].mean().reset_index()
average_transcripts['Length type'] = 'average'


#longest_transcripts2 = grouped['Transcript length (including UTRs and CDS)'].max().reset_index()
#longest_transcripts2['Length type'] = 'longest'

#shortest_transcripts = grouped['Transcript length (including UTRs and CDS)'].min().reset_index()
#shortest_transcripts['Length type'] = 'shortest'

# Concatenate the DataFrames together
# Print the first few rows of the final DataFrame

idx = gene_data_filtered.groupby([ 'Gene stable ID', 'gene_category'])['Transcript length (including UTRs and CDS)'].idxmax()
longest_transcripts = gene_data_filtered.loc[idx]

idx = gene_data_filtered.groupby([ 'Gene stable ID', 'gene_category'])['Transcript length (including UTRs and CDS)'].idxmin()
shortest_transcripts = gene_data_filtered.loc[idx]

longest_transcripts['Length type'] = 'longest'
shortest_transcripts['Length type'] = 'shortest'

geneindex = longest_transcripts.columns

average_merged = pd.merge( average_transcripts, longest_transcripts,  on='Gene stable ID')
average_merged = average_merged.drop( columns = [ "Transcript length (including UTRs and CDS)_y", "gene_category_y" ])
average_merged = average_merged.rename(columns={"gene_category_x": "gene_category", "Transcript length (including UTRs and CDS)_x": "Transcript length (including UTRs and CDS)"})
average_merged[ "Transcript stable ID"] = None
average_merged[ "Protein stable ID"] = None


final_data = pd.concat([longest_transcripts, shortest_transcripts, average_merged], ignore_index=True)



# check distribution of lengthes

plt.figure(figsize=(16,6))

# Longest transcript length
plt.subplot(141)
sns.barplot(x="gene_category", y="Transcript length (including UTRs and CDS)", hue="Length type", data=final_data[final_data["gene_category"]=="PCG"])
plt.title('Protein coding genes')
plt.ylim(0, 4000 )

# Shortest transcript length
plt.subplot(142)
sns.barplot(x="gene_category", y="Transcript length (including UTRs and CDS)", hue="Length type", data=final_data[final_data["gene_category"]=="lncRNA"])
plt.title('long non-cording RNA')
plt.ylim(0, 4000 )

# Shortest transcript length
plt.subplot(144)
sns.barplot(x="gene_category", y="Transcript length (including UTRs and CDS)", hue="Length type", data=final_data[final_data["gene_category"]=="ncRNA"])
plt.title('non-cording genes')
plt.ylim(0, 4000 )

# Shortest transcript length
plt.subplot(143)
sns.barplot(x="gene_category", y="Transcript length (including UTRs and CDS)", hue="Length type", data=final_data[final_data["gene_category"]=="pseudogene"])
plt.title('Pseudogene')
plt.ylim(0, 4000 )

plt.tight_layout()
plt.show()



# count # of unique gene # IDEA:
unique_gene_counts = final_data.groupby('gene_category')['Gene stable ID'].nunique()


# write to csv
final_data.to_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype.csv', index=False)



# genes tested for cis-eQTL in GTEx

import os
import pandas as pd
import glob

# a few genes have difference Stable gene ID for sthe same gene name.
# selection one recommended on ensembl

gene_biotype = pd.read_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_nonredun.csv' , low_memory = False )
gene_biotype_gid = gene_biotype[ "Gene stable ID"].tolist()


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
    matching_ids = df[df['gene_id'].isin( gene_biotype_gid )]

    # Append the gene_ids to the all_matching_ids list
    all_matching_ids.extend(matching_ids['gene_id'].tolist())



#
# Gene stable IDs included in all gene sets tested for cis-eQTL in GTEx
# for downstream analysis
gene_biotype_gtex = gene_biotype[ gene_biotype[ "Gene stable ID"].isin( all_matching_ids ) ]
gene_biotype_gtex.to_csv('/mnt/d/sandbox/0.gxe/data/gene_biotype_GTEx.csv', index=False)
