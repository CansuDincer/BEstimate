# BEstimate v7

BEstimate, a Python module that systematically analyses guide RNA (gRNA) targetable sites across given sequences for given Base Editors, and functional and clinical effects of the potential edits on the resulting proteins. It has the ability to provide in silico analysis of the sequences to identify positions that can be editable by Base Editors, and their features before starting experiments. 

## Requirements

Python 3.8
pandas 1.1.3
argparse 1.1
requests 2.24.0

## Inputs

*-gene* 
GENE NAME (Mandatory input): Hugo symbol of the gene of intereset 

*-assembly*
GENOME ASSEMBLY (Mandatory input - hg19/GRCh38): The genome assembly for interested genomic coordinates 

*-transcript*
ENSEMBL TRANSCRIPT ID (Optional input): The interested transcript for filtering the result, itherwise it uses canonical transcript from Ensembl which is basically the longest transcript. 

*-pamseq*
PAM SEQUENCE (Mandatory input - default = NGG): The sequence preference of the Cas9 protein 

*-pamwin*
PAM INDICES (Mandatory input - default = 21-23): The indices of the PAM sequence while counting 1 as the first nucleotide in the protospacer sequence. 

*-actwin*
ACTIVITY WINDOW INDICES (Mandatory input - default = 4-8): The indices of the activity window where the editable nucleotides will be searched.

*-protolen*
PROTOSPACER SEQUENCE LENGTH (Mandatory input - default = 20): The length of the protospacer sequence.

*-edit*
EDIT BASE (Mandatory input - default = C): The nucleotide which Base Editor can edit. 

*-edit_to*
EDITED BASE (Mandatory input - default = T): The nucleotide which Base Editor can alter the EDIT BASE into.

*-P*
PROTEIN ANALYSIS (Mandatory input - True/False - default = False): The boolean input for protein analysis. When it is True, the editable sites will be analysed by VEP API and Proteins API for their functional consequences. *Warning: Only one edit at a time can be analysed. Several changes in the activity window has not been implemented for further analysis yet.*

*-o*
OUTPUT PATH (Optional input - default = working directory): The interested output path where the files will be written. 

*-ofile*
OUTPUT INITIALS (Mandatory input): The initial name of the file before "_crispr_df.csv", "_edit_df.csv" or "_vep_df.csv".

## Examples

For evaluation with [CRISPR Finder](https://wge.stemcell.sanger.ac.uk//find_crisprs#Grch38/BRCA1), we tried to run BEstimate with GRCh38 genome assembly for *BRCA1* gene:

python3 BEstimate_v7.py -gene BRCA1 -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -ofile ./

The user also run the same analysis for different PAM only changing -pamseq NGN. (*Warning: Be caferul to write the PAM sequence to be in concordant with the length of the -pamwin. Here, NGN is in concordant with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.*) 

If you would like to run for a specific transcript and run the protein analysis:

python3 BEstimate_v7.py -gene BRAF -assembly GRCh38 **-transcript ENST00000646891** -edit C -edit_to T **-P** -ofile ./

## Contact

BEstimate is the product of Cansu Dincer, Dr Matthew Coelho and Dr Mathew Garnett from Garnett Group at the Wellcome Sanger Institute.

For any problems or feedback on BEstimate, you can contact [here](mailto:cd7@sanger.ac.uk).


