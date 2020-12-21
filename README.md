# BEstimate (beta)

BEstimate, a Python module that systematically analyses guide RNA (gRNA) targetable sites across given sequences for given Base Editors, and functional and clinical effects of the potential edits on the resulting proteins. It has the ability to provide in silico analysis of the sequences to identify positions that can be editable by Base Editors, and their features before starting experiments. 

*Warning: Protein level analysis is needed to be improved.* 
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
ENSEMBL TRANSCRIPT ID (Optional input): The interested transcript for filtering the result 

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

python3 BEstimate.py -gene BRCA1 -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -ofile ./

The user also run the same analysis for different PAM only changing -pamseq NGN. (*Warning: Be caferul to write the PAM sequence to be in concordant with the length of the -pamwin. Here, NGN is in concordant with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.*) 

## Contact

BEstimate is the product of Cansu Dincer, Dr Matthew Coelho and Dr Mathew Garnett from Garnett Group at the Wellcome Sanger Institute.

For any problems or feedback on BEstimate, you can contact [here](mailto:cd7@sanger.ac.uk).

## Terms and conditions

**BEstimate** is a Python tool that systematically analyses guide RNA (gRNA) targetable sites across given sequences for given Base Editors, and functional and clinical effects of the potential edits on the resulting proteins

Copyright © 2020 Wellcome Sanger Institute

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Neither the institution name nor the name BEstimate can be used to endorse or promote products derived from this software without prior written permission. 

Products derived from this software may not be called BEstimate nor may BEstimate appear in their names without prior written permission of the developers.

You should have received a copy of the GNU General Public License along with this program. If not, see [gnu.org/licenses/](http://www.gnu.org/licenses/).

You can also view the licence [here.](http://www.gnu.org/licenses/lgpl-3.0.txt)
