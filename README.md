# BEstimate

BEstimate, a Python module that systematically analyses guide RNA (gRNA) targetable sites across given sequences for given Base Editors, and functional and clinical effects of the potential edits on the resulting proteins. It has the ability to provide in silico analysis of the sequences to identify positions that can be editable by Base Editors, and their features before starting experiments. 

## Requirements

Python 3.8
pandas 1.1.3
argparse 1.1
requests 2.24.0

## Inputs

*-gene* 
GENE NAME (Mandatory input): Hugo symbol of the gene of interest 

*-assembly*
GENOME ASSEMBLY (Mandatory input - hg19/GRCh38): The genome assembly for interested genomic coordinates 

*-transcript*
ENSEMBL TRANSCRIPT ID (Optional input): The interested transcript for filtering the result, otherwise it uses canonical transcript from Ensembl which is the longest transcript (it is important to check before analysis). 

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
PROTEIN ANALYSIS (Mandatory input - True/False - default = False): The boolean input for protein analysis. When it is True, the editable sites will be analysed by [VEP API from Ensembl](https://rest.ensembl.org/) and [Proteins API from Uniprot](https://www.ebi.ac.uk/proteins/api/doc/) for their functional consequences. Additionally, the information if the resulting edit is on the interface region of the corresponding protein or not is also given by using [Interactome Insider](http://interactomeinsider.yulab.org/downloads.html). *Warning: Only one edit at a time can be analysed. Effect of several changes in the activity window as a whole has not been implemented.*

*-o*
OUTPUT PATH (Optional input - default = working directory): The interested output path where the files will be written. 

*-ofile*
OUTPUT INITIALS (Mandatory input): The initial name of the file before "_crispr_df.csv", "_edit_df.csv" or "_analysis_df.csv".

## Examples

For evaluation with [CRISPR Finder](https://wge.stemcell.sanger.ac.uk//find_crisprs#Grch38/BRCA1), we tried to run BEstimate with GRCh38 genome assembly for *BRCA1* gene:

python3 BEstimate.py -gene BRCA1 -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -ofile ./

The user also run the same analysis for different PAM only changing -pamseq NGN. (*Warning: Be caferul to write the PAM sequence to be in concordant with the length of the -pamwin. Here, NGN is in concordant with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.*) 

If you would like to run for a specific transcript and run the protein analysis:

python3 BEstimate_v7.py -gene BRAF -assembly GRCh38 **-transcript ENST00000646891** -edit C -edit_to T **-P** -ofile ./

## Results

BEstimate will produce two files if **-P** option is False, otherwise three files will be generated.

**CRISPR_df** represents all possible gRNA target sites from both dicrections.
- Hugo Symbol: Hugo Symbol of the corresponding gene location.
- CRISPR+PAM Sequence: gRNA target sequence with including PAM region.
- gRNA Target Sequence: gRNA target sequence.
- Location: Location of the gRNA target sequence (chromosome:start location:end location).
- Direction: Direction of the gRNA.
- Gene_ID: Ensembl Gene ID of the interested Hugo Symbol.
- Transcript_ID: Ensembl Transcript ID of the interested Hugo Symbol in the corresponding gene location.
- Exon_ID:  Ensembl Exon ID of the interested Hugo Symbol in the corresponding gene location.

**EDIT_df** adds information for all possible gRNA target sites if they can be edited with the given base editors and edit/edit_to parameters, if the edit is on the exon or not.
- Hugo Symbol: Hugo Symbol of the corresponding gene location.
- CRISPR+PAM Sequence: gRNA target sequence with including PAM region.
- gRNA Target Sequence: gRNA target sequence.
- Location: Location of the gRNA target sequence (chromosome:start location:end location).
- Edit Location: Location of the possible edit with the given Base Editor.
- Direction: Direction of the gRNA.
- Gene_ID: Ensembl Gene ID for the corresponding gRNA target site region.
- Transcript_ID: Ensembl Transcript ID for the corresponding gRNA target site region.
- Exon_ID:  Ensembl Exon ID for the corresponding gRNA target site region.
- Edit_in_Exon: Boolean output representing if the specified edit in the Edit Location happening on the exon or not.

**ANALYSIS_df** represents functional and clinical consequences of all possible edits in all possible gRNAs.
- Hugo Symbol: Hugo Symbol of the corresponding gene location.
- gRNA_target_sequence: gRNA target sequence.
- gRNA_target_location: Location of the gRNA target sequence (chromosome:start location:end location).
- Genomic_Position: Location (genomic position) of the possible edit with the given Base Editor.
- Direction: Direction of the gRNA.
- Transcript_ID: Ensembl Transcript ID for the corresponding gRNA target site region.
- Exon_ID:  Ensembl Exon ID for the corresponding gRNA target site region.
- cDNA_Change: cDNA changes resulting from the possible edits with the corresponding guides (position old nucleotide> new nucleotide)
- Edited_codon: Codon sequence before the corresponding edit.
- New_codon: Codon sequence after the corresponding edit.
- Protein_ID: Ensembl Protein ID for the corresponding gRNA target site region.
- Protein_Position: Location of the corresponding edit on the resulting protein (as Uniprot index).
- Protein_Change: Protein sequence change with the corresponding edit (old amino asid - position - new amino asid).
- Edited_aa: One letter representation of the amino asid before the corresponding edit.
- Edited_aa_prop: Chemical properties of the amino asid before the corresponding edit.
- New_aa: One letter representation of the amino asid after the corresponding edit.
- New_aa_prop: Chemical properties of the amino asid after the corresponding edit.
- is_synonymous: Boolean output representing if the resulting edit causes synonymous or non-synonymous mutations.
- is_stop: Boolean output representing if the resulting edit causes stop codon or not.
- polyphen_score: Polyphen Score of the corresponding edit.
- polyphen_prediction: Polyphen Prediction of the corresponding edit.
- sift_score: Sift Score of the corresponding edit.
- sift_prediction: Sift Prediction of the corresponding edit.
- cadd_phred: Curated CADD Score of the corresponding edit.
- cadd_raw: Raw CADD Score of the corresponding edit.
- lof: representing if the corresponding edit causes Loss of function with high confidence (HC)/Low confidence (LC) or not by implementing [LOFTEE](https://github.com/konradjk/loftee) through VEP.
- impact: Impact of the corresponding edit.
- blosum62: blosum62 Score of the corresponding edit.
- consequence_terms: representing the functional consequence of the corresponding edit such as splice region, missense, stop gain etc.
- is_ClinVar: representing if the corresponding edit has Clinical consequences in ClinVar or not. 
- clinical_allele: Nucleotide changes (old nucleotide > new nucleotide) from the corresponding edit with a Clinical consequence.
- clinical_id: dbSNP id of the clinical allele.
- clinical_significance: Clinical significance of the corresponding edit.
- Domain: Protein domain in which corresponding edit happening.
- Phosphorylation: Phosphorylation sites in which corresponding edit happening.
- Uniprot: Uniprot ID for the corresponding gRNA target site region.*(There might be several Uniprot IDs -SwissProt and Trembl- Each of them was represented separately.)*
- Reviewed: Boolean output representing ig the corresponding Uniprot ID is reviewed (SwissProt) or not (Trembl).
- is_disruptive_interface_EXP: Boolean output representing if the corresponding edit happening on the interface region and having an experimental evidence (PDB). *High confidence*
- is_disruptive_interface_MOD: Boolean output representing if the corresponding edit happening on the interface region and having evidence from a model (Interactome3D). *Medium confidence*
- is_disruptive_interface_PRED: Boolean output representing if the corresponding edit happening on the interface region and having evidence from prediction (Interactome Insider). *Low confidence*
- disrupted_PDB_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to experimental evidence - PDB. 
- disrupted_I3D_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to Interactome3D.
- disrupted_Eclair_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to Eclair - prediction algorithm of Interactome Insider. 

## Contact

BEstimate is the product of Cansu Dincer, Dr Matthew Coelho and Dr Mathew Garnett from Garnett Group at the Wellcome Sanger Institute.

For any problems or feedback on BEstimate, you can contact [here](mailto:cd7@sanger.ac.uk).


