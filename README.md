# BEstimate

BEstimate, a Python module that systematically identifies guide RNA (gRNA) targetable sites across given sequences for given Base Editors, functional and clinical effects of the potential edits on the resulting proteins and off target consequence of the found sequences. It has the ability to provide in silico analysis of the sequences to identify positions that can be editable by Base Editors, and their features before starting experiments. 

## Requirements

Python 3.8
pandas 1.1.3
argparse 1.4
biopython 1.78
requests 2.28.1

## Inputs

#### Sequence options:
*-gene* 
GENE NAME (Mandatory input): Hugo symbol of the gene of interest 

*-assembly*
GENOME ASSEMBLY (Mandatory input - hg19/GRCh38): The genome assembly for interested genomic coordinates 

*-transcript*
ENSEMBL TRANSCRIPT ID (Optional input): The interested transcript for filtering the result, otherwise it uses canonical transcript from Ensembl. BEstimate first tries to MANE selected transcript, if not the Ensembl canonical transcript is obtained. 

*-mutation*
MUTATION (Optional input - default= None): In the case that there is a mutation on the interested gene that you need to integrate into sequence to design gRNAs according to that. The mutation style should be in <chromosome:g.genomic_location edited_nucleotide>new_nucleotide> e.g. (3:g.179218303G>A).

*-mutation_file*
MUTATION FILE (Optional input - default= None): In the case that there are more than one mutations to be integrated into the sequence to design gRNAs according to them. The mutation style should be in <chromosome:g.genomic_location edited_nucleotide>new_nucleotide> e.g. (3:g.179218303G>A), and file should have one mutation in each row. 

*-flank*
FLANKING SEQUENCE (Optional input - True/False - default = False):  The boolean input specifying whether the user wants to retrieve the flanking sequences 3' and 5' of gRNA sequences.

*-flank3*  
3' FLANKING SEQUENCE LENGTH (Optional input - default= 7): If *-flank* input is provided, *flank3* input specifies the number of nucleotides in the 3' flanking region. 
  
 *-flank5*  
5' FLANKING SEQUENCE LENGTH (Optional input - default= 11): If *-flank* input is provided, *flank5* input specifies the number of nucleotides in the 5' flanking region. 


#### Base Editor options:
*-pamseq*
PAM SEQUENCE (Mandatory input - default = NGG): The sequence preference of the Cas9 protein.

*-pamwin*
PAM INDICES (Mandatory input - default = 21-23): The indices of the PAM sequence while counting 1 as the first nucleotide in the protospacer sequence. 

*-actwin*
ACTIVITY WINDOW INDICES (Mandatory input - default = 4-8): The indices of the activity window where the editable nucleotides will be searched on protospacer sequence.

*-protolen*
PROTOSPACER SEQUENCE LENGTH (Mandatory input - default = 20): The length of the protospacer sequence.

*-edit*
EDIT BASE (Mandatory input - default = C): The nucleotide which Base Editor can edit. 

*-edit_to*
EDITED BASE (Mandatory input - default = T): The nucleotide which Base Editor can alter the EDIT BASE into.

#### Annotation options:
*-vep*
ENSEMBL VEP and PROTEIN ANALYSIS (Optional input - True/False - default = False): The boolean input for VEP and protein analysis. When it is True, the editable sites will be analysed by [VEP API from Ensembl](https://rest.ensembl.org/) and [Proteins API from Uniprot](https://www.ebi.ac.uk/proteins/api/doc/) for their functional consequences on proteins. As well as the post translational modification and domain information, if the resulting edit is on the interface region of the corresponding protein is also given by using [Interactome Insider](http://interactomeinsider.yulab.org/downloads.html). 

*-ot*
OFF TARGET (Optional input - True/False - default= False): The boolean input for identification of off targets.

*-mm*
MISMATCH (Optional input - default= 4): **In the case that *-ot* provided**, Number of maximum mismatches allowed in off target analysis.  

*-genome*
GENOME (Optional input - default= Homo_sapiens_GRCh38_dna_sm_all_chromosomes): **In the case that *-ot* provided**, the name of the genome file in ./BEstimate/offtargets/genome/.

#### Output options:
*-o*
OUTPUT PATH (Optional input - default = working directory): The interested output path where the files will be written. 

*-ofile*
OUTPUT INITIALS (Mandatory input): The initial name of the file before "_crispr_df.csv", "_edit_df.csv" or "_hgvs_df.csv", "_vep_df.csv", "_protein_df.csv", "_summary_df.csv".

## Examples

`python3 BEstimate.py -gene BRCA1 -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -o ./output/ -ofile BRCA1_CBE_NGG`

The user also run the same analysis for different PAM only changing -pamseq NGN. 

*Warning: Be caferul to write the PAM sequence to be in concordant with the length of the -pamwin. Here, NGN is in concordant with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.* 

If you would like to run for a specific transcript and run the protein analysis:

`python3 BEstimate.py -gene BRAF -assembly GRCh38 -transcript ENST00000646891 -edit C -edit_to T -vep -o ./ -ofile BRAF_CBE_NGG`

If you would like to run with a specific point mutation, with NGN PAM and with VEP and protein analysis:

`python3 BEstimate.py -gene PIK3CA -assembly GRCh38 -pamseq NGN -pamwin 21-23 -actwin 4-8 -protolen 20 -mutation '3:g.179218303G>A' -edit A -edit_to G -vep -ofile ./PIK3CA_NGN_ABE_mE545K -o ./output/`

If you would like to see the off targets of WRN gene:
`python3 BEstimate.py -gene BRAF -assembly GRCh38 -pamseq NGN -edit A -edit_to G -vep -ot -mm 4 -o ./output/ -ofile BRAF_ABE_NGN`

## Interpretation of BEstimate results

BEstimate will produce different files.

**CRISPR_df** represents all possible gRNA target sites from both directions.
- Hugo_Symbol: Hugo Symbol of the corresponding gene location.
- CRISPR_PAM_Sequence gRNA target sequence with including PAM region.
- gRNA_Target_Sequence: gRNA target sequence without PAM region.
- Location: Location of the gRNA target sequence (chromosome:start location:end location).
- Direction: Direction of the gRNA.
- Gene_ID: Ensembl Gene ID of the interested Hugo Symbol.
- Transcript_ID: Ensembl Transcript ID of the interested Hugo Symbol in the corresponding gene location.
- Exon_ID:  Ensembl Exon ID of the interested Hugo Symbol in the corresponding gene location.
- guide_in_CDS: If gRNA has any nucleotide inside a coding sequence of the gene. 
- gRNA_flanking_sequences: In case that user has given *-flan* option. 

**EDIT_df** adds information for all possible gRNA target sites if they can be edited with the given base editors and edit/edit_to parameters, if the edit is on the exon or not. Additional columns are below:
- Edit_in_Exon: Boolean output representing if the specified edit in the Edit Location happening on the exon or not.
- Edit_in_CDS: Boolean output representing if the specified edit in the Edit Location happening on the coding sequence or not.
- \# Edits/guide: The number of editable nucleotide in gRNAs.
- Poly_T:  Boolean output representing if CRISPR_PAM_Sequence has consecutive 4 T nucleotides.
- guide_on_mutation: Boolean output representing if given mutation location inside gRNA (in the case that mutation information given by the user).
- guide_change_mutation: Boolean output representing if gRNA can change the mutated nucleotide (in the case that mutation information given by the user).

**HGVS_df** is a transient data frame including HGVS nomenclature of the potential edits for further VEP analysis. Additional columns are below:

- HGVS:  HGVS nomenclature of the potential edit.

**VEP_df** represents functional and clinical consequences of all possible edits in all possible gRNAs. Additional columns are below:

- Protein_ID:  Ensembl Protein ID of the potential edit.
- VEP_input: HGVS nomenclature of the potential edit.
- allele: Allelic change of the potential edit.
- variant_classification: VEP classification of the resulting variant of the potential edit. (Substitution, SNV etc.)
- most_severe_consequence: The most severe consequence of the resulting variant of the potential edit.
- consequence_terms: List of functional consequences of the resulting variant of the potential edit (splice region, missense, stop gain etc.)
- variant_biotype:  The biotype of the variant created by the potential edit.
- Regulatory_ID: Ensembl Regulatory ID of the potential edit.
- Motif_ID: Ensembl Model ID of the potential edit.
- TFs_on_motif: List of transcription factors on the Ensembl Regulatory ID of the potential edit.
- cDNA_Change: cDNA changes resulting from the possible edits with the corresponding guides (position old nucleotide> new nucleotide)
- Edited_codon: Codon sequence before the corresponding edit.
- New_codon: Codon sequence after the corresponding edit.
- CDS_Position: Position on the coding sequence of the potential edit.
- Protein_Position: Location of the corresponding edit on the resulting protein (as Uniprot index).
- Protein_Change: Protein sequence change with the corresponding edit (old amino asid - position - new amino asid).
- Edited_AA: One letter representation of the amino asid before the corresponding edit.
- Edited_AA_Prop: Chemical properties of the amino asid before the corresponding edit.
- New_AA: One letter representation of the amino asid after the corresponding edit.
- New_AA_Prop: Chemical properties of the amino asid after the corresponding edit.
- is_Synonymous: Boolean output representing if the resulting edit causes synonymous or non-synonymous mutations.
- is_Stop: Boolean output representing if the resulting edit causes stop codon or not.
- proline_addition: Boolean output representing if potential edit created a Proline amino acid or not. 
- swissprot: SwissProt ID of the corresponding gRNA target site region from Ensembl VEP.
- polyphen_score: Polyphen Score of the corresponding edit.
- polyphen_prediction: Polyphen Prediction of the corresponding edit.
- sift_score: Sift Score of the corresponding edit.
- sift_prediction: Sift Prediction of the corresponding edit.
- cadd_phred: Curated CADD Score of the corresponding edit.
- cadd_raw: Raw CADD Score of the corresponding edit.
- lof: representing if the corresponding edit causes Loss of function with high confidence (HC)/Low confidence (LC) or not by implementing [LOFTEE](https://github.com/konradjk/loftee) through VEP.
- impact: Impact of the corresponding edit.
- blosum62: blosum62 Score of the corresponding edit.
- is_clinical: Boolean output representing if there is a collocated clinical variant on the potential edit site.
- clinical_id: dbSNP id of the clinical allele.
- clinical_significance: Clinical significance of the corresponding edit.
- cosmic_id: COSMIC id of the clinical allele.
- clinvar_id: ClinVar id of the clinical allele.
- ancestral_populations: The conservation score from the Ensembl Compara databases  
for potential edit site (if any). 

**protein_df** represents functional consequences of all possible edits on corresponding protein. Additional columns are below:

- Domain: Functional domain name of the potential edit site.
- curated_Domain: Human friendly functional domain name of the potential edit site.
- PTM: Post translational modificaton on the potential edit site.
- is_disruptive_interface_EXP: Boolean output representing if the corresponding edit happening on the interface region and having an experimental evidence (PDB). *High confidence*
- is_disruptive_interface_MOD: Boolean output representing if the corresponding edit happening on the interface region and having evidence from a model (Interactome3D). *Medium confidence*
- is_disruptive_interface_PRED: Boolean output representing if the corresponding edit happening on the interface region and having evidence from prediction (Interactome Insider). *Low confidence*
- disrupted_PDB_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to experimental evidence - PDB. 
- disrupted_I3D_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to Interactome3D.
- disrupted_Eclair_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to Eclair - prediction algorithm of Interactome Insider. 
- disrupted_PDB_int_genes: List of partner genes (Hugo Symbols) whose interactions are disrupted by the corresponding edit according to experimental evidence - PDB. 
- disrupted_I3D_int_genes: List of partner genes (Hugo Symbols) whose interactions are disrupted by the corresponding edit according to Interactome3D.
- disrupted_Eclair_int_genes: List of partner genes (Hugo Symbols) whose interactions are disrupted by the corresponding edit according to Eclair - prediction algorithm of Interactome Insider. 

**summary_df** summarises all above information (each gRNA annotation has only one row). 

**ot_annotated_summary_df** adds off target information on summary file.  Additional columns are below:

- exact: The number of exact alignments of gRNA sequence.
- mm1: The number of alignments of gRNA sequence with one mismatch.
- mm2: The number of alignments of gRNA sequence with two mismatches.
- mm3: The number of alignments of gRNA sequence with three mismatches.
- mm4: The number of alignments of gRNA sequence with four mismatches.

## Contact

BEstimate is the product of Cansu Dincer, Dr Matthew Coelho and Dr Mathew Garnett from Garnett Group at the Wellcome Sanger Institute.

For any problems or feedback on BEstimate, you can contact [here](mailto:cd7@sanger.ac.uk).

## License

BEstimate, a Python module that systematically identifies guide RNA (gRNA) 
on and off target sites across given sequences for given Base Editors, and functional and clinical effects of the potential edits on the resulting proteins. 

Copyright (c) 2020-2023 Genome Research Ltd. 

Author: Cansu Dincer <cd7@sanger.ac.uk> 

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation; either version 3 of the License, or (at your option) any later 
version. 

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details. 

You should have received a copy of the GNU General Public License along with 
this program. If not, see <http://www.gnu.org/licenses/>. 

## Further Disclaimer

For policies regarding the underlying data, please also refer to:
- [Ensembl terms and conditions](https://www.ensembl.org/info/about/legal/code_licence.html#:~:text=Subject%20to%20the%20terms%20and,the%20Work%20and%20such%20Derivative)
- [Uniprot terms and conditions](https://www.uniprot.org/help/license#:~:text=We%20make%20no%20warranties%20regarding,by%20patents%20or%20other%20rights.)
- [Interactome Insider terms and conditions](http://interactomeinsider.yulab.org/)


