# BEstimate

BEstimate, a Python module that systematically identifies guide RNA (gRNA) targetable sites across given sequences for given Base Editors, functional and clinical effects of the potential edits on the resulting proteins and off target consequence of the found sequences. It has the ability to provide in silico analysis of the sequences to identify positions that can be editable by Base Editors, and their features before starting experiments. 

## Table of Contents
- [Quick start installation](#quick-start-installation)
- [Run BEstimate](#run-bestimate)
    - [Examples with BEstimate](#examples-with-bestimate)
    - [Off-Targets Analysis](#off-targets-analysis)
    - [Command line usage and options](#command-line-usage-and-options)
    - [Output Interpretation](#output-interpretation)
- [Contact](#contact)
- [License](#license)

## Quick start installation

You can directly use BEstimate environment if you have conda. Please follow below:

- `git clone https://github.com/CansuDincer/BEstimate.git`
- `cd BEstimate`
- `conda-env create -n bestimate -f=bestimate.yml`
- `conda activate bestimate`

If not, you should have python 3.13 and you can use requirements file:

- `pip3 install -r requirements.txt`

## Run BEstimate

### Examples with BEstimate

For example, if you would like to run for the *SRY* gene with NGG PAM sequence, with CBE (C to T editing) and without VEP and protein analysis:

```bash
BEstimate -gene SRY -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -o ../output/ -ofile SRY_CBE_NGG
```

The user also run the same analysis for different PAM only changing -pamseq NGN.

*Warning: Be careful to write the PAM sequence to be in concordance with the length of the -pamwin. Here, NGN is in concordance with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.*

If you would like to run for a specific transcript and run the protein analysis:

```bash
BEstimate -gene SRY -assembly GRCh38 -transcript ENST00000383070 -edit C -edit_to T -vep -o ../output/ -ofile SRY_CBE_NGG
```

If you would like to run with a specific point mutation, with NGN PAM and with VEP and protein analysis:
Prepare a `PIK3CA_mutation_file.txt` for example with 3:g.179218303G>A

```bash
BEstimate -gene PIK3CA -assembly GRCh38 -pamseq NGN -pamwin 21-23 -actwin 4-8 -protolen 20 -mutation_file PIK3CA_mutation_file.txt -edit A -edit_to G -vep -ofile PIK3CA_NGN_ABE_mE545K -o ../output/
```


### Off-Targets Analysis

To run the off-target analysis, first you need to have the [Ensembl](https://www.ensembl.org/) Genome indexed for the interested PAM sequence. 

The `x_genome.py` script will download the required files and index the genome for CRISPRs as follows.
- Download the specified FASTA genome assembly files from the Ensembl project,
- Gather CRISPRs from the FASTA files into CSV files detailing chromosome, position in chromosome, as well as PAM position,
- Generate a binary list of gRNA signatures (accounting for PAM position),
- Insert the CRISPRs into a SQLite database for cross-referencing the gRNAs found in the binary list.

to run the **x_genome.py** script:

The parameters are:
- *-p --pamseq* - the PAM sequence that will be used for indexing CRISPRs, e.g. "NGG" - *Required*,
- *-a --assembly* - the human reference genome targeted, "GRCh37" or "GRCh38" - *Required*,
- *-o --output_path* - the output directory, if not specified, will use the current directory - *Optional*,
- *-e --ensembl_version* - the Ensembl version of the genome being indexed (currently default is 113 for GRCh38, if using GRCh37, please use <=75) - *Required*,
- *-ot --offtargets_path* - path for the downloaded genome for the off-target analysis, defaults to "../offtargets" - *Optional*

For example:

```bash
python3 x_genome.py --pamseq NGN --assembly GRCh38 --ensembl_version 113
```

The gathering of CRISPRs from the genome assembly takes a while and requires a fair amount of disk storage. For example, using the GRCh38 genome assembly:

| Pam Sequence | Space (GB) | Run Time  |
| ------------ | ---------- | --------- |
| NGG          | 38         | ~3 Hours  |
| NGN          | 140        | ~9 Hours  |

Then, you can run the off-target analysis, see below for the *BRAF* gene:

```bash
python3 BEstimate.py -gene BRAF -assembly GRCh38 -pamseq NGN -edit A -edit_to G -vep -ot -o ../output -ot_path ../offtargets -ofile BRAF_ABE_NGN
```

### Command line usage and options

There are three programs when the package is installed available from the command line:
- `BEstimate` - the main program to find and analyse Base Editor sites
- `x_genome` - the program to download and index a genome for off-target analysis
- `x_crispranalyser` - the program to run off-target analysis on guides

<details>
<summary>Expand to see <strong>BEstimate</strong> command line options</summary>

```bash
BEstimate --help
usage: BEstimate [inputs]

********************************** Find and Analyse Base Editor sites **********************************

Mandatory Inputs:
  -h, --help            show this help message and exit
  --version             Show program's version number and exit.
  -gene GENE            The hugo symbol of the interested gene!
  -assembly ASSEMBLY    The genome assembly that will be used!
  -transcript TRANSCRIPT
                        The interested ensembl transcript id
  -uniprot UNIPROT      The interested Uniprot id
  -pamseq PAMSEQ        The PAM sequence in which features used for searching activity window and editable nucleotide.
  -pamwin PAMWINDOW     The index of the PAM sequence when starting from the first index of protospacer as 1.
  -actwin ACTWINDOW     The index of the activity window when starting from the first index of protospacer as 1.
  -protolen PROTOLEN    The total protospacer and PAM length.
  -vep                  The boolean option if user wants to analyse the edits through VEP and Uniprot.
  -mutation_file MUTATION_FILE
                        A file for the mutations on the interested gene that you need to integrate into guide and/or annotation analysis
  - library_file LIBRARY_FILE
                        Existing library file: should include:
			  1. gRNA + PAM Sequence (5>3'): column name -> CRISPR_PAM_Sequence
			  2. gRNA genomic location (smallest>largest): column name -> Location
			  3. gRNA directionality: column name -> Direction
  -flank                The boolean option if the user wants to add flanking sequences of the gRNAs
  -flank3 FLAN_3        The number of nucleotides in the 3' flanking region
  -flank5 FLAN_5        The number of nucleotides in the 5' flanking region
  -rs3 RS3              The boolean option if the user wants to add on target RuleSet3 scoring for the gRNAs
  -fc FCAST             The boolean option if the user wants to add on target ForeCast gRNAs efficiency info
  -edit {A,T,G,C}       The nucleotide which will be edited.
  -edit_to {A,T,G,C}    The nucleotide after edition.
  -o OUTPUT_PATH        The path for output. If not specified the current directory will be used!
  -ofile OUTPUT_PATH    The output file name, if not specified "position" will be used!
  -ot OFF_TARGET        The boolean option if off targets will be computed or not
  -genome GENOME        (If -ot provided) name of the genome file
  -v_ensembl VERSION    The ensembl version in which genome will be retrieved (if the assembly is GRCh37 then please use <=75)
  -ot_path OT_PATH      The path of the offtargets folder. default is os.getcwd() + "/../offtargets/
```
</details>

<details>
<summary>Expand to see <strong>x_genome</strong> command line options</summary>

```bash
usage: x_genome [inputs]

Script for indexing CRISPRs for finding off-targets

options:
  -h, --help            show this help message and exit
  --version             Show the version number and exit
  --pamseq PAMSEQ, -p PAMSEQ
                        The PAM sequence in which features used for searching activity window and editable nucleotide.
  --assembly {GRCh38,GRCh37}, -a {GRCh38,GRCh37}
                        The genome assembly that will be used!
  --output_path OUTPUT_PATH, -o OUTPUT_PATH
                        The path for output. If not specified the current directory will be used!
  --ensembl_version ENSEMBL_VERSION, -e ENSEMBL_VERSION
                        The ensembl version in which genome will be retrieved (if the assembly is GRCh37 then please use <=75)
  --offtargets_path OFFTARGETS_PATH, -ot OFFTARGETS_PATH
                        The path to the root offtargets output directory
```

</details>

<details>
<summary>Expand to see <strong>x_crispranalyser</strong> command line options</summary>

```bash
usage: x_crispranalyser [inputs]

Script for finding off-targets

options:
  -h, --help            show this help message and exit
  --version             Show the version number and exit
  --input_csv INPUT_CSV, -i INPUT_CSV
                        The input CSV file to be analysed
  --binary_index BINARY_INDEX, -b BINARY_INDEX
                        The CRISPR binary index file generated by x_genome.py
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        The output CSV generated
  --db_file DB_FILE, -d DB_FILE
                        The CRISPR DB file generated by x_genome.py
```

</details>

### Output Interpretation

BEstimate produces several files for different purposes:
- crispr_df: All gRNAs with and without editable sites 
<details>
<summary>Expand to see <strong>crispr_df</strong> column interpretation</summary>

```
columns:
    - Hugo_Symbol: Hugo Symbol of the corresponding gene location.
    - CRISPR_PAM_Sequence: gRNA target sequence with including PAM region.
    - gRNA_Target_Sequence: gRNA target sequence.
    - Location: Location of the gRNA target sequence (chromosome:start location:end location).
    - Direction: Direction of the gRNA.
    - Gene_ID: Ensembl Gene ID of the interested Hugo Symbol.
    - Transcript_ID: Ensembl Transcript ID of the interested Hugo Symbol in the corresponding gene location.
    - Exon_ID:  Ensembl Exon ID of the interested Hugo Symbol in the corresponding gene location.
    - guide_in_CDS: If gRNA has any nucleotide inside a coding sequence of the gene. 
    - gRNA_flanking_sequences: In case that user has given `-flan` option, then gRNA with the flanking sequences
    - Poly_T: Boolean output representing if the gRNA has a Poly-T region
    - GC%: GC content of the gRNA
```
</details>

- edit_df: gRNAs with editable sites within the targeted sequence

<details>
<summary>Expand to see <strong>edit_df</strong> column interpretation</summary>

```
columns:
    - Hugo_Symbol: Hugo Symbol of the corresponding gene location.
    - CRISPR_PAM_Sequence: gRNA target sequence with including PAM region.
    - gRNA_Target_Sequence: gRNA target sequence.
    - Location: Location of the gRNA target sequence (chromosome:start location:end location).
    - Edit_Location: Location of the possible edit with the given Base Editor.
    - Direction: Direction of the gRNA.
    - Strand: Strand of the interested gene on the genome (-1 or 1)
    - Gene_ID: Ensembl Gene ID for the corresponding gRNA target site region.
    - Transcript_ID: Ensembl Transcript ID for the corresponding gRNA target site region.
    - Exon_ID:  Ensembl Exon ID for the corresponding gRNA target site region.
    - guide_in_CDS: Boolean output representing if the any nucleotide on the gRNA is on the CDS region or not.
    - gRNA_flanking_sequences: In case that user has given `-flan` option, then gRNA with the flanking sequences
    - Edit_in_Exon: Boolean output representing if the specified edit in the Edit Location happening on the exon or not.
    - Edit_in_Exon: Boolean output representing if the specified edit in the Edit Location happening on the CDS or not.
    - GC%: GC content of the gRNA
    - # Edits/guide: Number of editable nuclotide within the activity window of the gRNA
    - Poly_T: Boolean output representing if CRISPR_PAM_Sequence has consecutive 4 T nucleotides.
    - mutation_on_guide: If any mutation is provided, boolean output representing if the mutation is included within the gRNA sequence. 
    - guide_change_mutation: If any mutation is provided, boolean output representing if the gRNA can make changed on specified mutation.
    - mutation_on_window: If any mutation is provided, boolean output representing if the mutation is included within the activity window of the gRNA sequence
    - mutation_on_PAM: If any mutation is provided, boolean output representing if the mutation is included within the PAM sequence of the gRNA sequence
```
</details>


- hgvs_df: If `-vep` provided, a  file with VEP API inputs for each gRNA
- vep_df: If `-vep` provided, an edit_df file with VEP API annotation for each gRNA
- protein_df: If `-vep` provided, a vep_df file with protein and structural annotation for each gRNA

<details>
<summary>Expand to see <strong>protein_df</strong> column interpretation</summary>

```
additional columns:
    - HGVS: HGVS nomenclature of the gRNA potential edits
    - Protein_ID:  Ensembl Protein ID of the potential edit.
    - VEP_input: HGVS nomenclature which was used for VEP API.
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
    - Protein_ID: Ensembl Protein ID for the corresponding gRNA target site region.
    - CDS_Position: Position on the coding sequence of the potential edit.
    - Protein_Position_ensembl: Location of the corresponding edit on the resulting protein (as Ensembl PID index).
    - Protein_Position: Location of the corresponding edit on the resulting protein (as Uniprot index).
    - Protein_Change: Protein sequence change with the corresponding edit (old amino asid - position - new amino asid).
    - Edited_AA: One letter representation of the amino asid before the corresponding edit.
    - Edited_AA_Prop: Chemical properties of the amino asid before the corresponding edit.
    - New_AA: One letter representation of the amino asid after the corresponding edit.
    - New_AA_Prop: Chemical properties of the amino asid after the corresponding edit.
    - is_Synonymous: Boolean output representing if the resulting edit causes synonymous or non-synonymous mutations.
    - is_Stop: Boolean output representing if the resulting edit causes stop codon or not.
    - proline_addition: Boolean output representing if potential edit created a Proline amino acid or not. 
    - swissprot_vep: SwissProt ID of the corresponding gRNA target site region from Ensembl VEP.
    - uniprot_provided: Uniprot ID from the user.
    - polyphen_score: Polyphen Score of the corresponding edit.
    - polyphen_prediction: Polyphen Prediction of the corresponding edit.
    - sift_score: Sift Score of the corresponding edit.
    - sift_prediction: Sift Prediction of the corresponding edit.
    - cadd_phred: CADD Prediction of the corresponding edit.
    - cadd_raw: CADD raw score of the corresponding edit.
    - lof: representing if the corresponding edit causes Loss of function with high confidence (HC)/Low confidence (LC) or not by implementing [LOFTEE](https://github.com/konradjk/loftee) through VEP.
    - impact: Impact of the corresponding edit.
    - blosum62: blosum62 Score of the corresponding edit.
    - is_clinical: representing if the corresponding edit has Clinical consequences in ClinVar or not. 
    - clinical_id: dbSNP id of the clinical allele.
    - clinical_significance: Clinical significance of the corresponding edit.
    - cosmic_id: COSMIC id of the clinical allele.
    - clinvar_id: ClinVar id of the clinical allele.
    - ancestral_populations: The conservation score from the Ensembl Compara databases
    - Domain: Protein domain in which corresponding edit happening.
    - curated_Domain: Curated name of the protein domain in which corresponding edit happening.
    - PTM: Post Translational Modification sites in which corresponding edit happening.
    - is_disruptive_interface_EXP: Boolean output representing if the corresponding edit happening on the interface region and having an experimental evidence (PDB). *High confidence*
    - is_disruptive_interface_MOD: Boolean output representing if the corresponding edit happening on the interface region and having evidence from a model (Interactome3D). *Medium confidence*
    - is_disruptive_interface_PRED: Boolean output representing if the corresponding edit happening on the interface region and having evidence from prediction (Interactome Insider). *Low confidence*
    - disrupted_PDB_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to experimental evidence - PDB. 
    - disrupted_I3D_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to Interactome3D.
    - disrupted_Eclair_int_partners: List of partner proteins whose interactions are disrupted by the corresponding edit according to Eclair - prediction algorithm of Interactome Insider. 
    - disrupted_PDB_int_genes: List of partner genes (Hugo Symbols) whose interactions are disrupted by the corresponding edit according to experimental evidence - PDB. 
    - disrupted_I3D_int_genes: List of partner genes (Hugo Symbols) whose interactions are disrupted by the corresponding edit according to Interactome3D.
    - disrupted_Eclair_int_genes: List of partner genes (Hugo Symbols) whose interactions are disrupted by the corresponding edit according to Eclair - prediction algorithm of Interactome Insider. 
```
</details>

- ot_annotated_df: If `-ot` provided, an edit_df/protein_df file with off-target annotation for each gRNA

<details>
<summary>Expand to see <strong>ot_annotated_df</strong> column interpretation</summary>

```
additional columns:
    - exact: The number of exact alignments of gRNA sequence.
    - mm1: The number of alignments of gRNA sequence with one mismatch.
    - mm2: The number of alignments of gRNA sequence with two mismatches.
    - mm3: The number of alignments of gRNA sequence with three mismatches.
    - mm4: The number of alignments of gRNA sequence with four mismatches.
```

</details>

- scored_df: If `-rs3` or `-fc` provided, an edit/protein_df/ot_annotated_df file with on-target scoring for each gRNA

<details>
<summary>Expand to see <strong>scored_df</strong> column interpretation</summary>

```
additional columns:
    - rs3_sequence: The flanking region and gRNA sequence for thr RS3 analysis.  
    - RuleSet3_Hsu2013: gRNA on-target activity score based on Hsu2013 tracRNA gRNA with Rule Set 3
    - RuleSet3_Chen2013: gRNA on-target activity score based on Chen2013 tracRNA gRNA with Rule Set 3
    - FORECasT-BE: gRNA predicted guide efficacy
```
</details>

## Contact

BEstimate is the product of Cansu Dinçer, Matthew Coelho and Mathew Garnett from Garnett Group at the Wellcome Sanger Institute. Off-target analysis has been adapted by Bo Fussing from the Cellular Informatics team within the Wellcome Sanger Institute.

If you have any problems or feedback regarding BEstimate, please contact [here](mailto:cd7@sanger.ac.uk).

## License

GNU AFFERO GENERAL PUBLIC LICENSE

BEstimate: A Python module to design and annotate base editor gRNAs

Copyright (C) 2025 Genome Research Limited

Authors: Cansu Dinçer (cd7@sanger.ac.uk), Bo Fussing (bf15@sanger.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

## Further Disclaimer
This tool is for research purposes and not for clinical use.
For policies regarding the underlying data, please also refer to:
- [Ensembl terms and conditions](https://www.ensembl.org/info/about/legal/code_licence.html#:~:text=Subject%20to%20the%20terms%20and,the%20Work%20and%20such%20Derivative)
- [Uniprot terms and conditions](https://www.uniprot.org/help/license#:~:text=We%20make%20no%20warranties%20regarding,by%20patents%20or%20other%20rights.)
- [Interactome Insider terms and conditions](http://interactomeinsider.yulab.org/)
- [Rule Set 3](https://gpp-rnd.github.io/rs3/)
- [FORECasT-BE](https://gpp-rnd.github.io/rs3/)
