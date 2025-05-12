# BEstimate

BEstimate, a Python module that systematically identifies guide RNA (gRNA) targetable sites across given sequences for given Base Editors, functional and clinical effects of the potential edits on the resulting proteins and off target consequence of the found sequences. It has the ability to provide in silico analysis of the sequences to identify positions that can be editable by Base Editors, and their features before starting experiments. 

You can directly use BEstimate environment if you have conda. Please follow below:

- `git clone https://github.com/CansuDincer/BEstimate.git`
- `cd BEstimate`
- `conda-env create -n bestimate -f=bestimate.yml`
- `conda activate bestimate`

If not, you should have python 3.13 and you can use requirements file:

- `pip3 install -r requirements.txt`

## Run BEstimate

`python3 BEstimate.py -gene BRCA1 -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -o ../output/ -ofile BRCA1_CBE_NGG`

The user also run the same analysis for different PAM only changing -pamseq NGN. 

*Warning: Be careful to write the PAM sequence to be in concordant with the length of the -pamwin. Here, NGN is in concordant with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.* 

If you would like to run for a specific transcript and run the protein analysis:

`python3 BEstimate.py -gene BRAF -assembly GRCh38 -transcript ENST00000646891 -edit C -edit_to T -vep -o ../output/ -ofile BRAF_CBE_NGG`

If you would like to run with a specific point mutation, with NGN PAM and with VEP and protein analysis:
Prepare a `PIK3CA_mutation_file.txt` for example with 3:g.179218303G>A

`python3 BEstimate.py -gene PIK3CA -assembly GRCh38 -pamseq NGN -pamwin 21-23 -actwin 4-8 -protolen 20 -mutation_file PIK3CA_mutation_file.txt -edit A -edit_to G -vep -ofile PIK3CA_NGN_ABE_mE545K -o ../output/`

To run off target analysis, first you need to have *Ensembl* Genome and its indexes for the interested PAM sequence. `x_genome.py` has been prepared for the user to download and index the genome.

What you will need:
- **pamseq** > PAM sequence as the genome will index accordingly
- **assembly** > The Ensembl genome assembly
- **v_ensembl** > Ensembl version (currently default is 113 for GRCh38, if the assembly is GRCh37 then please use <=75)
- **ot_path** > Path for the downloaded genome for the off-target analysis  

`python3 x_genome.py -pamseq NGN -assembly GRCh38 -v_ensembl 113 -ot_path ../`

Then, you can run the off target analysis, see below for *BRAF* gene:

`python3 BEstimate.py -gene BRAF -assembly GRCh38 -pamseq NGN -edit A -edit_to G -vep -ot -o ../output/ -ofile BRAF_ABE_NGN`

## Contact

BEstimate is the product of Cansu Dinçer, Matthew Coelho and Mathew Garnett from Garnett Group at the Wellcome Sanger Institute. Off-target analysis has been adapted by Bo Fussing from the Cellular Informatics team within the Wellcome Sanger Institute.

For any problems or feedback on BEstimate, you can contact [here](mailto:cd7@sanger.ac.uk).

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


