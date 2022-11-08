# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas, argparse
from BEstimate.BEstimate.BEstimate import find_pam_protospacer, add_genomic_location, take_input
from Bio import SeqIO


# -----------------------------------------------------------------------------------------#
# Take inputs

args = take_input()

output_path = ""
if args["OUTPUT_PATH"][-1] == "/":
	output_path = args["OUTPUT_PATH"]
else:
	output_path = args["OUTPUT_PATH"] + "/"

input_path = ""
if args["INPUT_PATH"][-1] == "/":
	input_path = args["INPUT_PATH"]
else:
	input_path = args["INPUT_PATH"] + "/"

# -----------------------------------------------------------------------------------------#

def get_genomic_crisprs(pam_sequence, pam_window, protospacer_length):

	nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

	chromosomes = [str(x) for x in list(range(1, 23))] + ["X", "Y"]

	chromosome_crisprs = pandas.DataFrame(columns = ["Chromosome", "CRISPR_PAM_Sequence",
													 "gRNA_Target_Sequence", "Location", "Direction"])

	for chr in chromosomes:
		if os.path.exists(input_path + "chromosomes/dna_chromosome_" + chr + ".fa"):
			crisprs_df = pandas.DataFrame(
				columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence", "Location", "Direction"])

			# Read fasta file
			file = input_path + "chromosomes/dna_chromosome_" + chr + ".fa"
			fasta_file = SeqIO.read(file, "fasta")

			# Extract chromosome sequence
			sequence = str(fasta_file.seq)

			# Reverse sequence for complementary DNA strand
			reverse_sequence = "".join([nucleotide_dict[n] for n in sequence[::-1]])

			# Sequence range
			start = fasta_file.description.split(" ")[2].split(":")[-2]
			end = fasta_file.description.split(" ")[2].split(":")[-3]
			chr_range = [int(start), int(end)]

			# PAM specific CRISPR location finding across chromosome
			right_chr_crispr = find_pam_protospacer(
				sequence=sequence, pam_sequence=pam_sequence,
				pam_window=pam_window, protospacer_length=protospacer_length)

			left_chr_crispr = find_pam_protospacer(
				sequence=reverse_sequence, pam_sequence=pam_sequence,
				pam_window=pam_window, protospacer_length=protospacer_length)

			crisprs_dict = {"left": right_chr_crispr, "right": left_chr_crispr}

			for direction, crispr in crisprs_dict.items():
				for cr in crispr:
					crispr_seq, genomic_location = add_genomic_location(sequence_range=chr_range, strand=1,
																		crispr_dict=cr, crispr_direction=direction)

					df = pandas.DataFrame([[crispr_seq, crispr_seq[:-len(pam_sequence)],
											chr + ":" + genomic_location, direction]],
										  columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
												   "Location", "Direction"])
					crisprs_df = pandas.concat([crisprs_df, df])

			crisprs_df["Chromosome"] = chr
			chromosome_crisprs = pandas.concat([chromosome_crisprs, crisprs_df])

	chromosome_crisprs.to_csv(output_path + args["PAMSEQ"] + "_genome_guides.csv", index=False)

	return chromosome_crisprs



###########################################################################################
# Execution


def main():
	print("Starting..")
	_ = get_genomic_crisprs(pam_sequence= args["PAMSEQ"], pam_window=args["PAMWINDOW"],
							protospacer_length= args["PROTOLEN"])

	return "Genomic guides were collected for GRCh38 genome assembly - %s PAM" % args["PAMSEQ"]


main()