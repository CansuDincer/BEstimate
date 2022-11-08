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
from BEstimate.BEstimate.BEstimate import find_pam_protospacer, add_genomic_location
from Bio import SeqIO


# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate",
									 usage="%(prog)s [inputs]",
									 description="""
                                     **********************************
                                     Find and Analyse Base Editor sites
                                     		 Off Target Analysis
                                     **********************************""")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# PAM AND PROTOSPACER INFORMATION

	# The NGG PAM will be used unless otherwise specified.
	parser.add_argument("-pamseq", dest="PAMSEQ", default="NGG",
						help="The PAM sequence in which features used "
							 "for searching activity window and editable nucleotide.")
	parser.add_argument("-pamwin", dest="PAMWINDOW", default="21-23",
						help="The index of the PAM sequence when starting "
							 "from the first index of protospacer as 1.")
	parser.add_argument("-protolen", dest="PROTOLEN", default="20",
						help="The total protospacer and PAM length.")


	# OUTPUT

	parser.add_argument("-o", dest="OUTPUT_PATH", default=os.getcwd() + "/",
						help="The path for output. If not specified the current directory will be used!")

	parser.add_argument("-ofile", dest="OUTPUT_FILE", default="output",
						help="The output file name, if not specified \"position\" will be used!")
	parser.add_argument("-ifile", dest="INPUT_FILE", default="output",
						help="The input file name, if not specified \"position\" will be used!")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


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

	chromosome_crisprs.to_csv(output_path + "_genome_guides.csv", index=False)

	return chromosome_crisprs



###########################################################################################
# Execution


def main():
	print("Starting..")
	_ = get_genomic_crisprs(pam_sequence= args["PAMSEQ"], pam_window=args["PAMWINDOW"],
							protospacer_length= args["PROTOLEN"])

	return "Genomic guides were collected for GRCh38 genome assembly - %s PAM" % args["PAMSEQ"]


main()