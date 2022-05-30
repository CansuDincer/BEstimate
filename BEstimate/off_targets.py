# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas, re, argparse
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

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


args = take_input()

path = ""
if args["OUTPUT_PATH"][-1] == "/":
	path = args["OUTPUT_PATH"]
else:
	path = args["OUTPUT_PATH"] + "/"

# -----------------------------------------------------------------------------------------#
# Functions


def find_pam_protospacer_genome(sequence, pam_sequence, pam_window, protospacer_length):
	"""
	Finding all possible PAM and protospacer regions on the genome
	:param sequence: The sequence of the interested gene.
	:param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
	:param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
	:param protospacer_length: The length of protospacer.
	:return crisprs: A list of dictionary having sequences and locations of the gRNA targeted
	gene parts. The indices are indicated when the first index of the gene is 1.
	"""
	# Since python index starts from 0, decrease the start position index given from the user
	pam_window = [int(pam_window.split("-")[0]), int(pam_window.split("-")[1])]

	# Using Regular Expressions, specify PAM pattern
	pam_pattern = ""
	for nuc in list(pam_sequence):
		if nuc != "N":
			pam_pattern += nuc + "{1}"
		else:
			pam_pattern += "[ATCG]{1}"

	# Search protospacer length of nucleatide sequence and add PAM pattern after that
	pattern = r'[ATCG]{%s}%s' % (protospacer_length, pam_pattern)
	print("Pattern is created!")

	crisprs = []
	print("Pattern is searching through the sequence...")
	for nuc_index in range(0, len(sequence)):

		# One by one in the given sequence
		if nuc_index + pam_window[1] <= len(sequence):

			# Add started nucleotide index total length of targeted base editing site (PAM index)
			sub_sequence = sequence[nuc_index:nuc_index + pam_window[1]]

			# Search regex pattern inside the sub sequence
			for match_sequence in re.finditer(pattern, sub_sequence):
				crisprs.append({"index": [nuc_index, nuc_index + pam_window[1]],
								"crispr": match_sequence.group()})

	if crisprs is not []: print("CRISPRs are found!")
	return crisprs


def add_genomic_location(sequence_range, crispr_dict, crispr_direction, strand):
	"""
	Adding genomic location info on crisprs
	:param sequence_range: The range of the sequence on the genome (from Ensembl)
	:param crispr_dict: The CRISPR dictionary created by extract_activity_window()
	:param crispr_direction: The direction of the created CRISPR (left or right)
	:param strand: The strand of the given gene (-1 or 1)
	:return crispr_seq: The sequence of the CRISPR
	:return genomic_location: The genomic coordinate of the above CRISPR on genome
	"""

	new_range = []
	if strand == 1:
		new_range = sequence_range
	elif strand == -1:
		new_range = [sequence_range[1], sequence_range[0]]

	# Look for both direction
	crispr_seq, genomic_location = crispr_dict["crispr"], ""

	if crispr_direction == "right":
		genomic_start = new_range[0] + crispr_dict["index"][0]
		genomic_end = (genomic_start + len(crispr_seq)) - 1
		genomic_location = str(genomic_start) + "-" + str(genomic_end)

	elif crispr_direction == "left":
		genomic_start = new_range[1] - crispr_dict["index"][0]
		genomic_end = (genomic_start - len(crispr_seq)) + 1
		genomic_location = str(genomic_end) + "-" + str(genomic_start)

	return crispr_seq, genomic_location


def get_genomic_crisprs(pam_sequence, pam_window, protospacer_length):

	nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

	chromosomes = [str(x) for x in list(range(1, 23))] + ["X", "Y"]

	chromosome_crisprs = pandas.DataFrame(columns = ["Chromosome", "CRISPR_PAM_Sequence",
													 "gRNA_Target_Sequence", "Location", "Direction"])

	for chr in chromosomes:
		if os.path.exists(os.getcwd() + "/../../extra_data/dna_chromosome_" + chr + ".fa"):
			crisprs_df = pandas.DataFrame(
				columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence", "Location", "Direction"])

			# Read fasta file
			file = os.getcwd() + "/../../extra_data/dna_chromosome_" + chr + ".fa"
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
			right_chr_crispr = find_pam_protospacer_genome(
				sequence=sequence, pam_sequence=pam_sequence,
				pam_window=pam_window, protospacer_length=protospacer_length)

			left_chr_crispr = find_pam_protospacer_genome(
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

	chromosome_crisprs.to_csv(path + args["OUTPUT_FILE"] + "_genome_guides.csv", index=False)

	return chromosome_crisprs


###########################################################################################
# Execution


def main():
	print("Starting..")
	_ = get_genomic_crisprs(pam_sequence= args["PAMSEQ"], pam_window=args["PAMWINDOW"],
							protospacer_length= args["PROTOLEN"])

	return "Genomic guides were collected for GRCh38 genome assembly - %s PAM" % args["PAMSEQ"]


main()