# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas, argparse, re
from Bio import SeqIO


# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate Off targets",
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
	parser.add_argument("-i", dest="INPUT_PATH", default=os.getcwd() + "/",
						help="The path for input. If not specified the current directory will be used!")
	parser.add_argument("-c", dest="CHROMOSOME", default=os.getcwd() + "/",
						help="The path for input. If not specified the current directory will be used!")

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
# Functions


def find_pam_protospacer(sequence, pam_sequence, pam_window, protospacer_length):
	"""
	Finding all possible PAM and protospacer regions on the sequence of the gene.
	:param sequence: The sequence of the interested gene.
	:param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
	:param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
	:return crisprs: A list of dictionary having sequences and locations of the gRNA targeted
	gene parts. The indices are indicated when the first index of the gene is 1.
	"""
	# Since python index starts from 0, decrease the start position index given from the user
	pam_window = [int(pam_window[0]) - 1, int(pam_window[1])]

	# Using Regular Expressions, specify PAM pattern
	pam_pattern = ""
	for nuc in list(pam_sequence):
		if nuc != "N":
			pam_pattern += nuc + "{1}"
		else:
			pam_pattern += "[ATCG]{1}"

	# Search protospacer length of nucleatide sequence and add PAM pattern after that
	pattern = r'[ATCG]{%s}%s' % (str(protospacer_length), pam_pattern)

	crisprs = []
	for nuc_index in range(len(sequence)):
		# One by one in the given sequence
		if nuc_index + pam_window[1] <= len(sequence):
			# Add started nucleotide index total length of targeted base editing site (PAM index)
			sub_sequence = sequence[nuc_index:nuc_index + pam_window[1]]

			# Search regex pattern inside the sub sequence
			for match_sequence in re.finditer(pattern, sub_sequence):
				crisprs.append({"index": [nuc_index, nuc_index + pam_window[1]],
								"crispr": match_sequence.group()})

	return crisprs

def add_genomic_location(sequence_range, crispr_dict, crispr_direction, strand, chromosome):
	"""
	Adding genomic location info on crisprs found by extract_activity_window() function.
	:param sequence_range: The range of the sequence on the genome (from Ensembl)
	:param crispr_dict: The CRISPR dictionary created by extract_activity_window()
	:param crispr_direction: The direction of the created CRISPR (left or right)
	:param strand: The strand of the given gene (-1 or 1)
	:return crispr_seq: The sequence of the CRISPR
	:return genomic_location: The genomic coordinate of the above CRISPR on genome
	"""
	# Prepare the sequence range according to the strand of the gene

	new_range = []
	if strand == 1:
		new_range = [int(sequence_range[0]), int(sequence_range[1])]
	elif strand == -1:
		new_range = [int(sequence_range[1]), int(sequence_range[0])]

	# Look for both direction
	crispr_seq, genomic_location = crispr_dict["crispr"], ""

	if crispr_direction == "right":
		genomic_start = new_range[0] + crispr_dict["index"][0]
		genomic_end = (genomic_start + len(crispr_seq)) - 1
		genomic_location = str(chromosome) + ":" + str(genomic_start) + "-" + str(genomic_end)

	elif crispr_direction == "left":
		genomic_start = new_range[1] - crispr_dict["index"][0]
		genomic_end = (genomic_start - len(crispr_seq)) + 1
		genomic_location = str(chromosome) + ":" + str(genomic_end) + "-" + str(genomic_start)

	return crispr_seq, genomic_location


def get_genomic_crisprs(pam_sequence, pam_window, protospacer_length, chromosome):

	nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

	chromosomes = [str(x) for x in list(range(1, 23))] + ["X", "Y"]

	chromosome_crisprs = pandas.DataFrame(columns = ["Chromosome", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
													 "Location", "Direction", "Strand"])

	crisprs_df = pandas.DataFrame(columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence", "Location", "Direction"])

	# Read fasta file
	file = input_path + "dna_chromosome_" + chromosome + ".fa"
	fasta_file = SeqIO.read(file, "fasta")

	# Extract chromosome sequence
	sequence = str(fasta_file.seq)

	# Reverse sequence for complementary DNA strand
	reverse_sequence = "".join([nucleotide_dict[n] for n in sequence[::-1]])

	# Sequence range
	start = fasta_file.description.split(" ")[2].split(":")[-3]
	end = fasta_file.description.split(" ")[2].split(":")[-2]
	chr_range = [int(start), int(end)]
	# PAM specific CRISPR location finding across chromosome
	right_chr_crispr = find_pam_protospacer(
		sequence=sequence, pam_sequence=pam_sequence,
		pam_window=[pam_window.split("-")[0], pam_window.split("-")[1]],
		protospacer_length=protospacer_length)

	left_chr_crispr = find_pam_protospacer(
		sequence=reverse_sequence, pam_sequence=pam_sequence,
		pam_window=[pam_window.split("-")[0], pam_window.split("-")[1]],
		protospacer_length=protospacer_length)

	crisprs_dict = {"left": right_chr_crispr, "right": left_chr_crispr}
	for direction, crispr in crisprs_dict.items():
		for cr in crispr:
			crispr_seq_forw, genomic_location = add_genomic_location(sequence_range=chr_range, strand=1, chromosome = chr,
																	 crispr_dict=cr, crispr_direction=direction)
			crispr_seq_backw, genomic_location = add_genomic_location(sequence_range=chr_range, strand=-1, chromosome = chr,
																	  crispr_dict=cr, crispr_direction=direction)

			forw_df = pandas.DataFrame([[crispr_seq_forw, crispr_seq_forw[:-len(pam_sequence)],
										 genomic_location, direction, 1]],
									   columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
												"Location", "Direction", "Strand"])
			crisprs_df = pandas.concat([crisprs_df, forw_df])
			backw_df = pandas.DataFrame([[crispr_seq_backw, crispr_seq_backw[:-len(pam_sequence)],
										  genomic_location, direction, -1]],
										columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
												 "Location", "Direction", "Strand"])
			crisprs_df = pandas.concat([crisprs_df, backw_df])

	crisprs_df["Chromosome"] = chromosome
	chromosome_crisprs = pandas.concat([chromosome_crisprs, crisprs_df])

	chromosome_crisprs.to_csv(output_path + args["PAMSEQ"] + "_" + chromosome + "_genome_guides.csv", index=False)

	return chromosome_crisprs



###########################################################################################
# Execution


def main():
	print("Starting..")
	_ = get_genomic_crisprs(pam_sequence= args["PAMSEQ"], pam_window=args["PAMWINDOW"],
							protospacer_length= args["PROTOLEN"], chromosome=args["CHROMOSOME"])

	return "Genomic guides of chromosome %s were collected for GRCh38 genome assembly - %s PAM" \
		   % (args["CHROMOSOME"], args["PAMSEQ"])


main()