# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas, argparse, re, pickle
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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

def find_pam_protospacer(sequence, pam_sequence, pam_window, protospacer_length, strand, chromosome):
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
				if strand == "1":
					genomic_start = nuc_index + 1
					genomic_end = (genomic_start + pam_window[1]) - 1
				elif strand == "-1":
					genomic_end = len(sequence) - nuc_index
					genomic_start = (genomic_end - pam_window[1]) + 1
				genomic_location = str(chromosome) + ":" + str(genomic_start) + "-" + str(genomic_end)

				crisprs.append({"crispr": match_sequence.group(),
								"grna": match_sequence.group()[:-len(pam_sequence)],
								"location": genomic_location, "strand": strand})

	crispr_df = pandas.DataFrame(crisprs)

	return crispr_df


def get_genomic_crisprs(pam_sequence, pam_window, protospacer_length, chromosome):
	if "guides_%s_%s.csv" % (chromosome, pam_sequence) not in os.listdir(output_path + "offtargets/"):
		nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

		# Read fasta file
		file = input_path + "chromosomes/dna_chromosome_" + chromosome + ".fa"
		fasta_file = SeqIO.read(file, "fasta")

		# Extract chromosome sequence
		sequence = str(fasta_file.seq)

		# Reverse sequence for complementary DNA strand
		reverse_sequence = "".join([nucleotide_dict[n] for n in sequence[::-1]])

		# PAM specific CRISPR location finding across chromosome
		positive_crispr_df = find_pam_protospacer(
			sequence=sequence, pam_sequence=pam_sequence,
			pam_window=[pam_window.split("-")[0], pam_window.split("-")[1]],
			protospacer_length=protospacer_length, strand="1", chromosome=chromosome)

		negative_crispr_df = find_pam_protospacer(
			sequence=reverse_sequence, pam_sequence=pam_sequence,
			pam_window=[pam_window.split("-")[0], pam_window.split("-")[1]],
			protospacer_length=protospacer_length, strand="-1", chromosome=chromosome)

		df = pandas.concat([positive_crispr_df, negative_crispr_df])
		df.columns = ["CRISPR_PAM_Sequence", "gRNA_Target_Sequence", "Location", "Strand"]
		df.to_csv(output_path + "offtargets/guides_%s_%s.csv" % (chromosome, pam_sequence), index=False)
	else:
		df = pandas.read_csv(output_path + "offtargets/guides_%s_%s.csv" % (chromosome, pam_sequence), index_col=0)

	return df


# -----------------------------------------------------------------------------------------#
# Object

class Guide:

	def __init__(self, guide):
		self.guide = guide
		self.guide_id = None
		self.locations = None
		self.one_mms = None
		self.two_mms = None
		self.three_mms = None
		self.four_mms = None

	def get_id(self, guide_index):
		self.guide_id = guide_index

	def get_locations(self, loc):
		if self.locations is None:
			self.locations = list()
			self.locations.append(loc)

		else:
			if loc not in self.locations:
				self.locations.append(loc)

	def mm_calculations(self, other_guide):
		# High gap opening and gap extension points to prevent indels
		# Exact match should equal to 40
		for alignm in pairwise2.align.globalms(self.guide, other_guide, 2, 0, -5, -5):
			if alignm[-2] == 0 and alignm[-1] == len(self.guide):
				score = alignm[2]
				if score >= 32:
					mm = (40 - score) / 2
					if mm == 1:
						if self.one_mms is None:
							self.one_mms = list()
							self.one_mms.append(other_guide)
						else:
							if other_guide not in self.one_mms:
								self.one_mms.append(other_guide)
					elif mm == 2:
						if self.two_mms is None:
							self.two_mms = list()
							self.two_mms.append(other_guide)
						else:
							if other_guide not in self.two_mms:
								self.two_mms.append(other_guide)
					elif mm == 3:
						if self.three_mms is None:
							self.three_mms = list()
							self.three_mms.append(other_guide)
						else:
							if other_guide not in self.three_mms:
								self.three_mms.append(other_guide)
					elif mm == 4:
						if self.four_mms is None:
							self.four_mms = list()
							self.four_mms.append(other_guide)
						else:
							if other_guide not in self.four_mms:
								self.four_mms.append(other_guide)


def serialise_object(genome_guide_df, pam):
	if pam == "NGG":
		genome_guide_df = genome_guide_df[genome_guide_df["CRISPR_PAM_Sequence"].str.endswith("GG")]

	guide_count = 0
	x = genome_guide_df.groupby(["gRNA_Target_Sequence"])
	t = len(x)
	for guide, guide_df in x:
		guide_count += 1
		if "offtargets/genome_%s_%s_object.p" % (pam, guide) not in os.listdir(output_path + "offtargets/guides/"):
			obj = Guide(guide=guide)
			obj.get_id(guide_index=guide_count)
			for group, _ in guide_df.groupby(["Location", "Strand"]):
				loc = group[0] + "|" + str(group[1])
				obj.get_locations(loc)

			for other_guide in genome_guide_df.gRNA_Target_Sequence.unique():
				if other_guide != guide:
					obj.mm_calculations(other_guide=other_guide)
			pickle.dump(guide, open(output_path + "offtargets/guides/genome_%s_%s_object.p" % (pam, guide), "wb"))
		print(guide_count * 100.0 / t)

def deserialise_viability_object(pam):
	guide_objs = pickle.load(open(output_path + "offtargets/genome_%s_guide_object.p" % pam, "rb"))
	return guide_objs


def summary_off_targets(pam):
	if "off_target_guide_summary_%s.csv" % pam not in os.listdir(output_path + "offtargets/"):
		guide_objs = deserialise_viability_object(pam)

		guide_df = pandas.DataFrame(columns=["guide_id", "full_match", "1mm", "2mm", "3mm", "4mm"],
									index=list(guide_objs.keys()))
		for guide, obj in guide_objs.items():
			guide_df.loc[guide, "guide_id"] = obj.guide_id
			guide_df.loc[guide, "full_match"] = len(obj.locations)
			guide_df.loc[guide, "1mm"] = len(obj.one_mms)
			guide_df.loc[guide, "2mm"] = len(obj.two_mms)
			guide_df.loc[guide, "3mm"] = len(obj.three_mms)
			guide_df.loc[guide, "4mm"] = len(obj.four_mms)

		guide_df.to_csv(output_path + "offtargets/off_target_guide_summary_%s.csv" % pam, index=True)

	else:
		guide_df = pandas.read_csv(output_path + "offtargets/off_target_guide_summary_%s.csv" % pam, index_col=0)

	return guide_df


###########################################################################################
# Execution


def main():
	if "genome_wide_%s.csv" % args["PAMSEQ"] not in os.listdir(output_path + "offtargets/"):
		all_dfs = list()
		for chromosome in [str(x) for x in list(range(1, 23))] + ["X", "Y"]:
			df = get_genomic_crisprs(pam_sequence=args["PAMSEQ"], pam_window=args["PAMWINDOW"],
									 protospacer_length=args["PROTOLEN"], chromosome=chromosome)
			all_dfs.append(df)
		genome_guides = pandas.concat(all_dfs)
		genome_guides.to_csv(output_path + "offtargets/genome_wide_%s.csv" % args["PAMSEQ"], index=True)
	else:
		print("Reading Guide File")
		genome_guides = pandas.read_csv(output_path + "offtargets/genome_wide_%s.csv" % args["PAMSEQ"], index_col=0)
		genome_guides.astype({col: int for col in desired_columns_to_convert})

	if "genome_%s_guide_object.p" % args["PAMSEQ"] not in os.listdir(output_path + "offtargets/"):
		serialise_object(genome_guide_df=genome_guides, pam=args["PAMSEQ"])

	summary_guide_df = summary_off_targets(pam=args["PAMSEQ"])

	return True

main()