# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                          CRISPR-Analyser | Genome Indexing                               #
#                        Author : Bo Fussing bf15@sanger.ac.uk                             #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

import collections, csv, getopt, re, sys, time

###########################################################################################
# Functions


def reverse_complement(sequence: str) -> str:
	"""Return the reverse complement of a DNA sequence.

	Args:
		sequence: The string DNA sequence to reverse complement.
	"""
	complement_map = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
	return "".join([complement_map[base] for base in reversed(sequence)])


def match_pam(
	dna_sequence: str,
	pam_sequence: str,
	pam_on_right: bool,
	legacy_mode: bool = False,
) -> bool:
	"""Check if the DNA sequence has a PAM sequence match.

	Args:
		dna_sequence: The string DNA sequence to check.
		pam_sequence: The string PAM sequence to match.
		pam_on_right: A boolean indicating if PAM sequence is on the right.
		legacy_mode: A boolean indicating if non-ACGT chars allowed in PAM
			region of the DNA sequence. Default is False.
	"""
	start = len(dna_sequence) - len(pam_sequence) if pam_on_right else 0
	for i in range(len(pam_sequence)):
		if legacy_mode is False and dna_sequence[start + i] not in "ACGT":
			return False
		if pam_sequence[i] == "N":
			continue
		if dna_sequence[start + i] != pam_sequence[i]:
			return False
	return True


def gather(inputfile: str, outputfile: str, pam: str, verbose: bool = False):
	"""Run the CRISPR gatherer.

	Args:
		inputfile: The input FASTA file.
		outputfile: The output CSV file generated.
		pam: The string PAM sequence to search for.
		verbose: A boolean indicating if verbose output is enabled.
			Default is False.
	"""
	if verbose:
		start = time.time()
	chromosome = ""
	crispr_count = 0
	position = 0
	buffer = collections.deque(maxlen=len(pam) + 20)

	with open(inputfile, "r") as infile:
		with open(outputfile, "w", newline="") as outfile:
			csvwriter = csv.writer(outfile)
			for line in infile:
				if len(line) == 0:
					continue
				if line[0] == ">":
					chromosome = re.search(
						r">(.*?) dna:chromosome", line
					).group(1)
					if verbose:
						print(f"Processing chromosome {chromosome}...")
					position = 0
					buffer.clear()
					continue
				else:
					line = line.strip()
					for base in line:
						buffer.append(base)
						if len(buffer) < len(pam) + 20:
							continue
						position += 1
						if match_pam(
							dna_sequence=buffer,
							pam_sequence=reverse_complement(pam),
							pam_on_right=False,
						):
							csvwriter.writerow(
								[chromosome, position, "".join(buffer), 0, 1]
							)
							crispr_count += 1
						if match_pam(
							dna_sequence=buffer,
							pam_sequence=pam,
							pam_on_right=True,
						):
							csvwriter.writerow(
								[chromosome, position, "".join(buffer), 1, 1]
							)
							crispr_count += 1
	if verbose:
		end = time.time()
		print(f"Gathered {crispr_count} CRISPRs in {end - start} seconds.")


def run(argv=sys.argv[1:]):
	"""Run the CRISPR gatherer."""
	inputfile = ""
	outputfile = ""
	pam = ""

	def usage():
		print("gather -i <inputfile> -o <CSV outputfile> -p <pam>")

	try:
		opts, args = getopt.getopt(
			argv, "hi:o:p:", ["ifile=", "ofile=", "pam="]
		)
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			usage()
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-p", "--pam"):
			pam = arg
	if inputfile == "" or outputfile == "" or pam == "":
		usage()
		sys.exit(2)

	gather(inputfile, outputfile, pam, verbose=True)


if __name__ == "__main__":
	run()
