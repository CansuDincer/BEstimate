# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                           Genome Retrieval and Indexing                                  #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                                                                                          #
# -----------------------------------------------------------------------------------------#


import argparse, pandas, os, subprocess, time
from crispr_analyser import index, gather
import x_index_db


# Extracting Humen Reference Genome

###########################################################################################
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate - Genome",
									 usage="%(prog)s [inputs]")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	parser.add_argument("-pamseq", dest="PAMSEQ", default="NGG",
						help="The PAM sequence in which features used "
							 "for searching activity window and editable nucleotide.")

	parser.add_argument("-assembly", dest="ASSEMBLY", default="GRCh38",
						help="The genome assembly that will be used!")

	parser.add_argument("-o", dest="OUTPUT_PATH", default=os.getcwd() + "/",
						help="The path for output. If not specified the current directory will be used!")

	parser.add_argument("-v_ensembl", dest="VERSION", default="113",
						help="The ensembl version in which genome will be retrieved "
							 "(if the assembly is GRCh37 then please use <=75)")

	parser.add_argument("-ot_path", dest="OT_PATH", default=os.getcwd() + "/../offtargets")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


###########################################################################################
# Functions

def check_genome_exist(assembly, ens_ver):
	global ot_path
	chromosomes = list(range(1, 23)) + ["X", "Y", "MT"]

	if assembly == "GRCh37":
		file_main_text = "Homo_sapiens.GRCh37.%s.dna.chromosome" % ens_ver
	elif assembly == "GRCh38":
		file_main_text = "Homo_sapiens.GRCh38.dna.chromosome"

	check_files = True
	for chromosome in chromosomes:
		if "%s.%s.fa.gz" % (file_main_text, chromosome) not in os.listdir(ot_path + "genome/"):
			check_files = False

	if check_files is False:

		print(
			"Genome is not found, BEstimate is downloading the %s Ensembl genome - version %s\n" % (assembly, ens_ver))

		if "chromosome_ftps.txt" not in os.listdir("%sgenome/" % ot_path):
			f = open("%sgenome/chromosome_ftps.txt" % ot_path, "w")
			for chromosome in chromosomes:
				f.writelines(
					"url=https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/dna/%s.%s.fa.gz\n" % (
						ens_ver, file_main_text, chromosome))
				f.writelines(
					"output=%sgenome/%s.%s.fa.gz\n" % (ot_path, file_main_text, chromosome))
			f.close()

		curl_command = "curl --parallel --parallel-immediate --parallel-max 25 --fail-with-body --retry 5 " \
					   "--config %sgenome/chromosome_ftps.txt -C -" % ot_path
		print("Collecting the genome files from Ensembl FTP..\n")

		_ = subprocess.Popen(curl_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							 text=True, shell=True)

		check_files = True
		for chromosome in chromosomes:
			if "%s.%s.fa.gz" % (file_main_text, chromosome) not in os.listdir(ot_path + "genome/"):
				check_files = False

		error_message = "Error in downloading genome, please manually downloading all chromosomes from:" \
						"https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/dna/ named as " \
						"Homo_sapiens.GRCh38.dna.chromosome.<chromosome>.fa.gz if the assembl is GRCh38, otherwise" \
						"Homo_sapiens.GRCh37.<version>.dna.chromosome.<chromosome>.fa.gz"

		if check_files:
			return True
		else:
			return error_message
	else:
		return True


def index_genome(assembly, ens_ver, pam_sequence):
	global ot_path
	chromosomes = list(range(1, 23)) + ["X", "Y", "MT"]

	if assembly == "GRCh37":
		file_main_text = "Homo_sapiens.GRCh37.%s.dna.chromosome" % ens_ver
	elif assembly == "GRCh38":
		file_main_text = "Homo_sapiens.GRCh38.dna.chromosome"

	if "%s.bin" % file_main_text not in os.listdir("%sgenome/" % ot_path):
		chromosome_input_text_list = list()
		for chromosome in chromosomes:
			chromosome_input_text_list.append("%sgenome/csv/c_%s.csv" % (ot_path, chromosome))


		# Gather all chromosome fasta files into csv files
		print("Gathering chromosomes..\n")
		for chromosome in chromosomes:
			file_name = "%s.%s.fa" % (file_main_text, chromosome)
			if "c1_%s.csv" % chromosome not in os.listdir("%sgenome/csv/" % ot_path):
				if file_name not in os.listdir("%sgenome/" % ot_path):
					os.system("gunzip --keep %sgenome/%s.gz" % (ot_path, file_name))

				print("Chromosome %s" % chromosome)
				gather.gather(inputfile=f"{ot_path}genome/{file_name}",
					outputfile=f"{ot_path}genome/csv/c_{chromosome}.csv", pam=pam_sequence)

				while "c_%s.csv" % chromosome not in os.listdir("%sgenome/csv/" % ot_path):
					print("Waiting chromosome %s.." % chromosome)
					time.sleep(10)
		print()
		# index the database with CRISPRs gathered in the CSV files
		x_index_db.index(chromosome_input_text_list)

		index.index(inputfiles=chromosome_input_text_list,
					outputfile=f"{ot_path}genome/{file_main_text}.bin",
					species="Human", assembly=ens_ver, offset=0, species_id="1")

	if "%s.bin" % file_main_text in os.listdir("%sgenome/" % ot_path):
		return True
	else:
		return False


def check_index_file(assembly, ens_ver):
	global ot_path

	if assembly == "GRCh37":
		file_main_text = "Homo_sapiens.GRCh37.%s.dna.chromosome" % ens_ver
	elif assembly == "GRCh38":
		file_main_text = "Homo_sapiens.GRCh38.dna.chromosome"

	if "%s.bin" % file_main_text not in os.listdir("%sgenome/" % ot_path):
		return False

	else:
		if "crisprs.db" in os.listdir(os.getcwd()):
			return True
		else:
			return False


###########################################################################################
# Execution


def get_genome():
	"""
	Run whole script with the input from terminal
	:return:
	"""

	try:
		os.mkdir(ot_path)
	except FileExistsError:
		pass

	try:
		os.mkdir(ot_path + "genome/")
	except FileExistsError:
		pass


	is_index = check_index_file(assembly=args["ASSEMBLY"], ens_ver=args["VERSION"])

	if is_index:
		is_genome = True
	else:
		is_genome = check_genome_exist(assembly=args["ASSEMBLY"], ens_ver=args["VERSION"])

	if is_genome:
		return True
	else:
		print("Error: Please download the Humen Reference Genome from Ensembl before continue!")
		return False


def get_index():
	try:
		os.mkdir(ot_path + "genome/csv/")
	except FileExistsError:
		pass

	is_index = check_index_file(assembly=args["ASSEMBLY"], ens_ver=args["VERSION"])

	if is_index:
		return True
	else:
		print("Indexing the genome..")
		res = index_genome(assembly=args["ASSEMBLY"], ens_ver=args["VERSION"], pam_sequence=args["PAMSEQ"])
		return res


if __name__ == '__main__':

	# -----------------------------------------------------------------------------------------#
	# Retrieve input

	args = take_input()
	# Output Path
	path = ""
	if args["OUTPUT_PATH"][-1] == "/":
		path = args["OUTPUT_PATH"]
	else:
		path = args["OUTPUT_PATH"] + "/"

	ot_path = ""
	if args["OT_PATH"][-1] == "/":
		ot_path = args["OT_PATH"]
	else:
		ot_path = args["OT_PATH"] + "/"

	g = get_genome()
	if g:
		get_index()

	else:
		print("no index")
