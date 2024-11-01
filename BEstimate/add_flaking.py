# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                 BEstimate Additional                                     #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, sys, pandas, argparse, requests
import time
import itertools

# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate Additional",
									 usage="%(prog)s [inputs]")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# BASIC INFORMATION

	parser.add_argument("-assembly", dest="ASSEMBLY", required=True,
						help="The genome assembly that will be used!")

	# gRNA FLANKING REGIONS

	parser.add_argument("-flank3", dest="FLAN_3", default="7",
						help="The number of nucleotides in the 3' flanking region")
	parser.add_argument("-flank5", dest="FLAN_5", default="11",
						help="The number of nucleotides in the 5' flanking region")

	# PATH

	parser.add_argument("-path", dest="PATH", default=os.getcwd() + "/",
						help="The path, if not specified the current directory will be used!")

	parser.add_argument("-file", dest="FILE", default="output",
						help="The file name")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict



# -----------------------------------------------------------------------------------------#
# Functions

def extract_gRNA_flan_sequence(location, direction, fivep, threep, assembly):
	"""
	Annotating flanking regions of the gRNAs
	:param location: Location of the gRNA target sequence
	:param direction: Direction of the gRNA
	:param fivep: Number of the nucleotides in the 5' flanking region
	:param threep: Number of the nucleotides in the 3' flanking region
	:return:
	"""

	if direction == "left":
		grna_strand = "-1"
	else:
		grna_strand = "1"
	
	server = "http://grch37.rest.ensembl.org" if assembly == "hg19" else "https://rest.ensembl.org"

	# X:154444204-154444226;X:154444207-154444229
	locations = location.split(";")
	resp = []
	for _loc in locations:

		grna_flan_ensembl = server + "/sequence/region/human/%s:%s?expand_3prime=%s;expand_5prime=%s;content-type=text/plain" \
							% (_loc, grna_strand, threep, fivep)
		grna_flan_request = requests.get(grna_flan_ensembl,
											headers={"Content-Type": "text/plain"})

		while grna_flan_request.status_code == 429:
			print("Received status 429. Waiting for 1 second before sending next request.")
			time.sleep(1)
			grna_flan_request = requests.get(grna_flan_ensembl, headers={"Content-Type": "text/plain"})

		if grna_flan_request.status_code != 200:
			print(f"No response from ensembl sequence! {grna_flan_request.status_code}\n")
			resp.append("API problem")
		else:
			resp.append(grna_flan_request.text)

	return ";".join(resp)


def main(path, args):

	df = pandas.read_csv(path + args["FILE"] + ".csv", index_col=0)
	flan_5 = args["FLAN_5"]
	flan_3 = args["FLAN_3"]
	assembly = args["ASSEMBLY"]
	df["gRNA_flanking_sequences"] = df.apply(
		lambda x: extract_gRNA_flan_sequence(location=x.CRISPR_PAM_Location, direction=x.Direction,
											 fivep=flan_5, threep=flan_3, assembly = assembly), axis=1)

	unique_responses = list(
		itertools.chain(
			*map(lambda x: x.split(";"), df.gRNA_flanking_sequences.unique())
		)
	)
	if "API problem" in unique_responses:
		print("API related problems")
		return False
	else:
		df.to_csv(path + args["FILE"] + "_flaking.csv", index=True)
		print("Flaking sequences have been retrieved.")
		return True


if __name__ == '__main__':

	args = take_input()
	# Path
	path = ""
	if args["PATH"][-1] == "/": path = args["PATH"]
	else: path = args["PATH"] + "/"

	_ = main(path, args)
