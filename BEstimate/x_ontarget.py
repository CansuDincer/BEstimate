# -----------------------------------------------------------------------------#
#                                                                              #
#                              B E s t i m a t e                               #
#                        	   On-target Scoring                               #
#                    Author : Cansu Dincer cd7@sanger.ac.uk                    #
#                    Copyright (C) 2025 Genome Research Ltd.                   #
#                                                                              #
# -----------------------------------------------------------------------------#

import argparse, os, requests, sys, pandas
from rs3.seq import predict_seq
from BEstimate import Ensembl

if f"{os.getcwd()}../FORECasT-BE/" not in sys.path:
	sys.path.append("{os.getcwd()}../FORECasT-BE/")

import forecast_be


###############################################################################
# capture command line arguments

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate On-target Scoring",
									 usage="%(prog)s [inputs]")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# BASIC INFORMATION
	# gRNA on target scoring
	parser.add_argument("-gene", dest="GENE", required=True,
						help="The hugo symbol of the interested gene!")
	parser.add_argument("-assembly", dest="ASSEMBLY", required=True, default="GRCh38",
						help="The genome assembly that will be used!")
	parser.add_argument("-mutation_file", dest="MUTATION_FILE", default=None, type=argparse.FileType('r'),
						help="A file for the mutations on the interested gene that you need to integrate "
							 "into guide and/or annotation analysis")
	parser.add_argument("-rs3", dest="RULESET3", action="store_true",
						help="The boolean option if the user wants to add on target RuleSet3 scoring for the gRNAs")
	parser.add_argument("-fc", dest="FCAST", action="store_true",
						help="The boolean option if the user wants to add on target ForeCast gRNAs efficiency info")
	parser.add_argument("-edit", dest="EDIT", help="The searched nuceleotide", required=True)
	parser.add_argument("-iname", dest="INPUT",
						help="The input BEstimate file extension (edit_df/ summary_df/ ot_annotated_df)",
						required=True)
	parser.add_argument("-ofile", dest="OUTPUT_FILE", default="output", required=True,
						help="The output file name, if not specified \"position\" will be used!")
	parser.add_argument("-o", dest="OUTPUT_PATH", default=os.getcwd() + "/", required=True,
						help="The path for output. If not specified the current directory will be used!")
	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


def run_ruleset3(final_df, ensembl_object, location_col):
	"""
	 You have to have 20 nuc in protospacer and 3 nuc of PAM
	"""

	final_df["rs3_sequence"] = final_df.apply(lambda x: ensembl_object.prep_rs3_seq(
		location=x[location_col], direction=x.Direction), axis=1)

	if len(final_df.rs3_sequence.unique()) == 1 and list(final_df.rs3_sequence.unique())[0] == None:
		print("No ensembl response > No Rule Set 3 analysis, please rerun!")
		return final_df
	else:
		final_df = final_df.reset_index(drop=True)
		final_df["RuleSet3_Hsu2013"], final_df["RuleSet3_Chen2013"] = None, None
		for grna, grna_df in final_df.groupby("rs3_sequence"):
			inds = list(grna_df.index)
			final_df.loc[inds, "RuleSet3_Hsu2013"] = predict_seq([grna], sequence_tracr='Hsu2013')[0]
			final_df.loc[inds, "RuleSet3_Chen2013"] = predict_seq([grna], sequence_tracr='Chen2013')[0]

		return final_df


def run_forecastbe(final_df, searched_nucleotide):
	"""
	You have to have 20 nuc in protospacer and used for ABE and CBE
	Z scores are calculated with 3-10th positions in the protospacer when 21-23 represents PAM
	"""

	if searched_nucleotide == "C":
		editor = "CBE"
	elif searched_nucleotide == "A":
		editor = "ABE"
	else:
		return None

	final_df = final_df.reset_index(drop=True)
	final_df["FORECasT-BE"] = None
	for grna, grna_df in final_df.groupby("gRNA_Target_Sequence"):
		inds = list(grna_df.index)
		final_df.loc[inds, "FORECasT-BE"] = forecast_be.predict_total(grna, editor=editor, mean=None, std=None)

	return final_df


def main():
	global args

	path = ""
	if args["OUTPUT_PATH"][-1] == "/":
		path = args["OUTPUT_PATH"]
	else:
		path = args["OUTPUT_PATH"] + "/"
	try:
		os.mkdir(path)
	except FileExistsError:
		pass

	if args["RULESET3"]:
		rs3_score = True
	else:
		rs3_score = False

	if args["FCAST"]:
		fc_score = True
	else:
		fc_score = False

	if args["MUTATION_FILE"]:
		mutations = list()
		for line in args["MUTATION_FILE"].readlines():
			mutations.append(line.strip())
	else:
		mutations = None

	if rs3_score or fc_score:
		print("""\n
--------------------------------------------------------------
		Annotation - On Target Annotation
--------------------------------------------------------------
					\n""")

		df = pandas.read_csv(f"{path}{args['OUTPUT_FILE']}_{args['INPUT']}.csv")

		# Create Ensembl object
		ensembl_obj = Ensembl(hugo_symbol=args["GENE"], assembly=args["ASSEMBLY"])
		ensembl_obj.extract_gene_id()

		if ensembl_obj.gene_id == '': sys.exit("No corresponding Ensembl Gene ID could be found!")

		ensembl_obj.extract_sequence(ensembl_obj.gene_id, mutations=mutations)


		if f"{args['OUTPUT_FILE']}_scored_{args['INPUT']}.csv" not in os.listdir(path):
			if rs3_score:
				if final_text != "edit_df":
					location_col = "CRISPR_PAM_Location"
				else:
					location_col = "Location"

				rs3_score_df = run_ruleset3(final_df=df, location_col=location_col, ensembl_object=ensembl_obj)
				df = rs3_score_df.copy()
				print("RuleSet3 on target scoring was added!")

			if fc_score:
				fc_score_df = run_forecastbe(final_df=df, searched_nucleotide=args["EDIT"])
				print("FORECast-BE gRNA efficiency was added!")
				df = fc_score_df.copy()

			df.to_csv(f"{path}{args['OUTPUT_FILE']}_scored_{args['INPUT']}.csv", index=False)
			print("On-target annotation was completed!")
			final_df = df.copy()
		else:
			df = pandas.read_csv(f"{path}{args['OUTPUT_FILE']}_scored_{args['INPUT']}.csv")
			print("On target annotation was read!")
			final_df = df.copy()

	return final_df


if __name__ == '__main__':
	# -----------------------------------------------------------------------------------------#
	# Retrieve input

	args = take_input()

	# -----------------------------------------------------------------------------------------#
	# Execution

	_ = main()

	print("""\n
--------------------------------------------------------------
	The BEstimate on-target analysis successfully finished!
--------------------------------------------------------------
	\n""")

