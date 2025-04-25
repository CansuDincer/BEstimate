# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                          CRISPR-Analyser | Genome Indexing                               #
#                        Author : Bo Fussing bf15@sanger.ac.uk                             #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

import csv, sys, argparse, sqlite3
import numpy as np
import pandas as pd
import multiprocessing as mp
import gc
from crispr_analyser import search, align, utils

###########################################################################################
# Functions

DNA_COMPLEMENT = {
	"A": "T",
	"C": "G",
	"G": "C",
	"T": "A",
}

def take_input():
	parser = argparse.ArgumentParser(prog="CRISPR-Analyser",
									 usage="%(prog)s [inputs]")

	parser.add_argument("-c", dest="I_CSV", required=True)
	parser.add_argument("-b", dest="BIN", required=True)
	parser.add_argument("-o", dest="O_CSV", required=True)
	parser.add_argument("-p", dest="PATH", required=True)

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


def find_crisprs(guide, wge_path, input_bin_file):
	"""Find the CRISPR IDs for a given guide sequence by running the CRISPR-Analyser search command"""
	exe = wge_path + "bin/crispr_analyser"
	search_process = subprocess.run([exe, 'search', '-i' , input_bin_file, '-s', guide, '-p', '1'], capture_output=True, text=True, timeout=6000)

	# get the response from subprocess and parse for the CRISPR IDs
	search_output = search_process.stdout.splitlines()
	print(f"search_output: {search_output}")
	return [int(crispr_id) for crispr_id in search_output[1:]]


def get_off_target_summaries(inputfile, chunk, wge_path):
	"""Get the off target summaries for the CRISPRs in the chunk by running the CRISPR-Analyser align command"""
	off_target_summaries = {}
	# get the idx from the dataframe, remove empty values
	crispr_ids = set(chunk["idx"].dropna().tolist())
	exe = wge_path + "bin/crispr_analyser"
	script_args = [exe, 'align', '-i', inputfile]

	for crispr_id in crispr_ids:
		script_args.append(str(crispr_id))

	align_process = subprocess.run(script_args, capture_output=True, text=True, timeout=1800)
	align_outputs = align_process.stdout.split("\n")
	for output in align_outputs:
		output = output.split("\t")
		if len(output) == 4:
			chunk.loc[chunk["idx"] == int(output[0]), "Off_target_summary"] = output[3]


def fetch_crispr_data(crispr_id: int, cur: sqlite3.Cursor):
	"""Get the CRISPR data for a given guide sequence"""
	sql = "SELECT seq, chr_name, chr_start, pam_right FROM crisprs WHERE id = %s" % crispr_id
	cur.execute(sql)
	(
		sequence,
		chr_name,
		chr_start,
		pam_right,
	) = cur.fetchone()
	return [sequence, chr_name, int(chr_start), "+" if pam_right == True else "-"]



def check_crispr_data(guide, crispr_seq):
	"""Check if the guide sequence is in the CRISPR sequence"""
	check = False
	reverse_guide = "".join(DNA_COMPLEMENT.get(base, base) for base in reversed(guide))
	if guide not in crispr_seq and reverse_guide not in crispr_seq:
		check = True
	return check


def usage() -> None:
	print("x_crispranalyser.py -c <input_csv_file> -b <input_bin_file> "
		  "-o <output_csv_file> -p <path of wge>")


def get_ots(chunk):
	global wge_path, input_bin_file

	COLUMN_NAMES = list(chunk.columns) + ["CRISPR_sequence", "Chromosome", "Start", "Strand", "Off_target_summary",
										  "idx"]
	output_chunk = pd.DataFrame(columns=COLUMN_NAMES)

	# connect to the database
	con = sqlite3.connect("crisprs.db", check_same_thread=False)
	cur = con.cursor()

	output_index = 0
	for _index, row in chunk.iterrows():

		# Get ot information
		crispr_ids = find_crisprs(row["gRNA_Target_Sequence"], wge_path, input_bin_file)

		for crispr_id in crispr_ids:
			crisp_data = fetch_crispr_data(crispr_id, cur)
			output_chunk.loc[output_index] = row
			output_chunk.loc[output_index, "CRISPR_sequence"] = (crisp_data[0])
			output_chunk.loc[output_index, "Chromosome"] = (crisp_data[1])
			output_chunk.loc[output_index, "Start"] = crisp_data[2]
			output_chunk.loc[output_index, "Strand"] = crisp_data[3]
			output_chunk.loc[output_index, "idx"] = crispr_ids[0]
			output_index += 1

	get_off_target_summaries(input_bin_file, output_chunk, wge_path)

	# drop the "idx" column and write the chunk to the output file
	chunk_df = output_chunk.reset_index()[[col for col in COLUMN_NAMES if col != "idx"]]
	return chunk_df


def main():
	global wge_path, input_bin_file, output_csv_file
	# read the input csv file in chunkss

	results_dfs = list()
	with mp.Pool(4) as pool:
		with pd.read_csv(input_csv_file, chunksize=100) as reader:
			for chunk in reader:
				results = pool.apply_async(get_ots, (chunk,))
				gc.collect()
				print(results.get())
				results_dfs.append(results.get())
	gc.collect()

	output_df = pd.concat(results_dfs, ignore_index=True)
	output_df.to_csv(output_csv_file + "wge_return.csv")

	if len(output_df.index) > 0:
		output_df = output_df[[col for col in output_df.columns if
							   col not in ["CRISPR_sequence", "Chromosome", "Start", "Strand"]]].drop_duplicates()
		output_df["exact"] = output_df.apply(lambda x: x.Off_target_summary.split(", ")[0].split(": ")[1], axis=1)
		output_df["mm1"] = output_df.apply(lambda x: x.Off_target_summary.split(", ")[1].split(": ")[1], axis=1)
		output_df["mm2"] = output_df.apply(lambda x: x.Off_target_summary.split(", ")[2].split(": ")[1], axis=1)
		output_df["mm3"] = output_df.apply(lambda x: x.Off_target_summary.split(", ")[3].split(": ")[1], axis=1)
		output_df["mm4"] = output_df.apply(lambda x: x.Off_target_summary.split(", ")[4].split(": ")[1].split("}")[0],
										   axis=1)
		output_df.to_csv(output_csv_file + "_ot_annotated_df.csv", index=False)

	return output_df


if __name__ == "__main__":

	args = take_input()
	# Output Path
	path = ""
	if args["PATH"][-1] == "/":
		path = args["PATH"]
	else:
		path = args["PATH"] + "/"

	input_csv_file = args["I_CSV"]
	input_bin_file = args["BIN"]
	output_csv_file = args["O_CSV"]
	wge_path = args["PATH"]

	if input_csv_file == "" or input_bin_file == "" or output_csv_file == "" or wge_path == "":
		sys.exit(2)

	main()