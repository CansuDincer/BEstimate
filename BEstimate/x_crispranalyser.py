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

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict

def format_summary(summary):
    """Format the summary of off-targets for printing"""
    return f"{{0: {summary[0]}, 1: {summary[1]}, 2: {summary[2]}, 3: {summary[3]}, 4: {summary[4]}}}"


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


def usage() -> None:
	print("x_crispranalyser.py -c <input_csv_file> -b <input_bin_file> "
		  "-o <output_csv_file>")


def get_ots(chunk, guides):
    global input_bin_file

    COLUMN_NAMES = list(chunk.columns) + ["CRISPR_sequence", "Chromosome", "Start", "Strand", "Off_target_summary",
    									  "idx"]
    output_chunk = pd.DataFrame(columns=COLUMN_NAMES)

    # connect to the database
    con = sqlite3.connect("crisprs.db", check_same_thread=False)
    cur = con.cursor()

    output_index = 0
    for _index, row in chunk.iterrows():
        sequence = row["gRNA_Target_Sequence"]
		# Find the CRISPR IDs for the given guide sequence
        crispr_ids = search.search(guides, sequence)
        # find the off-target summary for the CRISPR
        binary_sequence = utils.sequence_to_binary_encoding(sequence, 1)
        binary_reverse_sequence = align.reverse_complement_binary(binary_sequence, 20)
        summary, _off_target_ids = align.find_off_targets(guides, binary_sequence, binary_reverse_sequence)
        print(f"summary: {format_summary(summary)}")
        for crispr_id in crispr_ids:
                crispr_data = fetch_crispr_data(crispr_id, cur)
                output_chunk.loc[output_index] = row
                output_chunk.loc[output_index, "CRISPR_sequence"] = (crispr_data[0])
                output_chunk.loc[output_index, "Chromosome"] = (crispr_data[1])
                output_chunk.loc[output_index, "Start"] = crispr_data[2]
                output_chunk.loc[output_index, "Strand"] = crispr_data[3]
                output_chunk.loc[output_index, "Off_target_summary"] = format_summary(summary)
                output_chunk.loc[output_index, "idx"] = crispr_ids[0]
                output_index += 1

    # drop the "idx" column and write the chunk to the output file
    chunk_df = output_chunk.reset_index()[[col for col in COLUMN_NAMES if col != "idx"]]
    return chunk_df


def main():
    global input_bin_file, output_csv_file

    results_dfs = list()
    with open(input_bin_file, "rb") as guides_file:
        print(f"Reading guides from {input_bin_file}")
        # check that the binary guides file is valid
        utils.check_file_header(guides_file.read(utils.HEADER_SIZE))
        # fetch the guides as a numpy array
        guides = utils.get_guides(guides_file)
        with mp.Pool(4) as pool:
            # read the input csv file in chunkss
            with pd.read_csv(input_csv_file, chunksize=20) as reader:
                for chunk in reader:
                    results = pool.apply_async(get_ots, (chunk,guides,))
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

	if input_csv_file == "" or input_bin_file == "" or output_csv_file == "":
		sys.exit(2)

	main()
