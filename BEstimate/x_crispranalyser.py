# ----------------------------------------------------------------------------#
#                                                                             #
#                    CRISPR-Analyser | Off target Finding                     #
#                    Author : Bo Fussing bf15@sanger.ac.uk                    #
#                    Copyright (C) 2025 Genome Research Ltd.                  #
#                                                                             #
# ----------------------------------------------------------------------------#

import argparse
import csv
import sqlite3
import sys
import time
from crispr_analyser import search, align, utils
from multiprocessing import Pool, cpu_count

##############################################################################
# Functions


def take_input():
    parser = argparse.ArgumentParser(prog="CRISPR-Analyser", usage="%(prog)s [inputs]")

    parser.add_argument("-c", dest="I_CSV", required=True)
    parser.add_argument("-b", dest="BIN", required=True)
    parser.add_argument("-o", dest="O_CSV", required=True)

    parsed_input = parser.parse_args()
    input_dict = vars(parsed_input)

    return input_dict


def format_summary(summary: list[int]) -> str:
    """Format the summary of off-targets for printing"""
    return f"{{0: {summary[0]}, 1: {summary[1]}, 2: {summary[2]}, 3: {summary[3]}, 4: {summary[4]}}}"


def fetch_crispr_cursor(crispr_ids: list[int], cur: sqlite3.Cursor) -> sqlite3.Cursor:
    """Get the CRISPR data for a given guide sequence"""
    sql = f"SELECT seq, chr_name, chr_start, pam_right FROM crisprs WHERE id IN ({str(crispr_ids)[1:-1]})"
    cur.execute(sql)
    return cur


def usage() -> None:
    print(
        "x_crispranalyser.py -c <input_csv_file> -b <input_bin_file> "
        "-o <output_csv_file>"
    )


def get_ots_for_row(row, guides):
    output_details = list()
    # connect to the database
    con = sqlite3.connect("crisprs.db", check_same_thread=False)
    cur = con.cursor()

    sequence = row["gRNA_Target_Sequence"]
    # Find the CRISPR IDs for the given guide sequence
    crispr_ids = search.search(guides, sequence, True)
    # find the off-target summary for the CRISPR
    binary_sequence = utils.sequence_to_binary_encoding(sequence, 1)
    binary_reverse_sequence = align.reverse_complement_binary(binary_sequence, 20)
    summary, _off_target_ids = align.find_off_targets(
        guides, binary_sequence, binary_reverse_sequence
    )
    cursor = fetch_crispr_cursor(crispr_ids, cur)
    output_summary = {**row}
    output_summary["exact"] = summary[0]
    output_summary["mm1"] = summary[1]
    output_summary["mm2"] = summary[2]
    output_summary["mm3"] = summary[3]
    output_summary["mm4"] = summary[4]
    for record in cursor:
        output_row = {**row}
        output_row["CRISPR_sequence"] = record[0]
        output_row["Chromosome"] = record[1]
        output_row["Start"] = record[2]
        output_row["Strand"] = record[3]
        output_row["Off_target_summary"] = format_summary(summary)
        output_details.append(output_row)
    return output_details, output_summary


def get_off_targets(
    input_bin_file: str, input_csv_file: str, output_csv_file: str
) -> None:
    """Generate off-target summaries for a given input file and binary guides file"""
    with open(input_bin_file, "rb") as guides_file:
        # check that the binary guides file is valid
        utils.check_file_header(guides_file.read(utils.HEADER_SIZE))
        # fetch the guides as a numpy array
        guides = utils.get_guides(guides_file)
        details = list()
        summaries = list()
        details = list()
        results = list()
        with open(input_csv_file, "r") as input_csvfile:
            csvreader = csv.DictReader(input_csvfile)
            headers = csvreader.fieldnames
            details_headers = [
                "CRISPR_sequence",
                "Chromosome",
                "Start",
                "Strand",
                "Off_target_summary",
            ]
            summaries_headers = ["exact", "mm1", "mm2", "mm3", "mm4"]
            pool = Pool(processes=min(2, cpu_count()))
            for row in csvreader:
                result = pool.apply_async(
                    get_ots_for_row,
                    args=(
                        row,
                        guides,
                    ),
                )
                results.append(result)
            for res in results:
                try:
                    output_details, output_summary = res.get()
                    details += output_details
                    summaries.append(output_summary)
                except Exception as e:
                    print(f"Error processing row: {e}")
            pool.close()
            pool.join()
        with open(
            f"{output_csv_file}_ot_annotated_df.csv", "w"
        ) as output_csvfile_summaries:
            csv_summaries_writer = csv.DictWriter(
                output_csvfile_summaries, fieldnames=headers + summaries_headers
            )
            csv_summaries_writer.writeheader()
            csv_summaries_writer.writerows(summaries)
        with open(f"{output_csv_file}wge_return.csv", "w") as output_csvfile_details:
            csv_details_writer = csv.DictWriter(
                output_csvfile_details, fieldnames=headers + details_headers
            )
            csv_details_writer.writeheader()
            csv_details_writer.writerows(details)


def run():
    """Run the off-targets search from the command line."""
    args = take_input()
    input_csv_file = args["I_CSV"]
    input_bin_file = args["BIN"]
    output_csv_file = args["O_CSV"]

    if input_csv_file == "" or input_bin_file == "" or output_csv_file == "":
        sys.exit(2)

    start_time = time.time()
    get_off_targets(input_bin_file, input_csv_file, output_csv_file)
    print(f"Finished processing off-targets: {time.time() - start_time}")


if __name__ == "__main__":
    run()
