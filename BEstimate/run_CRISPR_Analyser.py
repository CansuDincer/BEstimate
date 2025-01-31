import csv
import getopt
import numpy as np
import pandas as pd
import psycopg
import subprocess
import sys

DNA_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


def find_crisprs(guide: str) -> list[int]:
    """Find the CRISPR IDs for a given guide sequence by running the CRISPR-Analyser search command"""
    search_process = subprocess.run(
        [
            "../CRISPR-Analyser/bin/crispr_analyser",
            "search",
            "-i",
            "grch38.bin",
            "-s",
            guide,
            "-p",
            "1",
        ],
        capture_output=True,
        text=True,
    )
    # get the response from subprocess and parse for the CRISPR IDs
    search_output = search_process.stdout.splitlines()
    print(f"search_output: {search_output}")
    return [int(crispr_id) for crispr_id in search_output[1:]]


def get_off_target_summaries(inputfile: str, chunk: pd.DataFrame) -> None:
    """Get the off target summaries for the CRISPRs in the chunk by running the CRISPR-Analyser align command"""
    off_target_summaries = {}
    # get the idx from the dataframe, remove empty values
    crispr_ids = set(chunk["idx"].dropna().tolist())
    script_args = [
        "../CRISPR-Analyser/bin/crispr_analyser",
        "align",
        "-i",
        inputfile,
    ]
    [script_args.append(str(crispr_id)) for crispr_id in crispr_ids]
    align_process = subprocess.run(
        script_args,
        capture_output=True,
        text=True,
    )
    align_outputs = align_process.stdout.split("\n")
    for output in align_outputs:
        output = output.split("\t")
        if len(output) == 4:
            chunk.loc[chunk["idx"] == int(output[0]), "Off_target_summary"] = output[3]


def fetch_crispr_data(
    crispr_id: int, cur: psycopg.Connection.cursor
) -> tuple[str, str, int, str]:
    """Get the CRISPR data for a given guide sequence"""
    cur.execute(
        "SELECT seq, chr_name, chr_start, pam_right FROM crisprs WHERE id = %s",
        (crispr_id,),
    )
    (
        sequence,
        chr_name,
        chr_start,
        pam_right,
    ) = cur.fetchone()
    return [sequence, chr_name, int(chr_start), "+" if pam_right == True else "-"]


def check_crispr_data(guide_sequence: str, crispr: str) -> bool:
    """Check if the guide sequence is in the CRISPR sequence"""
    check = False
    reverse_guide = "".join(DNA_COMPLEMENT.get(base, base) for base in reversed(guide))
    if guide not in crispr_seq and reverse_guide not in crispr_seq:
        check = True
    return check


def usage() -> None:
    print("process.py -c <input_csv_file> -b <input_bin_file> -o <output_csv_file>")


def main(argv) -> None:
    input_csv_file = ""
    input_bin_file = ""
    output_csv_file = ""

    try:
        opts, args = getopt.getopt(argv, "hc:o:b:", ["cfile=", "ofile=", "bfile="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt in ("-c", "--cfile"):
            input_csv_file = arg
        elif opt in ("-b", "--bfile"):
            input_bin_file = arg
        elif opt in ("-o", "--ofile"):
            output_csv_file = arg
    if input_csv_file == "" or input_bin_file == "" or output_csv_file == "":
        usage()
        sys.exit(2)

    input_df = pd.read_csv(input_csv_file)
    COLUMN_NAMES = input_df.columns + ["CRISPR_sequence", "Chromosome", "Start", "Strand", "Off_target_summary", "idx"]

    # connect to the database
    with psycopg.connect(dbname="off_targets", user="crispr", host="localhost") as con:
        processing_row = 0
        with con.cursor() as cur:
            # read the input csv file in chunks
            with pd.read_csv(input_csv_file, chunksize=20) as reader:
                with open(output_csv_file, "a", newline="") as writer:
                    writer_header = True
                    # as we are appending to the file, we need to truncate it first
                    writer.truncate(0)
                    for chunk in reader:
                        print(f"chunk: {chunk}")
                        output_chunk = pd.DataFrame(columns=COLUMN_NAMES)
                        output_index = 0
                        for _index, row in chunk.iterrows():
                            crispr_ids = find_crisprs(row["CRISPR_PAM_Sequence"])
                            for crispr_id in crispr_ids:
                                crisp_data = fetch_crispr_data(crispr_id, cur)
                                output_chunk.loc[output_index] = row
                                output_chunk.at[output_index, "CRISPR_sequence"] = (
                                    crisp_data[0]
                                )
                                output_chunk.at[output_index, "Chromosome"] = (
                                    crisp_data[1]
                                )
                                output_chunk.at[output_index, "Start"] = crisp_data[2]
                                output_chunk.at[output_index, "Strand"] = crisp_data[3]
                                output_chunk.at[output_index, "idx"] = crispr_ids[0]
                                output_index += 1
                        get_off_target_summaries(input_bin_file, output_chunk)
                        print(f"output_chunk: {output_chunk}")
                        # drop the "idx" column and write the chunk to the output file
                        writer.write(
                            output_chunk.to_csv(
                                columns=COLUMN_NAMES[:-1],
                                header=writer_header,
                                index=False,
                            )
                        )
                        writer_header = False


if __name__ == "__main__":
    main(sys.argv[1:])
