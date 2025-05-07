# Copyright (C) 2025 Genome Research Ltd.

import getopt
import numpy as np
import struct
import sys
import time

from .utils import (
    FILE_VERSION,
    sequence_to_binary_encoding,
)


def create_metadata(
    number_of_sequences: np.uint64,
    sequence_length: np.uint64,
    offset: np.uint64,
    species_id: np.uint8,
    species_name: str,
    assembly: str,
) -> bytes:
    """Create metadata for the version 3 of the binary output file

    Args:
        number_of_sequences: Number of sequences in the file
        sequence_length: Length of each sequence (guide + PAM)
        offset: Offset of the first sequence
        species_id: ID of the species e.g. 1
        species_name: Name of the species e.g. 'Human"
        assembly: Genome assembly used e.g. 'GRCh38'
    """
    format = "<QQQB30s30s"
    number_of_sequences = np.uint64(number_of_sequences)
    sequence_length = np.uint64(sequence_length)
    offset = np.uint64(offset)
    species_id = np.uint8(species_id)
    species_name = species_name.encode("utf-8")
    assembly = assembly.encode("utf-8")

    return struct.pack(
        format,
        number_of_sequences,
        sequence_length,
        offset,
        species_id,
        species_name,
        assembly,
    )


def parse_record(record: str, guide_length: int, pam_length: int) -> (str, int):
    """Parse a line from the input CSV file

    Args:
        record: A line from the input CSV file
        guide_length: The length of the guide sequence (CRISPR excluding PAM)
        pam_length: The length of the PAM sequence (CRISPR excluding guide)
    """
    records = record.split(",")
    if len(records) != 5:
        print(f"Record '{record}' contains {len(records)} columns, expected 5")
        sys.exit(2)
    pam_right = int(records[3])
    crispr_sequence = records[2]
    if len(crispr_sequence) != guide_length + pam_length:
        print(
            f"Record {record} has sequence_length {len(crispr_sequence)}, "
            "expected {guide_length + pam_length}"
        )
        sys.exit(2)
    if pam_right == 0:
        guide_sequence = crispr_sequence[pam_length:]
    else:
        guide_sequence = crispr_sequence[:guide_length]
    return guide_sequence, pam_right


def index(
    inputfiles: list[str],
    outputfile: str,
    species: str,
    assembly: str,
    offset: int,
    species_id: int,
    guide_length: int = 20,
    pam_length: int = 3,
    verbose: bool = False,
):
    """Run the CRISPR indexer.

    Args:
        inputfiles: The input CSV files e.g. ['input1.csv', 'input2.csv']
        outputfile: The name of the output binary file to be generated.
        species: The species name e.g. 'Human'.
        assembly: The assembly name e.g. 'GRCh38'.
        offset: The integer for offset after which to start numbering ID.
        species_id: The integer of the species ID e.g. 1.
        guide_length: The length of the guide sequence default is 20.
            (CRISPR excluding PAM)
        pam_length: The length of the PAM sequence default is 3.
            (CRISPR excluding guide)
        verbose: A boolean indicating if verbose output is enabled.
            Default is False.
    """
    if verbose:
        start = time.time()
    number_of_sequences = np.uint64(0)
    if verbose:
        print("Outfile:")
        print(f"\t{outputfile}")
        print("Inputfiles:")
        for inputfile in inputfiles:
            print(f"\t{inputfile}")
        print("writing metadata")
        print(f"Version: {FILE_VERSION}")

    with open(outputfile, "wb") as out_file:
        # write the file header
        out_file.write(struct.pack("<BL", np.uint8(1), np.uint(FILE_VERSION)))
        # write the metadata
        out_file.write(
            create_metadata(
                number_of_sequences,
                guide_length,
                offset,
                species_id,
                species,
                assembly,
            )
        )
        # put in a separator of 3 empty bytes before the vector of sequences
        out_file.write(struct.pack("<BBB", 0, 0, 0))
        for inputfile in inputfiles:
            if verbose:
                print(f"Processing {inputfile}")
            with open(inputfile, "r") as in_file:
                for line in in_file:
                    sequence, pam_right = parse_record(
                        line, guide_length, pam_length
                    )
                    record = sequence_to_binary_encoding(sequence, pam_right)
                    out_file.write(struct.pack("<Q", record))
                    number_of_sequences += 1
        # write the number of sequences in the correct position in the file
        out_file.seek(5)
        out_file.write(struct.pack("<Q", number_of_sequences))
        if verbose:
            total = time.time() - start
            print(
                f"Converted {number_of_sequences} sequences in {total} seconds"
            )


def run(argv=sys.argv[1:]):
    """Run the CRISPR indexer from the command line."""
    inputfiles = []
    outputfile = ""
    species = ""
    species_id = np.uint8(0)
    assembly = ""
    offset = np.uint64(0)
    guide_length = 20
    pam_length = 3

    def usage():
        print(
            """Usage: poetry run index [options...]
-a, --assembly <name>         The assembly name
-e, --species_id <integer>    The species ID
-f, --offset <integer>        The offset to start numbering from
-g, --guide_length <integer>  The length of the guide sequence
-h, --help                    Print this help message
-i, --ifile <file>            The input CSV file
-o, --ofile <file>            The ouput file
-p, --pam_length <integer>    The length of the PAM sequence
-s, --species <name>          The species name
"""
        )

    try:
        opts, args = getopt.getopt(
            argv,
            "hi:o:a:s:f:e:g:p:",
            [
                "help",
                "ifile=",
                "ofile=",
                "assembly=",
                "species=",
                "offset=",
                "species_id=",
                "guide_length=",
                "pam_length=",
            ],
        )
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfiles.append(arg)
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-a", "--assembly"):
            assembly = arg
        elif opt in ("-s", "--species"):
            species = arg
        elif opt in ("-f", "--offset"):
            offset = np.uint64(arg)
        elif opt in ("-e", "--species_id"):
            species_id = np.uint8(arg)
        elif opt in ("-g", "--guide_length"):
            guide_length = int(arg)
        elif opt in ("-p", "--pam_length"):
            pam_length = int(arg)
        else:
            print("Unhandled Option")
            usage()
            sys.exit(2)
    if inputfiles == [] or outputfile == "" or assembly == "" or species == "":
        usage()
        sys.exit(2)

    index(
        inputfiles,
        outputfile,
        species,
        assembly,
        offset,
        species_id,
        guide_length,
        pam_length,
        verbose=True,
    )


if __name__ == "__main__":
    run()
