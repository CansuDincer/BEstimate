# Copyright (C) 2025 Genome Research Ltd.

import getopt
import numpy as np
import sys

from .utils import (
    FILE_VERSION,
    HEADER_SIZE,
    METADATA_SIZE,
    check_file_header,
    get_guides,
    get_file_metadata,
    print_metadata,
    sequence_to_binary_encoding,
    reverse_complement,
)


def search(
    guides: np.ndarray,
    sequence: str,
    verbose: bool = False,
) -> list[int]:
    """Search for a sequence in an indexed binary file

    Args:
        guides: The numpy uint64 array of guides
        sequence: The query sequence to search for
        verbose: A boolean to print verbose output
    """
    reverse_sequence = reverse_complement(sequence)
    query_sequence = sequence_to_binary_encoding(sequence, 1)
    reverse_query_sequence = sequence_to_binary_encoding(reverse_sequence, 0)
    indices = np.where(
        (guides == query_sequence) | (guides == reverse_query_sequence)
    )
    # the binary index is 0-based,
    # so we add the offset and 1 to make it 1-based as per the db
    # this follows how we numbered the WGE index
    return [x + 1 for x in indices[0].tolist()]


def run(argv=sys.argv[1:]) -> None:
    """Run the search command from the command line."""
    inputfile = ""
    sequence = ""

    def usage() -> None:
        print(
            """Usage: poetry run search [options...]
-h, --help            Print this help message
-i, --ifile <file>    The input binary guides file
-s, --sequence <str>  The guide sequence to search for
"""
        )

    try:
        opts, args = getopt.getopt(
            argv,
            "hi:s:",
            [
                "help",
                "ifile=",
                "sequence=",
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
            inputfile = arg
        elif opt in ("-s", "--sequence"):
            sequence = arg
    if inputfile == "" or sequence == "":
        usage()
        sys.exit(2)

    with open(inputfile, "rb") as in_file:
        check_file_header(in_file.read(HEADER_SIZE))
        print(f"Version is {FILE_VERSION}", file=sys.stderr)
        metadata = get_file_metadata(in_file.read(METADATA_SIZE))
        print_metadata(metadata)
        guides = get_guides(in_file, verbose=True)
        print(f"Loaded {guides.size} sequences", file=sys.stderr)
        indices = search(
            guides=guides,
            sequence=sequence,
            verbose=True,
        )
        print(f"Found {len(indices)} exact matches", file=sys.stderr)
        print("Found the following matches:", file=sys.stderr)
        for idx in indices:
            print(f"\t{idx + metadata.offset}")


if __name__ == "__main__":
    run()
