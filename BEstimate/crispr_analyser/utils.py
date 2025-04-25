# Copyright (C) 2025 Genome Research Ltd.

from dataclasses import dataclass
import numpy as np
import struct
import sys
import time
import typing

COMPLEMENT_MAP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
ENCODING_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
ERROR_STR = np.uint64(0xFFFFFFFFFFFFFFFF)
FILE_HEADER_FORMAT = "<BL"
FILE_VERSION = np.uint16(3)
HEADER_SIZE = struct.calcsize(FILE_HEADER_FORMAT)
METADATA_FORMAT = "<QQQB30s30s"
METADATA_SIZE = struct.calcsize(METADATA_FORMAT)
PADDING_FORMAT = "<BBB"


@dataclass
class Metadata:
    number_of_sequences: np.uint64
    sequence_length: np.uint64
    offset: np.uint64
    species_id: np.uint8
    species_name: str
    assembly: str


def get_guides(
    guidesfile_handle: typing.BinaryIO, verbose: bool = False
) -> np.ndarray:
    """Get the array of guides from the binary guides file

    Args:
        guidesfile_handle: The file handle of the guides file
        verbose: A boolean to print verbose output
    """
    if verbose:
        start = time.time()
    guidesfile_handle.seek(HEADER_SIZE)
    # read the number of sequences in the header
    number_of_guides = struct.unpack("<Q", guidesfile_handle.read(8))[0]
    guidesfile_handle.seek(
        HEADER_SIZE + METADATA_SIZE + struct.calcsize(PADDING_FORMAT)
    )
    guides = np.fromfile(guidesfile_handle, dtype=np.uint64, count=-1)
    if number_of_guides != guides.size:
        raise ValueError("Invalid number of guides")
    if verbose:
        print(
            f"Loading took {time.time() - start:.2f} seconds", file=sys.stderr
        )
    return guides


def sequence_to_binary_encoding(sequence: str, pam_right: int) -> np.uint64:
    """Convert a string DNA sequence to bits accounting for pam right or left

    Args:
        sequence: The string DNA sequence to convert
        pam_right: An integer indicating if PAM is on the right
    """
    bits = np.uint64(pam_right)
    for character in sequence:
        if ENCODING_MAP[character] == 4:
            bits = ERROR_STR  # set the error flag if we encounter an 'N'
            break
        else:
            bits <<= 2  # shift left to make room for new char
            bits |= ENCODING_MAP[character]  # add our new char to the end
    return bits


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        sequence: The string DNA sequence to reverse complement.
    """
    return "".join([COMPLEMENT_MAP[base] for base in reversed(sequence)])


def check_file_header(bytes: bytes) -> None:
    """Check the header of the file from binary data

    Args:
        bytes: The binary data to check
    """
    if len(bytes) != HEADER_SIZE:
        raise ValueError("Invalid file header length")
    values = struct.unpack(FILE_HEADER_FORMAT, bytes)
    version = values[1]
    if version != FILE_VERSION:
        raise ValueError("Invalid file version")


def get_file_metadata(bytes: bytes) -> Metadata:
    """Parse the metadata from binary data

    Args:
        bytes: The binary data to parse
    """
    if len(bytes) != METADATA_SIZE:
        raise ValueError("Invalid metadata length")
    metadata = struct.unpack(METADATA_FORMAT, bytes)
    metadata = Metadata(
        number_of_sequences=metadata[0],
        sequence_length=metadata[1],
        offset=metadata[2],
        species_id=metadata[3],
        species_name=parse_c_string(metadata[4]),
        assembly=parse_c_string(metadata[5]),
    )
    return metadata


def parse_c_string(data: bytes) -> str:
    """Parse a null-terminated C string from a byte array

    Args:
        data: The byte array to parse
    """
    null = data.find(b"\x00")
    if null < 0:
        return ""
    cstr = data[:null].decode("utf-8")
    return cstr


def print_metadata(metadata: Metadata) -> None:
    """Print the metadata information to STDERR

    Args:
        metadata: The Metadata object
    """
    print(
        f"Assembly is {metadata.assembly} ({metadata.species_name})",
        file=sys.stderr,
    )
    print(f"File has {metadata.number_of_sequences} sequences", file=sys.stderr)
    print(f"Sequence length is {metadata.sequence_length}", file=sys.stderr)
    print(f"Offset is {metadata.offset}", file=sys.stderr)
    print(f"Species id is {metadata.species_id}", file=sys.stderr)
