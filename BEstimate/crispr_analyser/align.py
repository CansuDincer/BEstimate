# Copyright (C) 2025 Genome Research Ltd.

import getopt
import numpy as np
from numba import jit, cuda
import sys

from .utils import (
	ERROR_STR,
	FILE_VERSION,
	HEADER_SIZE,
	METADATA_SIZE,
	check_file_header,
	get_guides,
	get_file_metadata,
	print_metadata,
)

MAX_MISSMATCHES = 5
MAX_OFF_TARGETS = 2000
PAM_ON = np.left_shift(1, 40, dtype=np.uint64)
PAM_OFF = np.invert(PAM_ON, dtype=np.uint64)


@cuda.jit
def find_off_targets_kernel(
	guides: np.ndarray,
	query_sequence: np.uint64,
	reverse_query_sequence: np.uint64,
	summary: np.ndarray,
	off_target_ids_idx: np.ndarray,
	off_target_ids: np.ndarray,
	offset: np.uint64,
) -> None:
	"""Find off-targets for a given query sequence using CUDA

	Args:
		guides: The array of encoded gRNA sequences
		query_sequence: The query sequence
		reverse_query_sequence: The reverse complement of the query sequence
		summary: The array to store the results
		off_target_ids_idx: The index for the off-target_ids array
		off_target_ids: The array to store the off-target ids
		offset: The offset of the guides default 0
	"""
	index = cuda.grid(1)
	threads_per_grid = cuda.gridDim.x * cuda.blockDim.x

	for i in range(index, guides.size, threads_per_grid):
		if index < guides.size:
			if guides[i] == ERROR_STR:
				continue
			match = query_sequence ^ guides[i]
			if match & PAM_ON:
				match = reverse_query_sequence ^ guides[i]
			match = match & PAM_OFF
			match = (match | (match >> 1)) & 0x5555555555555555
			nos_off_targets = cuda.libdevice.popcll(match)
			if nos_off_targets < MAX_MISSMATCHES:
				cuda.atomic.add(summary, nos_off_targets, 1)
				idx = cuda.atomic.add(off_target_ids_idx, 0, 1)
				cuda.atomic.add(off_target_ids, idx, offset + i + 1)


@jit
def find_off_targets(
	guides: np.ndarray,
	query_sequence: np.uint64,
	reverse_query_sequence: np.uint64,
	offset: np.uint64 = 0):
	"""Find off-targets for a given query sequence using CPU

	Args:
		guides: The array of encoded gRNA sequences
		query_sequence: The query sequence
		reverse_query_sequence: The reverse complement of the query sequence
		offset: The offset of the guides default 0
	"""
	summary = [0] * MAX_MISSMATCHES
	off_target_ids = []
	for i, guide in enumerate(guides):
		if guide == ERROR_STR:
			continue
		match = query_sequence ^ guide
		if match & PAM_ON:
			match_r = reverse_query_sequence ^ guide
			match_count = _pop_count(match_r & PAM_OFF)
		else:
			match_count = _pop_count(match & PAM_OFF)
		if match_count < MAX_MISSMATCHES:
			summary[match_count] += 1
			off_target_ids.append(offset + i + 1)
	return summary, off_target_ids


@jit
def _pop_count(x: np.uint64):
	"""Count bits in integer accounting for encoding
	as everything is two bit we must convert them all to one bit,
	to do this we must turn off all MSBs, but before we can do that
	we need to ensure that when an MSB is set to 1, the LSB is also set.
	the 4 will change pam_right to 0 so it doesnt get counted.
	5 is 0101, 4 is 0100

	Args:
		x: The integer to count the bits of
	"""
	x = np.uint64((x | (x >> 1)) & 0x5555555555555555)
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
	return (0x0101010101010101 * x) >> 56


def reverse_complement_binary(sequence: np.uint64, size: int):
	"""Reverse complement a binary sequence

	Args:
		sequence: The binary sequence to reverse complement
		size: The size of the sequence
	"""
	mask = np.uint64(0xFFFFFFFFFFFFFFFF >> (63 - (size * 2)))
	sequence = ~sequence & mask
	reversed = sequence >> np.uint64(size * 2)
	shift = 0
	for i in range(size):
		reversed <<= np.uint64(2)
		reversed |= (sequence >> shift) & 0x3
		shift += 2
	return reversed


def print_off_targets(
	crispr_id: np.uint64,
	summary: np.ndarray,
	off_target_ids: np.ndarray,
	species_id: np.uint8,
) -> None:
	"""Print the off targets to the console

	Args:
		crispr_id: The id of the CRISPR
		summary: The summary of the off targets
		off_target_ids: The off target CRISPR ids
		species_id: The species id
	"""
	summary_output = ", ".join(
		[f'{i}: {summary[i]}' for i in range(MAX_MISSMATCHES)]
	)

	if len(off_target_ids) >= MAX_OFF_TARGETS:
		print(f'{crispr_id}\t{species_id}\t{{{summary_output}}}')
	else:
		print(f'{crispr_id}\t{species_id}\t'
			  f'{{{",".join(map(str, off_target_ids))}}}\t{{{summary_output}}}'
			 )


def run(argv=sys.argv[1:]) -> None:
    inputfile = ""
    use_cuda = True

    def usage() -> None:
        print(
            """Usage: poetry run align [options...] [ids...]
-h, --help            Print this help message
-i, --ifile <file>    The input binary guides file
--no-cuda             Do not use CUDA GPU acceleration
[ids...]              The ids of the CRISPRs to find off-targets for
"""
        )

    try:
        opts, args = getopt.getopt(
            argv,
            "hi:",
            [
                "help",
                "ifile=",
                "no-cuda",
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
        elif opt == "--no-cuda":
            use_cuda = False
    if inputfile == "" or len(args) == 0:
        usage()
        sys.exit(2)

    with open(inputfile, "rb") as in_file:
        check_file_header(in_file.read(HEADER_SIZE))
        print(f"Version is {FILE_VERSION}", file=sys.stderr)
        metadata = get_file_metadata(in_file.read(METADATA_SIZE))
        print_metadata(metadata)
        guides = get_guides(in_file, verbose=True)

        print("Searching for off targets", file=sys.stderr)
        if use_cuda & cuda.is_available():
            memory_required = guides.size * 8 / 1024 / 1024
            print(
                f"Requires {memory_required} MB of GPU memory",
                file=sys.stderr,
            )
            device_guides = cuda.to_device(guides)
            threads_per_block = 256
            blocks_per_grid = (
                guides.size + (threads_per_block - 1) // threads_per_block
            )

            for i in range(len(args)):
                query_sequence = guides[int(args[i]) - 1]
                reverse_query_sequence = reverse_complement_binary(
                    query_sequence, 20
                )
                summary = np.zeros(MAX_MISSMATCHES, dtype=np.uint32)
                off_target_ids_idx = np.zeros(1, dtype=np.uint32)
                off_target_ids = np.zeros(MAX_OFF_TARGETS, dtype=np.uint32)
                device_summary = cuda.to_device(summary)
                device_off_target_ids_idx = cuda.to_device(off_target_ids_idx)
                device_off_target_ids = cuda.to_device(off_target_ids)
                find_off_targets_kernel[blocks_per_grid, threads_per_block](
                    device_guides,
                    query_sequence,
                    reverse_query_sequence,
                    device_summary,
                    device_off_target_ids_idx,
                    device_off_target_ids,
                    metadata.offset,
                )
                host_summary = device_summary.copy_to_host()
                host_off_target_ids = np.trim_zeros(
                    device_off_target_ids.copy_to_host()
                )
                print_off_targets(
                    args[i],
                    host_summary,
                    np.sort(host_off_target_ids),
                    metadata.species_id,
                )
        else:
            for i in range(len(args)):
                query_sequence = guides[int(args[i]) - 1]
                reverse_query_sequence = reverse_complement_binary(
                    query_sequence, 20
                )
                summary, off_target_ids = find_off_targets(
                    guides,
                    query_sequence,
                    reverse_query_sequence,
                    metadata.offset,
                )
                print_off_targets(
                    args[i], summary, off_target_ids, metadata.species_id
                )


if __name__ == "__main__":
    run()
