# -----------------------------------------------------------------------------#
#                                                                              #
#                              B E s t i m a t e                               #
#                        Genome Retrieval and Indexing                         #
#                    Author : Cansu Dincer cd7@sanger.ac.uk                    #
#                    Copyright (C) 2025 Genome Research Ltd.                   #
#                                                                              #
# -----------------------------------------------------------------------------#

import argparse
import os
import requests
import shutil
import subprocess
import sys
from crispr_analyser import index, gather

CHROMOSOMES = list(range(1, 23)) + ["X", "Y", "MT"]

# Extracting CRISPRs from the Humen Reference Genome

###############################################################################
# capture command line arguments


def take_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="x_genome.py",
        description="Script for indexing CRISPRs for finding off-targets",
        usage="%(prog)s [inputs]",
    )
    parser.add_argument(
        "--pamseq",
        "-p",
        default="NGG",
        required=False,
        help="The PAM sequence in which features used "
        "for searching activity window and editable nucleotide.",
    )
    parser.add_argument(
        "--assembly",
        "-a",
        default="GRCh38",
        choices=["GRCh38", "GRCh37"],
        required=False,
        help="The genome assembly that will be used!",
    )
    parser.add_argument(
        "--output_path",
        "-o",
        default=f"{os.getcwd()}/",
        required=False,
        help="The path for output. If not specified the current directory "
        "will be used!",
    )
    parser.add_argument(
        "--ensembl_version",
        "-e",
        default="113",
        required=False,
        help="The ensembl version in which genome will be retrieved "
        "(if the assembly is GRCh37 then please use <=75)",
    )
    parser.add_argument(
        "--offtargets_path",
        "-ot",
        default=f"{os.getcwd()}/../offtargets",
        required=False,
        help="The path to the root offtargets output directory",
    )
    return parser.parse_args()


##############################################################################
# Functions


def base_file_name(assembly: str, ensembl_version: str) -> str:
    """Generate the base file name"""
    name = ""
    if assembly == "GRCh37":
        name = "Homo_sapiens.GRCh37.%s" % ensembl_version
    elif assembly == "GRCh38":
        name = "Homo_sapiens.GRCh38"
    return name


def check_genome_files_exist(
    assembly: str, ensembl_version: str, ot_path: str
) -> bool:
    """Check to see if the genome FASTA files have been dowloaded."""
    file_directory = f"{ot_path}/genome_files"
    if os.path.exists(file_directory) is False:
        return False
    files_exist = True
    base_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )
    for chromosome in CHROMOSOMES:
        file_name = "%s.dna.chromosome.%s.fa.gz" % (base_name, chromosome)
        if file_name not in os.listdir(file_directory):
            files_exist = False
    return files_exist


def fetch_genome_files(
    assembly: str, ensembl_version: str, ot_path: str
) -> None:
    """Download the genome FASTA files from Ensembl"""
    print(
        "Genome assembly is not found, BEstimate is downloading "
        "the %s Ensembl genome - version %s\n" % (assembly, ensembl_version)
    )
    try:
        os.mkdir(f"{ot_path}/genome_files")
    except FileExistsError:
        pass

    file_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )

    base_url = (
        "https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/dna"
        % ensembl_version
    )
    for chromosome in CHROMOSOMES:
        url = f"{base_url}/{file_name}.dna.chromosome.{chromosome}.fa.gz"
        response = requests.get(url)
        file_path = "%s/genome_files/%s.dna.chromosome.%s.fa.gz" % (
            ot_path,
            file_name,
            chromosome,
        )

        if response.status_code == 200:
            with open(file_path, "wb") as file:
                file.write(response.content)
        else:
            print(f"failed to download {url}")


def gather_crisprs_from_genome(
    assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str
) -> None:
    """Gather CRISPRs from the FASTA files and generate gRNA binary index"""
    file_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )
    try:
        os.mkdir(f"{ot_path}/crispr_csv")
    except FileExistsError:
        pass

    # Gather all chromosome fasta files into csv files
    print(
        "From genome assembly FASTA files gathering CRISPRs "
        "in chromosomes into CSVs.."
    )
    sequence_start = 0
    for chromosome in CHROMOSOMES:
        print("Chromosome %s" % chromosome)
        genome_file = "%s.dna.chromosome.%s.fa" % (file_name, chromosome)
        subprocess.run(
            [
                "gunzip",
                "--keep",
                f"{ot_path}/genome_files/{genome_file}.gz",
            ]
        )
        sequence_start = gather.gather(
            inputfile=f"{ot_path}/genome_files/{genome_file}",
            outputfile="%s/crispr_csv/%s.chromosome.%s.%s.csv"
            % (ot_path, file_name, chromosome, pam_sequence),
            pam=pam_sequence,
            sequence_start=sequence_start,
        )

    try:
        os.mkdir(f"{ot_path}/grna_bin")
    except FileExistsError:
        pass

    # create binary index file of CRISPRs
    print("\nCreating binary index of gRNAs from CSVs...")
    chromosome_input_text_list = list()
    for chromosome in CHROMOSOMES:
        chromosome_input_text_list.append(
            "%s/crispr_csv/%s.chromosome.%s.%s.csv"
            % (ot_path, file_name, chromosome, pam_sequence)
        )
    index.index(
        inputfiles=chromosome_input_text_list,
        outputfile=f"{ot_path}/grna_bin/{file_name}.{pam_sequence}.bin",
        species="Human",
        assembly=ensembl_version,
        offset=0,
        species_id=1,
    )
    print("CRISPRs indexed\n")


def init_db(db_file: str, ot_path: str) -> None:
    """Initialise the SQLite database with the crisprs table."""
    try:
        os.mkdir(f"{ot_path}/crispr_db")
    except FileExistsError:
        pass
    print("Creating CRISPR database...")
    subprocess.run(
        [
            "sqlite3",
            f"{ot_path}/crispr_db/{db_file}",
            "CREATE TABLE IF NOT EXISTS crisprs (id INT NOT NULL PRIMARY KEY, "
            "chr_name TEXT NOT NULL, chr_start INT NOT NULL, "
            "seq TEXT NOT NULL, "
            "pam_right INT NOT NULL CHECK (pam_right in (0, 1)))",
        ]
    )


def import_crisprs_to_db(
    assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str
) -> None:
    """Import CRISPRs in CSVs into the 'crisprs' table"""
    # Check to see if sqlite3 is available
    if shutil.which("sqlite3") is None:
        sys.exit("please install sqlite3 before proceeding")

    file_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )
    db_file = f"{file_name}.{pam_sequence}.db"

    # index the database with CRISPRs gathered in the CSV files
    init_db(db_file=db_file, ot_path=ot_path)
    print("Importing CRISPR data to db...")
    for chromosome in CHROMOSOMES:
        print(
            "importing %s.chromosome.%s.%s.csv"
            % (
                file_name,
                chromosome,
                pam_sequence,
            )
        )
        subprocess.run(
            [
                "sqlite3",
                "%s/crispr_db/%s" % (ot_path, db_file),
                ".import --csv %s/crispr_csv/%s.chromosome.%s.%s.csv crisprs"
                % (ot_path, file_name, chromosome, pam_sequence),
            ]
        )
    print("CRISPR data imported\n")


def check_grna_bin_exists(
    assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str
) -> bool:
    """Check to see if gRNA binary index file exists"""
    file_path = f"{ot_path}/grna_bin"
    if os.path.exists(file_path) is False:
        return False
    file_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )
    if "%s.%s.bin" % (file_name, pam_sequence) not in os.listdir(file_path):
        return False
    else:
        return True


def check_crispr_db_exists(
    assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str
) -> bool:
    """Check to see if CRISPR SQLite database file exists"""
    file_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )
    db_file = f"{file_name}.{pam_sequence}.db"
    return os.path.exists(f"{ot_path}/crispr_db/{db_file}")


def check_crispr_csvs_exist(
    assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str
) -> bool:
    """Check to see if all the indexed CRISPR csv files exist"""
    file_name = base_file_name(
        assembly=assembly, ensembl_version=ensembl_version
    )
    file_path = f"{ot_path}/crispr_csv"
    if os.path.exists(file_path) is False:
        return False
    files_exist = True
    for chromosome in CHROMOSOMES:
        if (
            f"{file_name}.chromosome.{chromosome}.{pam_sequence}.csv"
            not in os.listdir(file_path)
        ):
            files_exist = False
    return files_exist


def check_crispr_indexes_exist(
    assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str
) -> bool:
    """Check to see if both gRNA index and CRISPR CSVs exist"""
    return check_grna_bin_exists(
        assembly=assembly,
        ensembl_version=ensembl_version,
        pam_sequence=pam_sequence,
        ot_path=ot_path,
    ) & check_crispr_csvs_exist(
        assembly=assembly,
        ensembl_version=ensembl_version,
        pam_sequence=pam_sequence,
        ot_path=ot_path,
    )


################################################################################
# Execution


def run(assembly: str, ensembl_version: str, pam_sequence: str, ot_path: str):
    """Run all CRISPR and gRNA indexing from genome assembly FASTA files"""
    if os.path.exists(ot_path) is False:
        os.mkdir(ot_path)
    # Check if we have downloaded the genome files already
    if (
        check_genome_files_exist(
            assembly=assembly,
            ensembl_version=ensembl_version,
            ot_path=ot_path,
        )
        is False
    ):
        fetch_genome_files(
            assembly=assembly,
            ensembl_version=ensembl_version,
            ot_path=ot_path,
        )
        if (
            check_genome_files_exist(
                assembly=assembly,
                ensembl_version=ensembl_version,
                ot_path=ot_path,
            )
            is False
        ):
            sys.exit(
                "Error in downloading genome, please manually download all "
                "chromosomes from: "
                "https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/dna/"
                " named Homo_sapiens.GRCh38.dna.chromosome.<chromosome>.fa.gz "
                "if the assembly is GRCh38, otherwise "
                "Homo_sapiens.GRCh37.<version>"
                ".dna.chromosome.<chromosome>.fa.gz"
            )
    else:
        print("Genome FASTA files exist - nothing to be done")

    if (
        check_crispr_indexes_exist(
            assembly=assembly,
            ensembl_version=ensembl_version,
            pam_sequence=pam_sequence,
            ot_path=ot_path,
        )
        is False
    ):
        gather_crisprs_from_genome(
            assembly=assembly,
            ensembl_version=ensembl_version,
            pam_sequence=pam_sequence,
            ot_path=ot_path,
        )
    else:
        print("CRISPR CSV files exist - nothing to be done")
    if (
        check_crispr_db_exists(
            assembly=assembly,
            ensembl_version=ensembl_version,
            pam_sequence=pam_sequence,
            ot_path=ot_path,
        )
        is False
    ):
        import_crisprs_to_db(
            assembly=assembly,
            ensembl_version=ensembl_version,
            pam_sequence=pam_sequence,
            ot_path=ot_path,
        )
    else:
        print("CRISPR database exists - nothing to be done")


if __name__ == "__main__":
    # -------------------------------------------------------------------------#
    # Retrieve command line arguments
    args = take_input()

    # Off-targets Path (without trailing backslash)
    ot_path = (
        args.offtargets_path
        if args.offtargets_path[-1] != "/"
        else args.offtargets_path[:-1]
    )

    run(
        assembly=args.assembly,
        ensembl_version=args.ensembl_version,
        pam_sequence=args.pamseq,
        ot_path=ot_path,
    )
