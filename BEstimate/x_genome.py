# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                           Genome Retrieval and Indexing                                  #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                                                                                          #
# -----------------------------------------------------------------------------------------#


import argparse
import os
import requests
import shutil
import subprocess
import sys
import time
from crispr_analyser import index, gather


CHROMOSOMES = list(range(1, 23)) + ["X", "Y", "MT"]

# Extracting Humen Reference Genome

###########################################################################################
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
        help="The path for output. If not specified the current directory will be used!",
    )

    parser.add_argument(
        "--ensemble_version",
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


###########################################################################################
# Functions


def file_main_text(assembly: str, ensemble_version: str) -> str:
    """Generate the name of the saved FASTA files."""
    name = ""
    if assembly == "GRCh37":
        name = "Homo_sapiens.GRCh37.%s.dna.chromosome" % ensemble_version
    elif assembly == "GRCh38":
        name = "Homo_sapiens.GRCh38.dna.chromosome"
    return name


def check_genome_files_exist(
    assembly: str, ensemble_version: str, ot_path: str
) -> bool:
    """Check to see if the genome FASTA files have been dowloaded."""
    file_directory = f"{ot_path}/genome_files"
    if os.path.exists(file_directory) is False:
        return False
    files_exist = True
    file_name = file_main_text(assembly=assembly, ensemble_version=ensemble_version)
    for chromosome in CHROMOSOMES:
        if "%s.%s.fa.gz" % (file_name, chromosome) not in os.listdir(file_directory):
            files_exist = False
    return files_exist


def fetch_genome_files(assembly: str, ensemble_version: str, ot_path: str) -> None:
    """Download the genome FASTA files from Ensembl"""
    print(
        "Genome is not found, BEstimate is downloading the %s Ensembl genome - version %s\n"
        % (assembly, ensemble_version)
    )
    try:
        os.mkdir(f"{ot_path}/genome_files")
    except FileExistsError:
        pass

    file_name = file_main_text(assembly=assembly, ensemble_version=ensemble_version)

    for chromosome in CHROMOSOMES:
        url = (
            "https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/dna/%s.%s.fa.gz"
            % (ensemble_version, file_name, chromosome)
        )
        response = requests.get(url)
        file_path = "%s/genome_files/%s.%s.fa.gz" % (
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
    assembly: str, ensemble_version: str, pam_sequence: str, ot_path: str
) -> bool:
    file_name = file_main_text(assembly=assembly, ensemble_version=ensemble_version)

    if "%s.bin" % file_name not in os.listdir(f"{ot_path}/genome_files/"):
        chromosome_input_text_list = list()
        for chromosome in CHROMOSOMES:
            chromosome_input_text_list.append(
                "%s/crispr_csv/c_%s.csv" % (ot_path, chromosome)
            )

        # Gather all chromosome fasta files into csv files
        print("Gathering CRISPRs in chromosomes into CSVs..\n")
        sequence_start = 0
        for chromosome in CHROMOSOMES:
            genome_file = "%s.%s.fa" % (file_name, chromosome)
            if "c1_%s.csv" % chromosome not in os.listdir(f"{ot_path}/crispr_csv/"):
                if file_name not in os.listdir("%s/genome_files/" % ot_path):
                    subprocess.run(
                        ["gunzip", "--keep", f"{ot_path}/genome_files/{genome_file}.gz"]
                    )
                print("Chromosome %s" % chromosome)
                sequence_start = gather.gather(
                    inputfile=f"{ot_path}/genome_files/{genome_file}",
                    outputfile=f"{ot_path}/crispr_csv/c_{chromosome}.csv",
                    pam=pam_sequence,
                    sequence_start=sequence_start,
                )

                while "c_%s.csv" % chromosome not in os.listdir(
                    "%s/crispr_csv/" % ot_path
                ):
                    print("Waiting chromosome %s.." % chromosome)
                    time.sleep(10)

        try:
            os.mkdir(f"{ot_path}/crispr_bin")
        except FileExistsError:
            pass

        # create binary index file of CRISPRs
        print("Creating binary index of gRNAs...\n")
        index.index(
            inputfiles=chromosome_input_text_list,
            outputfile=f"{ot_path}/crispr_bin/{file_name}.bin",
            species="Human",
            assembly=ensemble_version,
            offset=0,
            species_id=1,
        )

    if "%s.bin" % file_name in os.listdir("%s/genome/" % ot_path):
        return True
    else:
        return False


def init_db(db_file: str, ot_path: str) -> None:
    """Initialise the SQLite database with the crisprs table."""
    try:
        os.mkdir(f"{ot_path}/crispr_db")
    except FileExistsError:
        pass
    subprocess.run(
        [
            "sqlite3",
            f"{ot_path}/crispr_db/{db_file}",
            "CREATE TABLE IF NOT EXISTS crisprs (id INT NOT NULL PRIMARY KEY, chr_name TEXT NOT NULL, chr_start INT NOT NULL, seq TEXT NOT NULL, pam_right INT NOT NULL CHECK (pam_right in (0, 1)))",
        ]
    )


def import_csvs_to_db(db_file: str, ot_path: str) -> None:
    # Check to see if sqlite3 is available
    if shutil.which("sqlite3") is None:
        sys.exit("please install sqlite3 before proceeding")

    # index the database with CRISPRs gathered in the CSV files
    init_db(db_file=db_file, ot_path=ot_path)

    for chromosome in CHROMOSOMES:
        subprocess.run(
            [
                "sqlite3",
                f"{ot_path}/crispr_db/{db_file}",
                f".import --csv {ot_path}/crispr_csv/csv_{chromosome}.csv crisprs",
            ]
        )


def check_crispr_bin_exists(assembly: str, ensemble_version: str, ot_path: str) -> bool:
    """Check to see if CRISPR binary index file exists"""
    file_name = file_main_text(assembly=assembly, ensemble_version=ensemble_version)
    if "%s.bin" % file_name not in os.listdir("%s/genome_files/" % ot_path):
        return False
    else:
        return True


def check_crispr_db_exists(ot_path: str) -> bool:
    """Check to see if CRISPR SQLite database file exists"""
    if "crisprs.db" in os.listdir(f"{ot_path}/crispr_db"):
        return True
    else:
        return False


def check_crispr_csvs_exist(ot_path: str) -> bool:
    """Check to see if all the indexed CRISPR csv files exist"""
    file_path = f"{ot_path}/crispr_csv"
    if os.path.exists(file_path) is False:
        return False
    files_exist = True
    for chromosome in CHROMOSOMES:
        if f"c_{chromosome}.csv" not in os.listdir(file_path):
            files_exist = False
    return files_exist


def check_crispr_indexes_exist(
    assembly: str, ensemble_version: str, ot_path: str
) -> bool:
    return check_crispr_bin_exists(
        assembly, ensemble_version, ot_path
    ) & check_crispr_csvs_exist(ot_path=ot_path)


def create_index(ot_path: str) -> bool:
    try:
        os.mkdir(f"{ot_path}/crispr_csv")
    except FileExistsError:
        pass

    is_index = check_crispr_bin_exists(
        assembly=args.assembly, ensemble_version=args.ensemble_version, ot_path=ot_path
    )

    if is_index:
        return True
    else:
        print("Indexing CRISPRS from the genome..")
        result = gather_crisprs_from_genome(
            assembly=args.assembly,
            ensemble_version=args.ensemble_version,
            pam_sequence=args.pamseq,
            ot_path=ot_path,
        )
        return result


###########################################################################################
# Execution


def run(assembly: str, ensemble_version: str, ot_path: str):
    """Run from command line."""
    if os.path.exists(ot_path) is False:
        os.mkdir(ot_path)
    # Check if we have downloaded the genome files already
    if (
        check_genome_files_exist(
            assembly=assembly, ensemble_version=ensemble_version, ot_path=ot_path
        )
        is False
    ):
        fetch_genome_files(
            assembly=assembly, ensemble_version=ensemble_version, ot_path=ot_path
        )
        # If we tried to download genome files and we still don't have the files with error
        if (
            check_genome_files_exist(
                assembly=assembly, ensemble_version=ensemble_version, ot_path=ot_path
            )
            is False
        ):
            sys.exit(
                "Error in downloading genome, please manually downloading all chromosomes from:"
                "https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/dna/ named as "
                "Homo_sapiens.GRCh38.dna.chromosome.<chromosome>.fa.gz if the assembl is GRCh38, otherwise"
                "Homo_sapiens.GRCh37.<version>.dna.chromosome.<chromosome>.fa.gz"
            )

    if (
        check_crispr_indexes_exist(
            assembly=assembly, ensemble_version=ensemble_version, ot_path=ot_path
        )
        is False
    ):
        # check that bin & CSVs exist (crispr_index) - generate the CSV & bin at same time
        # should also check that we have generated the CSV files - put index of bin

        # should we check that number of lines of CSV, bin and db match?
        create_index(ot_path)

    if check_crispr_db_exists(ot_path=ot_path) is False:
        import_csvs_to_db(db_file="crisprs.db", ot_path=ot_path)


if __name__ == "__main__":
    # -----------------------------------------------------------------------------------------#
    # Retrieve command line arguments
    args = take_input()

    # Off-targets Path (without trailing backslash)
    ot_path = (
        args.offtargets_path
        if args.offtargets_path[-1] != "/"
        else args.offtargets_path[:-1]
    )

    run(assembly=args.assembly, ensemble_version=args.ensemble_version, ot_path=ot_path)
