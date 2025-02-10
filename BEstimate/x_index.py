# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                          CRISPR-Analyser | Genome Indexing                               #
#                        Author : Bo Fussing bf15@sanger.ac.uk                             #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

import getopt, sqlite3, sys, time


###########################################################################################
# Functions

def usage():
    print("x_index.py -i <input CSV file> -d <database file> -f <offset>")


def main(argv):
    start = time.time()
    inputfiles = []
    dbfile = ""
    offset = 0
    number_of_sequences = 0
    sequence_id = 0

    try:
        opts, args = getopt.getopt(
            argv,
            "hi:d:f:",
            [
                "help",
                "ifile=",
                "dbfile=",
                "offset=",
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
        elif opt in ("-d", "--dbfile"):
            dbfile = arg
        elif opt in ("-f", "--offset"):
            offset = int(arg)
    if inputfiles == [] or dbfile == "":
        usage()
        sys.exit(2)
    if offset < 0:
        print("Offset must be a positive integer")
        sys.exit(2)
    if offset:
        sequence_id = offset

    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    cur.execute(
        "CREATE TABLE IF NOT EXISTS crisprs (id INT NOT NULL PRIMARY KEY, chr_name TEXT NOT NULL, chr_start INT NOT NULL, seq TEXT NOT NULL, pam_right INT NOT NULL CHECK (pam_right in (0, 1)))"
    )
    for inputfile in inputfiles:
        print(f"Processing {inputfile}")
        with open(inputfile, "r") as in_file:
            for line in in_file:
                line = line.strip()
                if line == "":
                    continue
                sequence_id += 1
                parts = line.split(",")
                if len(parts) != 5:
                    print(f"Invalid line: {line}")
                    continue
                cur.execute(
                    "INSERT INTO crisprs VALUES (?, ?, ?, ?, ?)",
                    (sequence_id, parts[0], parts[1], parts[2], parts[3]),
                )
                number_of_sequences += 1
            con.commit()


if __name__ == "__main__":
    main(sys.argv[1:])
