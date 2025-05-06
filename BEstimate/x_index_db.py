# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                          CRISPR-Analyser | Genome Indexing                               #
#                        Author : Bo Fussing bf15@sanger.ac.uk                             #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

import sqlite3
import pandas as pd

###########################################################################################
# Functions

def index(inputfiles: list[str]) -> None:
    con = sqlite3.connect("crispr_db", isolation_level=None)
    con.execute("PRAGMA temp_store=MEMORY")
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA cache_size=1000000")
    con.execute("PRAGMA locking_mode=EXCLUSIVE")
    con.execute("CREATE TABLE IF NOT EXISTS crisprs (id INT NOT NULL PRIMARY KEY, chr_name TEXT NOT NULL, chr_start INT NOT NULL, seq TEXT NOT NULL, pam_right INT NOT NULL CHECK (pam_right in (0, 1)))")
    cur = con.cursor()
    sequence_id = 0
    for inputfile in inputfiles:
        print(f"Importing CRISPRs from {inputfile}")
        for chunk in pd.read_csv(inputfile, chunksize=500_000, header=None):
            data = []
            for _, row in chunk.iterrows():
                sequence_id += 1
                if len(row) != 5:
                    print(f"Invalid line: {row}")
                    continue
                data.append((sequence_id, row[0], row[1], row[2], row[3]))
            cur.executemany("INSERT INTO crisprs VALUES (?, ?, ?, ?, ?)", data)

    cur.close()
    con.close()
