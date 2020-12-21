import os, pandas

###########################################################################################
# WGE Web

wge_path = os.getcwd() + "/website/"

crispr_files = os.listdir(wge_path)

dfs = []
for file in crispr_files:
    df = pandas.read_csv(wge_path + file, sep ="\t")[["location", "strand", "seq"]]
    df = df.set_index(["seq"])
    df["Location"] = df.apply(lambda x: x.location.split(":")[1].replace("-", ":"), axis=1)
    df["Strand"] = df.apply(lambda x: "right" if x.strand == "+" else "left", axis=1)
    df = df[["Location", "Strand"]]
    dfs.append(df)


CRISPR_Finder_df = pandas.concat(dfs)


script_df = pandas.read_csv(os.getcwd() + "/CBE_NGG_GRCh38_edit_df.csv")[[
    "CRISPR Sequence", "CRISPR Location", "Direction"]]

script_df = script_df.set_index(["CRISPR Sequence"])

script_df.columns = ["Location", "Strand"]

x = CRISPR_Finder_df.merge(script_df, indicator = True, how='left')

crispr_finder = x[x["_merge"] =="left_only"]
bestimate = x[x["_merge"] =="right_only"]

# Check if CRISPR Finder CRISPRs have editable nucleotide inside activity window
count = 0
for i, row in crispr_finder.iterrows():
    if row.Strand == "right":
        if "C" in left[3:8]:
            count +=1
            print(left)
    else:
        if "G" in left[-8:-3]:
            count +=1
            print(left)

###########################################################################################
