# -----------------------------------------------------------------------------------------#
# 						BEstimate VEP and Matt Coelho comaprison						   #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas, re
path = "/Users/cd7/Desktop/BEstimate_Matt/VEP_comparison/"

aa_3to1 = {"Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F", "Gly": "G",
		   "His": "H", "Ile": "I", "Lys": "K", "Leu": "L", "Met": "M", "Asn": "N",
		   "Pro": "P", "Gln": "Q", "Arg": "R", "Ser" : "S", "Thr": "T", "Val": "V",
		   "Trp": "W", "Tyr": "Y", "Asx": "B", "Glx": "Z", "Ter": "*"}
# -----------------------------------------------------------------------------------------#
# Matt File

matt = pandas.read_csv(path + "input.csv", index_col = 0)

# Only JAK1 NGG

jak1_ngg_matt = matt[(matt.Gene == "JAK1") & (matt.editor == "BE3-NGG (JAK1)")][[
	"guide", "Amino_Acid_Position", "Amino_Acid_Change", "Consequence"]]

matt_comp_df = pandas.DataFrame(columns = ["guide", "position", "aa", "new_aa", "consq"])

i, t = 0, len(jak1_ngg_matt.index)
for ind, row in jak1_ngg_matt.iterrows():
	matt_comp_df.loc[i, "guide"] = row.guide
	if pandas.isna(row.Amino_Acid_Position):
		matt_comp_df.loc[i, "position"] = None
	else:
		if len(row.Amino_Acid_Position.split(", ")) == 1:
			matt_comp_df.loc[i, "position"] = row.Amino_Acid_Position
		else:
			for pos in row.Amino_Acid_Position.split(", "):
				matt_comp_df.loc[i, "position"] = ";".join([k for k in row.Amino_Acid_Position.split("; ")])
	if pandas.isna(row.Amino_Acid_Change):
		matt_comp_df.loc[i, "aa"] = None
		matt_comp_df.loc[i, "new_aa"] = None
	else:
		if len(row.Amino_Acid_Change.split(", ")) == 1:
			if len(re.findall("delins", row.Amino_Acid_Change)) == 0:
				if len(re.findall("\?", row.Amino_Acid_Change)) == 0:
					matt_comp_df.loc[i, "aa"] = aa_3to1[row.Amino_Acid_Change[:3]]
					matt_comp_df.loc[i, "new_aa"] = aa_3to1[row.Amino_Acid_Change[-3:]]
				else:
					matt_comp_df.loc[i, "aa"] = aa_3to1[row.Amino_Acid_Change[:3]]
					matt_comp_df.loc[i, "new_aa"] = "None"
		else:
			aas, new_aas = list(), list()
			for ac in row.Amino_Acid_Change.split(", "):
				if len(re.findall("delins", ac)) == 0:
					if len(re.findall("\?", row.Amino_Acid_Change)) == 0:
						aas.append(aa_3to1[ac[:3]])
						new_aas.append(aa_3to1[ac[-3:]])
					else:
						aas.append(aa_3to1[ac[:3]])
						new_aas.append("None")
			matt_comp_df.loc[i, "aa"] = ";".join(aas)
			matt_comp_df.loc[i, "new_aa"] = ";".join(new_aas)
	matt_comp_df.loc[i, "consq"] = row.Consequence
	i += 1
	print(i * 100.0/t)

matt_comp_df_nan = matt_comp_df.dropna(subset=["position"])
matt_comp_df.to_csv(path+ "matt_df.csv")


# BEstimate file

def tidy_consequence(bestimate_cons):
	cons = list()
	if bestimate_cons is not None or type(bestimate_cons) != float:
		for c in bestimate_cons.split(";"):
			if len(re.findall("stop", c)) != 0:
				cons.append("stop codon")
			elif len(re.findall("splice", c)) != 0:
				cons.append("splice variant")
			elif len(re.findall("missense", c)) != 0:
				cons.append("missense")
			elif len(re.findall("UTR", c)) != 0:
				cons.append("UTR")
			elif len(re.findall("synonymous", c)) != 0:
				cons.append("synonymous")

	if "stop codon" in cons:
		return "stop codon"
	elif "missense" in cons:
		return "missense"
	elif "splice variant" in cons:
		return "splice variant"
	elif "UTR" in cons:
		return "UTR"
	elif "synonymous" in cons:
		return "synonymous"
	else:
		return None


bestimate = pandas.read_csv(path + "JAK1_NGG_CT_vep_df.csv", index_col = 0)
bestimate_comp_df = pandas.DataFrame(columns = ["guide", "position", "aa", "new_aa", "consq"])

"""
ORDER

stop codon
start lost
splice variant
missense
UTR
synonymous
"""

i, t = 0, len(bestimate.index)
for ind, row in bestimate.iterrows():
	bestimate_comp_df.loc[i, "guide"] = row.gRNA_Target_Sequence
	bestimate_comp_df.loc[i, "position"] = row.Protein_Position
	bestimate_comp_df.loc[i, "aa"] = row.Edited_AA
	bestimate_comp_df.loc[i, "new_aa"] = row.New_AA
	bestimate_comp_df.loc[i, "consq"] = tidy_consequence(row.consequence_terms)
	i += 1
	print(i*100.0/t)

bestimate_comp_df_nan = bestimate_comp_df.dropna(subset=["position"])
bestimate_comp_df.to_csv(path+ "bestimate_df.csv")


guides = list(set(bestimate_comp_df_nan.guide).intersection(set(matt_comp_df_nan.guide)))
matt_df = matt_comp_df_nan[matt_comp_df_nan.guide.isin(guides)]
bestimate_df = bestimate_comp_df_nan[bestimate_comp_df_nan.guide.isin(guides)]

for i,l in matt_df.groupby(["guide"]):
	check = False
	if i in list(bestimate_df.guide):
		x = bestimate_df[bestimate_df.guide == i]
		if len(l.index) == 1 and len(x.index) == 1:
			if l.position == x.position and l.position == x.position and l.new_aa == x.new_aa and l.consq == x.consq:
				check = True
		elif len(l.index) == 1:
			for ind, r in x.iterrows():


