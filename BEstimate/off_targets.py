# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas, argparse, re, pickle
import itertools


# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate Off targets",
									 usage="%(prog)s [inputs]",
									 description="""
                                     **********************************
                                     Find and Analyse Base Editor sites
                                     		 Off Target Analysis
                                     **********************************""")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# PAM AND PROTOSPACER INFORMATION

	# The NGG PAM will be used unless otherwise specified.
	parser.add_argument("-path", dest="PATH")
	parser.add_argument("-file", dest="FILE")
	parser.add_argument("-mm", dest="MM", default="4")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


args = take_input()

ot_path = "/Volumes/team215/Cansu/BEstimate/off_targets/"

# -----------------------------------------------------------------------------------------#
# Take gRNAs

def mm_combination_seq(mm, seq):
	"""
	Mismatched version of gRNAs
	:param mm:
	:param seq:
	:return:
	"""
	nuc_dict = {"A": ["T", "C", "G"],
				"T": ["A", "C", "G"],
				"C": ["A", "T", "G"],
				"G": ["A", "T", "C"]}
	pam = seq[-3:]
	seq = seq[:-3]
	if mm == 0: return seq
	else:
		ind_perm = list(itertools.permutations(list(range(len(seq))), mm))
		all_seqs = list()
		for pos in ind_perm:
			new_seq = []
			for ind in range(len(seq)):
				if ind not in pos:
					new_seq.append(list(seq[ind]))
				else:
					new_seq.append(list(nuc_dict[seq[ind]]))
			all_seqs.append(new_seq)

		all_mm_guides = list()
		for mm_seq in all_seqs:
			for s in list(itertools.product(*mm_seq)):
				all_mm_guides.append("".join(s) + pam)

		return list(set(all_mm_guides))


def grna_fasta(myc_guides, file, mm):
	#crisprs = pandas.read_csv("/Volumes/team215/Cansu/BEstimate/Matt/ABE_NGN_AKT1_39_summary_df.csv", index_col=0)["CRISPR_PAM_Sequence"].values[:100]
	#crisprs = pandas.read_csv(path + file, index_col=0)["CRISPR_PAM_Sequence"].values
	#crisprs = ["AAACGGGGCCATCTGTCACCAGG", "AAAGACGTTTTTGTGCTGTGGGC", "AAAGTTGCTTTTCAAATTTTTGG", "AAATTTGTTATTGTGTATTATGT","AAAAAAACGCCGTGGTGCAGCGG", "AAAAAACGCCGTGGTGCAGCGGC", "AAAAACCCCCAAAATGCATTTGA", "AAAAAGCTTCTCATGGTCCTGGT"]
	crisprs = myc_guides
	mm = int(mm)
	f = open(ot_path + "fasta/%s.fa" % "_".join(file.split(".")[0].split("_")[:-2]), "w")
	grna_dict = dict()
	count = 1
	for i in crisprs:
		seq = i[0]
		direction = i[1]
		f.write(">gRNA-%d-%s\n%s\n" %(count, direction, seq))
		if "gRNA-%d-%s" % (count, direction) not in grna_dict.keys():
			grna_dict["gRNA-%d-%s" %(count, direction)] = {"seq": seq}
		for k in range(1, mm + 1):
			mm_count = 1
			mm_seqs = mm_combination_seq(mm=k, seq=seq)
			for m_seq in mm_seqs:
				f.write(">gRNA-%d-%s:mm%d:count%s\n%s\n" % (count, direction, k, mm_count, m_seq))
				grna_dict["gRNA-%d-%s" % (count, direction)]["mm%d:mm_count%d" % (k, mm_count)] = m_seq
				mm_count += 1
		count += 1
		print(count)

	f.close()
	pickle.dump(grna_dict, open(ot_path + "fasta_dict/" + "_".join(file.split(".")[0].split("_")[:-2]) + "_dict.p", "wb"))
	return 1


def index_genome(genome, vcf):

	return 1

def run_mrsfast(fasta, threads_num):
	"""
	Example in cluster systems
	os.system("bsub -G team215-grp -R'select[mem>20000] rusage[mem=20000]' -M20000 -n 8 "
		  "-o /path/logs/X.o "
		  "-e /path/logs/X.e "
		  "mrsfast --search /path/genome/genome_38.fa --seq path/fasta/X.fa --threads 8 -e 0 "
		  "-o /path/sam_files/X")
	:return:
	"""
	os.system("mrsfast --search genome/genome_38.fa --seq fasta/%s.fa --threads %s -e 0"
			  "-o sam_files/%s_alignment.sam" % (fasta, threads_num, fasta))
	return 1


def read_sam(sam_file):
	ot_df = pandas.DataFrame(columns=["guide_flag", "guide_seq", "direction", "mm", "chr", "start", "mm_count"])
	f = open(ot_path + "sam_files/" + sam_file, "r")
	lines = f.readlines()
	n, t = 0, len(lines)
	ind_count = 0
	for line in lines:
		if line[:4] == "gRNA":
			print(line)
			l = line.split("\t")
			if len(l[0].split(":")) == 1:
				guide_flag = l[0]
				mm = "0"
				mm_count = None
			else:
				guide_flag = l[0].split(":")[0]
				mm = l[0].split(":")[1][2:]
				mm_count = l[0].split(":")[2][5:]

			ot_df.loc[ind_count, "guide_flag"] = guide_flag
			ot_df.loc[ind_count, "guide_seq"] = l[9]
			ot_df.loc[ind_count, "direction"] = guide_flag.split("-")[-1]
			ot_df.loc[ind_count, "mm"] = mm
			ot_df.loc[ind_count, "chr"] = l[2]
			ot_df.loc[ind_count, "start"] = l[3]
			ot_df.loc[ind_count, "mm_count"] = mm_count
			ind_count += 1
		print(n * 100.0 /t)
		n += 1
	f.close()

	chr_list = [str(i) for i in list(range(1, 23))] + ["X", "Y", "MT"]
	ot_df2 = ot_df[ot_df.chr.isin(chr_list)]

	summary_ot_df = pandas.DataFrame(columns = ["guide_seq", "direction", "exact", "mm1", "mm2", "mm3", "mm4"],
									 index = list(ot_df2.groupby(["guide_flag"]).groups.keys()))

	for g, g_df in ot_df2.groupby(["guide_flag"]):
		exact = g_df.groupby(["mm"]).size().loc["0"]
		summary_ot_df.loc[g, "guide_seq"] = g_df[g_df.mm == "0"]["guide_seq"].unique()[0]
		summary_ot_df.loc[g, "exact"] = exact
		for i in range(1, 5):
			if str(i) in g_df.mm.unique():
				summary_ot_df.loc[g, "mm%d" % i] = g_df.groupby(["mm"]).size().loc[str(i)]
			else:
				summary_ot_df.loc[g, "mm%d" % i] = 0

	summary_ot_df = summary_ot_df.reset_index()
	summary_ot_df.columns = ["guide_flag", "CRISPR_PAM_Sequence", "exact", "mm1", "mm2", "mm3", "mm4"]


	summary_ot_df2 = pandas.DataFrame(columns = ["guide_seq", "exact", "mm1", "mm2", "mm3", "mm4"],
									  index = list(ot_df.groupby(["guide_flag"]).groups.keys()))

	for g, g_df in ot_df.groupby(["guide_flag"]):
		exact = g_df.groupby(["mm"]).size().loc["0"]
		summary_ot_df2.loc[g, "guide_seq"] = g_df[g_df.mm == "0"]["guide_seq"].unique()[0]
		summary_ot_df2.loc[g, "exact"] = exact
		for i in range(1, 5):
			if str(i) in g_df.mm.unique():
				summary_ot_df2.loc[g, "mm%d" % i] = g_df.groupby(["mm"]).size().loc[str(i)]
			else:
				summary_ot_df2.loc[g, "mm%d" % i] = 0

	summary_ot_df2 = summary_ot_df2.reset_index()
	summary_ot_df2.columns = ["guide_flag", "CRISPR_PAM_Sequence", "exact", "mm1", "mm2", "mm3", "mm4"]

	crisprs = pandas.read_csv("/Volumes/team215/Cansu/BEstimate/Matt/ABE_NGN_AKT1_39_summary_df.csv", index_col=0)

	crisprs2 = pandas.merge(crisprs, summary_ot_df, how="inner", on=["CRISPR_PAM_Sequence"])



def read_unread(file):

	not_aligned = open(ot_path + "sam_files/AKT1_ABE_alignment.nohit", "r")
na_lines = not_aligned.readlines()
na_df = pandas.DataFrame(columns = ["guide_flag", "sequence"])
count = 0
for line in na_lines:
	if line[0] == ">":
		na_df.loc[count, "guide_flag"] = line.strip()[0]
	else:
		na_df.loc[count, "sequence"] = line.strip()[0]
	count += 1




# Benchmarking MYC gene with ENSE00003746860


def check_NGG(grna):
	if grna[-2:] == "GG":
		return True
	else:
		return False

all_myc = pandas.read_csv("/Volumes/team215/Cansu/BEstimate/off_targets/benchmarking/NGN_CBE_MYC_crispr_df.csv", index_col=0)
exon_myc = all_myc[all_myc.Exon_ID == "ENSE00003746860"]
exon_myc["NGG"] = exon_myc.apply(lambda x: check_NGG(x.CRISPR_PAM_Sequence), axis=1)
ngg_exon_myc = exon_myc[exon_myc.NGG]
myc_guides = list(ngg_exon_myc[["CRISPR_PAM_Sequence", "Direction"]].values)
summary_ot_df2

cf_myc = pandas.read_csv("/Volumes/team215/Cansu/BEstimate/off_targets/benchmarking/WGE-ENSE00003746860-crisprs.tsv", sep="\t")
cf_myc = cf_myc[["seq", "location", "strand", "off_target_summary"]]
cf_myc["exact"] = cf_myc.apply(lambda x: x.off_target_summary.split(", ")[0].split(":")[1], axis=1)
cf_myc["mm1"] = cf_myc.apply(lambda x: x.off_target_summary.split(", ")[1].split(":")[1], axis=1)
cf_myc["mm2"] = cf_myc.apply(lambda x: x.off_target_summary.split(", ")[2].split(":")[1], axis=1)
cf_myc["mm3"] = cf_myc.apply(lambda x: x.off_target_summary.split(", ")[3].split(":")[1], axis=1)
cf_myc["mm4"] = cf_myc.apply(lambda x: x.off_target_summary.split(", ")[4].split(":")[1], axis=1)















