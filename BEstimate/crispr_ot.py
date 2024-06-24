# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                 BEstimate Offtargets                                     #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, sys, pandas, argparse, requests, itertools


# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate Offtargets",
									 usage="%(prog)s [inputs]")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# BASIC INFORMATION

	parser.add_argument("-gene", dest="GENE", required=True,
						help="Gene name")

	parser.add_argument("-guide", dest="GRNA", required=True,
						help="gRNA sequence")

	parser.add_argument("-direction", dest="DIRECTION", required=True,
						help="The direction of the gRNA")

	# OFF TARGETS
	parser.add_argument("-mm", dest="MISMATCH", default=4,
						help="(If -ot provided) number of maximum mismatches allowed in off targets")
	parser.add_argument("-genome", dest="GENOME", default="Homo_sapiens_GRCh38_dna_sm_all_chromosomes",
						help="(If -ot provided) name of the genome file")
	parser.add_argument("-num", dest="NUM", default="1")

	# PATH

	parser.add_argument("-path", dest="PATH", default=os.getcwd() + "/",
						help="The path, if not specified the current directory will be used!")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict



# -----------------------------------------------------------------------------------------#
# Functions


def mm_combination_seq(mm, seq):

	nuc_dict = {"A": ["T", "C", "G"],
				"T": ["A", "C", "G"],
				"C": ["A", "T", "G"],
				"G": ["A", "T", "C"]}
	pam = seq[-3:]
	seq = seq[:-3]
	if mm == 0:
		return seq
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


def grna_fasta(gene, guide, direction, mm, path):

	try:
		os.mkdir("%sfasta_dict/" % path)
	except FileExistsError:
		pass

	if "%s_%s_%s_dict.p" % (gene, guide, direction) not in os.listdir(path + "fasta_dict/"):
		mm = int(mm)
		f = open(path + "fasta/%s_%s_%s_fasta.fa" % (gene, guide, direction), "w")
		grna_dict = dict()

		f.write(">gRNA-%s-%s\n%s\n" % (gene, direction, guide))
		if "gRNA-%s-%s" % (gene, direction) not in grna_dict.keys():
			grna_dict["gRNA-%s-%s" % (gene, direction)] = {"seq": guide}

		for k in range(1, mm + 1):
			mm_count = 1
			mm_seqs = mm_combination_seq(mm=k, seq=guide)
			for m_seq in mm_seqs:
				if len(m_seq) < 10:
					print(m_seq)
				f.write(">gRNA-%s-%s:mm%d:count%s\n%s\n" % (gene, direction, k, mm_count, m_seq))
				grna_dict["gRNA-%s-%s" % (gene, direction)]["mm%d:mm_count%d" % (k, mm_count)] = m_seq
				mm_count += 1

		f.close()
		pickle.dump(grna_dict, open(path + "fasta_dict/%s_%s_%s_dict.p" % (gene, guide, direction), "wb"))

	else:
		grna_dict = pickle.load(open(path + "fasta_dict/%s_%s_%s_dict.p" % (gene, guide, direction), "rb"))
	return grna_dict


def run_offtargets(genome, gene, guide, direction, mm, threads_num):

	try:
		os.mkdir("%ssam_files/" % path)
	except FileExistsError:
		pass

	try:
		os.mkdir("%sfasta/" % path)
	except FileExistsError:
		pass

	if "%s_%s_%s_alignment.sam" % (gene, guide, direction) not in os.listdir(path + "sam_files/"):
		# Create fasta file
		print("Writing Fasta file..")
		if "%s_%s_%s_fasta.fa" % (gene, guide, direction) not in os.listdir(path + "fasta/"):
			grna_fasta(gene= gene, guide=guide, direction=direction, mm=mm)
			print("Fasta file created.")

		if "%s.fa.index" % genome not in os.listdir("%sgenome/" % path):
			print("Indexing the genome..")
			os.system("mrsfast --index %sgenome/%s.fa" % (path, genome))

		print("\ngRNA alignment on the genome..\n")
		os.system(
			"mrsfast --search %sgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.all.final.fa --seq %sfasta/%s_%s_%s_fasta.fa --threads %s -e 0 -o %ssam_files/%s_%s_%s_alignment.sam --disable-nohits"
			% (path, path, gene, guide, direction, threads_num, path, gene, guide, direction))
		print("\ngRNA alignment on the genome was finished.\n")

	if "%s_%s_%s_alignment.sam" % (gene, guide, direction) in os.listdir(path + "sam_files/"):
		return True
	else:
		print("No alignment - off target")
		return False


def read_sam(gene, guide, direction, mm, path):

	try:
		os.mkdir("%sot_files/" % path)
	except FileExistsError:
		pass

	if "%s_%s_%s_ots.csv" % (gene, guide, direction) not in os.listdir(path + "ot_files/"):
		ot_df = pandas.DataFrame(columns=["gene_flag", "guide_seq", "direction", "mm", "chr", "start", "mm_count"])
		grna_dict = grna_fasta(gene=gene, guide=guide, direction=direction, mm=mm)
		f = open(path + "sam_files/%s_%s_%s_alignment.sam" % (gene, guide, direction), "r")
		lines = f.readlines()
		ind_count = 0
		for line in lines:
			if line[:4] == "gRNA":
				l = line.split("\t")
				if len(l[0].split(":")) == 1:
					gene_flag = l[0]
					mm = "0"
					mm_count = None
				else:
					gene_flag = l[0].split(":")[0]
					mm = l[0].split(":")[1][2:]
					mm_count = l[0].split(":")[2][5:]
				ot_df.loc[ind_count, "gene_flag"] = gene_flag
				ot_df.loc[ind_count, "guide_seq"] = l[9]
				ot_df.loc[ind_count, "direction"] = guide_flag.split("-")[-1]
				ot_df.loc[ind_count, "mm"] = mm
				ot_df.loc[ind_count, "chr"] = l[2]
				ot_df.loc[ind_count, "start"] = l[3]
				ot_df.loc[ind_count, "mm_count"] = mm_count
				ind_count += 1

		f.close()

		mm_col = ["mm%d" % (i + 1) for i in range(int(mm))]
		df_col = ["crispr_sequence", "exact"] + mm_col
		summary_ot_df = pandas.DataFrame(columns=df_col, index=list(ot_df.groupby(["gene_flag"]).groups.keys()))

		for g, g_df in ot_df.groupby(["gene_flag"]):
			exact = g_df.groupby(["mm"]).size().loc["0"]
			# Alignment sequence might be different strand but we want the sequence with PAM
			# Retrieve sequence of the actual gRNA
			grna_seq = grna_dict[g_df[g_df.mm == "0"]["gene_flag"].unique()[0]]["seq"]
			summary_ot_df.loc[g, "crispr_sequence"] = grna_seq
			summary_ot_df.loc[g, "exact"] = exact
			for i in range(1, int(mm) + 1):
				if str(i) in g_df.mm.unique():
					summary_ot_df.loc[g, "mm%d" % i] = g_df.groupby(["mm"]).size().loc[str(i)]
				else:
					summary_ot_df.loc[g, "mm%d" % i] = 0

		summary_ot_df = summary_ot_df.reset_index()
		changed_df_col = ["gene_flag", "crispr_sequence", "exact"] + mm_col
		summary_ot_df.columns = changed_df_col

		summary_ot_df.to_csv(path + "ot_files/%s_%s_%s_ots.csv" % (gene, guide, direction))
	else:
		summary_ot_df = pandas.read_csv(path + "ot_files/%s_%s_%s_ots.csv" % (gene, guide, direction))

	return summary_ot_df



def main(args):

	ot = run_offtargets(genome=args["GENOME"], gene=args["GENE"], guide=args["GRNA"],
						direction=args["DIRECTION"], mm=mm, threads_num=args["NUM"])

	ot_df = read_sam(gene=args["GENE"], guide=args["GRNA"], direction=args["DIRECTION"], mm=mm, path = args["PATH"])

	return True



if __name__ == '__main__':

	args = take_input()

	# Path
	path = ""
	if args["PATH"][-1] == "/": path = args["PATH"]
	else: path = args["PATH"] + "/"

	_ = main(args)


