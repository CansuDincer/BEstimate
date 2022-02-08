# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                            B E s t i m a t e - g M O D E L                               #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, sys, pandas, re, argparse, requests
import urllib.request

# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate - gMODEL",
									 usage="%(prog)s [inputs]",
									 description="""
                                     **********************************
                                     Find and Analyse Base Editor sites
                                           Cancer Model Specific
                                     **********************************""")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# CELL LINE SELECTION
	cell_line = parser.add_mutually_exclusive_group()
	cell_line.add_argument("-sanger", dest="SANGER", required=True,
						   help="The SANGER ID of the interested cancer cell line!")

	cell_line.add_argument("-cosmic", dest="COSMIC", required=True,
						   help="The COSMIC ID of the interested cancer cell line!")

	cell_line.add_argument("-name", dest="CELL_LINE", required=True,
						   help="The standard name of the interested cancer cell line!")

	cell_line.add_argument("-broad", dest="CELL_LINE", required=True,
						   help="The BROAD ID of the interested cancer cell line!")

	# OUTPUT

	parser.add_argument("-o", dest="OUTPUT_PATH", default=os.getcwd() + "/",
						help="The path for output. If not specified the current directory will be used!")

	parser.add_argument("-ofile", dest="OUTPUT_FILE", default="output",
						help="The output file name, if not specified \"position\" will be used!")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


args = take_input()


# -----------------------------------------------------------------------------------------#
# Functions

def find_sanger_model_id(model_name, model_source):

	sanger_id = None

	# Cell Line List from Cell Model Passports
	with urllib.request.urlopen("https://cog.sanger.ac.uk/cmp/download/model_list_20220124.csv") as f:
		cl_df = pandas.read_csv(f)[["model_id", "model_name", "mutational_burden",
									"COSMIC_ID", "BROAD_ID"]]

	# Searching Cell Model Passport for equivalent SANGER ID
	if model_source == "COSMIC":
		if model_name in list(cl_df["COSMIC_ID"]):
			sanger_id = cl_df[cl_df.COSMIC_ID == model_name]["model_id"]
		else:
			sys.exit("No SANGER Model ID equivalent")

	elif model_source == "BROAD":
		if model_name in list(cl_df["BROAD_ID"]):
			sanger_id = cl_df[cl_df.BROAD_ID == model_name]["model_id"]
		else:
			sys.exit("No SANGER Model ID equivalent")

	elif model_source == "NAME":
		if model_name in list(cl_df["model_name"]):
			sanger_id = cl_df[cl_df.model_name == model_name]["model_id"]
		elif model_name in ["".join(cl.split("-")) for cl in cl_df["model_name"]]:
			sanger_id = cl_df["".join(cl_df.model_name.split("-")) == model_name]["model_id"]
		else:
			sys.exit("No SANGER Model ID equivalent")

	return sanger_id


def find_gene_symbol(gene_id):

	gene_name = None

	# Gene List from Cell Model Passports
	with urllib.request.urlopen("https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20191101.csv") as f:
		gene_df = pandas.read_csv(f)

	if gene_id in list(gene_df["gene_id"]):
		gene_name = gene_df[gene_df.gene_id == gene_id]["hgnc_symbol"].values[0]
	return gene_name


def find_gene_id(gene_name):

	gene_id = None

	# Gene List from Cell Model Passports
	with urllib.request.urlopen("https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20191101.csv") as f:
		gene_df = pandas.read_csv(f)

	if gene_name in list(gene_df["hgnc_symbol"]):
		gene_id = gene_df[gene_df.hgnc_symbol == gene_name]["gene_id"].values[0]

	return gene_id


def extract_mutations(sanger_id):

	mut_request = requests.get("https://api.cellmodelpassports.sanger.ac.uk/models/%s/datasets/"
							   "mutations?page[size]=0&include=gene&fields[gene]=symbol" % sanger_id,
							   headers={"Content-Type": "application/json"})

	if mut_request.status_code == 200:
		dfs = list()
		k = 0
		for i in mut_request.json()["data"]:
			df = pandas.DataFrame(i["attributes"], index = [k])
			k += 1
			dfs.append(df)
		df = pandas.concat(dfs)

		return df

	else:
		return None


def mutate_gene(gene_name, model_name, model_source):

	if model_source != "SANGER":
		sanger_id = find_sanger_model_id(model_name=model_name, model_source=model_source)
	else:
		sanger_id = model_name

	if sanger_id is not None:
		mutation_df = extract_mutations(sanger_id = sanger_id)

		sanger_gene_id = find_gene_id(gene_name=gene_name)
		if sanger_gene_id is not None:
			if sanger_gene_id in list(mutation_df.gene_id):
				gene_mutation_df = mutation_df[mutation_df.gene_id == sanger_gene_id]


	else:
		return None




# -----------------------------------------------------------------------------------------#
# Execution

def main():




# -----------------------------------------------------------------------------------------#
# VCF related codes - complicated and computationally heavy + limitation --> need WGS data
"""
def extract_chromosome_info(assembly):

	if "chromosome_lengths.csv" not in os.listdir(os.getcwd() + "/BEstimate/data/"):
		server = "http://grch37.rest.ensembl.org" if assembly == "hg19" else \
			("https://rest.ensembl.org" if assembly == "GRCh38" else None)

		chr_request = requests.get(server + "/info/assembly/human?",
								   headers={"Content-Type": "application/json"})

		if chr_request.status_code == 200:
			dfs = []
			for i in chr_request.json()["top_level_region"]:
				if i["coord_system"] == "chromosome":
					df = pandas.DataFrame({"chromosome": [i["name"]], "length": [i["length"]]})
					dfs.append(df)
			chr_df = pandas.concat(dfs, ignore_index=True)
			chr_df.to_csv(os.getcwd() + "/BEstimate/data/chromosome_lengths.csv", index=False)
		else:
			chr_df = None

	else:
		chr_df = pandas.read_csv(os.getcwd() + "/BEstimate/data/chromosome_lengths.csv")

	return chr_df


def get_chromosome_number(chromosome, assembly):
	chr_df = extract_chromosome_info(assembly)
	return chr_df[chr_df.chromosome == chromosome]["length"].values[0]


def map_vcf_genome(assembly, vcf_path):

	# Unfinished  - complicated and computationally heavy

	chromosomes_map = dict()

	server = "http://grch37.rest.ensembl.org" if assembly == "hg19" else \
		("https://rest.ensembl.org" if assembly == "GRCh38" else None)

	# Read VCF file
	vcf = pandas.read_csv(vcf_path, comment='#', chunksize=10000,
						  delim_whitespace=True, header=None).read()[[0,1,3,4,6]]
	vcf.columns = ["chromosome", "position", "reference", "alteration", "quality"]
	vcf = vcf[vcf.quality == "PASS"]
	vcf["type"] = vcf.apply(lambda x: "insertion" if len(x.alteration) > len(x.reference) else (
		"deletion" if len(x.alteration) < len(x.reference) else None), axis=1)

	for chr, chr_df in vcf.groupby(["chromosome"]):
		chr_df = chr_df.sort_values(by=["position"], ascending=True)
		for ind, row in chr_df.iterrows():
			pos = int(row.position)
			which_1m = pos//1000000
			for mil in range(0, which_1m):
				region = "chr:%s..%s" % (str(mil*1000000), str(mil*1000000+ 1000000))





	# Retrieve reference genome

	r = requests.post(server + "/sequence/region/human",
					  data='{ "regions" : ["X:1000000..1000100:1", "ABBA01004489.1:1..100"] }',
					  headers= {"Content-Type": "application/json", "Accept": "application/json"})
"""



