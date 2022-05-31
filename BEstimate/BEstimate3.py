# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, sys, pandas, re, argparse, requests
from Bio import SeqIO


# -----------------------------------------------------------------------------------------#
# Take inputs

def take_input():
	parser = argparse.ArgumentParser(prog="BEstimate",
									 usage="%(prog)s [inputs]",
									 description="""
                                     **********************************
                                     Find and Analyse Base Editor sites
                                     **********************************""")

	for group in parser._action_groups:
		if group.title == "optional arguments":
			group.title = "Inputs"
		elif "positional arguments":
			group.title = "Mandatory Inputs"

	# BASIC INFORMATION

	parser.add_argument("-gene", dest="GENE", required=True,
						help="The hugo symbol of the interested gene!")

	parser.add_argument("-assembly", dest="ASSEMBLY", required=True,
						help="The genome assembly that will be used!")

	parser.add_argument("-transcript", dest="TRANSCRIPT", default=None,
						help="The interested ensembl transcript id")

	# PAM AND PROTOSPACER INFORMATION

	# The NGG PAM will be used unless otherwise specified.
	parser.add_argument("-pamseq", dest="PAMSEQ", default="NGG",
						help="The PAM sequence in which features used "
							 "for searching activity window and editable nucleotide.")
	parser.add_argument("-pamwin", dest="PAMWINDOW", default="21-23",
						help="The index of the PAM sequence when starting "
							 "from the first index of protospacer as 1.")
	parser.add_argument("-actwin", dest="ACTWINDOW", default="4-8",
						help="The index of the activity window when starting "
							 "from the first index of protospacer as 1.")
	parser.add_argument("-protolen", dest="PROTOLEN", default="20",
						help="The total protospacer and PAM length.")

	# VEP and PROTEIN LEVEL ANALYSIS
	parser.add_argument("-vep", dest="VEP", action="store_true",
						help="The boolean option if user wants to analyse the edits through VEP.")
	parser.add_argument("-P", dest="PROTEIN", action="store_true",
						help="The boolean option if user wants to analyse the edits.")

	# MUTATION INFORMATION
	parser.add_argument("-mutation", dest="MUTATION", default=None,
						help="The mutation on the interested gene that you need to integrate "
							 "into guide and/or annotation analysis")
	parser.add_argument("-mutation_file", dest="MUTATION_FILE", default=None, type=argparse.FileType('r'),
						help="If you have more than one mutations, a file for the mutations on the "
							 "interested gene that you need to integrate into guide and/or annotation analysis")

	# BE INFORMATION

	parser.add_argument("-edit", dest="EDIT", choices=["A", "T", "G", "C"],
						help="The nucleotide which will be edited.")

	parser.add_argument("-edit_to", dest="EDIT_TO", choices=["A", "T", "G", "C"],
						help="The nucleotide after edition.")

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
# Objects from APIs


class Uniprot:

	def __init__(self, uniprotid):
		self.uniprotid, self.reviewed = uniprotid, None
		self.sequence = None
		self.domains = dict()
		self.phosphorylation_sites = dict()
		self.ubiquitination_sites = dict()
		self.methylation_sites = dict()
		self.acetylation_sites = dict()
		self.mutagenesis = dict()
		self.server = "https://www.ebi.ac.uk/proteins/api/"

	def extract_uniprot_info(self):

		uniprot_api = "proteins?offset=0&size=-1&accession=%s" % self.uniprotid
		api_request = requests.get(self.server + uniprot_api,
								   headers={"Accept": "application/json"})

		# Check the response of the server for the request
		if api_request.status_code != 200:
			return "No response from UniProt!\n"

		else:
			for i in api_request.json():
				if len(i["accession"].split("-")) == 1:
					self.reviewed = False if i["info"]["type"] == "TrEMBL" else True
					self.sequence = i["sequence"]["sequence"]
					if "features" in i.keys() and i["features"] != []:
						for ftr in i["features"]:
							if ftr["type"] == "MOD_RES" and ftr["category"] == "PTM":
								if "description" in ftr.keys():
									pos, ptm = ftr["begin"], ftr["description"]
									ptm = ptm.split(";")[0]
									# Phosphorylation
									phos = ptm if re.search(r'Phospho', ptm) else None
									if phos is not None: self.phosphorylation_sites[pos] = phos
									# Methylation
									methy = ptm if re.search(r'Methyl', ptm) else None
									if methy is not None: self.methylation_sites[pos] = methy
									# Ubiquitination
									ubi = ptm if re.search(r'Ub', ptm) else None
									if ubi is not None: self.ubiquitination_sites[pos] = ubi
									# Acetylation
									acety = ptm if re.search(r'Ace', ptm) or re.search(r'N-ace', ptm) else None
									if acety is not None: self.acetylation_sites[pos] = acety

							if ftr["category"] == "DOMAINS_AND_SITES":
								if "description" in ftr.keys():
									domain = ftr["description"]
									domain_range = [int(ftr["begin"]), int(ftr["end"])]
									self.domains[domain] = domain_range

			if self.phosphorylation_sites == dict(): self.phosphorylation_sites = None
			if self.methylation_sites == dict(): self.methylation_sites = None
			if self.acetylation_sites == dict(): self.acetylation_sites = None
			if self.ubiquitination_sites == dict(): self.ubiquitination_sites = None
			if self.domains == dict(): self.domains = None

			return "UniProt API request is done."

	def find_domain(self, protein_edit_location, old_aa):

		edit_domain = None
		if self.domains != {} and self.domains is not None:
			for domain, domain_range in self.domains.items():
				if int(domain_range[0]) <= protein_edit_location <= int(domain_range[1]):
					if old_aa == self.sequence[protein_edit_location - 1]:
						edit_domain = domain
		return edit_domain

	def find_ptm_site(self, ptm_type, protein_edit_location, old_aa):

		edit_ptm_site = None
		if ptm_type == "phosphorylation": d = self.phosphorylation_sites
		if ptm_type == "methylation": d = self.methylation_sites
		if ptm_type == "ubiquitination": d = self.ubiquitination_sites
		if ptm_type == "acetylation": d = self.acetylation_sites

		if d != {} and d is not None:
			for ptm_pos, ptm in d.items():
				if int(ptm_pos) == int(protein_edit_location):
					if old_aa == self.sequence[int(protein_edit_location) - 1]:
						edit_ptm_site = ptm
		return edit_ptm_site

	def extract_mutagenesis(self):

		mut_api = "features/%s?types=MUTAGEN" % self.uniprotid
		mut_api_request = requests.get(self.server + mut_api,
									   headers={"Accept": "application/json"})

		# Check the response of the server for the request
		if mut_api_request.status_code != 200:
			return "No response from UniProt Mutagen!\n"

		else:
			if "features" in mut_api_request.json().keys():
				if mut_api_request.json()["features"]:
					for mut in mut_api_request.json()["features"]:
						if mut["category"] == "MUTAGENESIS":
							if mut["begin"] == mut["end"]:
								self.mutagenesis[mut["begin"]] = {mut["alternativeSequence"]: [mut["description"]]}
								if mut["alternativeSequence"] not in self.mutagenesis[mut["begin"]].keys():
									self.mutagenesis[mut["begin"]][mut["alternativeSequence"]] = [mut["description"]]
								else:
									if mut["description"] not in self.mutagenesis[mut["begin"]][
										mut["alternativeSequence"]]:
										self.mutagenesis[mut["begin"]][mut["alternativeSequence"]].append(
											mut["description"])

	def find_mutagenesis(self, protein_edit_location, new_aa):

		if str(protein_edit_location) in self.mutagenesis.keys():
			mut = self.mutagenesis[str(protein_edit_location)]
			if new_aa in mut.keys():
				mutagenesis = mut[new_aa]
				if len(mutagenesis) == 1:
					return mutagenesis[0]
				else:
					return ";".join(mutagenesis)


class Ensembl:

	def __init__(self, hugo_symbol, assembly):
		self.hugo_symbol = hugo_symbol
		self.assembly = assembly
		self.server = "http://grch37.rest.ensembl.org" \
			if self.assembly == "hg19" else "https://rest.ensembl.org"
		self.gene_id = ''
		self.info_dict = dict()
		self.sequence, self.flan_sequence = None, None
		self.right_sequence_analysis, self.flan_right_sequence_analysis = None, None
		self.left_sequence_analysis, self.flan_left_sequence_analysis = None, None
		self.chromosome, self.strand = None, None
		self.gene_range, self.flan_gene_range = list(), list()

	def extract_gene_id(self):

		hugo_ensembl = "/xrefs/symbol/homo_sapiens/%s?" % self.hugo_symbol

		print("Request to Ensembl REST API for Ensembl Gene ID:")
		gene_request = requests.get(self.server + hugo_ensembl,
									headers={"Content-Type": "application/json"})

		if gene_request.status_code != 200:
			print("No response from ensembl!\n")

		for x in gene_request.json():
			if x["id"][:4] == "ENSG":
				seq_ensembl = self.server + "/sequence/id/%s?" % x["id"]
				seq_request = requests.get(seq_ensembl,
										   headers={"Content-Type": "text/x-fasta"})
				chr = seq_request.text.split("\n")[0].split(":")[2].strip()
				try:
					if chr != "X" and chr != "Y":
						int(chr)
						self.gene_id = x["id"]
					elif chr == "X" or chr == "Y":
						self.gene_id = x["id"]
				except ValueError:
					print(" ")

		if self.gene_id != '':
			print("Ensembl Gene ID: %s\n" % self.gene_id)
			return 1
		else:
			return 0

	def extract_sequence(self, gene_id, mutations=None):

		seq_ensembl = self.server + "/sequence/id/%s?" % gene_id
		seq_flan_ensembl = self.server + "/sequence/id/%s?expand_3prime=23;expand_5prime=23" % gene_id

		print("Request to Ensembl REST API for sequence information:")
		seq_request = requests.get(seq_ensembl,
								   headers={"Content-Type": "text/x-fasta"})
		seq_flan_request = requests.get(seq_flan_ensembl,
										headers={"Content-Type": "text/x-fasta"})

		if seq_request.status_code != 200 and seq_flan_request.status_code != 200:
			print("No response from ensembl sequence!\n")

		# Sequence
		label_line = seq_request.text.split("\n")[0]
		print("The location of the interested gene: %s\n" % label_line.split(" ")[1])
		flan_label_line = seq_flan_request.text.split("\n")[0]
		self.sequence = "".join(seq_request.text.split("\n")[1:])
		self.flan_sequence = "".join(seq_flan_request.text.split("\n")[1:])
		self.gene_range = [int(label_line.split(":")[-3]), int(label_line.split(":")[-2])]
		self.flan_gene_range = [int(flan_label_line.split(":")[-3]),
								int(flan_label_line.split(":")[-2])]
		self.strand = int(label_line.split(":")[-1].strip())
		self.chromosome = label_line.split(":")[2].strip()

		# If strand is -1, the sequence has been reversed to be in 5'->3' direction
		# The genomic location should be reverse to match with the sequence too.
		if self.strand == -1: self.gene_range = [self.gene_range[1], self.gene_range[0]]
		if self.strand == -1: self.flan_gene_range = \
			[self.flan_gene_range[1], self.flan_gene_range[0]]

		# Preparation of the Ensembl sequence for analysis

		nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

		if mutations is None:

			if self.strand == 1:
				self.right_sequence_analysis = self.sequence
				self.flan_right_sequence_analysis = self.flan_sequence
				self.left_sequence_analysis = "".join([nucleotide_dict[n] for n in self.sequence[::-1]])
				self.flan_left_sequence_analysis = "".join([nucleotide_dict[n] for n in self.flan_sequence[::-1]])

			elif self.strand == -1:
				self.left_sequence_analysis = self.sequence
				self.flan_left_sequence_analysis = self.flan_sequence
				self.right_sequence_analysis = "".join([nucleotide_dict[n] for n in self.left_sequence_analysis[::-1]])
				self.flan_right_sequence_analysis = "".join(
					[nucleotide_dict[n] for n in self.flan_left_sequence_analysis[::-1]])

		else:
			for mutation in mutations:
				if mutation.split(":")[0] == self.chromosome:
					if mutation.split(":")[1].split(".")[0] == "g":
						alteration = mutation.split(":")[1].split(".")[1]
						mutation_location = int(re.match("([0-9]+)([a-z]+)", alteration.split(">")[0], re.I).groups()[0])
						altered_nuc = re.match("([0-9]+)([a-z]+)", alteration.split(">")[0], re.I).groups()[1]
						new_nuc = alteration.split(">")[1]

						# Check if altered nucleotide is in the given location
						if self.strand == 1:
							genomic_start = int(self.gene_range[0])
							genomic_flan_start = int(self.flan_gene_range[0])

							if self.sequence[mutation_location - genomic_start] == altered_nuc:
								# Altered nucleotide in the given mutation fits with the sequence
								s = list(self.sequence)
								s[mutation_location - genomic_start] = new_nuc
								self.sequence = "".join(s)
							else:
								print("\nGiven mutation location does not fit with the sequence."
									  "Nucleotides are different.\n")

							if self.flan_sequence[mutation_location - genomic_flan_start] == altered_nuc:
								# Altered nucleotide in the given mutation fits with the sequence
								s = list(self.flan_sequence)
								s[mutation_location - genomic_flan_start] = new_nuc
								self.flan_sequence = "".join(s)
							else:
								print("\nGiven mutation location does not fit with the sequence."
									  "Nucleotides are different.\n")

							self.right_sequence_analysis = self.sequence
							self.flan_right_sequence_analysis = self.flan_sequence
							self.left_sequence_analysis = "".join(
								[nucleotide_dict[n] for n in self.sequence[::-1]])
							self.flan_left_sequence_analysis = "".join(
								[nucleotide_dict[n] for n in self.flan_sequence[::-1]])

						elif self.strand == -1:
							genomic_end = int(self.gene_range[1])
							genomic_flan_end = int(self.flan_gene_range[1])

							if self.sequence[genomic_end - (mutation_location + 1)] == altered_nuc:
								# Altered nucleotide in the given mutation fits with the sequence
								s = list(self.sequence)
								s[genomic_end - (mutation_location + 1)] = new_nuc
								self.sequence = "".join(s)
							else:
								print("\nGiven mutation location does not fit with the sequence."
									  "\nNucleotides are different.\n")
							if self.flan_sequence[genomic_flan_end - (mutation_location + 1)] == altered_nuc:
								# Altered nucleotide in the given mutation fits with the sequence
								s = list(self.flan_sequence)
								s[genomic_flan_end - (mutation_location + 1)] = new_nuc
								self.flan_sequence = "".join(s)
							else:
								print("\nGiven mutation location does not fit with the sequence."
									  "\nNucleotides are different.\n")

							self.left_sequence_analysis = self.sequence
							self.flan_left_sequence_analysis = self.flan_sequence
							self.right_sequence_analysis = "".join(
								[nucleotide_dict[n] for n in self.left_sequence_analysis[::-1]])
							self.flan_right_sequence_analysis = "".join(
								[nucleotide_dict[n] for n in self.flan_left_sequence_analysis[::-1]])

	def extract_info(self, chromosome, loc_start, loc_end, transcript=None):

		ensembl = "/overlap/region/human/%s:%s-%s?feature=transcript;feature=exon;feature=mane;" \
				  "feature=cds" % (
			chromosome, int(loc_start), int(loc_end))

		request = requests.get(self.server + ensembl, headers={"Content-Type": "application/json"})

		if request.status_code != 200:
			print("No response from ensembl!")
		else:
			for output in request.json():
				if transcript is None:
					if output["feature_type"] == "transcript" and output["Parent"] == self.gene_id and output["source"] == "ensembl_havana":
						if output["id"] not in self.info_dict.keys():
							self.info_dict[output["id"]] = \
								[{"start": output["start"], "end": output["end"], "biotype": output["biotype"],
								  "canonical": output["is_canonical"], "source": output["source"]}]

						else:
							old_val = self.info_dict[output["id"]]
							if {"start": output["start"], "end": output["end"],
								"biotype": output["biotype"], "canonical": output["is_canonical"], "source": output["source"]} not in old_val:
								old_val.append(
									{"start": output["start"], "end": output["end"],
									 "biotype": output["biotype"], "canonical": output["is_canonical"], "source": output["source"]})
								self.info_dict[output["id"]] = old_val
				else:
					# Selected transcript
					if output["feature_type"] == "transcript" and output["Parent"] == self.gene_id:
						transcript_info = "/lookup/id/%s?expand=1;mane=1" % output["transcript_id"]
						transcript_request = requests.get(self.server + transcript_info,
														  headers={"Content-Type": "application/json"})
						transcript_output = transcript_request.json()
						if output["transcript_id"] not in self.info_dict.keys():
							self.info_dict[output["transcript_id"]] = \
								[{"start": transcript_output["start"], "end": transcript_output["end"],
								  "biotype": transcript_output["biotype"], "canonical": transcript_output["is_canonical"],
								  "source": output["source"]}]

						else:
							old_val = self.info_dict[output["transcript_id"]]
							if {"start": transcript_request["start"], "end": transcript_request["end"],
								"biotype": transcript_request["biotype"],
								"canonical": transcript_request["is_canonical"], "source": output["source"]} not in old_val:
								old_val.append(
									{"start": transcript_request["start"], "end": transcript_request["end"],
									 "biotype": transcript_request["biotype"],
									 "canonical": transcript_request["is_canonical"], "source": output["source"]})
								self.info_dict[output["transcript_id"]] = old_val

			for output in request.json():
				if output["feature_type"] == "exon" and self.info_dict != {} and \
						output["Parent"] in self.info_dict.keys() and output["source"] == "ensembl_havana":
					for d in self.info_dict[output["Parent"]]:
						if "exon" not in d.keys():
							d["exon"] = {output["exon_id"]: {"start": output["start"], "end": output["end"]}}
						else:
							if output["exon_id"] not in d["exon"].keys():
								d["exon"][output["exon_id"]] = {"start": output["start"], "end": output["end"]}

			for output in request.json():
				if output["feature_type"] == "cds" and self.info_dict != {} and \
						output["Parent"] in self.info_dict.keys() and output["source"] == "ensembl_havana":
					for d in self.info_dict[output["Parent"]]:
						coding_pos = list(range(output["start"], output["end"] + 1))
						if "cds" not in d.keys():
							d["cds"] = {output["protein_id"]: coding_pos}
						else:
							if output["protein_id"] not in d["cds"].keys():
								d["cds"][output["protein_id"]] = coding_pos
							else:
								t = d["cds"][output["protein_id"]]
								for i in coding_pos:
									if i not in t:
										t.append(i)
								d["cds"][output["protein_id"]] = t

		if self.info_dict != {}:
			return 1
		else:
			return 0

	def check_range_info(self, start, end):

		range_locations = list(range(start, end))
		transcripts_exons = dict()
		if self.info_dict != {}:
			for transcript, transcript_dict_list in self.info_dict.items():
				for transcript_dict in transcript_dict_list:
					in_transcript = False
					for loc in range_locations:
						if transcript_dict["start"] <= loc <= transcript_dict["end"]:
							in_transcript = True
					if in_transcript:
						if "exon" in transcript_dict.keys():
							for exon, exon_dict in transcript_dict["exon"].items():
								in_exon = False
								for loc in range_locations:
									if exon_dict["start"] <= loc <= exon_dict["end"]:
										in_exon = True
								if in_exon:
									if transcript not in transcripts_exons.keys():
										transcripts_exons[transcript] = [exon]
									else:
										if exon not in transcripts_exons[transcript]:
											transcripts_exons[transcript].append(exon)
						else:
							if transcript not in transcripts_exons.keys():
								transcripts_exons[transcript] = list()
		if transcripts_exons != {}:
			return transcripts_exons
		else:
			return None

	def check_cds(self, transcript_id, start, end):
		range_locations = list(range(start, end))
		in_cds = False
		if self.info_dict != {}:
			if transcript_id in self.info_dict.keys():
				for transcript_dict in self.info_dict[transcript_id]:
					if "cds" in transcript_dict.keys():
						for protein, protein_pos in transcript_dict["cds"].items():
							for loc in range_locations:
								if loc in protein_pos:
									in_cds = True

		return in_cds

	def extract_uniprot_info(self, ensembl_pid):

		protein_ensembl = "/xrefs/id/{0}?external_db=Uniprot/SWISSPROT%".format(ensembl_pid)
		protein_request = requests.get(self.server + protein_ensembl,
									   headers={"Content-Type": "application/json"})
		if protein_request.status_code != 200:
			print("No response from ensembl!")
			return 0
		else:
			"""
			Since some of the SwissProt object do not have xref mapping information 
			We are missing SwissProt
			Example : KRAS - ENSP00000308495 - 
			{'ensembl_start': 1, 'description': 'V-Ki-ras2 Kirsten rat sarcoma viral oncogene homolog, isoform CRA_b ', 
			'info_text': 'Generated via sequence_mapped', 'xref_end': 188, 'ensembl_identity': 100, 'info_type': 'SEQUENCE_MATCH',
			 'ensembl_end': 188, 'xref_identity': 100, 'display_id': 'A0A024RAV5.49', 'dbname': 'Uniprot/SPTREMBL', 
			 'synonyms': [], 'primary_id': 'A0A024RAV5', 'evalue': None, 'xref_start': 1, 'version': '49', 
			 'cigar_line': '188M', 'db_display_name': 'UniProtKB/TrEMBL', 'score': 972}

			{'description': 'GTPase KRas GTPase KRas, N-terminally processed', 'display_id': 'P01116.246', 
			'info_text': 'Generated via direct', 'dbname': 'Uniprot/SWISSPROT', 'primary_id': 'P01116', 'synonyms': [], 
			'info_type': 'DIRECT', 'version': '246', 'db_display_name': 'UniProtKB/Swiss-Prot'}

			{'info_type': 'DIRECT', 'db_display_name': 'UniProtKB isoform', 'version': '0', 'dbname': 'Uniprot_isoform', 
			'display_id': 'P01116-2', 'info_text': '', 'description': None, 'synonyms': [], 'primary_id': 'P01116-2'}

			We are loosing P01116 which was retrieved from VEP analysis too. 

			"""
			seq_mapping = dict()
			for i in protein_request.json():
				uniprot = i["primary_id"]
				if i["dbname"] == "Uniprot/SWISSPROT":
					if "ensembl_end" in i.keys() and "ensembl_start" in i.keys() and \
							"xref_end" in i.keys() and "xref_start" in i.keys():
						if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(i["xref_end"]) - int(i["xref_start"]):
							# Otherwise, there is an inconsistency --> Not take it
							seq_mapping[uniprot] = {i["ensembl_start"] + k: i["xref_start"] + k
													for k in range(
									int(i["ensembl_end"]) - int(i["ensembl_start"]) + 1)}
						break
				elif i["dbname"] == "Uniprot/SPTREMBL":
					if "ensembl_end" in i.keys() and "ensembl_start" in i.keys() and \
							"xref_end" in i.keys() and "xref_start" in i.keys():
						if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(i["xref_end"]) - int(i["xref_start"]):
							# Otherwise, there is an inconsistency --> Not take it
							seq_mapping[uniprot] = {i["ensembl_start"] + k: i["xref_start"] + k
													for k in range(
									int(i["ensembl_end"]) - int(i["ensembl_start"]) + 1)}

			if seq_mapping == dict():
				return None
			else:
				return seq_mapping


# -----------------------------------------------------------------------------------------#
# Data w/out API opportunity

yulab = pandas.read_table(os.getcwd() + "/../data/H_sapiens_interfaces.txt")


# -----------------------------------------------------------------------------------------#
# Functions


def find_pam_protospacer(sequence, pam_sequence, searched_nucleotide,
						 activity_window, pam_window, protospacer_length):
	"""
	Finding all possible PAM and protospacer regions on the sequence of the gene.
	:param sequence: The sequence of the interested gene.
	:param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
	:param searched_nucleotide: The interested nucleotide which will be changed with BE
	:param activity_window: The location of the activity window on the protospacer sequence.
	:param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
	:param protospacer_length: The length of protospacer.
	:return crisprs: A list of dictionary having sequences and locations of the gRNA targeted
	gene parts. The indices are indicated when the first index of the gene is 1.
	"""
	# Since python index starts from 0, decrease the start position index given from the user
	activity_window = [activity_window[0] - 1, activity_window[1]]
	pam_window = [pam_window[0] - 1, pam_window[1]]

	# Using Regular Expressions, specify PAM pattern
	pam_pattern = ""
	for nuc in list(pam_sequence):
		if nuc != "N":
			pam_pattern += nuc + "{1}"
		else:
			pam_pattern += "[ATCG]{1}"

	# Search protospacer length of nucleatide sequence and add PAM pattern after that
	pattern = r'[ATCG]{%s}%s' % (protospacer_length, pam_pattern)
	print("Pattern is created!")

	crisprs = []
	print("Pattern is searching through the sequence...")
	for nuc_index in range(0, len(sequence)):

		# One by one in the given sequence
		if nuc_index + pam_window[1] <= len(sequence):

			# Add started nucleotide index total length of targeted base editing site (PAM index)
			sub_sequence = sequence[nuc_index:nuc_index + pam_window[1]]

			# Search regex pattern inside the sub sequence
			for match_sequence in re.finditer(pattern, sub_sequence):

				# If there is match, then check if searched nucleotide inside the activity window
				if searched_nucleotide in list(match_sequence.group()[activity_window[0]:activity_window[1]]):
					activity_sequence = match_sequence.group()[activity_window[0]:activity_window[1]]

					# If searched nucleotide is also there, add the sequence inside crisprs!
					crisprs.append({"index": [nuc_index, nuc_index + pam_window[1]],
									"crispr": match_sequence.group(),
									"activity_seq": activity_sequence})
				else:
					crisprs.append({"index": [nuc_index, nuc_index + pam_window[1]],
									"crispr": match_sequence.group(),
									"activity_seq": "No window"})

	if crisprs is not []: print("CRISPRs are found!")
	return crisprs


def add_genomic_location(sequence_range, crispr_dict, crispr_direction, strand):
	"""
	Adding genomic location info on crisprs found by extract_activity_window() function.
	:param sequence_range: The range of the sequence on the genome (from Ensembl)
	:param crispr_dict: The CRISPR dictionary created by extract_activity_window()
	:param crispr_direction: The direction of the created CRISPR (left or right)
	:param strand: The strand of the given gene (-1 or 1)
	:return crispr_seq: The sequence of the CRISPR
	:return genomic_location: The genomic coordinate of the above CRISPR on genome
	"""
	# Prepare the sequence range according to the strand of the gene

	new_range = []
	if strand == 1:
		new_range = sequence_range
	elif strand == -1:
		new_range = [sequence_range[1], sequence_range[0]]

	# Look for both direction
	crispr_seq, genomic_location = crispr_dict["crispr"], ""

	if crispr_direction == "right":
		genomic_start = new_range[0] + crispr_dict["index"][0]
		genomic_end = (genomic_start + len(crispr_seq)) - 1
		genomic_location = str(genomic_start) + "-" + str(genomic_end)

	elif crispr_direction == "left":
		genomic_start = new_range[1] - crispr_dict["index"][0]
		genomic_end = (genomic_start - len(crispr_seq)) + 1
		genomic_location = str(genomic_end) + "-" + str(genomic_start)

	return crispr_seq, genomic_location


def extract_grna_sites(hugo_symbol, pam_sequence, searched_nucleotide,
					   activity_window, pam_window, protospacer_length,
					   ensembl_object):
	"""
	Extracting the gRNA targeted sites having editable nucleotide(s) on the interested genes
	:param hugo_symbol: The Hugo Symbol of the interested gene.
	:param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
	:param searched_nucleotide: The interested nucleotide which will be changed with BE
	:param activity_window: The location of the activity windiw on the protospacer sequence.
	:param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
	:param protospacer_length: The length of protospacer.
	:param ensembl_object: The Ensembl Object created with Ensembl().
	:return crispr_df: A data frame having sequence, location and direction information of the CRISPRs.
	"""

	print("Sequence is preparing...")

	left_sequence = ensembl_object.flan_left_sequence_analysis
	right_sequence = ensembl_object.flan_right_sequence_analysis
	strand, seq_range = ensembl_object.strand, ensembl_object.flan_gene_range
	chromosome = ensembl_object.chromosome

	# Right CRISPRs 5'-->3' : reversed and base changed sequence

	# Editted should be used
	print("Protospacer and PAM regions are searching for right direction...")
	right_crisprs = find_pam_protospacer(sequence=right_sequence,
										 pam_sequence=pam_sequence,
										 searched_nucleotide=searched_nucleotide,
										 activity_window=activity_window,
										 pam_window=pam_window,
										 protospacer_length=protospacer_length)

	# Left CRISPRs: raw sequence will be used
	print("\nProtospacer and PAM regions are searching for left direction...")
	left_crisprs = find_pam_protospacer(sequence=left_sequence,
										pam_sequence=pam_sequence,
										searched_nucleotide=searched_nucleotide,
										activity_window=activity_window,
										pam_window=pam_window,
										protospacer_length=protospacer_length)

	crisprs_dict = {"left": left_crisprs, "right": right_crisprs}

	print("\nCRISPR df is filling...")

	crisprs_df = pandas.DataFrame(columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
										   "Location", "Direction", "Gene_ID", "Transcript_ID", "Exon_ID",
										   "guide_in_CDS"])

	for direction, crispr in crisprs_dict.items():

		for cr in crispr:
			crispr_seq, genomic_location = add_genomic_location(sequence_range=seq_range, strand=strand,
																crispr_dict=cr, crispr_direction=direction)
			# Transcript & Exon Info
			transcript_exon = ensembl_object.check_range_info(int(genomic_location.split("-")[0]),
															  int(genomic_location.split("-")[1]))
			if transcript_exon is not None:
				for transcript, exon_list in transcript_exon.items():
					if exon_list:
						for exon in exon_list:
							df = pandas.DataFrame([[crispr_seq, crispr_seq[:-len(pam_sequence)],
													chromosome + ":" + genomic_location, direction,
													ensembl_object.gene_id, transcript, exon]],
												  columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
														   "Location", "Direction", "Gene_ID",
														   "Transcript_ID", "Exon_ID"])
							crisprs_df = pandas.concat([crisprs_df, df])
					else:
						df = pandas.DataFrame([[crispr_seq, crispr_seq[:-len(pam_sequence)],
												chromosome + ":" + genomic_location, direction,
												ensembl_object.gene_id, transcript, None]],
											  columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
													   "Location", "Direction", "Gene_ID",
													   "Transcript_ID", "Exon_ID"])
						crisprs_df = pandas.concat([crisprs_df, df])
			else:
				df = pandas.DataFrame([[crispr_seq, crispr_seq[:-len(pam_sequence)],
										chromosome + ":" + genomic_location, direction,
										ensembl_object.gene_id, None, None]],
									  columns=["CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
											   "Location", "Direction", "Gene_ID",
											   "Transcript_ID", "Exon_ID"])
				crisprs_df = pandas.concat([crisprs_df, df])

	crisprs_df["Hugo_Symbol"] = hugo_symbol
	crisprs_df["guide_in_CDS"] = crisprs_df.apply(
		lambda x: ensembl_object.check_cds(x["Transcript_ID"], int(x["Location"].split(":")[1].split("-")[0]),
										   int(x["Location"].split(":")[1].split("-")[1]) + 1)
		if int(x["Location"].split(":")[1].split("-")[0]) < int(x["Location"].split(":")[1].split("-")[1])
		else ensembl_object.check_cds(x["Transcript_ID"], int(x["Location"].split(":")[1].split("-")[1]),
									  int(x["Location"].split(":")[1].split("-")[0]) + 1), axis=1)

	return crisprs_df


def collect_mutation_location(mutations):
	locations = list()
	for mutation in mutations:
		alteration = mutation.split(":")[1].split(".")[1]
		mutation_location = re.match("([0-9]+)([a-z]+)", alteration.split(">")[0], re.I).groups()[0]
		if int(mutation_location) not in locations:
			locations.append(int(mutation_location))
	if locations:
		return locations
	else: return None


def check_genome_for_mutation(genomic_range, direction, mutations):

	yes_mutation = False

	if direction == "left":
		end = int(genomic_range.split("-")[0])
		start = int(genomic_range.split("-")[1])

	elif direction == "right":
		end = int(genomic_range.split("-")[1])
		start = int(genomic_range.split("-")[0])

	for loc in mutations:
		if start <= loc <= end:
			yes_mutation = True
	return yes_mutation


def find_editable_nucleotide(crispr_df, searched_nucleotide, activity_window,
							 ensembl_object, mutations):
	"""
	Finding editable nucleotides and their genomic coordinates
	:param crispr_df: A data frame having sequence, location and direction information of
	the CRISPRs from extract_crisprs().
	:param searched_nucleotide: The interested nucleotide which will be changed with BE
	:param activity_window: The location of the activity window on the protospacer sequence.
	:param ensembl_object: The Ensembl Object created with Ensembl().
	:param mutations: Given mutation list from the user
	:return edit_df: A data frame having sequence, edit_location, location and direction
	information of the CRISPRs.
	"""

	actual_seq_range = ensembl_object.gene_range
	if actual_seq_range[0] > actual_seq_range[1]:
		actual_seq_range = [actual_seq_range[1], actual_seq_range[0]]

	actual_locations = list(range(actual_seq_range[0], actual_seq_range[1]))

	activity_window = [activity_window[0] - 1, activity_window[1]]

	print("Edit df is filling...")
	edit_df = pandas.DataFrame(columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence", "Location",
										"Edit_Location", "Direction", "Strand", "Gene_ID", "Transcript_ID", "Exon_ID",
										"guide_in_CDS", "Edit_in_Exon", "Edit_in_CDS", "guide_on_mutation",
										"guide_change_mutation"])

	for ind, row in crispr_df.iterrows():

		# Check only with the sequence having PAM since it only has the searched nucleotide!
		searched_ind = [nuc_ind for nuc_ind in range(0, len(row["gRNA_Target_Sequence"]))
						if nuc_ind in list(range(activity_window[0], activity_window[1])) and
						row["gRNA_Target_Sequence"][nuc_ind] == searched_nucleotide]

		if searched_ind is not []:
			# If there is an editable nucleotide in the activity sites
			actual_inds = []
			if row["Direction"] == "left":
				for nuc_ind in searched_ind:
					if int(row["Location"].split(":")[1].split("-")[1]) - nuc_ind in actual_locations:
						actual_inds.append(int(row["Location"].split(":")[1].split("-")[1]) - nuc_ind)

			elif row["Direction"] == "right":
				for nuc_ind in searched_ind:
					if int(row["Location"].split(":")[1].split("-")[0]) + nuc_ind in actual_locations:
						actual_inds.append(int(row["Location"].split(":")[1].split("-")[0]) + nuc_ind)

			for actual_ind in actual_inds:

				transcript_exon = ensembl_object.check_range_info(actual_ind, actual_ind + 1)

				if transcript_exon is not None:
					for transcript, exon_list in transcript_exon.items():
						if row["Exon_ID"] is not None and row["Exon_ID"] in exon_list:
							edit_in_exon = True
						else:
							edit_in_exon = False
				else:
					edit_in_exon = False

				edit_in_cds = ensembl_object.check_cds(row["Transcript_ID"], actual_ind, actual_ind + 1)

				df = pandas.DataFrame([[row["Hugo_Symbol"], row["CRISPR_PAM_Sequence"],
										row["gRNA_Target_Sequence"], row["Location"], actual_ind,
										row["Direction"], ensembl_object.strand, ensembl_object.gene_id,
										row["Transcript_ID"], row["Exon_ID"], row["guide_in_CDS"],
										edit_in_exon, edit_in_cds]],
									  columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
											   "Location", "Edit_Location", "Direction", "Strand",
											   "Gene_ID", "Transcript_ID", "Exon_ID", "guide_in_CDS",
											   "Edit_in_Exon", "Edit_in_CDS"])
				edit_df = pandas.concat([edit_df, df])

		else:
			# If not --> no edit
			df = pandas.DataFrame([[row["Hugo_Symbol"], row["CRISPR_PAM_Sequence"],
									row["gRNA_Target_Sequence"], row["Location"], "No edit",
									row["Direction"], ensembl_object.strand, ensembl_object.gene_id,
									row["Transcript_ID"], row["Exon_ID"], row["guide_in_CDS"], False, False]],
								  columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
										   "Location", "Edit_Location", "Direction", "Strand", "Gene_ID",
										   "Transcript_ID", "Exon_ID", "guide_in_CDS", "Edit_in_Exon", "Edit_in_CDS"])
			edit_df = pandas.concat([edit_df, df])

	edit_df["# Edits/guide"] = 0
	for guide, g_df in edit_df.groupby(["gRNA_Target_Sequence"]):
		unique_edits_per_guide = len(set(list(g_df["Edit_Location"])))
		edit_df.loc[edit_df.gRNA_Target_Sequence == guide, "# Edits/guide"] = unique_edits_per_guide

	edit_df["Poly_T"] = edit_df.apply(
		lambda x: True if re.search("T{4,}", x.CRISPR_PAM_Sequence) is not None else False, axis=1)

	mutation_locations = collect_mutation_location(mutations=mutations)

	edit_df["guide_on_mutation"] = edit_df.apply(
		lambda x: check_genome_for_mutation(genomic_range=x.Location.split(":")[1],
											direction=x.Direction, mutations=mutation_locations), axis=1)

	edit_df["guide_change_mutation"] = edit_df.apply(
		lambda x: True if int(x.Edit_Location) in mutation_locations else False, axis=1)

	return edit_df


def extract_hgvs(edit_df, ensembl_object, transcript_id, edited_nucleotide,
				 new_nucleotide, activity_window, mutations):
	"""
	Collect Ensembl VEP information for given edits
	:param edit_df: Edit data frame created with find_editable_nucleotide()
	:param ensembl_object: The Ensembl Object created with Ensembl().
	:param transcript_id: Ensembl Transcript id for the filtration
	:param edited_nucleotide: The interested nucleotide which will be changed with BE.
	:param new_nucleotide: The new nucleotide which will be changed to with BE.
	:param activity_window: The location of the activity window on the protospacer sequence.
	:param mutations: Given mutation list from the user
	:return hgvs_df: The HGVS notations of all possible variants
	"""
	# For (-) direction crisprs, base reversion should be done.
	nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

	# Collect chromosome
	chromosome, strand = ensembl_object.chromosome, ensembl_object.strand
	activity_window = [activity_window[0] - 1, activity_window[1]]

	# Collect mutations
	mutation_locations = collect_mutation_location(mutations)

	# Transcript filtration
	if transcript_id is not None:
		loc_edit_df = edit_df[edit_df.Transcript_ID == transcript_id]
	else:
		for transcript, transcript_dict in ensembl_object.info_dict.items():
			for d in transcript_dict:
				if d["source"] == "ensembl_havana":
					loc_edit_df = edit_df[edit_df.Transcript_ID == transcript]

	# Each gRNA at a time
	row_dicts = list()
	for direction, direction_df in loc_edit_df.groupby(["Direction"]):
		if direction == "left":
			# Base reversion of the (-) direction crisprs
			rev_edited_nucleotide, rev_new_nucleotide = \
				nucleotide_dict[edited_nucleotide], nucleotide_dict[new_nucleotide]

			for grna, grna_df in direction_df.groupby(["gRNA_Target_Sequence"]):

				total_edit = len(set(list(grna_df["Edit_Location"].values)))

				# For individual edits
				for edit_loc, grna_edit_df in grna_df.groupby(["Edit_Location"]):

					if True not in grna_edit_df.guide_on_mutation.unique():

						hgvs = "%s:g.%s%s>%s" \
							   % (str(chromosome), str(edit_loc), rev_edited_nucleotide, rev_new_nucleotide)

					else:
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]

						guide_change_mutation = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								guide_change_mutation.append(mutation)

						first_mut = max(guide_change_mutation)
						last_mut = min(guide_change_mutation)
						if first_mut > start:
							start_ind = activity_window[0] - (first_mut - start)
							start = first_mut
						else:
							if last_mut < end:
								end_ind = activity_window[1] + (end - last_mut)
								end = last_mut

						position = str(start) + "_" + str(end)

						activity_sites = grna[start_ind: end_ind]
						activity_sites = "".join([nucleotide_dict[n] for n in activity_sites[::-1]])
						# Mutation included with the sequence of guide
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)

						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					d = {"Hugo_Symbol": list(grna_edit_df["Hugo_Symbol"].values)[0], "Edit_Type": "individual",
						 "CRISPR_PAM_Sequence": grna_edit_df["CRISPR_PAM_Sequence"].values[0],
						 "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
						 "gRNA_Target_Sequence": grna,
						 "gRNA_Target_Location": grna_edit_df["Location"].values[0].split(":")[0] + ":" +
												 str(int(
													 grna_edit_df["Location"].values[0].split(":")[1].split("-")[
														 0]) - 3) + "-" +
												 grna_edit_df["Location"].values[0].split(":")[1].split("-")[1],
						 "Total_Edit": total_edit, "Edit_Location": edit_loc,
						 "Direction": direction,
						 "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
						 "Exon_ID": grna_edit_df["Exon_ID"].values[0],
						 "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
						 "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
						 "guide_on_mutation": grna_edit_df["guide_on_mutation"].values[0],
						 "guide_change_mutation": grna_edit_df["guide_change_mutation"].values[0],
						 "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
						 "Poly_T": grna_edit_df["Poly_T"].values[0],
						 "HGVS": hgvs}

					row_dicts.append(d)

				if total_edit > 1:

					if True not in grna_df.guide_on_mutation.unique():
						# For multiple edits
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]
						position = str(start) + "_" + str(end)

						guide_change_mutation = False
						for mutation in mutation_locations:
							if start <= mutation <= end:
								guide_change_mutation = True

						activity_sites = grna[activity_window[0]: activity_window[1]]
						activity_sites = "".join([nucleotide_dict[n] for n in activity_sites[::-1]])
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)
						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					else:
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]

						guide_change_mutation = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								guide_change_mutation.append(mutation)

						first_mut = max(guide_change_mutation)
						last_mut = min(guide_change_mutation)
						if first_mut > start:
							start_ind = activity_window[0] - (first_mut - start)
							start = first_mut
						else:
							if last_mut < end:
								end_ind = activity_window[1] + (end - last_mut)
								end = last_mut

						position = str(start) + "_" + str(end)

						activity_sites = grna[start_ind: end_ind]
						activity_sites = "".join([nucleotide_dict[n] for n in activity_sites[::-1]])
						# Mutation included with the sequence of guide
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)

						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					d = {"Hugo_Symbol": list(grna_df["Hugo_Symbol"].values)[0], "Edit_Type": "multiple",
						 "CRISPR_PAM_Sequence": grna_df["CRISPR_PAM_Sequence"].values[0],
						 "CRISPR_PAM_Location": grna_df["Location"].values[0],
						 "gRNA_Target_Sequence": grna,
						 "gRNA_Target_Location": grna_df["Location"].values[0].split(":")[0] + ":" +
												 str(int(grna_df["Location"].values[0].split(":")[1].split("-")[
															 0]) - 3) + "-" +
												 grna_df["Location"].values[0].split(":")[1].split("-")[1],
						 "Total_Edit": total_edit,
						 "Edit_Location": position.split("_")[0] + "-" + position.split("_")[1],
						 "Direction": direction,
						 "Transcript_ID": grna_df["Transcript_ID"].values[0],
						 "Exon_ID": grna_df["Exon_ID"].values[0],
						 "guide_in_CDS": grna_df["guide_in_CDS"].values[0],
						 "Edit_in_Exon": grna_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_df["Edit_in_CDS"].values[0],
						 "guide_on_mutation": grna_df["guide_on_mutation"].values[0],
						 "guide_change_mutation": guide_change_mutation,
						 "# Edits/guide": grna_df["# Edits/guide"].values[0],
						 "Poly_T": grna_df["Poly_T"].values[0],
						 "HGVS": hgvs}
					row_dicts.append(d)

		elif direction == "right":

			for grna, grna_df in direction_df.groupby(["gRNA_Target_Sequence"]):

				total_edit = len(set(list(grna_df["Edit_Location"].values)))

				# For individual edits

				for edit_loc, grna_edit_df in grna_df.groupby(["Edit_Location"]):

					if True not in grna_edit_df.guide_on_mutation.unique():

						hgvs = "%s:g.%s%s>%s" % (str(chromosome), str(edit_loc), edited_nucleotide, new_nucleotide)

					else:
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]

						guide_change_mutation = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								guide_change_mutation.append(mutation)

						first_mut = min(guide_change_mutation)
						last_mut = max(guide_change_mutation)
						if first_mut < start:
							start_ind = activity_window[0] - (start - first_mut)
							start = first_mut
						else:
							if last_mut > end:
								end_ind = activity_window[1] + (last_mut - end)
								end = last_mut

						position = str(start) + "_" + str(end)

						activity_sites = grna[start_ind: end_ind]
						# Mutation included with the sequence of guide
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)

						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					d = {"Hugo_Symbol": list(grna_edit_df["Hugo_Symbol"].values)[0], "Edit_Type": "individual",
						 "CRISPR_PAM_Sequence": grna_edit_df["CRISPR_PAM_Sequence"].values[0],
						 "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
						 "gRNA_Target_Sequence": grna,
						 "gRNA_Target_Location": grna_edit_df["Location"].values[0].split(":")[0] + ":" +
												 grna_edit_df["Location"].values[0].split(":")[1].split("-")[0] + "-" + \
												 str(int(grna_edit_df["Location"].values[0].split(":")[1].split("-")[
															 1]) - 3),
						 "Total_Edit": total_edit, "Edit_Location": edit_loc,
						 "Direction": direction,
						 "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
						 "Exon_ID": grna_edit_df["Exon_ID"].values[0],
						 "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
						 "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
						 "guide_on_mutation": grna_edit_df["guide_on_mutation"].values[0],
						 "guide_change_mutation": grna_edit_df["guide_change_mutation"].values[0],
						 "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
						 "Poly_T": grna_edit_df["Poly_T"].values[0],
						 "HGVS": hgvs}

					row_dicts.append(d)

				if total_edit > 1:

					if True not in grna_df.guide_on_mutation.unique():
						# For multiple edits
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]
						position = str(start) + "_" + str(end)

						guide_change_mutation = False
						for mutation in mutation_locations:
							if start <= mutation <= end:
								guide_change_mutation = True

						activity_sites = grna[activity_window[0]: activity_window[1]]
						edited_activity_sites = activity_sites.replace(edited_nucleotide, new_nucleotide)
						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					else:
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]

						guide_change_mutation = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								guide_change_mutation.append(mutation)

						first_mut = min(guide_change_mutation)
						last_mut = max(guide_change_mutation)
						if first_mut < start:
							start_ind = activity_window[0] - (start - first_mut)
							start = first_mut
						else:
							if last_mut > end:
								end_ind = activity_window[1] + (last_mut - end)
								end = last_mut

						position = str(start) + "_" + str(end)

						activity_sites = grna[start_ind: end_ind]
						# Mutation included with the sequence of guide
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)

						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					d = {"Hugo_Symbol": grna_df["Hugo_Symbol"].values[0], "Edit_Type": "multiple",
						 "CRISPR_PAM_Sequence": grna_df["CRISPR_PAM_Sequence"].values[0],
						 "CRISPR_PAM_Location": grna_df["Location"].values[0],
						 "gRNA_Target_Sequence": grna,
						 "gRNA_Target_Location": grna_df["Location"].values[0].split(":")[0] + ":" +
												 grna_df["Location"].values[0].split(":")[1].split("-")[0] + "-" + \
												 str(int(grna_df["Location"].values[0].split(":")[1].split("-")[
															 1]) - 3),
						 "Total_Edit": total_edit,
						 "Edit_Location": position.split("_")[0] + "-" + position.split("_")[1],
						 "Direction": direction,
						 "Transcript_ID": grna_df["Transcript_ID"].values[0],
						 "Exon_ID": grna_df["Exon_ID"].values[0],
						 "guide_in_CDS": grna_df["guide_in_CDS"].values[0],
						 "Edit_in_Exon": grna_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_df["Edit_in_CDS"].values[0],
						 "guide_on_mutation": grna_df["guide_on_mutation"].values[0],
						 "guide_change_mutation": guide_change_mutation,
						 "# Edits/guide": grna_df["# Edits/guide"].values[0],
						 "Poly_T": grna_df["Poly_T"].values[0],
						 "HGVS": hgvs}
					row_dicts.append(d)

	hgvs_df = pandas.DataFrame(row_dicts)

	return hgvs_df


def extract_hgvsp(hgvsp, which):
	aa_3to1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q",
			   "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
			   "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
	if hgvsp is not None:
		ensembl_protein = hgvsp.split(":p.")[0].split(".")[0]
		protein_change = hgvsp.split("p.")[1]
		if len(protein_change.split("delins")) == 1:
			# SNP
			if len(protein_change.split("=")) == 1:
				if which == "old_aa":
					return aa_3to1[protein_change[:3]]
				if which == "new_aa":
					return aa_3to1[protein_change[-3:]]
				if which == "position":
					return protein_change[3:-3]
			else:
				if which == "old_aa":
					# Synonymous variant
					aa1 = list()
					aa_string = re.match("([a-z]+)([0-9]+)", protein_change.split("=")[0], re.I).groups()[0]
					for i in [aa_string[x: x + 3] for x in range(0, len(aa_string), 3)]:
						aa1.append(aa_3to1[i])
					return ";".join(aa1)
				if which == "new_aa":
					# Synonymous variant
					aa1 = list()
					aa_string = re.match("([a-z]+)([0-9]+)", protein_change.split("=")[0], re.I).groups()[0]
					for i in [aa_string[x: x + 3] for x in range(0, len(aa_string), 3)]:
						aa1.append(aa_3to1[i])
					return ";".join(aa1)
				if which == "position":
					return re.match("([a-z]+)([0-9]+)", protein_change.split("=")[0], re.I).groups()[1]

		elif len(protein_change.split("delins")) > 1:
			# Substitution
			if which == "old_aa":
				aa1 = list()
				for i in protein_change.split("delins")[0].split("_"):
					aa1.append(aa_3to1[i[:3]])
				return ";".join(aa1)

			if which == "new_aa":
				aa2 = list()
				for i in [protein_change.split("delins")[1][x: x + 3] for x in range(0, len(protein_change.split("delins")[1]), 3)]:
					aa2.append(aa_3to1[i])
				return ";".join(aa2)

			if which == "position":
				pos = list()
				for i in protein_change.split("delins")[0].split("_"):
					pos.append(re.match("([a-z]+)([0-9]+)", i, re.I).groups()[1])
				return ";".join(pos)
	else:
		return None


def retrieve_vep_info(hgvs_df, ensembl_object, transcript_id=None):
	"""
	Collect Ensembl VEP information for given edits
	:param hgvs_df: The HGVS notations of all possible variants
	:param ensembl_object: The Ensembl Object created with Ensembl().
	:param transcript_id: The interested Ensembl transcript id
	:return uniprot_results: The Uniprot IDs in which edit occurs (swissprot or trembl)
	"""

	# Dictionary to find the chemical properperty change due to the edit
	aa_chem = {"G": "Non-Polar", "A": "Non-Polar", "V": "Non-Polar", "C": "Polar", "P": "Non-Polar",
			   "L": "Non-Polar", "I": "Non-Polar", "M": "Non-Polar", "W": "Non-Polar", "F": "Non-Polar",
			   "S": "Polar", "T": "Polar", "Y": "Polar", "N": "Polar", "Q": "Polar", "K": "Charged",
			   "R": "Charged", "H": "Charged", "D": "Charged", "E": "Charged", "*": "-"}

	chromosome, strand = ensembl_object.chromosome, ensembl_object.strand
	"""
	if transcript_id is None:
		for transcript, transcript_dict in ensembl_object.info_dict.items():
			for d in transcript_dict:
				if d["canonical"]:
					transcript_id = transcript
	"""
	if transcript_id is None:
		for transcript, transcript_dict in ensembl_object.info_dict.items():
			for d in transcript_dict:
				if d["source"] == "ensembl_havana":
					transcript_id = transcript


	print("VEP df is filling...")

	# Decide the server
	server = "http://grch37.rest.ensembl.org" if ensembl_object.assembly == "hg19" \
		else "https://rest.ensembl.org"

	all_vep_dfs = list()
	for ind, row in hgvs_df.iterrows():
		# VEP API request
		vep = "/vep/human/hgvs/%s?Blosum62=1;Conservation=1;Mastermind=1;" \
			  "LoF=1;CADD=1;protein=1;variant_class=1;hgvs=1;uniprot=1;transcript_id=%s" \
			  % (row.HGVS, transcript_id)
		vep_request = requests.get(server + vep, headers={"Content-Type": "application/json"})

		vep_dfs = list()
		# Check the response of the server for the request
		if vep_request.status_code != 200:
			print("No response from VEP for %s" % row.HGVS)
			na_df = pandas.DataFrame(
				[[row.Hugo_Symbol, row.CRISPR_PAM_Sequence, row.CRISPR_PAM_Location, row.gRNA_Target_Sequence,
				  row.gRNA_Target_Location, row.Edit_Location, row.Edit_Type, row.Direction, row.Transcript_ID,
				  row.Exon_ID, row.guide_on_mutation, row.guide_change_mutation, None, None, None, None, None, None,
				  None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
				  None, None, None, None, None, None, None, None, None]],
				columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "CRISPR_PAM_Location", "gRNA_Target_Sequence",
						 "gRNA_Target_Location", "Edit_Position", "Edit_Type", "Direction", "Transcript_ID",
						 "Exon_ID", "guide_on_mutation", "guide_change_mutation", "Protein_ID", "VEP_input",
						 "variant_classification", "cDNA_Change", "Edited_Codon", "New_Codon", "Protein_Position",
						 "Protein_Change", "Edited_AA", "Edited_AA_Prop", "New_AA", "New_AA_Prop", "is_Synonymous",
						 "is_Stop", "proline_addition", "swissprot", "polyphen_score", "polyphen_prediction",
						 "sift_score", "sift_prediction", "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
						 "consequence_terms", "is_clinical", "clinical_allele", "clinical_id", "clinical_significance",
						 "var_synonyms"])
			vep_dfs.append(na_df)

		else:
			vep_initials, vep_transcript_results, vep_clinical_results, vep_most_severe_results = \
				list(), list(), list(), list()
			vep_transcript_interested = ["hgvsc", "hgvsp", "protein_id", "transcript_id",
										 "amino_acids", "codons", "polyphen_score",
										 "polyphen_prediction", "sift_score", "sift_prediction",
										 "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
										 "protein_start", "consequence_terms", "swissprot", "protein_id"]

			for x in vep_request.json():
				vep_initial_result = {}
				if "variant_class" in x.keys():
					vep_initial_result["variant_classification"] = x["variant_class"]
					vep_initial_result["VEP_input"] = x["id"]
					vep_initials.append(vep_initial_result)

				vep_most_severe = {}
				if "most_severe_consequence" in x.keys():
					vep_most_severe["most_severe_consequence"] = x["most_severe_consequence"]
					vep_most_severe_results.append(vep_most_severe)

				for t in x["transcript_consequences"]:
					vep_transcript_result = {}
					if t["strand"] == strand and t["gene_id"] == ensembl_object.gene_id \
							and t["transcript_id"] == transcript_id:
						for opt in vep_transcript_interested:
							if opt in t.keys():
								vep_transcript_result[opt] = t[opt]
							else:
								vep_transcript_result[opt] = None

							vep_transcript_results.append(vep_transcript_result)

				if "colocated_variants" in x.keys():
					for c in x["colocated_variants"]:
						vep_clinical_result = {}
						if c["start"] == c["end"]:
							vep_clinical_result["clinical_allele"] = c["allele_string"] \
								if "allele_string" in c.keys() else None
							vep_clinical_result["clinical_id"] = c["id"] if "id" in c.keys() else None

							clinical_sig = ''
							if "clin_sig_allele" in c.keys():
								clinical_sig = c["clin_sig_allele"]
							else:
								clinical_sig = None
							vep_clinical_result["clinical_significance"] = clinical_sig
							var_syns = ''
							if "var_synonyms" in c.keys():
								if type(c["var_synonyms"]) == str:
									for clnv in c["var_synonyms"]["ClinVar"]:
										clnv += ", "
										var_syns += clnv
									if var_syns != '' and var_syns[-2:] == ", ":
										var_syns = var_syns[:-2]
									elif var_syns == '':
										var_syns = None
							vep_clinical_result["var_synonyms"] = var_syns

							vep_clinical_results.append(vep_clinical_result)

				if vep_initials:
					vep_i_df = pandas.DataFrame(vep_initials)
				else:
					vep_i_df = pandas.DataFrame(
						None, index=[0], columns=["VEP_input", "variant_classification"])
				if vep_transcript_results:
					vep_t_df = pandas.DataFrame(vep_transcript_results)
				else:
					vep_t_df = pandas.DataFrame(
						None, index=[0], columns=["hgvsc", "hgvsp", "protein_id", "transcript_id",
												  "amino_acids", "codons", "polyphen_score",
												  "polyphen_prediction", "sift_score", "sift_prediction",
												  "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
												  "protein_start", "consequence_terms"])
				if vep_clinical_results:
					vep_c_df = pandas.DataFrame(vep_clinical_results)
				else:
					vep_c_df = pandas.DataFrame(
						None, index=[0], columns=["clinical_allele", "clinical_id",
												  "clinical_significance", "var_synonyms"])
				if vep_most_severe_results:
					vep_m_df = pandas.DataFrame(vep_most_severe_results)
				else:
					vep_m_df = pandas.DataFrame(None, index=[0], columns = ["most_severe_consequence"])

				vep_t_df["merging"], vep_c_df["merging"], vep_i_df["merging"], vep_m_df["merging"] = 1, 1, 1, 1
				vep_df1 = pandas.merge(vep_t_df, vep_c_df, on="merging")
				vep_df2 = pandas.merge(vep_df1, vep_m_df, on="merging")
				vep_df = pandas.merge(vep_df2, vep_i_df, on="merging").drop(columns=["merging"])

				if vep_df.empty is False: vep_dfs.append(vep_df)

			VEP_df = pandas.concat(vep_dfs)

			if VEP_df.empty is False:
				VEP_df["Hugo_Symbol"] = row["Hugo_Symbol"]
				VEP_df["CRISPR_PAM_Sequence"] = row["CRISPR_PAM_Sequence"]
				VEP_df["CRISPR_PAM_Location"] = row["CRISPR_PAM_Location"]
				VEP_df["gRNA_Target_Sequence"] = row["gRNA_Target_Sequence"]
				VEP_df["gRNA_Target_Location"] = row["gRNA_Target_Location"]
				VEP_df["Edit_Position"] = row["Edit_Location"]
				VEP_df["Edit_Type"] = row["Edit_Type"]
				VEP_df["Direction"] = row.Direction
				VEP_df["Transcript_ID"] = VEP_df.apply(
					lambda x: x.transcript_id
					if x.transcript_id is not None and pandas.isna(x.transcript_id) is False and
					   type(x.transcript_id) != float else None, axis=1)

				VEP_df["Exon_ID"] = row["Exon_ID"]
				VEP_df["guide_in_CDS"] = row["guide_in_CDS"]
				VEP_df["Edit_in_Exon"] = row["Edit_in_Exon"]
				VEP_df["Edit_in_CDS"] = row["Edit_in_CDS"]
				VEP_df["guide_on_mutation"] = row["guide_on_mutation"]
				VEP_df["guide_change_mutation"] = row["guide_change_mutation"]
				VEP_df["# Edits/guide"] = row["# Edits/guide"]
				VEP_df["Poly_T"] = row["Poly_T"]
				VEP_df["cDNA_Change"] = VEP_df.apply(
					lambda x: x.hgvsc.split(":")[1].split(".")[1]
					if x.hgvsc is not None and pandas.isna(x.hgvsc) is False and
					   type(x.hgvsc) != float else None, axis=1)

				VEP_df["Protein_ID"] = VEP_df.apply(
					lambda x: x.protein_id if x.protein_id is not None else None, axis=1)

				VEP_df["Protein_Change"] = VEP_df.apply(
					lambda x: x.hgvsp.split(":")[1].split(".")[1]
					if x.hgvsp is not None and pandas.isna(x.hgvsp) is False and type(x.hgvsp) != float else None,
					axis=1)

				VEP_df["Protein_Position"] = VEP_df.apply(
					lambda x: extract_hgvsp(hgvsp = x.hgvsp, which="position"), axis=1)

				VEP_df["Edited_AA"] = VEP_df.apply(
					lambda x: extract_hgvsp(hgvsp = x.hgvsp, which="old_aa")
					if x.amino_acids is not None and pandas.isna(x.amino_acids) is False and
					   type(x.amino_acids) != float else None, axis=1)

				VEP_df["New_AA"] = VEP_df.apply(
					lambda x: extract_hgvsp(hgvsp = x.hgvsp, which="new_aa")
					if x.amino_acids is not None and pandas.isna(x.amino_acids) is False and
					   type(x.amino_acids) != float else None, axis=1)

				VEP_df["Edited_AA_Prop"] = VEP_df.apply(
					lambda x: aa_chem[x.Edited_AA]
					if x.Edited_AA is not None and x.Edited_AA in aa_chem.keys()
					   and len(x.Edited_AA) == 1 else (
						";".join([aa_chem[i] for i in x.Edited_AA.split(";") if i in aa_chem.keys()])
						if x.Edited_AA is not None and len(x.Edited_AA) > 1 else None), axis=1)

				VEP_df["New_AA_Prop"] = VEP_df.apply(
					lambda x: aa_chem[x.New_AA]
					if x.New_AA is not None and x.New_AA in aa_chem.keys() and len(x.New_AA) == 1 else (
						";".join([aa_chem[i] for i in x.New_AA.split(";") if i in aa_chem.keys()])
						if x.New_AA is not None and len(x.New_AA) > 1 else None), axis=1)

				VEP_df["Edited_Codon"] = VEP_df.apply(
					lambda x: x.codons.split("/")[0]
					if x.codons is not None and pandas.isna(x.codons) is False and
					   type(x.codons) != float else None, axis=1)

				VEP_df["New_Codon"] = VEP_df.apply(
					lambda x: x.codons.split("/")[1]
					if x.codons is not None and pandas.isna(x.codons) is False and
					   type(x.codons) != float else None, axis=1)

				VEP_df["is_Synonymous"] = VEP_df.apply(
					lambda x: True if x.Edited_Codon is not None and x.New_Codon is not None and
									  x.Edited_AA is not None and x.New_AA is not None and
									  x.Edited_Codon != x.New_Codon and
									  x.Edited_AA == x.New_AA
					else (None if x.Edited_Codon is None and x.New_Codon is None or
								  x.Edited_AA is None and x.New_AA is None else False), axis=1)

				VEP_df["proline_addition"] = VEP_df.apply(
					lambda x: True if x.is_Synonymous is not None and x.is_Synonymous == False and
									  x.New_AA is not None and "P" in x.New_AA.split(";") else False, axis=1)

				VEP_df["is_Stop"] = VEP_df.apply(
					lambda x: True if x.New_AA is not None and
									  x.New_AA == "*" and len(x.New_AA) == 1
					else (True if x.New_AA is not None and "*" in x.New_AA and len(x.New_AA) > 1
						  else (None if x.New_AA is None else False)), axis=1)

				VEP_df["swissprot"] = VEP_df.apply(
					lambda x: x.swissprot[0].split(".")[0] if x.swissprot is not None else False, axis=1)

				VEP_df["is_clinical"] = VEP_df.apply(
					lambda x: True if x.clinical_allele is not None and type(x.clinical_allele) != float and
									  pandas.isna(x.clinical_allele) is False else False, axis=1)

				VEP_df["consequence_terms"] = VEP_df.apply(
					lambda x: ";".join(x.consequence_terms) if type(x.consequence_terms) == list
					else x.consequence_terms, axis=1)

				VEP_df = VEP_df[["Hugo_Symbol", "CRISPR_PAM_Sequence", "CRISPR_PAM_Location",
								 "gRNA_Target_Sequence", "gRNA_Target_Location",
								 "Edit_Position", "Edit_Type", "Direction", "Transcript_ID",
								 "Exon_ID", "Protein_ID", "guide_in_CDS", "Edit_in_Exon", "Edit_in_CDS",
								 "guide_on_mutation", "guide_change_mutation", "# Edits/guide", "Poly_T",
								 "VEP_input", "hgvsp", "variant_classification", "cDNA_Change", "Edited_Codon",
								 "New_Codon", "Protein_Position", "Protein_Change", "Edited_AA", "Edited_AA_Prop",
								 "New_AA", "New_AA_Prop", "is_Synonymous", "is_Stop", "proline_addition",
								 "swissprot", "polyphen_score", "polyphen_prediction", "sift_score",
								 "sift_prediction", "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
								 "consequence_terms", "is_clinical", "clinical_allele", "clinical_id",
								 "clinical_significance", "var_synonyms", "most_severe_consequence"]]
				VEP_df = VEP_df.drop_duplicates()
				all_vep_dfs.append(VEP_df)

	whole_VEP_df = pandas.concat(all_vep_dfs)
	return whole_VEP_df


def annotate_edits(ensembl_object, vep_df):
	"""
	Adding Uniprot API Information on VEP DF
	:param ensembl_object: The object of the Ensembl from Ensembl API
	:param vep_df: The data frame filled with the information from VEP API
	Ensembl Protein ID to Uniprot IDs (SwissProt/Reviewed)
	:return: analysis_df: The data frame enriched with the information from Uniprot API
	"""

	analysis_dfs = list()

	"""
	Since some of the SwissProt object do not have xref mapping information 
	We are missing SwissProt - No sequence mapping
	"""

	ensembl_seq_mapping = {}
	for p in set(list(vep_df["Protein_ID"])):
		if p is not None:
			seq_mapping = ensembl_object.extract_uniprot_info(p)
			if seq_mapping is not None:
				ensembl_seq_mapping[p] = seq_mapping

	if ensembl_seq_mapping != {}:
		for ind, row in vep_df.iterrows():
			uniprot, domain, ptm, reviewed = None, None, None, None
			if row.Protein_ID is not None and row.Protein_ID in ensembl_seq_mapping.keys() and \
					ensembl_seq_mapping[row.Protein_ID] is not None and ensembl_seq_mapping[row.Protein_ID] != []:
				seq_mapping = ensembl_seq_mapping[row.Protein_ID]
				reviewed = list()
				uniprot_tbc = list()
				if row["swissprot"] is not None and row["swissprot"] in seq_mapping.keys():
					uniprot_tbc = [row["swissprot"]]
				else:
					for uniprot in seq_mapping.keys():
						uniprot_object = Uniprot(uniprotid=uniprot)
						if uniprot_object.reviewed:
							reviewed.append(uniprot)

					uniprot_tbc = list()
					if len(reviewed) > 0:
						uniprot_tbc.extend(reviewed)

				# else:
				# uniprot_tbc.extend(seq_mapping.keys())

				if uniprot_tbc:
					for uniprot in uniprot_tbc:
						# Only one SwissProt
						uniprot_object = Uniprot(uniprotid=uniprot)
						reviewed = uniprot_object.reviewed
						uniprot_object.extract_uniprot_info()
						print(row["Protein_Position"])
						if row["Protein_Position"] is not None and pandas.isna(row["Protein_Position"]) is False:
							ptms, domains = list(), list()
							for position in str(row["Protein_Position"]).split(";"):
								if position is not None and position != "None" and type(position) != float:
									if int(position) in seq_mapping[uniprot].keys():
										dom = uniprot_object.find_domain(
											seq_mapping[uniprot][int(position)], row["Edited_AA"])
										phos = uniprot_object.find_ptm_site(
											"phosphorylation", seq_mapping[uniprot][int(position)],
											row["Edited_AA"])
										meth = uniprot_object.find_ptm_site(
											"methylation", seq_mapping[uniprot][int(position)],
											row["Edited_AA"])
										ubi = uniprot_object.find_ptm_site(
											"ubiquitination", seq_mapping[uniprot][int(position)],
											row["Edited_AA"])
										acet = uniprot_object.find_ptm_site(
											"acetylation", seq_mapping[uniprot][int(position)],
											row["Edited_AA"])

										if dom is not None: domains.append(dom + "-" + position)
										if phos is not None: ptms.append(phos + "-" + position)
										if meth is not None: ptms.append(meth + "-" + position)
										if ubi is not None: ptms.append(ubi + "-" + position)
										if acet is not None: ptms.append(acet + "-" + position)
							if ptms: ptm = ";".join([i for i in ptms])
							if domains: domain = ";".join([i for i in domains])

			df_d = {"Hugo_Symbol": [row["Hugo_Symbol"]], "CRISPR_PAM_Sequence": [row["CRISPR_PAM_Sequence"]],
					"CRISPR_PAM_Location": [row["CRISPR_PAM_Location"]],
					"gRNA_Target_Sequence": [row["gRNA_Target_Sequence"]],
					"gRNA_Target_Location": [row["gRNA_Target_Location"]], "Edit_Position": [row["Edit_Position"]],
					"Edit_Type": [row["Edit_Type"]], "Direction": [row["Direction"]],
					"Transcript_ID": [row["Transcript_ID"]], "Exon_ID": [row["Exon_ID"]],
					"Protein_ID": [row["Protein_ID"]], "guide_in_CDS": [row["guide_in_CDS"]],
					"Edit_in_Exon": [row["Edit_in_Exon"]], "Edit_in_CDS": [row["Edit_in_CDS"]],
					"guide_on_mutation" : [row["guide_on_mutation"]],
					"guide_change_mutation" : [row["guide_change_mutation"]],
					"# Edits/guide": [row["# Edits/guide"]], "Poly_T": [row["Poly_T"]],
					"VEP_input": [row["VEP_input"]], "hgvsp": [row["hgvsp"]],
					"variant_classification": [row["variant_classification"]], "cDNA_Change": [row["cDNA_Change"]],
					"Edited_Codon": [row["Edited_Codon"]], "New_Codon": [row["New_Codon"]],
					"Protein_Position": [row["Protein_Position"]], "Protein_Change": [row["Protein_Change"]],
					"Edited_AA": [row["Edited_AA"]], "Edited_AA_Prop": [row["Edited_AA_Prop"]],
					"New_AA": [row["New_AA"]], "New_AA_Prop": [row["New_AA_Prop"]],
					"is_Synonymous": [row["is_Synonymous"]], "is_Stop": [row["is_Stop"]],
					"proline_addition": [row["proline_addition"]], "polyphen_score": [row["polyphen_score"]],
					"polyphen_prediction": [row["polyphen_prediction"]], "sift_score": [row["sift_score"]],
					"sift_prediction": [row["sift_prediction"]], "cadd_phred": [row["cadd_phred"]],
					"cadd_raw": [row["cadd_raw"]], "lof": [row["lof"]], "impact": [row["impact"]],
					"blosum62": [row["blosum62"]], "consequence_terms": [row["consequence_terms"]],
					"is_clinical": [row["is_clinical"]], "clinical_allele": [row["clinical_allele"]],
					"clinical_id": [row["clinical_id"]], "clinical_significance": [row["clinical_significance"]],
					"var_synonyms" : [row["var_synonyms"]], "most_severe_consequence" : [row["most_severe_consequence"]],
					"Domain": [domain], "PTM": [ptm], "SwissProt_VEP": [row["swissprot"]], "Uniprot": [uniprot]}

			df = pandas.DataFrame.from_dict(df_d)
			analysis_dfs.append(df)

	if analysis_dfs:
		analysis_df = pandas.concat(analysis_dfs, ignore_index=True)
	else:
		analysis_df = None

	return analysis_df


def extract_pis(pis):
	"""
	Curating the protein interaction sites from the YULab data
	:param pis: Interaction string from YUlab data (row --> P_IRES)
	:return: List of interacting uniprot indices
	"""
	sites = list()
	if pis != "[]":
		for site in pis.split(","):
			if site[0] == "[" and site[-1] != "]":
				s_first = site[1:]
				if len(s_first.split("-")) > 1:
					for s in list(range(int(s_first.split("-")[0]), int(s_first.split("-")[1]) + 1)):
						sites.append(int(s))
				else:
					sites.append(int(s_first))
			elif site[-1] == "]" and site[0] != "[":
				s_last = site[:-1]
				if len(s_last.split("-")) > 1:
					for s in list(range(int(s_last.split("-")[0]), int(s_last.split("-")[1]) + 1)):
						sites.append(int(s))
				else:
					sites.append(int(s_last))
			elif site[-1] == "]" and site[0] == "[":
				s_only = site[1:-1]
				if len(s_only.split("-")) > 1:
					for s in list(range(int(s_only.split("-")[0]), int(s_only.split("-")[1]) + 1)):
						sites.append(int(s))
				else:
					sites.append(int(s_only))
			else:
				if len(site.split("-")) > 1:
					for s in list(range(int(site.split("-")[0]), int(site.split("-")[1]) + 1)):
						sites.append(int(s))
				else:
					sites.append(int(site))
		sites.sort()
		return sites
	else:
		return None


def collect_pis(uniprot):
	"""
	Collecting protein interaction position for a given uniprot id
	:param uniprot: Uniprot ID
	:return: positional dictionary specify the position and their source and partner (PDB/I3D/ECLAIR)
	"""
	pis_dict = dict()
	df1 = yulab[yulab.P1 == uniprot]
	df2 = yulab[yulab.P2 == uniprot]

	if len(df1.index) != 0:
		for partner, partner_df in df1.groupby(["P2"]):
			for p_ind, p_row in partner_df.iterrows():
				source = p_row["Source"]
				interface_indices = extract_pis(p_row["P1_IRES"])
				if interface_indices is not None:
					for ind in interface_indices:
						if ind not in pis_dict.keys():
							pis_dict[ind] = [{"partner": partner, "source": source}]
						else:
							t = pis_dict[ind]
							if {"partner": partner, "source": source} not in t:
								t.append({"partner": partner, "source": source})
							pis_dict[ind] = t
	if len(df2.index) != 0:
		for partner, partner_df in df1.groupby(["P1"]):
			for p_ind, p_row in partner_df.iterrows():
				source = p_row["Source"]
				interface_indices = extract_pis(p_row["P2_IRES"])
				if interface_indices is not None:
					for ind in interface_indices:
						if ind not in pis_dict.keys():
							pis_dict[ind] = [{"partner": partner, "source": source}]
						else:
							t = pis_dict[ind]
							if {"partner": partner, "source": source} not in t:
								t.append({"partner": partner, "source": source})
							pis_dict[ind] = t
	return pis_dict


def disrupt_interface(uniprot, pos):
	"""
	Checking if the given position disrupts the interfaces in the given uniprot
	:param uniprot: Uniprot ID
	:param pos: Uniprot index
	:return: True/False, effected PDB partners, effected I3D partners, effected ECLAIR partners
	"""
	d = collect_pis(uniprot)
	if pos in d.keys():
		pdb_partner_list = list()
		i3d_partner_list = list()
		eclair_partner_list = list()
		for k in d[pos]:
			if k["source"] == "PDB":
				if k["partner"] not in pdb_partner_list: pdb_partner_list.append(k["partner"])
			elif k["source"] == "I3D":
				if k["partner"] not in i3d_partner_list: i3d_partner_list.append(k["partner"])
			elif k["source"] == "ECLAIR":
				if k["partner"] not in eclair_partner_list: eclair_partner_list.append(k["partner"])
		if len(pdb_partner_list) == 0:
			pdb = None
		else:
			pdb = ",".join(pdb_partner_list)
		if len(i3d_partner_list) == 0:
			i3d = None
		else:
			i3d = ",".join(i3d_partner_list)
		if len(eclair_partner_list) == 0:
			eclair = None
		else:
			eclair = ",".join(eclair_partner_list)
		return pdb, i3d, eclair
	else:
		return None, None, None


def annotate_interface(annotated_edit_df):
	"""
	Add Interactome Insider protein interface information for edgetic perturbation.
	:param annotated_edit_df: Data frame created with annotate_edits
	:return: Added interface annotation on edit table
	"""
	global yulab
	server_url = "https://www.ebi.ac.uk/proteins/api/proteins?"
	df = annotated_edit_df.copy()
	df["is_disruptive_interface_EXP"] = None
	df["is_disruptive_interface_MOD"] = None
	df["is_disruptive_interface_PRED"] = None
	df["disrupted_PDB_int_partners"] = None
	df["disrupted_I3D_int_partners"] = None
	df["disrupted_Eclair_int_partners"] = None
	df["disrupted_PDB_int_genes"] = None
	df["disrupted_I3D_int_genes"] = None
	df["disrupted_Eclair_int_genes"] = None
	for group, group_df in df.groupby(["Uniprot", "Protein_Position"]):
		if group[1] is not None and group[1] != "None" and pandas.isna(group[1]) == False:
			if group[0] in list(yulab.P1) or group[0] in list(yulab.P2):
				all_pdb_partners, all_i3d_partners, all_eclair_partners = list(), list(), list()
				if len(group[1].split(";")) == 1:
					pdb_partners, i3d_partners, eclair_partners = disrupt_interface(
						uniprot=group[0], pos=int(group[1]))
					if pdb_partners is not None:
						all_pdb_partners.append(pdb_partners)
					if i3d_partners is not None:
						all_i3d_partners.append(i3d_partners)
					if eclair_partners is not None:
						all_eclair_partners.append(eclair_partners)
				else:
					for pos in group[1].split(";"):
						pdb_partners, i3d_partners, eclair_partners = disrupt_interface(
							uniprot=group[0], pos=int(pos))
						if pdb_partners is not None:
							all_pdb_partners.append(pdb_partners)
						if i3d_partners is not None:
							all_i3d_partners.append(i3d_partners)
						if eclair_partners is not None:
							all_eclair_partners.append(eclair_partners)

				if all_pdb_partners:
					df.loc[list(group_df.index), "is_disruptive_interface_EXP"] = True
					df.loc[list(group_df.index), "disrupted_PDB_int_partners"] = ";".join(all_pdb_partners)
					for uniprot in all_pdb_partners:
						api_url = "offset=0&size=-1&accession=%s&reviewed=true&isoform=0" %uniprot

						api_request = requests.get(server_url + api_url,
												   headers={"Accept": "application/json"})

						# Check the response of the server for the request
						genes = list()
						if api_request.status_code == 200:
							for i in api_request.json():
								for k in i["gene"]: gene_name = k["name"]["value"]
								if gene_name not in genes:
									genes.append(gene_name)
							genes = ";".join(genes)
						else:
							genes = None
					df.loc[list(group_df.index), "disrupted_PDB_int_genes"] = genes
				else:
					df.loc[list(group_df.index), "is_disruptive_interface_EXP"] = False
					df.loc[list(group_df.index), "disrupted_PDB_int_partners"] = None
					df.loc[list(group_df.index), "disrupted_PDB_int_genes"] = None
				if all_i3d_partners:
					df.loc[list(group_df.index), "is_disruptive_interface_MOD"] = True
					df.loc[list(group_df.index), "disrupted_I3D_int_partners"] = ";".join(all_i3d_partners)
					for uniprot in all_i3d_partners:
						api_url = "offset=0&size=-1&accession=%s&reviewed=true&isoform=0" %uniprot

						api_request = requests.get(server_url + api_url,
												   headers={"Accept": "application/json"})

						# Check the response of the server for the request
						genes = list()
						if api_request.status_code == 200:
							for i in api_request.json():
								for k in i["gene"]: gene_name = k["name"]["value"]
								if gene_name not in genes:
									genes.append(gene_name)
							genes = ";".join(genes)
						else:
							genes = None
					df.loc[list(group_df.index), "disrupted_I3D_int_genes"] = genes
				else:
					df.loc[list(group_df.index), "is_disruptive_interface_MOD"] = False
					df.loc[list(group_df.index), "disrupted_I3D_int_partners"] = None
					df.loc[list(group_df.index), "disrupted_I3D_int_genes"] = None
				if all_eclair_partners:
					df.loc[list(group_df.index), "is_disruptive_interface_PRED"] = True
					df.loc[list(group_df.index), "disrupted_Eclair_int_partners"] = ";".join(all_eclair_partners)
					for uniprot in all_eclair_partners:
						api_url = "offset=0&size=-1&accession=%s&reviewed=true&isoform=0" %uniprot

						api_request = requests.get(server_url + api_url,
												   headers={"Accept": "application/json"})

						# Check the response of the server for the request
						genes = list()
						if api_request.status_code == 200:
							for i in api_request.json():
								for k in i["gene"]: gene_name = k["name"]["value"]
								if gene_name not in genes:
									genes.append(gene_name)
							genes = ";".join(genes)
						else:
							genes = None
					df.loc[list(group_df.index), "disrupted_Eclair_int_genes"] = genes
				else:
					df.loc[list(group_df.index), "is_disruptive_interface_PRED"] = False
					df.loc[list(group_df.index), "disrupted_Eclair_int_partners"] = None
					df.loc[list(group_df.index), "disrupted_Eclair_int_genes"] = None
	return df


def rename_mutational_consequences(mutation_consequence):
	consequence = list()
	for consq in mutation_consequence.split(";"):
		if consq == "missense_variant":
			consequence.append("missense")
		if consq == "missense_mutation":
			consequence.append("missense")
		elif consq == "missense_variant_splice_region_variant":
			consequence.append("splice variant")
		elif consq == "splice_region_variant":
			consequence.append("splice variant")
		elif consq == "stop_retained_variant":
			consequence.append("synonymous")
		elif consq == "synonymous_variant":
			consequence.append("synonymous")
		elif consq == "splice_region_variant_synonymous_variant":
			consequence.append("splice variant")
		elif consq == "splice_acceptor_variant":
			consequence.append("splice variant")
		elif consq == "splice_donor_variant":
			consequence.append("splice variant")
		elif consq == "splice_region_variant_intron_variant":
			consequence.append("splice variant")
		elif consq == "splice_region_variant,intron_variant":
			consequence.append("splice variant")
		elif consq == "splice_donor_region_variant_intron_variant":
			consequence.append("splice variant")
		elif consq == "splice_polypyrimidine_tract_variant_intron_variant":
			consequence.append("splice variant")
		elif consq == "splice_polypyrimidine_tract_variant_splice_region_variant_intron_variant":
			consequence.append("splice variant")
		elif consq == "splice_donor_5th_base_variant_intron_variant":
			consequence.append("splice variant")
		elif consq == "downstream_gene_variant":
			consequence.append("UTR")
		elif consq == "stop_gained_splice_region_variant":
			consequence.append("stop codon")
		elif consq == "stop_gained,splice_region_variant":
			consequence.append("stop codon")
		elif consq == "start_lost":
			consequence.append("start lost")
		elif consq == "stop_gained_start_lost":
			consequence.append("stop codon")
		elif consq == "upstream_gene_variant":
			consequence.append("promoter")
		elif consq == "intron_variant":
			consequence.append("intron")
		elif consq == "5_prime_UTR_variant":
			consequence.append("5'UTR")
		elif consq == "stop_gained":
			consequence.append("stop codon")
	consequence = ";".join(consequence)
	return consequence


def select_severe_effects(mutation_consequence):
	if mutation_consequence is None or pandas.isna(mutation_consequence) or mutation_consequence == "":
		return ""
	elif "stop codon" in mutation_consequence.split(";"):
		return "stop codon"
	elif "start lost" in mutation_consequence.split(";"):
		return "start lost"
	elif "splice variant" in mutation_consequence.split(";"):
		return "splice variant"
	elif "missense" in mutation_consequence.split(";"):
		return "missense"
	elif "UTR" in mutation_consequence.split(";"):
		return "UTR"
	elif "synonymous" in mutation_consequence.split(";"):
		return "synonymous"


def summarise_3di(list_of_partners):
	all_partners = list()
	if list_of_partners is not None:
		for partner_list in list_of_partners:
			if partner_list is not None and pandas.isna(partner_list) == False and partner_list != []:
				for partner in partner_list.split(";"):
					if partner not in all_partners:
						all_partners.append(partner)
	if len(all_partners) > 0:
		return ";".join(all_partners)
	else: return None


def summarise_guides(last_df):

	summary_df = pandas.DataFrame(index = list(range(0, len(last_df.groupby(["CRISPR_PAM_Sequence"])))),
		columns = ["Hugo_Symbol", "CRISPR_PAM_Sequence", "CRISPR_PAM_Location", "gRNA_Target_Sequence",
				   "gRNA_Target_Location", "Edit_Position", "Direction", "Transcript_ID", "Exon_ID", "Protein_ID",
				   "guide_in_CDS", "Edit_in_Exon", "Edit_in_CDS", "guide_on_mutation", "guide_change_mutation",
				   "# Edits/guide", "Poly_T", "cDNA_Change", "Protein_Position", "Protein_Change", "is_stop",
				   "is_synonymous", "proline_addition", "SwissProt_VEP",  "consequence", "most_severe_consequence",
				   "is_clinical", "clinical_allele", "clinical_id", "clinical_significance", "Domain", "PTM",
				   "SwissProt_VEP", "Uniprot", "is_disruptive_interface_EXP", "disrupted_PDB_int_partners",
				   "disrupted_PDB_int_genes", "is_disruptive_interface_MOD", "disrupted_I3D_int_partners",
				   "disrupted_I3D_int_genes", "is_disruptive_interface_PRED", "disrupted_Eclair_int_partners",
				   "disrupted_Eclair_int_genes"])

	i = 0
	for guide, guide_df in last_df.groupby(["CRISPR_PAM_Sequence"]):
		summary_df.loc[i, "Hugo_Symbol"] = ";".join([x for x in list(guide_df.Hugo_Symbol.unique()) if x is not None])
		summary_df.loc[i, "CRISPR_PAM_Sequence"] = ";".join([x for x in list(guide_df.CRISPR_PAM_Sequence.unique()) if x is not None])
		summary_df.loc[i, "CRISPR_PAM_Location"] = ";".join([x for x in list(guide_df.CRISPR_PAM_Location.unique()) if x is not None])
		summary_df.loc[i, "gRNA_Target_Sequence"] = ";".join([x for x in list(guide_df.gRNA_Target_Sequence.unique()) if x is not None])
		summary_df.loc[i, "gRNA_Target_Location"] = ";".join([x for x in list(guide_df.gRNA_Target_Location.unique()) if x is not None])
		summary_df.loc[i, "Edit_Position"] = ";".join([str(x) for x in list(guide_df.Edit_Position.unique()) if x is not None])
		summary_df.loc[i, "Direction"] = ";".join([x for x in list(guide_df.Direction.unique()) if x is not None])
		summary_df.loc[i, "Transcript_ID"] = ";".join([x for x in list(guide_df.Transcript_ID.unique()) if x is not None])
		summary_df.loc[i, "Exon_ID"] = ";".join([x for x in list(guide_df.Exon_ID.unique()) if x is not None])
		summary_df.loc[i, "Protein_ID"] = ";".join([x for x in list(guide_df.Protein_ID.unique()) if x is not None])
		summary_df.loc[i, "guide_in_CDS"] = True if True in guide_df.guide_in_CDS.unique() else False
		summary_df.loc[i, "Edit_in_Exon"] = True if True in guide_df.Edit_in_Exon.unique() else False
		summary_df.loc[i, "Edit_in_CDS"] = True if True in guide_df.Edit_in_CDS.unique() else False
		summary_df.loc[i, "guide_on_mutation"] = True if True in guide_df.guide_on_mutation.unique() else False
		summary_df.loc[i, "guide_change_mutation"] = True if True in guide_df.guide_change_mutation.unique() else False
		summary_df.loc[i, "# Edits/guide"] = guide_df["# Edits/guide"].unique()[0]
		summary_df.loc[i, "Poly_T"] = True if True in guide_df.Poly_T.unique() else False
		summary_df.loc[i, "cDNA_Change"] = ";".join([x for x in list(guide_df.cDNA_Change.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "Protein_Position"] = ";".join([str(x) for x in list(guide_df.Protein_Position.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "Protein_Change"] = ";".join([x for x in list(guide_df.Protein_Change.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "SwissProt_VEP"] = ";".join([x for x in list(guide_df.SwissProt_VEP.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "consequence"] = ";".join([
			select_severe_effects(x)
			for x in [select_severe_effects(i)
					  for i in [rename_mutational_consequences(c)
								for c in [x for x in list(guide_df.consequence_terms.unique()) if x is not None and pandas.isna(x) == False]]]])
		summary_df.loc[i, "most_severe_consequence"] = \
			select_severe_effects(";".join([rename_mutational_consequences(x)
											for x in guide_df.most_severe_consequence.unique()
											if x is not None and pandas.isna(x) == False]))
		summary_df.loc[i, "is_stop"] = True if "stop codon" in guide_df.most_severe_consequence.unique() else False
		summary_df.loc[i, "is_synonymous"] = True if "synonymous" in guide_df.most_severe_consequence.unique() else False
		summary_df.loc[i, "proline_addition"] = True if True in guide_df.proline_addition.unique() else False
		summary_df.loc[i, "is_clinical"] = True if True in guide_df.is_clinical.unique() else False
		summary_df.loc[i, "clinical_allele"] = ";".join([x for x in list(guide_df.clinical_allele.unique())
														 if x is not None and type(x) != float])
		summary_df.loc[i, "clinical_id"] = ";".join([x for x in list(guide_df.clinical_id.unique())
													 if x is not None and type(x) != float])
		summary_df.loc[i, "clinical_significance"] = ";".join([x for x in list(guide_df.clinical_significance.unique())
															   if x is not None and type(x) != float])
		summary_df.loc[i, "Domain"] = ";".join([x for x in list(guide_df.Domain.unique())
												if x is not None and type(x) != float])
		summary_df.loc[i, "PTM"] = ";".join([x for x in list(guide_df.PTM.unique())
											 if x is not None and type(x) != float])
		summary_df.loc[i, "SwissProt_VEP"] = ";".join([x for x in list(guide_df.SwissProt_VEP.unique())
													   if x is not None and type(x) != float])
		summary_df.loc[i, "Uniprot"] = ";".join([x for x in list(guide_df.Uniprot.unique())
												 if x is not None and type(x) != float])
		summary_df.loc[i, "is_disruptive_interface_EXP"] = True if True in guide_df.is_disruptive_interface_EXP.unique() else False
		summary_df.loc[i, "disrupted_PDB_int_partners"] = summarise_3di(list(guide_df.disrupted_PDB_int_partners.unique()))
		summary_df.loc[i, "disrupted_PDB_int_genes"] = summarise_3di(list(guide_df.disrupted_PDB_int_genes.unique()))
		summary_df.loc[i, "is_disruptive_interface_MOD"] = True if True in guide_df.is_disruptive_interface_MOD.unique() else False
		summary_df.loc[i, "disrupted_I3D_int_partners"] = summarise_3di(list(guide_df.disrupted_I3D_int_partners.unique()))
		summary_df.loc[i, "disrupted_I3D_int_genes"] = summarise_3di(list(guide_df.disrupted_I3D_int_genes.unique()))
		summary_df.loc[i, "is_disruptive_interface_PRED"] = True if True in guide_df.is_disruptive_interface_PRED.unique() else False
		summary_df.loc[i, "disrupted_Eclair_int_partners"] = summarise_3di(list(guide_df.disrupted_Eclair_int_partners.unique()))
		summary_df.loc[i, "disrupted_Eclair_int_genes"] = summarise_3di(list(guide_df.disrupted_Eclair_int_genes.unique()))
		i += 1

	return summary_df

"""
def extract_offtargets(guide_df):

	guide = "GTCCTTCCCCCAATCCCCTCAGG"
	m = re.findall(guide, chr1_n[0])
	m = re.findall("(%s){s<=1}" %guide, chr1_n)
"""


###########################################################################################
# Execution


def main():
	"""
	Run whole script with the input from terminal
	:return:
	"""
	print("""
------------------------------------------------------                                                                                         
                B E s t i m a t e                                      

            Wellcome Sanger Institute                                  
------------------------------------------------------
    """)
	if args["PROTEIN"]:
		protein = True
	else:
		protein = False

	if args["MUTATION"]:
		mutations = [args["MUTATION"]]
	else:
		if args["MUTATION_FILE"]:
			mutations = list()
			for line in args["MUTATION_FILE"].readlines():
				mutations.append(line.strip())
		else:
			mutations = None


	print("""
The given arguments are:\nGene: %s\nAssembl: %s\nEnsembl transcript ID: %s\nPAM sequence: %s\nPAM window: %s
Protospacer length: %s\nActivity window: %s\nEdited nucleotide: %s\nNew nucleotide: %s\nVEP and Uniprot analysis: %s\nMutation on genome: %s """
		  % (args["GENE"], args["ASSEMBLY"], args["TRANSCRIPT"], args["PAMSEQ"], args["PAMWINDOW"], args["PROTOLEN"],
			 args["ACTWINDOW"], args["EDIT"], args["EDIT_TO"], protein, mutations))

	print("""\n
------------------------------------------------------ 
              Ensembl Gene Information
------------------------------------------------------ 
    \n""")

	ensembl_obj = Ensembl(hugo_symbol=args["GENE"], assembly=args["ASSEMBLY"])
	ensembl_obj.extract_gene_id()

	if ensembl_obj.gene_id == '': sys.exit("No corresponding Ensembl Gene ID could be found!")

	if mutations:
		ensembl_obj.extract_sequence(ensembl_obj.gene_id, mutations=mutations)
	else:
		ensembl_obj.extract_sequence(ensembl_obj.gene_id)

	if ensembl_obj.gene_range[0] < ensembl_obj.gene_range[1]:
		ensembl_obj.extract_info(chromosome=ensembl_obj.chromosome,
								 loc_start=ensembl_obj.gene_range[0],
								 loc_end=ensembl_obj.gene_range[1])
	else:
		ensembl_obj.extract_info(chromosome=ensembl_obj.chromosome,
								 loc_start=ensembl_obj.gene_range[1],
								 loc_end=ensembl_obj.gene_range[0])

	print("""\n
------------------------------------------------------
                gRNAs - Targetable Sites
------------------------------------------------------
    \n""")
	path = ""
	if args["OUTPUT_PATH"][-1] == "/":
		path = args["OUTPUT_PATH"]
	else:
		path = args["OUTPUT_PATH"] + "/"

	if args["OUTPUT_FILE"] + "_crispr_df.csv" not in os.listdir(path):
		crispr_df = extract_grna_sites(hugo_symbol=args["GENE"], searched_nucleotide=args["EDIT"],
									   pam_window=[int(args["PAMWINDOW"].split("-")[0]),
												   int(args["PAMWINDOW"].split("-")[1])],
									   activity_window=[int(args["ACTWINDOW"].split("-")[0]),
														int(args["ACTWINDOW"].split("-")[1])],
									   pam_sequence=args["PAMSEQ"], protospacer_length=args["PROTOLEN"],
									   ensembl_object=ensembl_obj)

		if len(crispr_df.index) != 0: print("CRISPR Data Frame was created!")
		crispr_df.to_csv(path + args["OUTPUT_FILE"] + "_crispr_df.csv", index=False)

		print("CRISPR Data Frame was written in %s as %s\n" % (path, args["OUTPUT_FILE"] + "_crispr_df.csv"))

	else:
		crispr_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_crispr_df.csv")
	print("""\n
------------------------------------------------------
               gRNAs - Editable Sites
------------------------------------------------------
    \n""")
	if args["OUTPUT_FILE"] + "_edit_df.csv" not in os.listdir(path):
		edit_df = find_editable_nucleotide(crispr_df=crispr_df, searched_nucleotide=args["EDIT"],
										   activity_window=[int(args["ACTWINDOW"].split("-")[0]),
															int(args["ACTWINDOW"].split("-")[1])],
										   ensembl_object=ensembl_obj, mutations=mutations)

		if len(edit_df.index) != 0: print("Edit Data Frame was created!")

		edit_df.to_csv(path + args["OUTPUT_FILE"] + "_edit_df.csv", index=False)

		print("Edit Data Frame was written in %s as %s" % (path, args["OUTPUT_FILE"] + "_edit_df.csv\n"))

	else:
		edit_df =pandas.read_csv(path + args["OUTPUT_FILE"] + "_edit_df.csv")

	if args["VEP"]:
		print("""\n
------------------------------------------------------
               Edits - VEP Annotation
------------------------------------------------------
        \n""")
		if args["OUTPUT_FILE"] + "_vep_df.csv" not in os.listdir(path):
			hgvs_df = extract_hgvs(edit_df=edit_df, ensembl_object=ensembl_obj,
								   transcript_id=args["TRANSCRIPT"],
								   edited_nucleotide=args["EDIT"], new_nucleotide=args["EDIT_TO"],
								   activity_window=[int(args["ACTWINDOW"].split("-")[0]),
													int(args["ACTWINDOW"].split("-")[1])],
								   mutations=mutations)

			if hgvs_df is not None and len(hgvs_df.index) != 0:

				whole_vep_df = retrieve_vep_info(hgvs_df=hgvs_df, ensembl_object=ensembl_obj,
												 transcript_id=args["TRANSCRIPT"])
				if len(whole_vep_df.index) != 0:
					print("VEP Data Frame was created!")
					whole_vep_df.to_csv(path + args["OUTPUT_FILE"] + "_vep_df.csv")
					print("VEP Data Frame was written in %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_vep_df.csv"))
				else:
					print("VEP Data Frame cannot be created because it is empty!")
		else:
			whole_vep_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_vep_df.csv")

		if args["PROTEIN"]:
			print("""\n
------------------------------------------------------
		   Edits - Uniprot Annotation
------------------------------------------------------
			\n""")

			if args["OUTPUT_FILE"] + "_protein_df.csv" not in os.listdir(path):
				print("Adding Uniprot ID, corresponding Domain and PTM information..\n")
				if len(whole_vep_df.index) != 0:
					uniprot_df = annotate_edits(ensembl_object=ensembl_obj, vep_df=whole_vep_df)
					print("Uniprot information was added!\n")
				print("""\n
	------------------------------------------------------
		   Edits - 3D Interface Annotation
	------------------------------------------------------
				\n""")

				if uniprot_df is not None and len(uniprot_df.index) != 0:

					print("Adding affected interface and interacting partners..")

					interaction_df = annotate_interface(annotated_edit_df=uniprot_df)

					if interaction_df is not None and len(interaction_df.index) != 0:
						print("Annotation Data Frame was created!")

						interaction_df.to_csv(path + args["OUTPUT_FILE"] + "_protein_df.csv", index=False)

						print("Protein Data Frame was written in %s as %s\n\n" % (path, args["OUTPUT_FILE"] +
																				  "_protein_df.csv"))
					else:
						print("Protein Data Frame cannot be created because it is empty.")
			else:
				interaction_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_protein_df.csv")

				print("""\n
------------------------------------------------------
			Summary Data Frame 
------------------------------------------------------
				\n""")

			if interaction_df is not None and len(interaction_df.index) != 0:
				print("Summarising information..")

				summary_df = summarise_guides(last_df=interaction_df)

				if summary_df is not None and len(summary_df.index) != 0:
					print("Summary Data Frame was created!")

					summary_df.to_csv(path + args["OUTPUT_FILE"] + "_summary_df.csv", index=False)

					print("Summary Data Frame was written in %s as %s\n\n" % (path, args["OUTPUT_FILE"] +
																			  "_summary_df.csv"))
				else:
					print("Summary Data Frame cannot be created because it is empty.")


	return True


###########################################################################################
# Execution

main()

print("""\n
------------------------------------------------------
           The BEstimate analysis finished!
------------------------------------------------------
\n""")





