# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#

# Import necessary packages
import os, sys, pandas, re, argparse, requests, json, itertools, pickle
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


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
						help="The boolean option if user wants to analyse the edits through VEP and Uniprot.")

	# MUTATION INFORMATION
	parser.add_argument("-mutation", dest="MUTATION", default=None,
						help="The mutation on the interested gene that you need to integrate "
							 "into guide and/or annotation analysis")
	parser.add_argument("-mutation_file", dest="MUTATION_FILE", default=None, type=argparse.FileType('r'),
						help="If you have more than one mutations, a file for the mutations on the "
							 "interested gene that you need to integrate into guide and/or annotation analysis")

	# gRNA FLANKING REGIONS
	parser.add_argument("-flank", dest="FLAN", action="store_true",
						help="The boolean option if the user wants to add flanking sequences of the gRNAs")
	parser.add_argument("-flank3", dest="FLAN_3", default="7",
						help="The number of nucleotides in the 3' flanking region")
	parser.add_argument("-flank5", dest="FLAN_5", default="11",
						help="The number of nucleotides in the 5' flanking region")

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

	# OFF TARGETS
	parser.add_argument("-ot", dest="OFF_TARGET", action="store_true",
						help="Whether off targets will be computed or not")
	parser.add_argument("-mm", dest="MISMATCH", default=4,
						help="(If -ot provided) number of maximum mismatches allowed in off targets")
	parser.add_argument("-genome", dest="GENOME", default="Homo_sapiens_GRCh38_dna_sm_all_chromosomes",
						help="(If -ot provided) name of the genome file")

	parsed_input = parser.parse_args()
	input_dict = vars(parsed_input)

	return input_dict


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

	def extract_uniprot(self):

		uniprot_api = "proteins?offset=0&size=-1&accession=%s" % self.uniprotid
		api_request = requests.get(self.server + uniprot_api,
								   headers={"Accept": "application/json"})

		# Check the response of the server for the request
		if api_request.status_code != 200:
			return "No response from UniProt!\n"

		else:
			for i in api_request.json():
				self.reviewed = False if i["info"]["type"] == "TrEMBL" else True
				self.sequence = i["sequence"]["sequence"]
				if "features" in i.keys() and i["features"] != []:
					for ftr in i["features"]:
						if ftr["type"] == "MOD_RES" and ftr["category"] == "PTM":
							if "description" in ftr.keys():
								pos, ptm = ftr["begin"], ftr["description"]
								ptm = ptm.split(";")[0]
								# Phosphorylation
								phos = ptm if re.search(r'Phospho', ptm) or re.search(r'phospho', ptm) else None
								if phos is not None: self.phosphorylation_sites[pos] = phos
								# Methylation
								methy = ptm if re.search(r'Methyl', ptm) or re.search(r'methyl', ptm) else None
								if methy is not None: self.methylation_sites[pos] = methy
								# Ubiquitination
								ubi = ptm if re.search(r'Ub', ptm) or re.search(r'ub', ptm) else None
								if ubi is not None: self.ubiquitination_sites[pos] = ubi
								# Acetylation
								acety = ptm if re.search(r'Ace', ptm) or re.search(r'acetyl', ptm) else None
								if acety is not None: self.acetylation_sites[pos] = acety

						if ftr["category"] == "DOMAINS_AND_SITES":
							if ftr["type"] == "BINDING":
								domain = ftr["ligand"]["name"] + " binding site"
							else:
								if "description" in ftr.keys():
									domain = ftr["description"]
							domain_range = list(range(int(ftr["begin"]), int(ftr["end"])))
							if domain not in self.domains.keys():
								self.domains[domain] = domain_range
							else:
								t = self.domains[domain]
								for p in domain_range:
									if p not in t:
										t.append(p)
								self.domains[domain] = t

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
				if int(protein_edit_location) in domain_range:
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
		self.p_sequence = None

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
				print(seq_request.text.split("\n")[0])
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

		print("Request to Ensembl REST API for sequence information")
		seq_request = requests.get(seq_ensembl,
								   headers={"Content-Type": "text/x-fasta"})
		seq_flan_request = requests.get(seq_flan_ensembl,
										headers={"Content-Type": "text/x-fasta"})

		if seq_request.status_code != 200 and seq_flan_request.status_code != 200:
			print("No response from ensembl sequence!\n")

		# Sequence
		label_line = seq_request.text.split("\n")[0]
		if label_line[0] != "{":
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

			nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

			if mutations is None:

				if self.strand == 1:
					self.right_sequence_analysis = self.sequence
					self.flan_right_sequence_analysis = self.flan_sequence
					self.left_sequence_analysis = "".join([nucleotide_dict[n] for n in self.sequence[::-1]])
					self.flan_left_sequence_analysis = "".join([nucleotide_dict[n] for n in self.flan_sequence[::-1]])

				elif self.strand == -1:
					self.left_sequence_analysis = self.sequence
					self.flan_left_sequence_analysis = self.flan_sequence
					self.right_sequence_analysis = "".join(
						[nucleotide_dict[n] for n in self.left_sequence_analysis[::-1]])
					self.flan_right_sequence_analysis = "".join(
						[nucleotide_dict[n] for n in self.flan_left_sequence_analysis[::-1]])

			else:
				for mutation in mutations:
					if mutation.split(":")[0] == self.chromosome:
						if mutation.split(":")[1].split(".")[0] == "g":
							alteration = mutation.split(":")[1].split(".")[1]
							mutation_location = int(
								re.match("([0-9]+)([a-z]+)", alteration.split(">")[0], re.I).groups()[0])
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
		if self.sequence is not None:
			print("Sequence was retrieved.")

	def extract_gRNA_flan_sequence(self, location, direction, fivep, threep):
		"""
		Annotating flanking regions of the gRNAs
		:param grna: gRNA target sequence
		:param location: Location of the gRNA target sequence
		:param fivep: Number of the nucleotides in the 5' flanking region
		:param threep: Number of the nucleotides in the 3' flanking region
		:return:
		"""

		if direction == "left":
			grna_strand = "-1"
		else:
			grna_strand = "1"

		grna_flan_ensembl = self.server + "/sequence/region/human/%s:%s?expand_3prime=%s;expand_5prime=%s;content-type=text/plain" \
							% (location, grna_strand, threep, fivep)
		grna_flan_request = requests.get(grna_flan_ensembl,
										 headers={"Content-Type": "text/plain"})

		if grna_flan_request.status_code != 200:
			print("No response from ensembl sequence!\n")

		return grna_flan_request.text

	def extract_info(self, chromosome, loc_start, loc_end, transcript=None):

		if loc_start < loc_end:
			ensembl = "/overlap/region/human/%s:%s-%s?feature=transcript;feature=exon;feature=mane;" \
					  "feature=cds" % (
						  chromosome, str(loc_start), str(loc_end))
		else:
			ensembl = "/overlap/region/human/%s:%s-%s?feature=transcript;feature=exon;feature=mane;" \
					  "feature=cds" % (
						  chromosome, str(loc_end), str(loc_start))

		request = requests.get(self.server + ensembl, headers={"Content-Type": "application/json"})

		info_dict = dict()

		if request.status_code != 200:
			print("No response from ensembl!")
		else:
			canonicals = list()
			for output in request.json():
				if transcript is None:
					if output["feature_type"] == "mane" and output["Parent"] == self.gene_id:
						if "refseq_match" in output.keys():
							if output["type"] != "MANE_Plus_Clinical":
								if output["id"].split(".")[0] not in canonicals:
									canonicals.append(output["id"].split(".")[0])

								if output["id"].split(".")[0] in canonicals:
									if output["id"] not in info_dict.keys():
										info_dict[output["id"]] = \
											[{"start": output["start"], "end": output["end"]}]

									else:
										old_val = info_dict[output["id"]]
										if {"start": output["start"], "end": output["end"]} not in old_val:
											old_val.append(
												{"start": output["start"], "end": output["end"]})
											info_dict[output["id"]] = old_val

					elif output["feature_type"] == "transcript" and output["Parent"] == self.gene_id:
						if "is_canonical" in output.keys():
							if output["is_canonical"] == 1:
								if output["id"].split(".")[0] not in canonicals:
									canonicals.append(output["id"].split(".")[0])
						if "source" in output.keys():
							if output["source"] == "ensembl_havana":
								if output["id"].split(".")[0] not in canonicals:
									canonicals.append(output["id"].split(".")[0])
						if output["id"].split(".")[0] in canonicals:
							if output["id"] not in info_dict.keys():
								info_dict[output["id"]] = \
									[{"start": output["start"], "end": output["end"]}]

							else:
								old_val = info_dict[output["id"]]
								if {"start": output["start"], "end": output["end"]} not in old_val:
									old_val.append(
										{"start": output["start"], "end": output["end"]})
									info_dict[output["id"]] = old_val

				else:
					# Selected transcript
					if output["feature_type"] == "transcript" and output["Parent"] == self.gene_id:
						transcript_info = "/lookup/id/%s?expand=1;mane=1" % output["transcript_id"]
						transcript_request = requests.get(self.server + transcript_info,
														  headers={"Content-Type": "application/json"})
						transcript_output = transcript_request.json()
						if output["transcript_id"] not in info_dict.keys():
							info_dict[output["transcript_id"]] = \
								[{"start": transcript_output["start"], "end": transcript_output["end"]}]

						else:
							old_val = info_dict[output["transcript_id"]]
							if {"start": transcript_request["start"], "end": transcript_request["end"]} not in old_val:
								old_val.append(
									{"start": transcript_request["start"], "end": transcript_request["end"]})
								info_dict[output["transcript_id"]] = old_val

			for output in request.json():
				if output["feature_type"] == "cds" and info_dict != {} and output["Parent"] in info_dict.keys():
					for k in range(len(info_dict[output["Parent"]])):
						d = info_dict[output["Parent"]][k]
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
						if d not in info_dict[output["Parent"]]:
							info_dict[output["Parent"]][k] = d

			protein_ids = list()
			swiss_protein_ids = list()
			selected_transcript = list()
			for ids in info_dict.keys():
				for d in info_dict[ids]:
					if "cds" in d.keys():
						protein_ids.extend(d["cds"].keys())
			if protein_ids:
				for p in protein_ids:
					protein_ensembl = "/xrefs/id/{0}?external_db=Uniprot/SWISSPROT%".format(p)
					protein_request = requests.get(self.server + protein_ensembl,
												   headers={"Content-Type": "application/json"})
					for i in protein_request.json():
						if i["dbname"] == "Uniprot/SWISSPROT":
							swiss_protein_ids.append(p)
				if swiss_protein_ids:
					for p in swiss_protein_ids:
						for ids in info_dict.keys():
							for d in info_dict[ids]:
								if "cds" in d.keys():
									if p in d["cds"].keys():
										if ids not in selected_transcript:
											selected_transcript.append(ids)
					else:
						selected_transcript = list(info_dict.keys())
			else:
				selected_transcript = list(info_dict.keys())

			if selected_transcript:
				info_dict2 = {key: info_dict[key] for key in info_dict.keys() if key in selected_transcript}

				for output in request.json():
					if output["feature_type"] == "exon" and info_dict2 != {} and output["Parent"] in info_dict2.keys():
						for d in info_dict2[output["Parent"]]:
							if "exon" not in d.keys():
								d["exon"] = {output["exon_id"]: {"start": output["start"], "end": output["end"]}}
							else:
								if output["exon_id"] not in d["exon"].keys():
									d["exon"][output["exon_id"]] = {"start": output["start"], "end": output["end"]}
				if info_dict2 != {}:
					self.info_dict = info_dict2
					return 1
			else:
				self.info_dict = None
				return 0

	def check_range_info(self, start, end):

		range_locations = list(range(start, end))
		transcripts_exons = dict()
		if self.info_dict != {} and self.info_dict is not None:
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
		if self.info_dict != {} and self.info_dict is not None:
			if transcript_id in self.info_dict.keys():
				for transcript_dict in self.info_dict[transcript_id]:
					if "cds" in transcript_dict.keys():
						for protein, protein_pos in transcript_dict["cds"].items():
							for loc in range_locations:
								if loc in protein_pos:
									in_cds = True

		return in_cds

	def extract_uniprot_info(self, ensembl_pid, uniprot):
		global path
		uniprot = uniprot.split(".")[0]
		protein_ensembl = "/xrefs/id/{0}?external_db=Uniprot/SWISSPROT%".format(ensembl_pid)
		protein_request = requests.get(self.server + protein_ensembl,
									   headers={"Content-Type": "application/json"})
		if protein_request.status_code != 200:
			print("No response from ensembl!")
			return 0
		else:
			seq_mapping = dict()
			for i in protein_request.json():
				if i["primary_id"] == uniprot:
					if i["dbname"] == "Uniprot/SWISSPROT":
						if "ensembl_end" in i.keys() and "ensembl_start" in i.keys() and \
								"xref_end" in i.keys() and "xref_start" in i.keys():
							if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(i["xref_end"]) - int(
									i["xref_start"]):
								# Otherwise, there is an inconsistency --> Not take it
								seq_mapping[uniprot] = {i["ensembl_start"] + k: i["xref_start"] + k
														for k in range(
										int(i["ensembl_end"]) - int(i["ensembl_start"]) + 1)}
							break

					elif i["dbname"] == "Uniprot/SPTREMBL":
						if "ensembl_end" in i.keys() and "ensembl_start" in i.keys() and \
								"xref_end" in i.keys() and "xref_start" in i.keys():
							if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(i["xref_end"]) - int(
									i["xref_start"]):
								# Otherwise, there is an inconsistency --> Not take it
								seq_mapping[uniprot] = {i["ensembl_start"] + k: i["xref_start"] + k
														for k in range(
										int(i["ensembl_end"]) - int(i["ensembl_start"]) + 1)}

			if seq_mapping == dict():
				# Alignment
				uniprot_obj = Uniprot(uniprotid=uniprot)
				uniprot_obj.extract_uniprot()
				uniprot_seq = uniprot_obj.sequence

				if self.p_sequence is not None:
					ensembl_seq = self.p_sequence
				else:
					seq_ensembl = self.server + "/sequence/id/%s?" % ensembl_pid
					seq_request = requests.get(seq_ensembl,
											   headers={"Content-Type": "text/x-fasta"})

					if seq_request.status_code != 200 and seq_flan_request.status_code != 200:
						print("No response from ensembl protein sequence!\n")

					else:
						label_line = seq_request.text.split("\n")[0]
					if label_line[0] != "{":
						self.p_sequence = "".join(seq_request.text.split("\n")[1:])
						ensembl_seq = self.p_sequence

				alignment_f = open(path + args["OUTPUT_FILE"] + "_%s_alignment.txt" % uniprot, "w")
				alignments = list()
				for a in pairwise2.align.globalms(ensembl_seq, uniprot_seq, 2, -1, -1, -0.1):
					alignment_f.write(format_alignment(*a, full_sequences=True))
					alignments.append(format_alignment(*a, full_sequences=True))
				alignment_f.close()

				ens = alignments[0].split("\n")[:3]
				alignment_df = pandas.DataFrame(0, index=list(range(len(ens[0]))),
												columns=["e_aa", "e_index"])
				alignment_df["e_aa"] = [aa for aa in ens[0]]

				alignment_dict = {i: alignments[i] for i in range(len(alignments))}
				for a_num, align in alignment_dict.items():
					alignment_df["u_aa_%d" % a_num] = None
					alignment_df["u_index_%d" % a_num] = None
					alignment_df["mismatch_%d" % a_num] = None
					align_list = align.split("\n")[:3]
					uni_index = 0
					ens_index = 0
					for i in range(len(align_list[0])):
						ens = align_list[0][i]
						a_style = align_list[1][i]
						uni = align_list[2][i]

						if a_style == "|":
							# Match
							alignment_df.loc[i, "e_index"] = ens_index
							alignment_df.loc[i, "u_aa_%d" % a_num] = uni
							alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
							alignment_df.loc[i, "mismatch_%d" % a_num] = False
							uni_index += 1
							ens_index += 1
						elif a_style == ".":
							# Mismatch
							alignment_df.loc[i, "e_index"] = ens_index
							alignment_df.loc[i, "u_aa_%d" % a_num] = uni
							alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
							alignment_df.loc[i, "mismatch_%d" % a_num] = True
							uni_index += 1
							ens_index += 1
						elif a_style == " ":
							# Gap
							if ens == "-":
								alignment_df.loc[i, "e_index"] = ens_index
								alignment_df.loc[i, "u_aa_%d" % a_num] = uni
								alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
								alignment_df.loc[i, "mismatch_%d" % a_num] = False
								uni_index += 1
							elif uni == "-":
								alignment_df.loc[i, "e_index"] = ens_index
								alignment_df.loc[i, "u_aa_%d" % a_num] = "-"
								alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
								alignment_df.loc[i, "mismatch_%d" % a_num] = False
								ens_index += 1

				alignment_df["inconsistent"] = False
				mismatch_indices = list()
				for i in range(len(alignments)):
					mismatch_indices.extend(list(alignment_df[alignment_df["mismatch_%d" % a_num]].index))
				mismatch_indices = list(set(mismatch_indices))
				different_indices = list()
				missing_indices = list()
				uni_ind_cols = [col for col in alignment_df if col.startswith("u_index")]
				uni_aa_cols = [col for col in alignment_df if col.startswith("u_aa")]
				for ind, row in alignment_df.iterrows():
					indices = list()
					aas = list()
					for col in uni_ind_cols:
						indices.append(row[col])
					for col in uni_aa_cols:
						aas.append(row[col])
						if row[col] == "-":
							if ind not in missing_indices:
								missing_indices.append(ind)
					if len(list(set(indices))) > 1:
						different_indices.append(ind)
					if len(list(set(aas))) > 1:
						if ind not in different_indices:
							different_indices.append(ind)

				flagged_ind = list(set(mismatch_indices).union(set(different_indices).union(set(missing_indices))))
				alignment_df.loc[flagged_ind, "inconsistent"] = True

				alignment_df.to_csv(path + args["OUTPUT_FILE"] + "_%s_alignment_df.csv" % uniprot, index=True)

				alignment_df2 = alignment_df[alignment_df.inconsistent == False]
				seq_mapping[uniprot] = {r["e_index"]: r["u_index_0"] for i, r in alignment_df2.iterrows()}

			if seq_mapping:
				return seq_mapping
			else:
				return None


class Variant:
	def __init__(self, hgvs, gene, strand, transcript):
		self.hgvs, self.hgvsc, self.hgvsp = hgvs, None, None
		self.vep = None
		self.gene, self.strand = gene, strand
		self.transcript = transcript
		self.allele = None
		self.regulatory = None
		self.motif, self.motif_TFs = None, None
		self.variant_class, self.consequence_terms, self.biotype = None, None, None
		self.most_severe_consequence = None
		self.cdna_change, self.cds_position = None, None
		self.old_codon, self.new_codon = None, None
		self.protein_change, self.protein_position = None, None
		self.old_aa, self.new_aa = None, None
		self.old_aa_chem, self.new_aa_chem = None, None
		self.synonymous, self.stop, self.proline = None, None, None
		self.protein, self.swissprot = None, None
		self.polyphen_score, self.polyphen_prediction = None, None
		self.sift_score, self.sift_prediction = None, None
		self.cadd_phred, self.cadd_raw, self.lof = None, None, None
		self.impact, self.blosum62 = None, None
		self.clinical, self.clinical_id, self.clinical_sig = None, None, None
		self.clinvar_id = None
		self.cosmic, self.cosmic_id = None, None
		self.ancestral_populations = None

	def extract_vep_obj(self, vep_json):
		for vep in vep_json:
			if vep["input"] == self.hgvs:
				self.vep = vep

	def extract_hgvsp(self, hgvsp, which):
		aa_3to1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q",
				   "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
				   "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
		if hgvsp is not None:
			protein_change = hgvsp.split("p.")[1]
			if len(protein_change.split("delins")) == 1:
				# SNP
				if len(protein_change.split("=")) == 1:
					if len(protein_change.split("?")) == 1:
						if len(protein_change.split("ext")) == 1:
							if which == "old_aa":
								return aa_3to1[protein_change[:3]]
							if which == "new_aa":
								return aa_3to1[protein_change[-3:]]
							if which == "position":
								return protein_change[3:-3]
						else:
							# Extension for termination or start Ter629GlnextTer1 | Met1ext-5
							if protein_change[:3] == "Ter":
								alteration = protein_change.split("ext")[0]
								extension_amount = int(protein_change.split("ext")[1][3:]) - 1
								if which == "old_aa":
									return aa_3to1[alteration[:3]]
								if which == "new_aa":
									return aa_3to1[alteration[-3:]] + "X%s" % extension_amount + "*"
								if which == "position":
									return alteration[3:-3]
							else:
								if which == "old_aa":
									return aa_3to1[protein_change[:3]]
								if which == "new_aa":
									extension_amount = abs(int(protein_change.split("ext")[1])) - 1
									return aa_3to1[protein_change[:3]] + "X-%s" % extension_amount + aa_3to1[
										protein_change[:3]]
								if which == "position":
									return protein_change.split("ext")[0][3:]

					else:
						# Start codon lost - Met1? | MetAla1_?2
						if which == "old_aa":
							aa1 = list()
							aa_string = re.match("([a-z]+)([0-9]+)", protein_change.split("?")[0], re.I).groups()[0]
							for i in [aa_string[x: x + 3] for x in range(0, len(aa_string), 3)]:
								aa1.append(aa_3to1[i])
							return ";".join(aa1)

						if which == "new_aa":
							if protein_change[-1] == "?" or protein_change[-2] == "?":
								return "-"
						if which == "position":
							return re.match("([a-z]+)([0-9]+)", protein_change.split("?")[0], re.I).groups()[1]

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
					for i in [protein_change.split("delins")[1][x: x + 3] for x in
							  range(0, len(protein_change.split("delins")[1]), 3)]:
						aa2.append(aa_3to1[i])
					return ";".join(aa2)

				if which == "position":
					pos = list()
					for i in protein_change.split("delins")[0].split("_"):
						pos.append(re.match("([a-z]+)([0-9]+)", i, re.I).groups()[1])
					return ";".join(pos)

		else:
			return None

	def extract_consequences(self):
		consequence_terms = list()
		ancestral_populations = list()
		# Dictionary to find the chemical properperty change due to the edit
		aa_chem = {"G": "Non-Polar", "A": "Non-Polar", "V": "Non-Polar", "C": "Polar", "P": "Non-Polar",
				   "L": "Non-Polar", "I": "Non-Polar", "M": "Non-Polar", "W": "Non-Polar", "F": "Non-Polar",
				   "S": "Polar", "T": "Polar", "Y": "Polar", "N": "Polar", "Q": "Polar", "K": "Charged",
				   "R": "Charged", "H": "Charged", "D": "Charged", "E": "Charged", "*": "-"}

		if "allele_string" in self.vep.keys():
			self.allele = self.vep["allele_string"]

		if "most_severe_consequence" in self.vep.keys():
			self.most_severe_consequence = self.vep["most_severe_consequence"]

		if "variant_class" in self.vep.keys():
			self.variant_class = self.vep["variant_class"]

		if "regulatory_feature_consequences" in self.vep.keys():
			for r in self.vep["regulatory_feature_consequences"]:
				if "strand" in r.keys():
					if r["strand"] == self.strand:
						if "regulatory_feature_id" in r.keys():
							self.regulatory = r["regulatory_feature_id"]
						if "consequence_terms" in r.keys():
							for cons_term in r["consequence_terms"]:
								if cons_term not in consequence_terms:
									consequence_terms.append(cons_term)

		if "motif_feature_consequences" in self.vep.keys():
			for m in self.vep["motif_feature_consequences"]:
				if "motif_feature_id" in m.keys():
					self.motif = m["motif_feature_id"]
				if "transcription_factors" in m.keys():
					self.motif_TFs = ", ".join([tf for tf in m["transcription_factors"]])
				if "consequence_terms" in m.keys():
					for cons_term in m["consequence_terms"]:
						if cons_term not in consequence_terms:
							consequence_terms.append(cons_term)

		if "transcript_consequences" in self.vep.keys():
			for t in self.vep["transcript_consequences"]:
				if t["gene_symbol"] == self.gene and t["transcript_id"] == self.transcript:
					if "hgvsc" in t.keys():
						self.hgvsc = t["hgvsc"]
					if "biotype" in t.keys():
						self.biotype = t["biotype"]
					if "hgvsp" in t.keys():
						self.hgvsp = t["hgvsp"]
						self.protein_position = self.extract_hgvsp(hgvsp=self.hgvsp, which="position")
						self.old_aa = self.extract_hgvsp(hgvsp=self.hgvsp, which="old_aa")
						self.new_aa = self.extract_hgvsp(hgvsp=self.hgvsp, which="new_aa")
						self.old_aa_chem = \
							aa_chem[self.old_aa] if self.old_aa is not None and self.old_aa in aa_chem.keys() and \
													len(self.old_aa) == 1 \
								else (";".join([aa_chem[i] for i in self.old_aa.split(";") if i in aa_chem.keys()])
									  if self.old_aa is not None and len(self.old_aa) > 1 else None)
						self.new_aa_chem = \
							aa_chem[self.new_aa] if self.new_aa is not None and self.new_aa in aa_chem.keys() and \
													len(self.new_aa) == 1 \
								else (";".join([aa_chem[i] for i in self.new_aa.split(";") if i in aa_chem.keys()])
									  if self.new_aa is not None and len(self.new_aa) > 1 else None)
					if "protein_id" in t.keys():
						self.protein = t["protein_id"]
					if "amino_acids" in t.keys():
						self.protein_change = t["amino_acids"]
					if "codons" in t.keys():
						self.cdna_change = t["codons"]
						self.old_codon = self.cdna_change.split("/")[0] \
							if self.cdna_change is not None and pandas.isna(self.cdna_change) is False \
							   and type(self.cdna_change) != float else None
						self.new_codon = self.cdna_change.split("/")[1] \
							if self.cdna_change is not None and pandas.isna(self.cdna_change) is False \
							   and type(self.cdna_change) != float else None
					if "cds_start" in t.keys() and "cds_end" in t.keys():
						self.cds_position = str(t["cds_start"]) + "-" + str(t["cds_end"])

					self.synonymous = True if self.old_codon is not None and self.new_codon is not None and \
											  self.old_aa is not None and self.new_aa is not None and \
											  self.old_codon != self.new_codon and self.old_aa == self.new_aa else (
						None if self.old_codon is None and self.new_codon is None or self.old_aa is None and
								self.new_aa is None else False)
					self.proline = True if self.synonymous is not None and self.synonymous == False and \
										   self.new_aa is not None and "P" in self.new_aa.split(";") else False
					self.stop = True if self.new_aa is not None and self.new_aa == "*" and len(self.new_aa) == 1 else (
						True if self.new_aa is not None and "*" in self.new_aa and len(self.new_aa) > 1 else (
							None if self.new_aa is None else False))

					if "swissprot" in t.keys():
						self.swissprot = t["swissprot"][0]
					if "polyphen_score" in t.keys():
						self.polyphen_score = t["polyphen_score"]
					if "polyphen_prediction" in t.keys():
						self.polyphen_prediction = t["polyphen_prediction"]
					if "sift_score" in t.keys():
						self.sift_score = t["sift_score"]
					if "sift_prediction" in t.keys():
						self.sift_prediction = t["sift_prediction"]
					if "cadd_phred" in t.keys():
						self.cadd_phred = t["cadd_phred"]
					if "cadd_raw" in t.keys():
						self.cadd_raw = t["cadd_raw"]
					if "lof" in t.keys():
						self.lof = t["lof"]
					if "impact" in t.keys():
						self.impact = t["impact"]
					if "blosum62" in t.keys():
						self.blosum62 = t["blosum62"]
					if "consequence_terms" in t.keys():
						for cons_term in t["consequence_terms"]:
							if cons_term not in consequence_terms:
								consequence_terms.append(cons_term)

		if consequence_terms:
			self.consequence_terms = ", ".join(consequence_terms)

		if "colocated_variants" in self.vep.keys():
			self.clinical = True
			cosmic_id = list()
			clinvar_id = list()
			for c in self.vep["colocated_variants"]:
				if "allele_string" in c.keys():
					if c["allele_string"] == "COSMIC_MUTATION":
						self.cosmic = True
						if "id" in c.keys():
							if c["id"] not in cosmic_id:
								cosmic_id.append(c["id"])

				if "clin_sig" in c.keys():
					self.clinical_sig = ", ".join([i for i in c["clin_sig"]])
				if "id" in c.keys():
					self.clinical_id = c["id"]

				if "var_synonyms" in c.keys():
					if type(c["var_synonyms"]) == str:
						for clnv in c["var_synonyms"]["ClinVar"]:
							for cl_id in clnv:
								if cl_id not in clinvar_id:
									clinvar_id.append(cl_id)
				if "frequencies" in c.keys():
					for alele, freq_dict in c["frequencies"].items():
						for pop, val in freq_dict.items():
							if val >= 0.01:
								if pop not in ancestral_populations:
									ancestral_populations.append(pop)

			if cosmic_id:
				self.cosmic_id = ", ".join(cosmic_id)
			if clinvar_id:
				self.clinvar_id = ", ".join(clinvar_id)
			if ancestral_populations:
				self.ancestral_populations = ", ".join(ancestral_populations)


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
					   flan, flan_3, flan_5, ensembl_object):
	"""
	Extracting the gRNA targeted sites having editable nucleotide(s) on the interested genes
	:param hugo_symbol: The Hugo Symbol of the interested gene.
	:param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
	:param searched_nucleotide: The interested nucleotide which will be changed with BE
	:param activity_window: The location of the activity windiw on the protospacer sequence.
	:param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
	:param protospacer_length: The length of protospacer.
	:param flan: The boolean parameter if the user wants to add flanking regions of the gRNAs
	:param flan_3: The number of nucleotides which will be added 3' flanking region
	:param flan_5: The number of nucleotides which will be added 5' flanking region
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

	if flan:
		crisprs_df["gRNA_flanking_sequences"] = crisprs_df.apply(
			lambda x: ensembl_object.extract_gRNA_flan_sequence(
				location=x.Location, direction=x.Direction, fivep=flan_5, threep=flan_3), axis=1)
	else:
		crisprs_df["gRNA_flanking_sequences"] = None

	return crisprs_df


def collect_mutation_location(mutations):
	if mutations:
		locations = list()
		for mutation in mutations:
			alteration = mutation.split(":")[1].split(".")[1]
			mutation_location = re.match("([0-9]+)([a-z]+)", alteration.split(">")[0], re.I).groups()[0]
			if int(mutation_location) not in locations:
				locations.append(int(mutation_location))
		if locations:
			return locations
		else:
			return None
	else:
		return None


def check_genome_for_mutation(genomic_range, direction, mutations, window_type, window):
	yes_mutation = False

	if direction == "left":
		end = int(genomic_range.split("-")[0])
		start = int(genomic_range.split("-")[1])

	elif direction == "right":
		end = int(genomic_range.split("-")[1])
		start = int(genomic_range.split("-")[0])

	if window_type == "gRNA":
		if mutations:
			for loc in mutations:
				if start <= loc <= end:
					yes_mutation = True

	elif window_type == "activity":
		if direction == "right":
			act_start = start + window[0]
			act_end = start + window[1]
			activity_sites = list(range(act_start, act_end))
		elif direction == "left":
			act_start = start - window[0]
			act_end = start - window[1]
			activity_sites = list(range(act_end, act_start))
		if mutations:
			for loc in mutations:
				if loc in activity_sites:
					yes_mutation = True

	elif window_type == "PAM":
		if direction == "right":
			pam_start = start + window[0]
			pam_end = start + window[1]
			pam_sites = list(range(pam_start, pam_end))
		elif direction == "left":
			pam_start = start - window[0]
			pam_end = start - window[1]
			pam_end = list(range(pam_end, pam_start))
		if mutations:
			for loc in mutations:
				if loc in pam_sites:
					yes_mutation = True

	return yes_mutation


def find_editable_nucleotide(crispr_df, searched_nucleotide, activity_window, pam_window,
							 ensembl_object, mutations):
	"""
	Finding editable nucleotides and their genomic coordinates
	:param crispr_df: A data frame having sequence, location and direction information of
	the CRISPRs from extract_crisprs().
	:param searched_nucleotide: The interested nucleotide which will be changed with BE
	:param activity_window: The location of the activity window on the protospacer sequence.
	:param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
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
	pam_window = [pam_window[0] - 1, pam_window[1]]

	print("Edit Data Frame is filling...")
	edit_df = pandas.DataFrame(columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence", "Location",
										"Edit_Location", "Direction", "Strand", "Gene_ID", "Transcript_ID", "Exon_ID",
										"guide_in_CDS", "gRNA_flanking_sequences", "Edit_in_Exon", "Edit_in_CDS", "GC%",
										"# Edits/guide", "Poly_T", "mutation_on_guide", "guide_change_mutation",
										"mutation_on_window", "mutation_on_PAM"])

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
										row["gRNA_flanking_sequences"], edit_in_exon, edit_in_cds]],
									  columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
											   "Location", "Edit_Location", "Direction", "Strand",
											   "Gene_ID", "Transcript_ID", "Exon_ID", "guide_in_CDS",
											   "gRNA_flanking_sequences", "Edit_in_Exon", "Edit_in_CDS"])
				edit_df = pandas.concat([edit_df, df])

		else:
			# If not --> no edit
			df = pandas.DataFrame([[row["Hugo_Symbol"], row["CRISPR_PAM_Sequence"],
									row["gRNA_Target_Sequence"], row["Location"], "No edit",
									row["Direction"], ensembl_object.strand, ensembl_object.gene_id,
									row["Transcript_ID"], row["Exon_ID"], row["guide_in_CDS"],
									row["gRNA_flanking_sequences"], False, False]],
								  columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "gRNA_Target_Sequence",
										   "Location", "Edit_Location", "Direction", "Strand", "Gene_ID",
										   "Transcript_ID", "Exon_ID", "guide_in_CDS", "gRNA_flanking_sequences",
										   "Edit_in_Exon", "Edit_in_CDS"])
			edit_df = pandas.concat([edit_df, df])

	edit_df["# Edits/guide"] = 0
	for guide, g_df in edit_df.groupby(["gRNA_Target_Sequence"]):
		unique_edits_per_guide = len(set(list(g_df["Edit_Location"])))
		edit_df.loc[edit_df.gRNA_Target_Sequence == guide, "# Edits/guide"] = unique_edits_per_guide

	edit_df["Poly_T"] = edit_df.apply(
		lambda x: True if re.search("T{4,}", x.CRISPR_PAM_Sequence) is not None else False, axis=1)

	edit_df["GC%"] = edit_df.apply(
		lambda x: (x.CRISPR_PAM_Sequence.count("C") + x.CRISPR_PAM_Sequence.count("G")) * 100.0 / len(
			x.CRISPR_PAM_Sequence), axis=1)

	mutation_locations = collect_mutation_location(mutations=mutations)

	edit_df["mutation_on_guide"] = edit_df.apply(
		lambda x: check_genome_for_mutation(genomic_range=x.Location.split(":")[1], direction=x.Direction,
											mutations=mutation_locations, window_type="gRNA", window=None), axis=1)
	edit_df["guide_change_mutation"] = edit_df.apply(
		lambda x: True if mutation_locations is not None and int(x.Edit_Location) in mutation_locations else False,
		axis=1)

	edit_df["mutation_on_window"] = edit_df.apply(
		lambda x: check_genome_for_mutation(genomic_range=x.Location.split(":")[1], direction=x.Direction,
											mutations=mutation_locations, window_type="activity",
											window=activity_window), axis=1)

	edit_df["mutation_on_PAM"] = edit_df.apply(
		lambda x: check_genome_for_mutation(genomic_range=x.Location.split(":")[1], direction=x.Direction,
											mutations=mutation_locations, window_type="PAM", window=pam_window), axis=1)

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
	loc_edit_df = None
	if transcript_id is not None:
		loc_edit_df = edit_df[edit_df.Transcript_ID == transcript_id]
	else:
		print("Ensembl object and info dict")
		print(ensembl_object.info_dict)
		for transcript, transcript_dict in ensembl_object.info_dict.items():
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

					if True not in grna_edit_df.mutation_on_window.unique():

						hgvs = "%s:g.%s%s>%s" \
							   % (str(chromosome), str(edit_loc), rev_edited_nucleotide, rev_new_nucleotide)

					elif len(list(grna_edit_df.mutation_on_window.unique())) == 1 and \
							list(grna_edit_df.mutation_on_window.unique())[0] is None:

						hgvs = "%s:g.%s%s>%s" \
							   % (str(chromosome), str(edit_loc), rev_edited_nucleotide, rev_new_nucleotide)

					else:
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]

						mutations_on_window = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								mutations_on_window.append(mutation)

						first_mut = max(mutations_on_window)
						last_mut = min(mutations_on_window)
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
						 "gRNA_flanking_sequences": grna_edit_df["gRNA_flanking_sequences"].values[0],
						 "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
						 "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[0],
						 "guide_change_mutation": grna_edit_df["guide_change_mutation"].values[0],
						 "mutation_on_window": grna_edit_df["mutation_on_window"].values[0],
						 "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
						 "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
						 "Poly_T": grna_edit_df["Poly_T"].values[0],
						 "GC%": grna_edit_df["GC%"].values[0],
						 "HGVS": hgvs}

					row_dicts.append(d)

				if total_edit > 1:

					if True not in grna_df.mutation_on_window.unique():
						# For multiple edits
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]
						position = str(start) + "_" + str(end)

						activity_sites = grna[activity_window[0]: activity_window[1]]
						activity_sites = "".join([nucleotide_dict[n] for n in activity_sites[::-1]])
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)
						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					elif len(list(grna_edit_df.mutation_on_window.unique())) == 1 and \
							list(grna_edit_df.mutation_on_window.unique())[0] is None:
						# For multiple edits
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]
						position = str(start) + "_" + str(end)

						activity_sites = grna[activity_window[0]: activity_window[1]]
						activity_sites = "".join([nucleotide_dict[n] for n in activity_sites[::-1]])
						edited_activity_sites = activity_sites.replace(rev_edited_nucleotide, rev_new_nucleotide)
						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					else:
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) - \
								activity_window[1] + 1
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[1]) - \
							  activity_window[0]

						mutations_on_window = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								mutations_on_window.append(mutation)

						first_mut = max(mutations_on_window)
						last_mut = min(mutations_on_window)
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

					d = {"Hugo_Symbol": list(grna_edit_df["Hugo_Symbol"].values)[0], "Edit_Type": "multiple",
						 "CRISPR_PAM_Sequence": grna_edit_df["CRISPR_PAM_Sequence"].values[0],
						 "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
						 "gRNA_Target_Sequence": grna,
						 "gRNA_Target_Location": grna_edit_df["Location"].values[0].split(":")[0] + ":" +
												 str(int(grna_edit_df["Location"].values[0].split(":")[1].split("-")[
															 0]) - 3) + "-" +
												 grna_edit_df["Location"].values[0].split(":")[1].split("-")[1],
						 "Total_Edit": total_edit,
						 "Edit_Location": position.split("_")[0] + "-" + position.split("_")[1],
						 "Direction": direction,
						 "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
						 "Exon_ID": grna_edit_df["Exon_ID"].values[0],
						 "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
						 "gRNA_flanking_sequences": grna_edit_df["gRNA_flanking_sequences"].values[0],
						 "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
						 "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[0],
						 "guide_change_mutation": grna_edit_df["guide_change_mutation"].values[0],
						 "mutation_on_window": grna_edit_df["mutation_on_window"].values[0],
						 "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
						 "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
						 "Poly_T": grna_edit_df["Poly_T"].values[0],
						 "GC%": grna_edit_df["GC%"].values[0],
						 "HGVS": hgvs}
					row_dicts.append(d)

		elif direction == "right":

			for grna, grna_df in direction_df.groupby(["gRNA_Target_Sequence"]):

				total_edit = len(set(list(grna_df["Edit_Location"].values)))

				# For individual edits

				for edit_loc, grna_edit_df in grna_df.groupby(["Edit_Location"]):

					if True not in grna_edit_df.mutation_on_window.unique():

						hgvs = "%s:g.%s%s>%s" % (str(chromosome), str(edit_loc), edited_nucleotide, new_nucleotide)

					elif len(list(grna_edit_df.mutation_on_window.unique())) == 1 and \
							list(grna_edit_df.mutation_on_window.unique())[0] is None:

						hgvs = "%s:g.%s%s>%s" % (str(chromosome), str(edit_loc), edited_nucleotide, new_nucleotide)

					else:
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]

						mutations_on_window = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								mutations_on_window.append(mutation)

						first_mut = min(mutations_on_window)
						last_mut = max(mutations_on_window)
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
						edited_activity_sites = activity_sites.replace(edited_nucleotide, new_nucleotide)

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
						 "gRNA_flanking_sequences": grna_edit_df["gRNA_flanking_sequences"].values[0],
						 "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
						 "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[0],
						 "guide_change_mutation": grna_edit_df["guide_change_mutation"].values[0],
						 "mutation_on_window": grna_edit_df["mutation_on_window"].values[0],
						 "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
						 "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
						 "Poly_T": grna_edit_df["Poly_T"].values[0],
						 "GC%": grna_edit_df["GC%"].values[0],
						 "HGVS": hgvs}

					row_dicts.append(d)

				if total_edit > 1:

					if True not in grna_df.mutation_on_window.unique():
						# For multiple edits
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]
						position = str(start) + "_" + str(end)

						activity_sites = grna[activity_window[0]: activity_window[1]]
						edited_activity_sites = activity_sites.replace(edited_nucleotide, new_nucleotide)
						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					elif len(list(grna_edit_df.mutation_on_window.unique())) == 1 and \
							list(grna_edit_df.mutation_on_window.unique())[0] is None:

						# For multiple edits
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]
						position = str(start) + "_" + str(end)

						activity_sites = grna[activity_window[0]: activity_window[1]]
						edited_activity_sites = activity_sites.replace(edited_nucleotide, new_nucleotide)
						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					else:
						end = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
							  activity_window[1] - 1
						start = int(list(grna_df["Location"].values)[0].split(":")[1].split("-")[0]) + \
								activity_window[0]

						mutations_on_window = list()
						for mutation in mutation_locations:
							if start <= mutation <= end:
								mutations_on_window.append(mutation)

						first_mut = min(mutations_on_window)
						last_mut = max(mutations_on_window)
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
						edited_activity_sites = activity_sites.replace(edited_nucleotide, new_nucleotide)

						hgvs = "%s:g.%sdelins%s" % (str(chromosome), position, edited_activity_sites)

					d = {"Hugo_Symbol": grna_edit_df["Hugo_Symbol"].values[0], "Edit_Type": "multiple",
						 "CRISPR_PAM_Sequence": grna_edit_df["CRISPR_PAM_Sequence"].values[0],
						 "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
						 "gRNA_Target_Sequence": grna,
						 "gRNA_Target_Location": grna_edit_df["Location"].values[0].split(":")[0] + ":" +
												 grna_edit_df["Location"].values[0].split(":")[1].split("-")[0] + "-" + \
												 str(int(grna_edit_df["Location"].values[0].split(":")[1].split("-")[
															 1]) - 3),
						 "Total_Edit": total_edit,
						 "Edit_Location": position.split("_")[0] + "-" + position.split("_")[1],
						 "Direction": direction,
						 "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
						 "Exon_ID": grna_edit_df["Exon_ID"].values[0],
						 "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
						 "gRNA_flanking_sequences": grna_edit_df["gRNA_flanking_sequences"].values[0],
						 "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
						 "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
						 "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[0],
						 "guide_change_mutation": grna_edit_df["guide_change_mutation"].values[0],
						 "mutation_on_window": grna_edit_df["mutation_on_window"].values[0],
						 "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
						 "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
						 "Poly_T": grna_edit_df["Poly_T"].values[0],
						 "GC%": grna_edit_df["GC%"].values[0],
						 "HGVS": hgvs}
					row_dicts.append(d)

	hgvs_df = pandas.DataFrame(row_dicts)

	return hgvs_df


def retrieve_vep_info(hgvs_df, ensembl_object, transcript_id=None):
	"""
	Collect Ensembl VEP information for given edits
	:param hgvs_df: The HGVS notations of all possible variants
	:param ensembl_object: The Ensembl Object created with Ensembl().
	:param transcript_id: The interested Ensembl transcript id
	:return uniprot_results: The Uniprot IDs in which edit occurs (swissprot or trembl)
	"""

	chromosome, strand = ensembl_object.chromosome, ensembl_object.strand

	if transcript_id is None:
		transcript_id = list(ensembl_object.info_dict.keys())[0]

	vep_columns = ["Protein_ID", "VEP_input", "allele", "variant_classification", "most_severe_consequence",
				   "consequence_terms", "variant_biotype", "Regulatory_ID", "Motif_ID", "TFs_on_motif",
				   "cDNA_Change", "Edited_Codon", "New_Codon", "CDS_Position", "Protein_Position",
				   "Protein_Change", "Edited_AA", "Edited_AA_Prop", "New_AA", "New_AA_Prop", "is_Synonymous",
				   "is_Stop", "proline_addition", "swissprot",
				   "polyphen_score", "polyphen_prediction", "sift_score", "sift_prediction", "cadd_phred",
				   "cadd_raw", "lof", "impact", "blosum62", "is_clinical", "clinical_id",
				   "clinical_significance", "cosmic_id", "clinvar_id", "ancestral_populations"]

	for c in vep_columns:
		hgvs_df[c] = None

	print("VEP Data Frame is filling with VEP API.")
	# Decide the server
	server = "http://grch37.rest.ensembl.org" if ensembl_object.assembly == "hg19" \
		else "https://rest.ensembl.org"
	ext = "/vep/human/hgvs"
	headers = {"Content-Type": "application/json", "Accept": "application/json"}
	params = {"AncestralAllele": 1, "Blosum62": 1, "Conservation": 1,
			  "LoF": 1, "CADD": 1, "protein": 1, "variant_class": 1, "hgvs": 1, "uniprot": 1,
			  "transcript_id": transcript_id}

	hgvs_index = pandas.DataFrame(columns=["HGVS"], index=list(range(len(list(hgvs_df["HGVS"].unique())))))
	count = 0
	for hgvs in list(hgvs_df["HGVS"].unique()):
		hgvs_index.loc[count, "HGVS"] = hgvs
		count += 1

	t = count // 200
	r = count - (200 * t)
	hgvs_obj = dict()
	for i in range(t):
		x = 200 * i
		hgvs_list = list(hgvs_index.loc[x: x + 199]["HGVS"].values)
		hgvs_json = json.dumps(hgvs_list)
		vep_request = requests.post(server + ext, headers=headers, params=params,
									data='{ "hgvs_notations" : %s }' % hgvs_json)
		print(vep_request.ok)
		if not vep_request.ok:
			print("No response from VEP %d - %d" % (x, x + 199))

		else:
			whole_vep = vep_request.json()
			for hgvs in hgvs_list:
				obj = Variant(hgvs=hgvs, gene=ensembl_object.hugo_symbol,
							  transcript=transcript_id, strand=strand)
				obj.extract_vep_obj(vep_json=whole_vep)
				obj.extract_consequences()
				hgvs_obj[hgvs] = obj

	hgvs_list = list(hgvs_index.loc[x + 200: x + 200 + r]["HGVS"].values)
	hgvs_json = json.dumps(hgvs_list)
	vep_request = requests.post(server + ext, headers=headers, params=params,
								data='{ "hgvs_notations" : %s }' % hgvs_json)
	if not vep_request.ok:
		print("No response from VEP %d - %d" % (x + 200, x + 200 + r))

	else:
		whole_vep = vep_request.json()
		for hgvs in hgvs_list:
			obj = Variant(hgvs=hgvs, gene=ensembl_object.hugo_symbol,
						  transcript=transcript_id, strand=strand)
			obj.extract_vep_obj(vep_json=whole_vep)
			obj.extract_consequences()
			hgvs_obj[hgvs] = obj

	for hgvs, obj in hgvs_obj.items():
		ind = list(hgvs_df[hgvs_df.HGVS == hgvs].index)
		hgvs_df.loc[ind, "Protein_ID"] = obj.protein
		hgvs_df.loc[ind, "VEP_input"] = obj.hgvs
		hgvs_df.loc[ind, "allele"] = obj.allele
		hgvs_df.loc[ind, "variant_classification"] = obj.variant_class
		hgvs_df.loc[ind, "most_severe_consequence"] = obj.most_severe_consequence
		hgvs_df.loc[ind, "consequence_terms"] = obj.consequence_terms
		hgvs_df.loc[ind, "variant_biotype"] = obj.biotype
		hgvs_df.loc[ind, "Regulatory_ID"] = obj.regulatory
		hgvs_df.loc[ind, "Motif_ID"] = obj.motif
		hgvs_df.loc[ind, "TFs_on_motif"] = obj.motif_TFs
		hgvs_df.loc[ind, "cDNA_Change"] = obj.cdna_change
		hgvs_df.loc[ind, "Edited_Codon"] = obj.old_codon
		hgvs_df.loc[ind, "New_Codon"] = obj.new_codon
		hgvs_df.loc[ind, "CDS_Position"] = obj.cds_position
		hgvs_df.loc[ind, "Protein_Position"] = obj.protein_position
		hgvs_df.loc[ind, "Protein_Change"] = obj.protein_change
		hgvs_df.loc[ind, "Edited_AA"] = obj.old_aa
		hgvs_df.loc[ind, "Edited_AA_Prop"] = obj.old_aa_chem
		hgvs_df.loc[ind, "New_AA"] = obj.new_aa
		hgvs_df.loc[ind, "New_AA_Prop"] = obj.new_aa_chem
		hgvs_df.loc[ind, "is_Synonymous"] = obj.synonymous
		hgvs_df.loc[ind, "is_Stop"] = obj.stop
		hgvs_df.loc[ind, "proline_addition"] = obj.proline
		hgvs_df.loc[ind, "swissprot"] = obj.swissprot
		hgvs_df.loc[ind, "polyphen_score"] = obj.polyphen_score
		hgvs_df.loc[ind, "polyphen_prediction"] = obj.polyphen_prediction
		hgvs_df.loc[ind, "sift_score"] = obj.sift_score
		hgvs_df.loc[ind, "sift_prediction"] = obj.sift_prediction
		hgvs_df.loc[ind, "cadd_phred"] = obj.cadd_phred
		hgvs_df.loc[ind, "cadd_raw"] = obj.cadd_raw
		hgvs_df.loc[ind, "lof"] = obj.lof
		hgvs_df.loc[ind, "impact"] = obj.impact
		hgvs_df.loc[ind, "blosum62"] = obj.blosum62
		hgvs_df.loc[ind, "is_clinical"] = obj.clinical
		hgvs_df.loc[ind, "clinical_id"] = obj.clinical_id
		hgvs_df.loc[ind, "clinical_significance"] = obj.clinical_sig
		hgvs_df.loc[ind, "cosmic_id"] = obj.cosmic_id
		hgvs_df.loc[ind, "clinvar_id"] = obj.clinvar_id
		hgvs_df.loc[ind, "ancestral_populations"] = obj.ancestral_populations

	vep_df = hgvs_df.copy()
	vep_df = vep_df.drop_duplicates()

	"""
	vep_df["cosmic_freq"] = vep_df.apply(
		lambda x: cosmic_freq.loc[x.cosmic_id]["frequency"] if x.cosmic_id is not None and x.cosmic_id in cosmic_freq.index else None, axis=1)
	"""
	return vep_df


def annotate_edits(ensembl_object, vep_df):
	"""
	Adding Uniprot API Information on VEP DF
	:param ensembl_object: The object of the Ensembl from Ensembl API
	:param vep_df: The data frame filled with the information from VEP API
	Ensembl Protein ID to Uniprot IDs (SwissProt/Reviewed)
	:return: analysis_df: The data frame enriched with the information from Uniprot API
	"""

	uniprot_df = vep_df.copy()
	uniprot_df["Domain"] = None
	uniprot_df["curated_Domain"] = None
	uniprot_df["PTM"] = None
	print(vep_df["swissprot"].unique())
	uniprot = list(vep_df["swissprot"].unique())[0]
	ensembl_p = list(vep_df["Protein_ID"].unique())[0]
	seq_mapping = ensembl_object.extract_uniprot_info(ensembl_pid=ensembl_p, uniprot=uniprot)
	if seq_mapping:
		uniprot = uniprot.split(".")[0].split("-")[0]
		print("Uniprot ID: %s" % uniprot)
		smap = seq_mapping[uniprot]
		obj = Uniprot(uniprotid=uniprot)
		obj.extract_uniprot()
		for ind, row in uniprot_df.iterrows():
			ptm, domain, c_domain = None, None, None
			if row["Protein_Position"] is not None and pandas.isna(row["Protein_Position"]) is False:
				ptms, domains = list(), list()
				for position in str(row["Protein_Position"]).split(";"):
					if position is not None and position != "None" and type(position) != float:
						if int(position) in smap.keys():
							dom = obj.find_domain(smap[int(position)], row["Edited_AA"])
							phos = obj.find_ptm_site("phosphorylation", smap[int(position)], row["Edited_AA"])
							meth = obj.find_ptm_site("methylation", smap[int(position)], row["Edited_AA"])
							ubi = obj.find_ptm_site("ubiquitination", smap[int(position)], row["Edited_AA"])
							acet = obj.find_ptm_site("acetylation", smap[int(position)], row["Edited_AA"])
							if dom is not None: domains.append(dom + "-" + position)
							if phos is not None: ptms.append(phos + "-" + position)
							if meth is not None: ptms.append(meth + "-" + position)
							if ubi is not None: ptms.append(ubi + "-" + position)
							if acet is not None: ptms.append(acet + "-" + position)
				if ptms: ptm = ";".join([i for i in ptms])
				if domains:
					domain = ";".join([i for i in domains])
					c_domain = ";".join(["-".join(i.split("-")[:-1]) for i in domains])
			uniprot_df.loc[ind, "Domain"] = domain
			uniprot_df.loc[ind, "curated_Domain"] = c_domain
			uniprot_df.loc[ind, "PTM"] = ptm
	return uniprot_df


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
	for group, group_df in df.groupby(["swissprot", "Protein_Position"]):
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
						api_url = "offset=0&size=-1&accession=%s&reviewed=true&isoform=0" % uniprot

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
						api_url = "offset=0&size=-1&accession=%s&reviewed=true&isoform=0" % uniprot

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
						api_url = "offset=0&size=-1&accession=%s&reviewed=true&isoform=0" % uniprot

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
	elif "intron" in mutation_consequence.split(";"):
		return "intron"
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
	else:
		return None


def summarise_guides(last_df):
	summary_df = pandas.DataFrame(index=list(range(0, len(last_df.groupby(["CRISPR_PAM_Sequence"])))),
								  columns=["Hugo_Symbol", "CRISPR_PAM_Sequence", "CRISPR_PAM_Location",
										   "gRNA_Target_Sequence",
										   "gRNA_Target_Location", "Edit_Location", "Direction", "Transcript_ID",
										   "Exon_ID", "Protein_ID",
										   "guide_in_CDS", "Edit_in_Exon", "Edit_in_CDS", "mutation_on_guide",
										   "guide_change_mutation",
										   "mutation_on_window", "mutation_on_PAM", "# Edits/guide", "Poly_T", "GC%",
										   "allele", "cDNA_Change",
										   "CDS_Position", "Protein_Position", "Protein_Change", "Edited_AA",
										   "Edited_AA_Prop",
										   "New_AA", "New_AA_Prop", "is_stop", "is_synonymous", "proline_addition",
										   "variant_classification",
										   "consequence_terms", "most_severe_consequence", "variant_biotype",
										   "Regulatory_ID",
										   "Motif_ID", "TFs_on_motif", "polyphen_prediction", "sift_prediction",
										   "impact",
										   "is_clinical", "clinical_id", "clinical_significance", "cosmic_id",
										   "clinvar_id",
										   "ancestral_populations", "swissprot", "Domain", "curated_Domain", "PTM",
										   "is_disruptive_interface_EXP", "disrupted_PDB_int_partners",
										   "disrupted_PDB_int_genes", "is_disruptive_interface_MOD",
										   "disrupted_I3D_int_partners",
										   "disrupted_I3D_int_genes", "is_disruptive_interface_PRED",
										   "disrupted_Eclair_int_partners",
										   "disrupted_Eclair_int_genes"])
	# cosmic_freq

	i = 0
	for guide, guide_df in last_df.groupby(["CRISPR_PAM_Sequence"]):
		summary_df.loc[i, "Hugo_Symbol"] = ";".join([x for x in list(guide_df.Hugo_Symbol.unique()) if x is not None])
		summary_df.loc[i, "CRISPR_PAM_Sequence"] = ";".join(
			[x for x in list(guide_df.CRISPR_PAM_Sequence.unique()) if x is not None])
		summary_df.loc[i, "CRISPR_PAM_Location"] = ";".join(
			[x for x in list(guide_df.CRISPR_PAM_Location.unique()) if x is not None])
		summary_df.loc[i, "gRNA_Target_Sequence"] = ";".join(
			[x for x in list(guide_df.gRNA_Target_Sequence.unique()) if x is not None])
		summary_df.loc[i, "gRNA_Target_Location"] = ";".join(
			[x for x in list(guide_df.gRNA_Target_Location.unique()) if x is not None])
		summary_df.loc[i, "Edit_Location"] = ";".join(
			[str(x) for x in list(guide_df.Edit_Location.unique()) if x is not None])
		summary_df.loc[i, "Direction"] = ";".join([x for x in list(guide_df.Direction.unique()) if x is not None])
		summary_df.loc[i, "Transcript_ID"] = ";".join(
			[x for x in list(guide_df.Transcript_ID.unique()) if x is not None])
		summary_df.loc[i, "Exon_ID"] = ";".join([x for x in list(guide_df.Exon_ID.unique()) if x is not None])
		summary_df.loc[i, "Protein_ID"] = ";".join([x for x in list(guide_df.Protein_ID.unique()) if x is not None])
		if guide_df[~pandas.isna(guide_df.Regulatory_ID)].Regulatory_ID.unique() is not None and type(
				guide_df.Regulatory_ID) != float and list(
				guide_df[~pandas.isna(guide_df.Regulatory_ID)].Regulatory_ID.unique()):
			summary_df.loc[i, "Regulatory_ID"] = ";".join(
				[x for x in list(guide_df.Regulatory_ID.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "Regulatory_ID"] = None
		if guide_df[~pandas.isna(guide_df.Motif_ID)].Motif_ID.unique() is not None and type(
				guide_df.Motif_ID) != float and list(guide_df[~pandas.isna(guide_df.Motif_ID)].Motif_ID.unique()):
			summary_df.loc[i, "Motif_ID"] = ";".join(
				[x for x in list(guide_df.Motif_ID.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "Motif_ID"] = None
		if guide_df[~pandas.isna(guide_df.TFs_on_motif)].TFs_on_motif.unique() is not None and type(
				guide_df.TFs_on_motif) != float and list(
				guide_df[~pandas.isna(guide_df.TFs_on_motif)].TFs_on_motif.unique()):
			summary_df.loc[i, "TFs_on_motif"] = ";".join(
				[x for x in list(guide_df.TFs_on_motif.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "TFs_on_motif"] = None
		summary_df.loc[i, "guide_in_CDS"] = True if True in guide_df.guide_in_CDS.unique() else False
		summary_df.loc[i, "Edit_in_Exon"] = True if True in guide_df.Edit_in_Exon.unique() else False
		summary_df.loc[i, "Edit_in_CDS"] = True if True in guide_df.Edit_in_CDS.unique() else False
		summary_df.loc[i, "mutation_on_guide"] = True if True in guide_df.mutation_on_guide.unique() else False
		summary_df.loc[i, "guide_change_mutation"] = True if True in guide_df.guide_change_mutation.unique() else False
		summary_df.loc[i, "mutation_on_window"] = True if True in guide_df.mutation_on_window.unique() else False
		summary_df.loc[i, "mutation_on_PAM"] = True if True in guide_df.mutation_on_PAM.unique() else False
		summary_df.loc[i, "# Edits/guide"] = guide_df["# Edits/guide"].unique()[0]
		summary_df.loc[i, "Poly_T"] = True if True in guide_df.Poly_T.unique() else False
		summary_df.loc[i, "GC%"] = True if True in guide_df["GC%"].unique() else False
		summary_df.loc[i, "allele"] = ";".join(
			[x for x in list(guide_df.allele.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "cDNA_Change"] = ";".join(
			[x for x in list(guide_df.cDNA_Change.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "CDS_Position"] = ";".join(
			[x for x in list(guide_df.CDS_Position.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "Protein_Position"] = ";".join(
			[str(x) for x in list(guide_df.Protein_Position.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "Protein_Change"] = ";".join(
			[x for x in list(guide_df.Protein_Change.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "Edited_AA"] = ";".join(
			[x for x in list(guide_df.Edited_AA.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "New_AA"] = ";".join(
			[x for x in list(guide_df.New_AA.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "Edited_AA_Prop"] = ";".join(
			[x for x in list(guide_df.Edited_AA_Prop.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "New_AA_Prop"] = ";".join(
			[x for x in list(guide_df.New_AA_Prop.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "swissprot"] = ";".join(
			[x for x in list(guide_df.swissprot.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "variant_classification"] = ";".join(
			[x for x in list(guide_df.variant_classification.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "variant_biotype"] = ";".join(
			[x for x in list(guide_df.variant_biotype.unique()) if x is not None and type(x) != float])
		summary_df.loc[i, "consequence_terms"] = ";".join([
			select_severe_effects(x)
			for x in [select_severe_effects(i)
					  for i in [rename_mutational_consequences(c)
								for c in [x for x in list(guide_df.consequence_terms.unique()) if
										  x is not None and pandas.isna(x) == False]]]])
		summary_df.loc[i, "most_severe_consequence"] = \
			select_severe_effects(";".join([rename_mutational_consequences(x)
											for x in guide_df.most_severe_consequence.unique()
											if x is not None and pandas.isna(x) == False]))
		summary_df.loc[i, "is_stop"] = True if "stop codon" in guide_df.most_severe_consequence.unique() else False
		summary_df.loc[
			i, "is_synonymous"] = True if "synonymous" in guide_df.most_severe_consequence.unique() else False
		summary_df.loc[i, "proline_addition"] = True if True in guide_df.proline_addition.unique() else False
		if guide_df[~pandas.isna(guide_df.polyphen_prediction)].polyphen_prediction.unique() is not None and type(
				guide_df.polyphen_prediction) != float and \
				list(guide_df[~pandas.isna(guide_df.polyphen_prediction)].polyphen_prediction.unique()):
			summary_df.loc[i, "polyphen_prediction"] = ";".join(
				[str(x) for x in list(guide_df.polyphen_prediction.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "polyphen_prediction"] = None
		summary_df.loc[i, "sift_prediction"] = ";".join([x for x in list(guide_df.sift_prediction.unique()) if
														 x is not None and pandas.isna(x) == False and type(
															 x) != float])
		summary_df.loc[i, "impact"] = ";".join([x for x in list(guide_df.impact.unique()) if
												x is not None and pandas.isna(x) == False and type(x) != float])
		summary_df.loc[i, "is_clinical"] = True if True in guide_df.is_clinical.unique() else False
		if guide_df[~pandas.isna(guide_df.clinical_id)].clinical_id.unique() is not None and type(
				guide_df.clinical_id) != float and \
				list(guide_df[~pandas.isna(guide_df.clinical_id)].clinical_id.unique()):
			summary_df.loc[i, "clinical_id"] = ";".join([x for x in list(guide_df.clinical_id.unique()) if
														 x is not None and pandas.isna(x) == False and type(
															 x) != float])
		else:
			summary_df.loc[i, "clinical_id"] = None
		if guide_df[~pandas.isna(guide_df.clinical_significance)].clinical_significance.unique() is not None and type(
				guide_df.clinical_significance) != float and \
				list(guide_df[~pandas.isna(guide_df.clinical_significance)].clinical_significance.unique()):
			summary_df.loc[i, "clinical_significance"] = ";".join(
				[x for x in list(guide_df.clinical_significance.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "clinical_significance"] = None
		if guide_df[~pandas.isna(guide_df.cosmic_id)].cosmic_id.unique() is not None and type(
				guide_df.cosmic_id) != float and \
				list(guide_df[~pandas.isna(guide_df.cosmic_id)].cosmic_id.unique()):
			summary_df.loc[i, "cosmic_id"] = ";".join(
				[str(x) for x in list(guide_df.cosmic_id.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "cosmic_id"] = None
		if guide_df[~pandas.isna(guide_df.clinvar_id)].clinvar_id.unique() is not None and type(
				guide_df.clinvar_id) != float and \
				list(guide_df[~pandas.isna(guide_df.clinvar_id)].clinvar_id.unique()):
			summary_df.loc[i, "clinvar_id"] = ";".join(
				[str(x) for x in list(guide_df.clinvar_id.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "clinvar_id"] = None
		"""
		if guide_df[~pandas.isna(guide_df.cosmic_freq)].cosmic_freq.unique() is not None and type(guide_df.cosmic_freq) != float and \
				list(guide_df[~pandas.isna(guide_df.cosmic_freq)].cosmic_freq.unique()):
			summary_df.loc[i, "cosmic_freq"] = ";".join([str(x) for x in list(guide_df.cosmic_freq.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "cosmic_freq"] = None
		"""
		if guide_df[~pandas.isna(guide_df.ancestral_populations)].ancestral_populations.unique() is not None and type(
				guide_df.ancestral_populations) != float and \
				list(guide_df[~pandas.isna(guide_df.ancestral_populations)].ancestral_populations.unique()):
			summary_df.loc[i, "ancestral_populations"] = ";".join(
				[x for x in list(guide_df.ancestral_populations.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "ancestral_populations"] = None
		if guide_df[~pandas.isna(guide_df.Domain)].Domain.unique() is not None and type(guide_df.Domain) != float and \
				list(guide_df[~pandas.isna(guide_df.Domain)].Domain.unique()):
			summary_df.loc[i, "Domain"] = ";".join(
				[x for x in list(guide_df.Domain.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "Domain"] = None
		if guide_df[~pandas.isna(guide_df.curated_Domain)].curated_Domain.unique() is not None and type(
				guide_df.curated_Domain) != float and \
				list(guide_df[~pandas.isna(guide_df.curated_Domain)].curated_Domain.unique()):
			summary_df.loc[i, "curated_Domain"] = ";".join(
				[x for x in list(guide_df.curated_Domain.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "curated_Domain"] = None
		if guide_df[~pandas.isna(guide_df.PTM)].PTM.unique() is not None and type(guide_df.PTM) != float and list(
				guide_df[~pandas.isna(guide_df.PTM)].PTM.unique()):
			summary_df.loc[i, "PTM"] = ";".join(
				[x for x in list(guide_df.PTM.unique()) if x is not None and type(x) != float])
		else:
			summary_df.loc[i, "PTM"] = None

		summary_df.loc[
			i, "is_disruptive_interface_EXP"] = True if True in guide_df.is_disruptive_interface_EXP.unique() else False
		summary_df.loc[i, "disrupted_PDB_int_partners"] = summarise_3di(
			list(guide_df.disrupted_PDB_int_partners.unique()))
		summary_df.loc[i, "disrupted_PDB_int_genes"] = summarise_3di(list(guide_df.disrupted_PDB_int_genes.unique()))
		summary_df.loc[
			i, "is_disruptive_interface_MOD"] = True if True in guide_df.is_disruptive_interface_MOD.unique() else False
		summary_df.loc[i, "disrupted_I3D_int_partners"] = summarise_3di(
			list(guide_df.disrupted_I3D_int_partners.unique()))
		summary_df.loc[i, "disrupted_I3D_int_genes"] = summarise_3di(list(guide_df.disrupted_I3D_int_genes.unique()))
		summary_df.loc[
			i, "is_disruptive_interface_PRED"] = True if True in guide_df.is_disruptive_interface_PRED.unique() else False
		summary_df.loc[i, "disrupted_Eclair_int_partners"] = summarise_3di(
			list(guide_df.disrupted_Eclair_int_partners.unique()))
		summary_df.loc[i, "disrupted_Eclair_int_genes"] = summarise_3di(
			list(guide_df.disrupted_Eclair_int_genes.unique()))
		i += 1

	return summary_df


def mm_combination_seq(mm, seq):
	"""
	Mismatched version of gRNAs
	:param mm: The max number of mismatches wanted
	:param seq: gRNA sequencd
	:return: list of potentially mutated sequences of the related gRNA
	"""
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


def grna_fasta(crisprs, output_name, mm):
	"""
	Create fats file
	:param crisprs: the list of gRNAs
	:param output_name: the result file
	:param mm: The max number of mismatches wanted
	:return:
	"""
	global ot_path

	if "%s_dict.p" % output_name not in os.listdir(ot_path + "/fasta_dict/"):
		mm = int(mm)
		f = open(ot_path + "/fasta/%s_fasta.fa" % output_name, "w")
		grna_dict = dict()
		count = 1
		for i in crisprs:
			seq = i[0]
			direction = i[1]
			f.write(">gRNA-%d-%s\n%s\n" % (count, direction, seq))
			if "gRNA-%d-%s" % (count, direction) not in grna_dict.keys():
				grna_dict["gRNA-%d-%s" % (count, direction)] = {"seq": seq}
			for k in range(1, mm + 1):
				mm_count = 1
				mm_seqs = mm_combination_seq(mm=k, seq=seq)
				for m_seq in mm_seqs:
					if len(m_seq) < 10:
						print(m_seq)
					f.write(">gRNA-%d-%s:mm%d:count%s\n%s\n" % (count, direction, k, mm_count, m_seq))
					grna_dict["gRNA-%d-%s" % (count, direction)]["mm%d:mm_count%d" % (k, mm_count)] = m_seq
					mm_count += 1
			count += 1

		f.close()
		pickle.dump(grna_dict, open(ot_path + "/fasta_dict/%s_dict.p" % output_name, "wb"))

	else:
		grna_dict = pickle.load(open(ot_path + "/fasta_dict/%s_dict.p" % output_name, "rb"))
	return grna_dict


def run_offtargets(genome, crisprs, output_name, mm, ot_vcf, threads_num):
	"""
	Run mrsfast from python with gRNAs created with BEstimate
	:param genome: fasta file of the genome
	:param crisprs: List of gRNA sequence with their directionality
	:param output_name: Name of the output file initial
	:param mm: Number of max mismatch allowed
	:param mm: Number of max mismatch allowed
	:return:
	"""
	global ot_path, mrsfast_path
	if "%s_alignment.sam" % output_name not in os.listdir(ot_path + "/sam_files/"):
		# Create fasta file
		print("Writing Fasta file..")
		if "%s_fasta.fa" % output_name not in os.listdir(ot_path + "/fasta/"):
			grna_fasta(crisprs=crisprs, output_name=output_name, mm=mm)
			print("Fasta file created.")

		print(ot_vcf)
		if ot_vcf is False:
			if "%s.fa.index" % genome not in os.listdir("%s/genome/" % ot_path):
				print("Indexing the genome..")
				os.system("mrsfast --index %s/genome/%s.fa" % (ot_path, genome))

			print("\ngRNA alignment on the genome..\n")
			print(genome)
			os.system(
				"mrsfast --search %s/genome/Homo_sapiens.GRCh38.dna_sm.chromosome.all.final.fa --seq %s/fasta/%s_fasta.fa --threads %s -e 0 -o %s/sam_files/%s_alignment.sam --disable-nohits"
				% (ot_path, ot_path, output_name, threads_num, ot_path, output_name))
			print("\ngRNA alignment on the genome was finished.\n")
	if "%s_alignment.sam" % output_name in os.listdir(ot_path + "/sam_files/"):
		return True
	else:
		print("No alignment - off target")
		return False


def read_sam(output_name, crisprs, mm):
	"""
	Reading the output of msrfast algorithm
	:param output_name: Name of the output file initial
	:param crisprs: List of gRNA sequence with their directionality
	:param mm: Number of max mismatch allowed
	:return:
	"""
	if "%s_ots.csv" % output_name not in os.listdir(ot_path + "/ot_files/"):
		ot_df = pandas.DataFrame(columns=["guide_flag", "guide_seq", "direction", "mm", "chr", "start", "mm_count"])
		grna_dict = grna_fasta(crisprs, output_name, mm)
		f = open(ot_path + "/sam_files/" + output_name + "_alignment.sam", "r")
		lines = f.readlines()
		ind_count = 0
		for line in lines:
			if line[:4] == "gRNA":
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

		f.close()
		mm_col = ["mm%d" % (i + 1) for i in range(int(mm))]
		df_col = ["CRISPR_PAM_Sequence", "exact"] + mm_col
		summary_ot_df = pandas.DataFrame(columns=df_col, index=list(ot_df.groupby(["guide_flag"]).groups.keys()))

		for g, g_df in ot_df.groupby(["guide_flag"]):
			exact = g_df.groupby(["mm"]).size().loc["0"]
			# Alignment sequence might be different strand but we want the sequence with PAM
			# Retrieve sequence of the actual gRNA
			grna_seq = grna_dict[g_df[g_df.mm == "0"]["guide_flag"].unique()[0]]["seq"]
			summary_ot_df.loc[g, "CRISPR_PAM_Sequence"] = grna_seq
			summary_ot_df.loc[g, "exact"] = exact
			for i in range(1, int(mm) + 1):
				if str(i) in g_df.mm.unique():
					summary_ot_df.loc[g, "mm%d" % i] = g_df.groupby(["mm"]).size().loc[str(i)]
				else:
					summary_ot_df.loc[g, "mm%d" % i] = 0

		summary_ot_df = summary_ot_df.reset_index()
		changed_df_col = ["guide_flag", "CRISPR_PAM_Sequence", "exact"] + mm_col
		summary_ot_df.columns = changed_df_col

		summary_ot_df.to_csv(ot_path + "/ot_files/%s_ots.csv" % output_name)
	else:
		summary_ot_df = pandas.read_csv(ot_path + "/ot_files/%s_ots.csv" % output_name)

	return summary_ot_df


def add_offtargets(genome, output_name, df, mm, threads_num):
	"""
	Adding off targets information on BEstimate results
	:param genome: fasta file of the genome
	:param output_name: Name of the output file initial
	:param df: BEstimate final data frame
	:param mm: Number of max mismatch allowed
	:param threads_num: Number of threads will be used for alignment
	:return:
	"""
	global ot_path
	crisprs = df[["CRISPR_PAM_Sequence", "Direction"]].values
	ot = run_offtargets(genome=genome, crisprs=crisprs, output_name=output_name, mm=mm, ot_vcf=False,
						threads_num=threads_num)
	if ot:
		ot_df = read_sam(output_name=output_name, crisprs=crisprs, mm=mm)
		if len(ot_df.index) > 0:
			df_final = pandas.merge(df, ot_df, how="outer", on=["CRISPR_PAM_Sequence"])
			return df_final
		else:
			return False
	else:
		return False


###########################################################################################
# Execution


def main():
	"""
	Run whole script with the input from terminal
	:return:
	"""

	global args

	print("""
--------------------------------------------------------------                                                                                         
		   B E s t i m a t e                                      

	       Wellcome Sanger Institute          

--------------------------------------------------------------
    """)
	if args["VEP"]:
		vep = True
	else:
		vep = False

	if args["OFF_TARGET"]:
		ot_analysis = True
		mm = int(args["MISMATCH"])
	else:
		ot_analysis = False

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
Protospacer length: %s\nActivity window: %s\nNucleotide change: %s>%s\nVEP and Uniprot analysis: %s\nMutation on genome: %s
Off target analysis: %s"""
		  % (args["GENE"], args["ASSEMBLY"], args["TRANSCRIPT"], args["PAMSEQ"], args["PAMWINDOW"], args["PROTOLEN"],
			 args["ACTWINDOW"], args["EDIT"], args["EDIT_TO"], vep, ", ".join("" if mutations is None else mutations),
			 ot_analysis))

	print("""\n
-------------------------------------------------------------- 
		Ensembl Gene Information
-------------------------------------------------------------- 
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
--------------------------------------------------------------
		gRNAs - Targetable Sites
--------------------------------------------------------------
    \n""")
	path = ""
	if args["OUTPUT_PATH"][-1] == "/":
		path = args["OUTPUT_PATH"]
	else:
		path = args["OUTPUT_PATH"] + "/"

	final_df = None

	if args["OUTPUT_FILE"] + "_crispr_df.csv" not in os.listdir(path):
		crispr_df = extract_grna_sites(hugo_symbol=args["GENE"], searched_nucleotide=args["EDIT"],
									   pam_window=[int(args["PAMWINDOW"].split("-")[0]),
												   int(args["PAMWINDOW"].split("-")[1])],
									   activity_window=[int(args["ACTWINDOW"].split("-")[0]),
														int(args["ACTWINDOW"].split("-")[1])],
									   pam_sequence=args["PAMSEQ"], protospacer_length=args["PROTOLEN"],
									   flan=args["FLAN"], flan_3=args["FLAN_3"], flan_5=args["FLAN_5"],
									   ensembl_object=ensembl_obj)

		if len(crispr_df.index) != 0: print("CRISPR Data Frame was created!")
		crispr_df.to_csv(path + args["OUTPUT_FILE"] + "_crispr_df.csv", index=False)

		print("CRISPR Data Frame was written in %s as %s\n" % (path, args["OUTPUT_FILE"] + "_crispr_df.csv"))

	else:
		print("CRISPR Data Frame was readed from %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_crispr_df.csv"))
		crispr_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_crispr_df.csv")
	print("""\n
--------------------------------------------------------------
		gRNAs - Editable Sites
--------------------------------------------------------------
    \n""")
	if args["OUTPUT_FILE"] + "_edit_df.csv" not in os.listdir(path):
		edit_df = find_editable_nucleotide(crispr_df=crispr_df, searched_nucleotide=args["EDIT"],
										   activity_window=[int(args["ACTWINDOW"].split("-")[0]),
															int(args["ACTWINDOW"].split("-")[1])],
										   pam_window=[int(args["PAMWINDOW"].split("-")[0]),
													   int(args["PAMWINDOW"].split("-")[1])],
										   ensembl_object=ensembl_obj, mutations=mutations)

		if len(edit_df.index) != 0: print("Edit Data Frame was created!")

		edit_df.to_csv(path + args["OUTPUT_FILE"] + "_edit_df.csv", index=False)

		print("Edit Data Frame was written in %s as %s" % (path, args["OUTPUT_FILE"] + "_edit_df.csv\n"))

	else:
		print("Edit Data Frame was readed from %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_edit_df.csv"))
		edit_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_edit_df.csv")

	if args["VEP"]:
		print("""\n
--------------------------------------------------------------
		Annotation - VEP Annotation
--------------------------------------------------------------
        \n""")
		if args["OUTPUT_FILE"] + "_vep_df.csv" not in os.listdir(path):
			if args["OUTPUT_FILE"] + "_hgvs_df.csv" not in os.listdir(path):
				hgvs_df = extract_hgvs(edit_df=edit_df, ensembl_object=ensembl_obj,
									   transcript_id=args["TRANSCRIPT"],
									   edited_nucleotide=args["EDIT"], new_nucleotide=args["EDIT_TO"],
									   activity_window=[int(args["ACTWINDOW"].split("-")[0]),
														int(args["ACTWINDOW"].split("-")[1])],
									   mutations=mutations)

				if hgvs_df is not None and len(hgvs_df.index) != 0:
					hgvs_df.to_csv(path + args["OUTPUT_FILE"] + "_hgvs_df.csv")
					print("HGVS nomenclatures were collected.\n")
			else:
				hgvs_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_hgvs_df.csv")
				print("HGVS nomenclatures were collected.\n")

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
			print("VEP Data Frame was readed from %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_vep_df.csv"))
			whole_vep_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_vep_df.csv")
		print(whole_vep_df)
		print("""\n
--------------------------------------------------------------
		Annotation - Protein Annotation
--------------------------------------------------------------
		\n""")
		if args["OUTPUT_FILE"] + "_protein_df.csv" not in os.listdir(path):
			print("Adding Uniprot ID, corresponding Domain and PTM information..")
			if len(whole_vep_df.index) != 0:
				uniprot_df = annotate_edits(ensembl_object=ensembl_obj, vep_df=whole_vep_df)
				if uniprot_df is not None and len(uniprot_df.index) != 0:
					print("Adding affected interface and interacting partners..")
					protein_df = annotate_interface(annotated_edit_df=uniprot_df)

					if protein_df is not None and len(protein_df.index) != 0:
						print("Protein Data Frame was created!")
						protein_df.to_csv(path + args["OUTPUT_FILE"] + "_protein_df.csv", index=False)
						print("Protein Data Frame was written in %s as %s\n" % (
						path, args["OUTPUT_FILE"] + "_protein_df.csv\n"))
					else:
						print("Protein Data Frame cannot be created because it is empty.")
				else:
					print("Protein Data Frame cannot be created because it is empty.")
			else:
				print("Protein Data Frame cannot be created because it is empty.")
		else:
			print("Protein Data Frame was readed from %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_protein_df.csv"))
			protein_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_protein_df.csv")

		if len(protein_df.index) > 0:
			if args["OUTPUT_FILE"] + "_summary_df.csv" not in os.listdir(path):
				print("Summarising information..")
				summary_df = summarise_guides(last_df=protein_df)

				if summary_df is not None and len(summary_df.index) != 0:
					print("Summary Data Frame was created!")
					summary_df.to_csv(path + args["OUTPUT_FILE"] + "_summary_df.csv", index=False)
					final_df = summary_df.copy()
					print("Summary Data Frame was written in %s as %s\n\n" % (
					path, args["OUTPUT_FILE"] + "_summary_df.csv"))
				else:
					print("Summary Data Frame cannot be created because it is empty.")

			else:
				print(
					"Summary Data Frame was readed from %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_summary_df.csv"))
				summary_df = pandas.read_csv(path + args["OUTPUT_FILE"] + "_summary_df.csv")
				final_df = summary_df.copy()
		else:
			print("Protein Data Frame cannot be created because it is empty.")

	else:
		final_df = edit_df.copy()

	if args["OFF_TARGET"]:
		print("""\n
--------------------------------------------------------------
		Annotation - Off Target Annotation
--------------------------------------------------------------
				\n""")

		bestimate_ot_df = add_offtargets(genome=args["GENOME"], output_name=args["OUTPUT_FILE"],
										 df=final_df, mm=mm, threads_num=4)
		if bestimate_ot_df is not False:
			bestimate_ot_df.to_csv(path + args["OUTPUT_FILE"] + "_ot_annotated_summary_df.csv", index=False)
		else:
			print("Off target information cannot be added.")
	return True



if __name__ == '__main__':

	# -----------------------------------------------------------------------------------------#
	# Retrieve input

	args = take_input()
	# Output Path
	path = ""
	if args["OUTPUT_PATH"][-1] == "/":
		path = args["OUTPUT_PATH"]
	else:
		path = args["OUTPUT_PATH"] + "/"

	ot_path = os.getcwd() + "/../offtargets"
	mrsfast_path = os.getcwd() + "/../bin/mrsfast"

	# -----------------------------------------------------------------------------------------#
	# Data w/out API opportunity

	yulab = pandas.read_table(os.getcwd() + "/../data/H_sapiens_interfaces.txt")

	# -----------------------------------------------------------------------------------------#

	###########################################################################################
	# Execution

	main()

	print("""\n
	--------------------------------------------------------------
		The BEstimate analysis successfully finished!
	--------------------------------------------------------------
	\n""")