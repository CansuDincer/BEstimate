###########################################################################################
##### Find editable nucleotides of a gene from a given PAM sequence with Base Editors #####
###########################################################################################

# Import necessary packages
import os, pandas, re, argparse, requests

###########################################################################################
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

    # FILTERING OPTION

    parser.add_argument("-transcript", dest="TRANSCRIPT", default = None,
                        help="If the user wants to filter the results according to given transcript ensembl id")

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

    # PROTEIN LEVEL ANALYSIS

    parser.add_argument("-P", dest="PROTEIN", action = "store_true",
                        help="The boolean option if user wants to analyse the edits.")

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


###########################################################################################
# Objects from APIs

class Uniprot:

    def __init__(self, uniprotid):
        self.uniprotid, self.reviewed = uniprotid, None
        self.sequence = None
        self.domains = dict()
        self.phosphorylation_sites = dict()
        self.server = "https://www.ebi.ac.uk/proteins/api/"

    def extract_uniprot_info(self):

        uniprot_api = "proteins?offset=0&size=-1&accession=%s" % self.uniprotid
        api_request = requests.get(self.server + uniprot_api,
                                   headers={"Accept": "application/json"})

        # Check the response of the server for the request
        if api_request.status_code != 200:
            return "No response from UniProt!"

        else:
            for i in api_request.json():
                if len(i["accession"].split("-")) == 1:
                    self.reviewed = False if i["info"]["type"] == "TrEMBL" else True
                    self.sequence = i["sequence"]["sequence"]
                    ptms = {}
                    if "features" in i.keys() and i["features"] != []:
                        for ftr in i["features"]:
                            if ftr["type"] == "MOD_RES" and ftr["category"] == "PTM":
                                if "description" in ftr.keys():
                                    pos, ptm = ftr["begin"], ftr["description"]
                                    ptm = ptm.split(";")[0]
                                    phos = ptm if re.search(r'Phospho', ptm) else None
                                    if phos is not None: self.phosphorylation_sites[pos] = phos

                            if ftr["category"] == "DOMAINS_AND_SITES":
                                if "description" in ftr.keys():
                                    domain = ftr["description"]
                                    domain_range = list(int(ftr["begin"]), int(ftr["end"]))
                                    self.domains[domain] = domain_range

            if self.phosphorylation_sites == dict(): self.phosphorylation_sites = None
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

    def find_psite(self, protein_edit_location, old_aa):

        edit_psite = None
        if self.phosphorylation_sites != {} and self.phosphorylation_sites is not None:
            for phos, phos_pos in self.phosphorylation_sites.items():
                if phos_pos == protein_edit_location:
                    if old_aa == self.sequence[protein_edit_location - 1]:
                        edit_psite = phos
        return edit_psite


class Ensembl:

    def __init__(self, hugo_symbol, assembly):
        self.hugo_symbol = hugo_symbol
        self.assembly = assembly
        self.server = "http://grch37.rest.ensembl.org" \
            if self.assembly == "hg19" else "https://rest.ensembl.org"
        self.gene_id = ''
        self.info_dict = dict()
        self.sequence, self.flan_sequence = None, None
        self.sequence_analysis, self.flan_sequence_analysis = None, None
        self.chromosome, self.strand = None, None
        self.gene_range, self.flan_gene_range = list(), list()

    def extract_gene_id(self):

        hugo_ensembl = "/xrefs/symbol/homo_sapiens/%s?" % self.hugo_symbol

        print("Request is done to Ensembl REST API for Ensembl Gene information:")
        gene_request = requests.get(self.server + hugo_ensembl,
                                    headers={"Content-Type": "application/json"})

        if gene_request.status_code != 200:
            print("No response from ensembl!")
        else:
            print("Successful response!")

        for x in gene_request.json():
            if x["id"][:4] == "ENSG": self.gene_id = x["id"]

        if self.gene_id != '':
            print("Corresponding Ensembl Gene ID: %s" % self.gene_id)
            return 1
        else:
            return 0

    def extract_sequence(self, gene_id):

        seq_ensembl = self.server + "/sequence/id/%s?" % gene_id
        seq_flan_ensembl = self.server + "/sequence/id/%s?expand_3prime=23;expand_5prime=23" % gene_id

        print("Request is done to Ensembl REST API for sequence information:")
        seq_request = requests.get(seq_ensembl,
                                   headers={"Content-Type": "text/x-fasta"})
        seq_flan_request = requests.get(seq_flan_ensembl,
                                        headers={"Content-Type": "text/x-fasta"})

        if seq_request.status_code != 200 and seq_flan_request.status_code != 200:
            print("No response from ensembl sequence!")
        else: print("Successful response!")

        # Sequence
        label_line = seq_request.text.split("\n")[0]
        flan_label_line = seq_flan_request.text.split("\n")[0]
        self.sequence = "".join(seq_request.text.split("\n")[1:])
        self.flan_sequence = "".join(seq_flan_request.text.split("\n")[1:])
        self.gene_range = [int(label_line.split(":")[-3]), int(label_line.split(":")[-2])]
        self.flan_gene_range = [int(flan_label_line.split(":")[-3]),
                                int(flan_label_line.split(":")[-2])]
        self.strand, self.chromosome = int(label_line.split(":")[-1].strip()),\
                                       int(label_line.split(":")[2].strip())

        # If strand is -1, the sequence has been reversed to be in 5'->3' direction
        # The genomic location should be reverse to match with the sequence too.
        if self.strand == -1: self.gene_range = [self.gene_range[1], self.gene_range[0]]
        if self.strand == -1: self.flan_gene_range = \
            [self.flan_gene_range[1], self.flan_gene_range[0]]

        # Preparation of the Ensembl sequence for analysis

        nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

        if self.strand == 1:
            self.sequence_analysis = self.sequence
            self.flan_sequence_analysis = self.flan_sequence
        elif self.strand == -1:
            pre_sequence = self.sequence[::-1]
            pre_flan_sequence = self.flan_sequence[::-1]
            self.sequence_analysis = "".join([nucleotide_dict[n] for n in pre_sequence])
            self.flan_sequence_analysis = "".join([nucleotide_dict[n] for n in pre_flan_sequence])

    def extract_info(self, chromosome, loc_start, loc_end):

        ensembl = "/overlap/region/human/%s:%s-%s?feature=transcript;feature=exon" % (
            chromosome, int(loc_start), int(loc_end))

        request = requests.get(self.server + ensembl, headers={"Content-Type": "application/json"})

        if request.status_code != 200:
            print("No response from ensembl!")
        else:
            for output in request.json():
                if output["feature_type"] == "transcript" and output["Parent"] == self.gene_id:
                    if output["transcript_id"] not in self.info_dict.keys():
                        self.info_dict[output["transcript_id"]] = \
                            [{"start": output["start"], "end": output["end"], "biotype": output["biotype"]}]
                    else:
                        old_val = self.info_dict[output["transcript_id"]]
                        if {"start": output["start"], "end": output["end"],
                            "biotype": output["biotype"]} not in old_val:
                            old_val.append(
                                {"start": output["start"], "end": output["end"], "biotype": output["biotype"]})
                            self.info_dict[output["transcript_id"]] = old_val

            for output in request.json():
                if output["feature_type"] == "exon" and self.info_dict != {} and \
                        output["Parent"] in self.info_dict.keys():
                    for d in self.info_dict[output["Parent"]]:
                        if "exon" not in d.keys():
                            d["exon"] = {output["exon_id"]: {"start": output["start"], "end": output["end"]}}
                        else:
                            if output["exon_id"] not in d["exon"].keys():
                                d["exon"][output["exon_id"]] = {"start": output["start"], "end": output["end"]}

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

    def extract_uniprot_info(self, ensembl_pid):

        protein_ensembl = "/xrefs/id/{0}?external_db=Uniprot/SPTREMBL".format(ensembl_pid)
        print("\nRequest is done to Ensembl REST API for Ensembl Protein information:\n")
        protein_request = requests.get(self.server + protein_ensembl,
                                       headers={"Content-Type": "application/json"})
        if protein_request.status_code != 200:
            print("No response from ensembl!")
            return 0
        else:
            for i in protein_request.json():
                uniprot = i["primary_id"]
                if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(i["xref_end"]) - int(i["xref_start"]):
                    # Otherwise, there is an inconsistency --> Not take it
                    seq_mapping[uniprot] = {i["ensembl_start"] + k: i["xref_start"] + k
                                            for k in range(
                            int(i["ensembl_end"]) - int(i["ensembl_start"]) + 1)}

            if seq_mapping == dict():
                return None
            else:
                return seq_mapping

###########################################################################################
# Functions


def find_pam_protospacer(sequence, pam_sequence, searched_nucleotide,
                         activity_window, pam_window, protospacer_length):
    """
    Finding all possible PAM and protospacer regions on the sequence of the gene.
    :param sequence: The sequence of the interested gene.
    :param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
    :param searched_nucleotide: The interested nucleotide which will be changed with BE
    :param activity_window: The location of the activity windiw on the protospacer sequence.
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

    crisprs, nuc_crisprs = [], []
    print("Pattern is searching through the sequence...")
    for nuc_index in range(0, len(sequence)):

        # One by one in the given sequence
        if nuc_index + pam_window[1] <= len(sequence):

            # Add started nucleotide index total length of targeted base editing site (PAM index)
            sub_sequence = sequence[nuc_index:nuc_index + pam_window[1]]

            # Search regex pattern inside the sub sequence
            for match_sequence in re.finditer(pattern, sub_sequence):

                # If there is match, then check if searched nucleotide inside the activity windiw
                if searched_nucleotide in list(match_sequence.group()[activity_window[0]:activity_window[1]]):
                    activity_sequence = match_sequence.group()[activity_window[0]:activity_window[1]]

                    # If searched nucleotide is also there, add the sequence inside crisprs!
                    nuc_crisprs.append({"index": [nuc_index, nuc_index + pam_window[1]],
                                        "crispr": match_sequence.group(),
                                        "activity_seq": activity_sequence})

                # For this one no need to check searched nucleotide.
                crisprs.append({"index": [nuc_index, nuc_index + pam_window[1]],
                                "crispr": match_sequence.group()})
    if crisprs != [] and nuc_crisprs != []: print("CRISPRs are found!")
    return crisprs, nuc_crisprs


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
        genomic_location = str(genomic_start) + ":" + str(genomic_end)

    elif crispr_direction == "left":
        genomic_start = new_range[1] - crispr_dict["index"][0]
        genomic_end = (genomic_start - len(crispr_seq)) + 1
        genomic_location = str(genomic_end) + ":" + str(genomic_start)

    return crispr_seq, genomic_location


def extract_grna_sites(pam_sequence, searched_nucleotide,
                       activity_window, pam_window, protospacer_length,
                       ensembl_object):
    """
    Extracting the gRNA targeted sites having editable nucleotide(s) on the interested genes
    :param pam_sequence: The sequence pattern of the PAM region (NGG/NG etc)
    :param searched_nucleotide: The interested nucleotide which will be changed with BE
    :param activity_window: The location of the activity windiw on the protospacer sequence.
    :param pam_window: The location of the PAM sequence when 1st index of the protospacer is 1.
    :param protospacer_length: The length of protospacer.
    :param ensembl_object: The Ensembl Object created with Ensembl().
    :return crispr_df: A data frame having sequence, location and direction information of the CRISPRs.
    """

    print("Inside extract_grna_sites function:")

    # If strand -1 : raw and edit sequences are different, otherwise they are same
    print("Sequence is preparing...")

    raw_sequence = ensembl_object.flan_sequence
    edit_sequence = ensembl_object.flan_sequence_analysis
    strand, seq_range = ensembl_object.strand, ensembl_object.flan_gene_range
    chromosome = ensembl_object.chromosome

    # Right CRISPRs 5'-->3' : reversed and base changed sequence

    # Editted should be used
    print("Protospacer and PAM regions are searching for right direction...")
    right_crisprs, right_nuc_crisprs = find_pam_protospacer(sequence=edit_sequence,
                                                            pam_sequence=pam_sequence,
                                                            searched_nucleotide=searched_nucleotide,
                                                            activity_window=activity_window,
                                                            pam_window=pam_window,
                                                            protospacer_length=protospacer_length)

    # Left CRISPRs: raw sequence will be used
    print("Protospacer and PAM regions are searching for left direction...")
    left_crisprs, left_nuc_crisprs = find_pam_protospacer(sequence=raw_sequence,
                                                          pam_sequence=pam_sequence,
                                                          searched_nucleotide=searched_nucleotide,
                                                          activity_window=activity_window,
                                                          pam_window=pam_window,
                                                          protospacer_length=protospacer_length)

    crisprs_dict = {"left": left_crisprs, "right": right_crisprs}
    crisprs_nuc_dict = {"left": left_nuc_crisprs, "right": right_nuc_crisprs}

    crispr_df = pandas.DataFrame(columns=["gRNA Target Sequence", "Location", "Direction",
                                          "Gene_ID"])
    print("CRISPR df is filling...")
    for direction, crispr in crisprs_dict.items():

        for cr in crispr:
            crispr_seq, genomic_location = add_genomic_location(sequence_range=seq_range, strand=strand,
                                                                crispr_dict=cr, crispr_direction=direction)
            """
            transcript_exon = ensembl_object.check_range_info(int(genomic_location.split(":")[0]),
                                                              int(genomic_location.split(":")[1]))
            if transcript_exon is not None:
                for transcript, exon_list in transcript_exon.items():
                    if exon_list:
                        for exon in exon_list:
                            df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                                    ensembl_object.gene_id, transcript, exon]],
                                                  columns=["gRNA Target Sequence", "Location",
                                                           "Direction", "Gene_ID", "Transcript_ID",
                                                           "Exon_ID"])
                            crispr_df = pandas.concat([crispr_df, df])
                    else:
                        df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                                ensembl_object.gene_id, transcript, None]],
                                              columns=["gRNA Target Sequence", "Location",
                                                       "Direction", "Gene_ID", "Transcript_ID",
                                                       "Exon_ID"])
                        crispr_df = pandas.concat([crispr_df, df])
            else:
                df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                        ensembl_object.gene_id, None, None]],
                                      columns=["gRNA Target Sequence", "Location",
                                               "Direction", "Gene_ID", "Transcript_ID",
                                               "Exon_ID"])
                crispr_df = pandas.concat([crispr_df, df])
            """
            df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                    ensembl_object.gene_id]],
                                  columns=["gRNA Target Sequence", "Location",
                                           "Direction", "Gene_ID"])
            crispr_df = pandas.concat([crispr_df, df])
    crispr_nuc_df = pandas.DataFrame(columns=["gRNA Target Sequence", "Location", "Direction",
                                              "Gene_ID"])

    for direction, crispr_nuc in crisprs_nuc_dict.items():

        for cr in crispr_nuc:
            crispr_seq, genomic_location = add_genomic_location(sequence_range=seq_range, strand=strand,
                                                                crispr_dict=cr, crispr_direction=direction)
            """
            transcript_exon = ensembl_object.check_range_info(int(genomic_location.split(":")[0]),
                                                              int(genomic_location.split(":")[1]))
            if transcript_exon is not None:
                for transcript, exon_list in transcript_exon.items():
                    if exon_list:
                        for exon in exon_list:
                            df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                                    ensembl_object.gene_id, transcript, exon]],
                                                  columns=["gRNA Target Sequence", "Location",
                                                           "Direction", "Gene_ID", "Transcript_ID",
                                                           "Exon_ID"])
                            crispr_nuc_df = pandas.concat([crispr_nuc_df, df])
                    else:
                        df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                                ensembl_object.gene_id, transcript, None]],
                                              columns=["gRNA Target Sequence", "Location",
                                                       "Direction", "Gene_ID", "Transcript_ID",
                                                       "Exon_ID"])
                        crispr_nuc_df = pandas.concat([crispr_nuc_df, df])
            else:
                df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                        ensembl_object.gene_id, None, None]],
                                      columns=["gRNA Target Sequence", "Location",
                                               "Direction", "Gene_ID", "Transcript_ID",
                                               "Exon_ID"])
                crispr_nuc_df = pandas.concat([crispr_nuc_df, df])
            """
            df = pandas.DataFrame([[crispr_seq, genomic_location, direction,
                                    ensembl_object.gene_id]],
                                  columns=["gRNA Target Sequence", "Location",
                                           "Direction", "Gene_ID"])

            crispr_nuc_df = pandas.concat([crispr_nuc_df, df])

    return crispr_df, crispr_nuc_df


def find_editable_nucleotide(crispr_df, searched_nucleotide, activity_window,
                             ensembl_object, transcript_id):
    """
    Finding editable nucleotides and their genomic coordinates
    :param crispr_df: A data frame having sequence, location and direction information of
    the CRISPRs from extract_crisprs().
    :param searched_nucleotide: The interested nucleotide which will be changed with BE
    :param activity_window: The location of the activity windiw on the protospacer sequence.
    :param ensembl_object: The Ensembl Object created with Ensembl().
    :param transcript_id: The user input for transcript id or defualt is all
    :return edit_df: A data frame having sequence, edit_location, location and direction
    information of the CRISPRs.
    """

    nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    actual_seq_range = ensembl_object.gene_range
    if actual_seq_range[0] > actual_seq_range[1]:
        actual_seq_range = [actual_seq_range[1], actual_seq_range[0]]

    actual_locations = list(range(actual_seq_range[0], actual_seq_range[1]))

    activity_window = [activity_window[0] - 1, activity_window[1]]

    print("EDIT df is filling...")
    edit_df = pandas.DataFrame(columns=["CRISPR Sequence", "gRNA Target Sequence", "Edit Location",
                                        "CRISPR Location", "Direction", "Gene ID (edit)",
                                        "Transcript ID (edit)"])

    for ind, row in crispr_df.iterrows():
        # Check inly with the sequence having PAM since it only has the searched nucleotide!
        searched_ind = [nuc_ind for nuc_ind in range(0, len(row["gRNA Target Sequence"]))
                        if nuc_ind in list(range(activity_window[0], activity_window[1])) and
                        row["gRNA Target Sequence"][nuc_ind] == searched_nucleotide]

        actual_inds = []
        if row["Direction"] == "left":
            for nuc_ind in searched_ind:
                if int(row["Location"].split(":")[1]) - nuc_ind in actual_locations:
                    actual_inds.append(int(row["Location"].split(":")[1]) - nuc_ind)

        elif row["Direction"] == "right":
            for nuc_ind in searched_ind:
                if int(row["Location"].split(":")[0]) + nuc_ind in actual_locations:
                    actual_inds.append(int(row["Location"].split(":")[0]) + nuc_ind)
        for actual_ind in actual_inds:

            if row["Direction"] == "left":

                seq = row["gRNA Target Sequence"][::-1]
                sequence = "".join([nucleotide_dict[n] for n in seq])

            else:
                sequence = row["gRNA Target Sequence"]

            transcript_exon = ensembl_object.check_range_info(actual_ind, actual_ind + 1)

            if transcript_id is not None:
                if transcript_id in transcript_exon.keys():
                    transcript_exon = transcript_exon[transcript_id]
                else:
                    transcript_exon = None

            if transcript_exon is not None:
                for transcript in transcript_exon.keys():
                    df = pandas.DataFrame([[sequence, row["gRNA Target Sequence"],
                                            actual_ind, row["Location"], row["Direction"],
                                            ensembl_object.gene_id, transcript]],
                                          columns=["CRISPR Sequence", "gRNA Target Sequence",
                                                   "Edit Location", "CRISPR Location",
                                                   "Direction", "Gene ID (edit)",
                                                   "Transcript ID (edit)"])
                    edit_df = pandas.concat([edit_df, df])
            else:
                df = pandas.DataFrame([[sequence, row["gRNA Target Sequence"], actual_ind,
                                        row["Location"], row["Direction"], ensembl_object.gene_id,None]],
                                      columns=["CRISPR Sequence", "gRNA Target Sequence", "Edit Location",
                                               "CRISPR Location", "Direction", "Gene ID (edit)",
                                               "Transcript ID (edit)"])

                edit_df = pandas.concat([edit_df, df])

    return edit_df


def extract_vep_info(hugo_symbol, edit_location, direction,ensembl_object,
                     edited_nucleotide, new_nucleotide, transcript_id):
    """
    Collect Emsembl VEP infomation for given edits
    :param hugo_symbol: The hugo symbol of the interested gene.
    :param direction: The direction of the CRISPR that can edit the location.
    :param edit_location: The genomic location of the editable nucleotide in activity window.
    :param edited_nucleotide: The interested nucleotide which will be changed with BE.
    :param new_nucleotide: The new nucleotide which will be changed to with BE.
    :param ensembl_object: The Ensembl Object created with Ensembl().
    :param transcript_id: The interested ensembl transcript id
    :return vep_df: The information from VEP API due to the edit.
    :return uniprot_results: The Uniprot IDs in which edit occurs (swissprot or trembl)
    """
    # For (-) direction crisprs, base reversion should be done.
    nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    # Collect chromosome
    chromosome, strand = ensembl_object.chromosome, ensembl_object.strand

    # Dictionary to find the chemical properperty change due to the edit
    aa_chem = {"G": "Non-Polar", "A": "Non-Polar", "V": "Non-Polar", "C": "Polar", "P": "Non-Polar",
               "L": "Non-Polar", "I": "Non-Polar", "M": "Non-Polar", "W": "Non-Polar", "F": "Non-Polar",
               "S": "Polar", "T": "Polar", "Y": "Polar", "N": "Polar", "Q": "Polar", "K": "Charged",
               "R": "Charged", "H": "Charged", "D": "Charged", "E": "Charged"}

    # Decide the server
    server = "http://grch37.rest.ensembl.org" if ensembl_object.assembly == "hg19" \
        else "https://rest.ensembl.org"

    # Base reversion of the (-) direction crisprs
    if direction == "left":
        edited_nucleotide, new_nucleotide = nucleotide_dict[edited_nucleotide], \
                                            nucleotide_dict[new_nucleotide]

    # VEP API request

    vep = "/vep/human/hgvs/%s:g.%s%s>%s?Blosum62=1;Conservation=1;" \
          "LoF=1;CADD=1;protein=1;variant_class=1;canonical=1;hgvs=1;uniprot=1" \
          % (str(chromosome), str(edit_location), edited_nucleotide, new_nucleotide)
    vep_request = requests.get(server + vep, headers={"Content-Type": "application/json"})

    # Check the response of the server for the request
    if vep_request.status_code != 200:
        print("No response from VEP!")
        return None, None

    else:
        uniprot_dict = dict()
        vep_transcript_results, vep_clinical_results = [], []
        vep_transcript_interested = ["hgvsc", "hgvsp", "protein_id", "transcript_id",
                                     "amino_acids", "codons", "polyphen_score",
                                     "polyphen_prediction", "sift_score", "sift_prediction",
                                     "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                     "protein_start", "consequence_terms"]

        VEP_df = pandas.DataFrame(columns = ["hgvsc", "hgvsp", "protein_id", "transcript_id",
                                             "amino_acids", "codons", "polyphen_score",
                                             "polyphen_prediction", "sift_score", "sift_prediction",
                                             "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                             "protein_start", "consequence_terms", "clinical_allele",
                                             "clinical_id","clinical_significance", "ClinVar"])
        vep_dfs = list()
        for x in vep_request.json():
            for t in x["transcript_consequences"]:
                vep_transcript_result = {}
                if t["strand"] == strand:
                    if transcript_id is None:
                        if "canonical" in t.keys():
                            if t["gene_symbol_source"] == "HGNC" and t["gene_symbol"] == hugo_symbol:
                                for opt in vep_transcript_interested:
                                    if opt in t.keys():
                                        vep_transcript_result[opt] = t[opt]
                                    else:
                                        vep_transcript_result[opt] = None
                                if "swissprot" in t.keys():
                                    if t["protein_id"] not in uniprot_dict.keys():
                                        uniprot_dict[t["protein_id"]] = t["swissprot"]
                                    else:
                                        for sw in t["swissprot"]:
                                            if sw not in uniprot_dict[t["protein_id"]]:
                                                uniprot_dict[t["protein_id"]].append(sw)

                                elif "trembl" in t.keys():
                                    if t["protein_id"] not in uniprot_dict.keys():
                                        uniprot_dict[t["protein_id"]] = t["trembl"]
                                    else:
                                        for tr in t["trembl"]:
                                            if tr not in uniprot_dict[t["protein_id"]]:
                                                uniprot_dict[t["protein_id"]].append(tr)

                    else:
                        if t["transcript_id"] == transcript_id:
                            if t["gene_symbol_source"] == "HGNC" and t["gene_symbol"] == hugo_symbol:
                                for opt in vep_transcript_interested:
                                    if opt in t.keys():
                                        vep_transcript_result[opt] = t[opt]
                                    else:
                                        vep_transcript_result[opt] = None
                                if "swissprot" in t.keys():
                                    if t["protein_id"] not in uniprot_dict.keys():
                                        uniprot_dict[t["protein_id"]] = t["swissprot"]
                                    else:
                                        for sw in t["swissprot"]:
                                            if sw not in uniprot_dict[t["protein_id"]]:
                                                uniprot_dict[t["protein_id"]].append(sw)

                                elif "trembl" in t.keys():
                                    if t["protein_id"] not in uniprot_dict.keys():
                                        uniprot_dict[t["protein_id"]] = t["trembl"]
                                    else:
                                        for tr in t["trembl"]:
                                            if tr not in uniprot_dict[t["protein_id"]]:
                                                uniprot_dict[t["protein_id"]].append(tr)

                if vep_transcript_result != {}:
                    trascript_df = pandas.DataFrame.from_dict(vep_transcript_result)
                    vep_transcript_results.append(trascript_df)

            if "colocated_variants" in x.keys():
                for c in x["colocated_variants"]:
                    vep_clinical_result = {}
                    if c["start"] == c["end"] and c["end"] == edit_location:
                        vep_clinical_result["clinical_allele"] = c["allele_string"] \
                            if "allele_string" in c.keys() else None
                        vep_clinical_result["clinical_id"] = c["id"] if "id" in c.keys() else None
                        clinical_sig = ''
                        if "clin_sig" in c.keys():
                            for sig in c["clin_sig"]:
                                sig += ", "
                                clinical_sig += sig
                            if clinical_sig[-2:] == ", ": clinical_sig = clinical_sig[:-2]
                        else:
                            clinical_sig = None
                        vep_clinical_result["clinical_significance"] = clinical_sig
                        clinvars = ''
                        if "var_synonyms" in c.keys():
                            print(c)
                            if type(c["var_synonyms"]) == str:
                            #if "ClinVar" in c["var_synonyms"].keys():
                                for clnv in c["var_synonyms"]:#["ClinVar"]:
                                    clnv += ", "
                                    clinvars += clnv
                                if clinvars != '' and clinvars[-2:] == ", ":
                                    clinvars = clinvars[:-2]
                                elif clinvars == '':
                                    clinvars = None
                        vep_clinical_result["ClinVar"] = clinvars

                    if vep_clinical_result != {}:
                        clinical_df = pandas.DataFrame([vep_clinical_result])
                        vep_clinical_results.append(clinical_df)

            vep_t_df = pandas.DataFrame(columns = ["hgvsc", "hgvsp", "protein_id", "transcript_id",
                                                   "amino_acids", "codons", "polyphen_score",
                                                   "polyphen_prediction", "sift_score", "sift_prediction",
                                                   "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                                   "protein_start", "consequence_terms"])
            vep_c_df = pandas.DataFrame(columns = ["clinical_allele","clinical_id",
                                                   "clinical_significance", "ClinVar"])

            if vep_transcript_results: vep_t_df = pandas.concat(vep_transcript_results)
            if vep_clinical_results: vep_c_df = pandas.concat(vep_clinical_results)

            if vep_t_df.empty:
                vep_t_df = pandas.DataFrame(
                    None, index=[0], columns=["hgvsc", "hgvsp", "protein_id", "transcript_id",
                                              "amino_acids", "codons", "polyphen_score",
                                              "polyphen_prediction", "sift_score", "sift_prediction",
                                              "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                              "protein_start", "consequence_terms"])
            if vep_c_df.empty:
                vep_c_df = pandas.DataFrame(None, index = [0],
                                            columns = ["clinical_allele","clinical_id",
                                                       "clinical_significance", "ClinVar"])

            vep_t_df["merging"], vep_c_df["merging"] = 1, 1
            vep_df = pandas.merge(vep_t_df, vep_c_df, on ="merging").drop(columns=["merging"])

            if vep_df.empty is False:
                vep_dfs.append(vep_df)

        VEP_df = pandas.concat(vep_dfs)

        if VEP_df.empty is False:
            VEP_df["Genomic_Position"] = edit_location
            VEP_df["Direction"] = direction
            VEP_df["Transcrip_ID"] = VEP_df.apply(
                lambda x: x.transcript_id
                if x.transcript_id is not None and pandas.isna(x.transcript_id) is False and
                   type(x.transcript_id) != float else None, axis=1)

            VEP_df["cDNA_Change"] = VEP_df.apply(
                lambda x: x.hgvsc.split(":")[1].split(".")[1]
                if x.hgvsc is not None and pandas.isna(x.hgvsc) is False and
                   type(x.hgvsc) != float else None, axis=1)

            VEP_df["Protein_ID"] = VEP_df.apply(
                lambda x: x.protein_id if x.protein_id is not None else None, axis=1)

            VEP_df["Protein_Change"] = VEP_df.apply(
                lambda x: x.hgvsp.split(":")[1].split(".")[1]
                if x.hgvsp is not None and pandas.isna(x.hgvsp) is False and type(x.hgvsp) != float else None, axis=1)

            VEP_df["Protein_Position"] = VEP_df.apply(
                lambda x: x.protein_start if x.protein_start is not None else None, axis=1)

            VEP_df["Edited_aa"] = VEP_df.apply(
                lambda x: x.amino_acids.split("/")[0]
                if x.amino_acids is not None and pandas.isna(x.amino_acids) is False and
                   type(x.amino_acids) != float and len(x.amino_acids.split("/")) > 1 else None, axis=1)

            VEP_df["New_aa"] = VEP_df.apply(
                lambda x: x.amino_acids.split("/")[1]
                if x.amino_acids is not None and pandas.isna(x.amino_acids) is False and
                   type(x.amino_acids) != float and len(x.amino_acids.split("/")) > 1 else (
                    x.amino_acids.split("/")[0] if x.amino_acids is not None and pandas.isna(x.amino_acids) is False and
                                                   type(x.amino_acids) != float and len(x.amino_acids.split("/")) == 1
                    else None), axis=1)

            VEP_df["Edited_aa_prop"] = VEP_df.apply(
                lambda x: aa_chem[x.Edited_aa] if x.Edited_aa is not None else None, axis=1)
            VEP_df["New_aa_prop"] = VEP_df.apply(
                lambda x: aa_chem[x.New_aa] if x.New_aa is not None and x.New_aa != "*" else None, axis=1)

            VEP_df["Edited_codon"] = VEP_df.apply(
                lambda x: x.codons.split("/")[0]
                if x.codons is not None and pandas.isna(x.codons) is False and
                   type(x.codons) != float else None, axis=1)

            VEP_df["New_codon"] = VEP_df.apply(
                lambda x: x.codons.split("/")[1]
                if x.codons is not None and pandas.isna(x.codons) is False and
                   type(x.codons) != float else None, axis=1)

            VEP_df["is_synonymous"] = VEP_df.apply(
                lambda x: True if x.Edited_codon is not None and x.New_codon is not None and
                                  x.Edited_aa is not None and x.New_aa is not None and
                                  x.Edited_codon != x.New_codon and
                                  x.Edited_aa == x.New_aa else False, axis=1)
            VEP_df["is_stop"] = VEP_df.apply(
                lambda x: True if x.New_codon is not None and
                                  x.New_codon in ["UAA", "UAG", "UGA"] else False, axis=1)
            VEP_df["is_ClinVar"] = VEP_df.apply(
                lambda x: True if x.ClinVar is not None else False, axis=1)

            VEP_df = VEP_df[["Genomic_Position", "Direction", "Transcrip_ID", "cDNA_Change",
                             "Edited_codon", "New_codon", "Protein_ID", "Protein_Position",
                             "Protein_Change", "Edited_aa", "Edited_aa_prop", "New_aa", "New_aa_prop",
                             "is_synonymous", "is_stop", "polyphen_score", "polyphen_prediction",
                             "sift_score", "sift_prediction", "cadd_phred", "cadd_raw", "lof",
                             "impact", "blosum62", "consequence_terms", "is_ClinVar", "clinical_allele",
                             "clinical_id", "clinical_significance", "ClinVar"]]

        return vep_request, VEP_df, uniprot_dict


def after_vep_analysis(ensembl_object, vep_df):
    """
    Adding Uniprot API Information on VEP DF
    :param ensembl_object: The object of the Ensembl from Ensembl API
    :param vep_df: The data frame filled with the information from VEP API
    :return: analysis_df: The data frame ensriched with the information from Uniprot API
    """

    analysis_df = pandas.DataFrame(columns=["Genomic_Position", "Direction", "Transcrip_ID",
                                            "cDNA_Change", "Edited_codon", "New_codon",
                                            "Protein_ID", "Protein_Position", "Protein_Change",
                                            "Edited_aa", "Edited_aa_prop", "New_aa", "New_aa_prop",
                                            "is_synonymous", "is_stop", "polyphen_score",
                                            "polyphen_prediction", "sift_score", "sift_prediction",
                                            "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                            "consequence_terms", "is_ClinVar", "clinical_allele",
                                            "clinical_id", "clinical_significance", "ClinVar",
                                            "Domain", "Phosphorylation", "Uniprot", "Reviewed"])
    analysis_dfs = list()

    for ind, row in vep_df.iterrows():
        seq_mapping = ensembl_object.extract_uniprot_info(row.Protein_ID)
        uniprot_obj_dict = dict()
        if seq_mapping is not None:
            doms, phoss = {}, {}
            no_dom_phos_uniprot = []

            for uniprot in seq_mapping.keys():
                uniprot_object = Uniprot(uniprotid=uniprot)
                uniprot_object.extract_uniprot_info()
                if row["Protein_Position"] in seq_mapping[uniprot].keys():
                    if uniprot in uniprot_object.keys():
                        dom = uniprot_object.find_domain(
                            seq_mapping[uniprot][row["Protein_Position"]], row["Edited_aa"])
                        phos = uniprot_object.find_psite(
                            seq_mapping[uniprot][row["Protein_Position"]], row["Edited_aa"])

                        if dom is not None:
                            if uniprot not in doms.keys():
                                doms[uniprot] = [dom]
                            else:
                                if dom not in doms[uniprot]:
                                    doms[uniprot].append(dom)

                        if phos is not None:
                            if uniprot not in phoss.keys():
                                phoss[uniprot] = phos
                            else:
                                if phos != phoss[uniprot]:
                                    # Cannot be a different phosphosite
                                    del phoss[uniprot]

                        if dom is None and phos is None:
                            no_dom_phos_uniprot.append(uniprot)
                else:
                    no_dom_phos_uniprot.append(uniprot)

            if doms != {}:
                for u, dom_list in doms.items():
                    for d in dom_list:
                        dom = d
                        if phoss != {} and u in phoss.keys():
                            phos = phoss[u]
                        else:
                            phos = None

                        df = pandas.DataFrame([[row["Genomic_Position"], row["Direction"],
                                                row["Transcrip_ID"], row["cDNA_Change"],
                                                row["Edited_codon"], row["New_codon"],
                                                row["Protein_ID"], row["Protein_Position"],
                                                row["Protein_Change"], row["Edited_aa"],
                                                row["Edited_aa_prop"], row["New_aa"],
                                                row["New_aa_prop"], row["is_synonymous"],
                                                row["is_stop"], row["polyphen_score"],
                                                row["polyphen_prediction"], row["sift_score"],
                                                row["sift_prediction"], row["cadd_phred"],
                                                row["cadd_raw"], row["lof"], row["impact"],
                                                row["blosum62"], row["consequence_terms"],
                                                row["is_ClinVar"], row["clinical_allele"],
                                                row["clinical_id"], row["clinical_significance"],
                                                row["ClinVar"], dom, phos, u,
                                                uniprot_obj_dict[u].reviewed]],
                                              columns=["Genomic_Position", "Direction", "Transcrip_ID",
                                                       "cDNA_Change", "Edited_codon", "New_codon",
                                                       "Protein_ID", "Protein_Position", "Protein_Change",
                                                       "Edited_aa", "Edited_aa_prop", "New_aa", "New_aa_prop",
                                                       "is_synonymous", "is_stop", "polyphen_score",
                                                       "polyphen_prediction", "sift_score", "sift_prediction",
                                                       "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                                       "consequence_terms", "is_ClinVar", "clinical_allele",
                                                       "clinical_id", "clinical_significance", "ClinVar",
                                                       "Domain", "Phosphorylation", "Uniprot", "Reviewed"])
                        analysis_dfs.append(df)

            elif phoss != {}:
                for u, p in phoss.items():
                    dom, phos = None, p
                    df = pandas.DataFrame([[row["Genomic_Position"], row["Direction"],
                                            row["Transcrip_ID"], row["cDNA_Change"],
                                            row["Edited_codon"], row["New_codon"],
                                            row["Protein_ID"], row["Protein_Position"],
                                            row["Protein_Change"], row["Edited_aa"],
                                            row["Edited_aa_prop"], row["New_aa"],
                                            row["New_aa_prop"], row["is_synonymous"],
                                            row["is_stop"], row["polyphen_score"],
                                            row["polyphen_prediction"], row["sift_score"],
                                            row["sift_prediction"], row["cadd_phred"],
                                            row["cadd_raw"], row["lof"], row["impact"],
                                            row["blosum62"], row["consequence_terms"],
                                            row["is_ClinVar"], row["clinical_allele"],
                                            row["clinical_id"], row["clinical_significance"],
                                            row["ClinVar"], dom, phos, u,
                                            uniprot_obj_dict[u].reviewed]],
                                          columns=["Genomic_Position", "Direction", "Transcrip_ID",
                                                   "cDNA_Change", "Edited_codon", "New_codon",
                                                   "Protein_ID", "Protein_Position", "Protein_Change",
                                                   "Edited_aa", "Edited_aa_prop", "New_aa", "New_aa_prop",
                                                   "is_synonymous", "is_stop", "polyphen_score",
                                                   "polyphen_prediction", "sift_score", "sift_prediction",
                                                   "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                                   "consequence_terms", "is_ClinVar", "clinical_allele",
                                                   "clinical_id", "clinical_significance", "ClinVar",
                                                   "Domain", "Phosphorylation", "Uniprot", "Reviewed"])
                    analysis_dfs.append(df)

            for u in no_dom_phos_uniprot:
                dom, phos = None, None
                df = pandas.DataFrame([[row["Genomic_Position"], row["Direction"],
                                        row["Transcrip_ID"], row["cDNA_Change"],
                                        row["Edited_codon"], row["New_codon"],
                                        row["Protein_ID"], row["Protein_Position"],
                                        row["Protein_Change"], row["Edited_aa"],
                                        row["Edited_aa_prop"], row["New_aa"],
                                        row["New_aa_prop"], row["is_synonymous"],
                                        row["is_stop"], row["polyphen_score"],
                                        row["polyphen_prediction"], row["sift_score"],
                                        row["sift_prediction"], row["cadd_phred"],
                                        row["cadd_raw"], row["lof"], row["impact"],
                                        row["blosum62"], row["consequence_terms"],
                                        row["is_ClinVar"], row["clinical_allele"],
                                        row["clinical_id"], row["clinical_significance"],
                                        row["ClinVar"], dom, phos, u,
                                        uniprot_obj_dict[u].reviewed]],
                                      columns=["Genomic_Position", "Direction", "Transcrip_ID",
                                               "cDNA_Change", "Edited_codon", "New_codon",
                                               "Protein_ID", "Protein_Position", "Protein_Change",
                                               "Edited_aa", "Edited_aa_prop", "New_aa", "New_aa_prop",
                                               "is_synonymous", "is_stop", "polyphen_score",
                                               "polyphen_prediction", "sift_score", "sift_prediction",
                                               "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                               "consequence_terms", "is_ClinVar", "clinical_allele",
                                               "clinical_id", "clinical_significance", "ClinVar",
                                               "Domain", "Phosphorylation", "Uniprot", "Reviewed"])
                analysis_dfs.append(df)

    if analysis_dfs:
        analysis_df = pandas.concat(analysis_dfs)
        return analysis_df
    else: return None


###########################################################################################
# Execution

def main():
    """
    Run whole script with the input from terminal
    :return:
    """

    print("************************************************************")
    print("\n\nRetrieving Ensembl Gene Information...")

    ensembl_obj = Ensembl(hugo_symbol=args["GENE"], assembly=args["ASSEMBLY"])
    ensembl_obj.extract_gene_id()

    if ensembl_obj.gene_id == '': sys.exit("No corresponding Ensembl Gene ID could be found!")

    ensembl_obj.extract_sequence(ensembl_obj.gene_id)

    if ensembl_obj.gene_range[0] < ensembl_obj.gene_range[1]:
        ensembl_obj.extract_info(chromosome=ensembl_obj.chromosome,
                                 loc_start=ensembl_obj.gene_range[0],
                                 loc_end=ensembl_obj.gene_range[1])
    else:
        ensembl_obj.extract_info(chromosome=ensembl_obj.chromosome,
                                 loc_start=ensembl_obj.gene_range[1],
                                 loc_end=ensembl_obj.gene_range[0])

    print("\nEntering extract_gRNA_sites function...\n")
    crispr_df, crispr_nuc_df = extract_grna_sites(searched_nucleotide=args["EDIT"],
                                                  pam_window=[int(args["PAMWINDOW"].split("-")[0]),
                                                              int(args["PAMWINDOW"].split("-")[1])],
                                                  activity_window=[int(args["ACTWINDOW"].split("-")[0]),
                                                                   int(args["ACTWINDOW"].split("-")[1])],
                                                  pam_sequence=args["PAMSEQ"],
                                                  protospacer_length=args["PROTOLEN"],
                                                  ensembl_object=ensembl_obj)

    path = ""
    if args["OUTPUT_PATH"][-1] == "/":
        path = args["OUTPUT_PATH"]
    else:
        path = args["OUTPUT_PATH"] + "/"

    if len(crispr_df.index) != 0: print("CRISPR Data Frame Created!")
    crispr_df.to_csv(path + args["OUTPUT_FILE"] + "_crispr_df.csv")

    print("\nCRISPR Data Frame wrote in %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_crispr_df.csv"))

    print("Entering find_editable_nucleotide function...\n")
    edit_df = find_editable_nucleotide(crispr_df=crispr_nuc_df,
                                       searched_nucleotide=args["EDIT"],
                                       activity_window=[int(args["ACTWINDOW"].split("-")[0]),
                                                        int(args["ACTWINDOW"].split("-")[1])],
                                       ensembl_object=ensembl_obj, transcript_id=args["TRANSCRIPT"])

    if len(edit_df.index) != 0: print("EDIT Data Frame Created!")

    edit_df.to_csv(path + args["OUTPUT_FILE"] + "_edit_df.csv")

    print("\nEDIT Data Frame wrote in %s as %s" % (path, args["OUTPUT_FILE"] + "_edit_df.csv\n\n"))
    
    if args["PROTEIN"]:
        print("\nExtracting VEP Information...")

        loc_edit_df = edit_df[["Edit Location", "Direction"]]
        loc_edit_df = loc_edit_df.drop_duplicates()

        whole_vep_df = pandas.DataFrame(columns=["Genomic_Position", "Direction", "Transcrip_ID",
                                                 "cDNA_Change", "Edited_codon", "New_codon",
                                                 "Protein_ID", "Protein_Position", "Protein_Change",
                                                 "Edited_aa", "Edited_aa_prop", "New_aa", "New_aa_prop",
                                                 "is_synonymous", "is_stop", "polyphen_score",
                                                 "polyphen_prediction", "sift_score", "sift_prediction",
                                                 "cadd_phred", "cadd_raw", "lof", "impact", "blosum62",
                                                 "consequence_terms", "is_ClinVar", "clinical_allele",
                                                 "clinical_id", "clinical_significance", "ClinVar"])

        for ind, row in loc_edit_df.iterrows():

            api_response, vep_df, uniprots = extract_vep_info(hugo_symbol=args["GENE"],
                                                              edit_location=row["Edit Location"],
                                                              edited_nucleotide=args["EDIT"],
                                                              new_nucleotide=args["EDIT_TO"],
                                                              transcript_id=args["TRANSCRIPT"],
                                                              direction=row["Direction"],
                                                              ensembl_object=ensembl_obj)

            if len(vep_df.index) != 0:
                whole_vep_df = pandas.concat([whole_vep_df, vep_df])

        if len(whole_vep_df.index) != 0:
            print("\nVEP Data Frame Created!")
            whole_vep_df.to_csv(path + args["OUTPUT_FILE"] + "_vep_df.csv")
            print("\nVEP Data Frame wrote in %s as %s\n\n" % (path, args["OUTPUT_FILE"] + "_vep_df.csv"))
        else:
            print("\nVEP Data Frame cannot be created because it is empty.")

        print("\nExtracting Uniprot Information...")
        analysis_df = after_vep_analysis(ensembl_object= ensembl_obj, vep_df= whole_vep_df)

        if len(analysis_df.index) != 0:
            print("\nAnalysis Data Frame Created!")
            analysis_df.to_csv(path + args["OUTPUT_FILE"] + "_analysis_df.csv")
            print("\nAnalysis Data Frame wrote in %s as %s\n\n" % (path, args["OUTPUT_FILE"] +
                                                                   "_analysis_df.csv"))
        else:
            print("\nAnalysis Data Frame cannot be created because it is empty.")
        print("************************************************************\n\n")
        return True

    else: return True

main()





