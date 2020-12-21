###########################################################################################
#################### Down sampling of the SGE data to compare with BEs ####################
###########################################################################################
# Import necessary packages
import os, pandas, numpy, matplotlib, requests, re, warnings
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.stats as stats
from scipy.stats import chi2_contingency
import scikit_posthocs
from statannot import add_stat_annotation
warnings.simplefilter(action='ignore', category=FutureWarning)

font = {'family': 'normal','weight': 'normal','size': 8}
matplotlib.rc('font', **font)


###########################################################################################
# Classes

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
                                    if phos is not None:
                                        self.phosphorylation_sites[pos] = phos

                            if ftr["category"] == "DOMAINS_AND_SITES":
                                if "description" in ftr.keys():
                                    domain = ftr["description"]
                                    domain_range = list(range(int(ftr["begin"]), int(ftr["end"])+1))
                                    self.domains[domain] = domain_range

            if self.phosphorylation_sites == dict(): self.phosphorylation_sites = None
            if self.domains == dict(): self.domains = None

            return "UniProt API request is done."

###########################################################################################
# Function


def sge_file(file_path, wanted_columns, new_columns):

    # Read the SGE file
    sge_df = pandas.read_csv(file_path, index_col = 0)

    # Take only the interested information
    sge_df = sge_df[wanted_columns]

    # Change the column name as interested
    sge_df.columns = new_columns

    return sge_df


def uniprot_seq_checking(uniprot_obj, sge_df, aa_pos_column, aa_ref_column):
    # Take only the coding region information from SGE data
    aa_df = sge_df[[aa_pos_column, aa_ref_column]]
    aa_df = aa_df[~pandas.isna(aa_df[aa_pos_column])]
    aa_df = aa_df.drop_duplicates()
    aa_df = aa_df.set_index([aa_pos_column])
    aa_df["uniprot_aa"] = aa_df.apply(lambda x: uniprot_obj.sequence[int(x.name) - 1], axis=1)
    aa_df["check"] = aa_df.apply(
        lambda x: True if x[aa_ref_column] == x["uniprot_aa"] else False, axis=1)
    if set(aa_df["check"]) == {True}:
        return True
    else:
        return False


def find_domain(pos, dom_dict):
    d = []
    for i,l in dom_dict.items():
        if pos in l:
            if i not in d:
                d.append(i)
    if len(d) == 1:
        return d[0]
    else: return "No Domain"


def add_domain_info(sge_df, aa_pos_column, uniprot_obj):
    sge_df["domain"] = sge_df.apply(
        lambda x: find_domain(x[aa_pos_column], uniprot_obj.domains), axis=1)

    return sge_df


def be_positions(file_path, be, pam):

    # Read the edit_df from PAM_Extractor and take the editable positions
    bes = list(set(pandas.read_csv(file_path, index_col=0)["Edit Location"]))

    return bes


def add_be_annotation(sge_df, ref_columns_name, alt_column_name,
                      position_column_name,  be_file, be, pam, ref, alt):

    # Take editable positions according to given BE information
    be_pos = be_positions(file_path=be_file, be=be, pam=pam)

    # Prepare the column name for BE
    be_annot = be + "_" + pam

    # Only take the edits from SGE both reference and altered nucleotides are given BE target
    sge_filt_df = sge_df[(sge_df[ref_columns_name] == ref) & (sge_df[alt_column_name] == alt)]

    # Label the positions if BE can potentially edit
    print(sge_filt_df)
    sge_filt_df[be_annot] = sge_filt_df.apply(
        lambda x: True if x[position_column_name] in be_pos else False, axis=1)

    # Merge the all SGE data with the labelled data
    sge_enrich_df = pandas.merge(sge_df, sge_filt_df, how='left',
                                 on=list(sge_filt_df.columns).remove(be_annot))

    # Label the positions if BE can potentially edit
    # We make it again since we did not want to label True for uneditable
    # changes because of their positions
    sge_enrich_df[be_annot] = sge_enrich_df.apply(
        lambda x: False if pandas.isna(x[be_annot]) else (True if x[be_annot] else False), axis=1)
    return sge_enrich_df


def percentage(total, sample):
    if sample != numpy.nan:
        return (sample * 100.0)/total
    else: return numpy.nan


def func_class_percentages(row, df):
    d = {"LOF": 0, "INT": 0, "FUNC": 0}
    for i in df.loc[row]:
        if pandas.isna(i) is False:
            if i > 0.99:
                d["LOF"] += 1
            elif 0.01 < i < 0.99:
                d["INT"] += 1
            elif i < 0.01:
                d["FUNC"] += 1

    p = {i: ((l * 100.0) / sum(d.values())) / 100.0 for i, l in d.items()}
    return p


def pivot_sge(df, index_col, column_col, value_col, aa_position_col, coding):
    if coding:
        df = df[~pandas.isna(df[aa_position_col])]
        df[aa_position_col] = df.apply(lambda x: int(x[aa_position_col]),axis=1)
    piv_df = pandas.pivot_table(df, index= index_col, columns=column_col, values=value_col,
                                fill_value=numpy.nan)

    return piv_df


def rgb_rgbtples(rgb):
    t = list()
    for i in rgb:
        rbg_ = i/255.0
        t.append(rbg_)
    return tuple(t)


def normalisation(val, all):
    n = (val - min(all)) / float((max(all) - min(all)))
    return n


# according to NGG and NG
def check_pam_positions(pam_df, index):
    if len(set(pam_df.loc[index])) == 1 and set(pam_df.loc[index]) == {1}:
        return 2
    elif len(set(pam_df.loc[index])) == 1 and set(pam_df.loc[index]) == {0}:
        return 0
    elif len(set(pam_df.loc[index])) == 2: return 1


def figure_heatmap(sge_df, aa_position_col, aa_ref_col,pam_list, be_dict,func_class_col,
                   aa_alt_col, func_score_col, nonfunctional_p_col, uniprot_object):

    coding_score_piv_df = pivot_sge(df=sge_df, index_col=[aa_position_col, aa_ref_col],
                                    column_col=[aa_alt_col], value_col=func_score_col,
                                    aa_position_col=aa_position_col, coding=True)

    # Domain Info
    domain_df = pandas.DataFrame(None, index=coding_score_piv_df.index, columns=["Domain"])
    domain_df = domain_df.reset_index()
    domain_df["domain"] = domain_df.apply(
        lambda x: find_domain(x[aa_position_col], uniprot_object.domains), axis=1)
    domain_df = domain_df.set_index([aa_position_col, aa_ref_col])
    #dom_value_to_int = {j: i for i, j in enumerate(pandas.unique(domain_df.values.ravel()))}
    domain_df["Domain"] = domain_df.apply(
        lambda x: 0 if x.domain == "No Domain" else (
            1 if x.domain == "BRCT 1" else (
                2 if x.domain == "BRCT 2" else (
                    3 if x.domain == "RING-type" else None))), axis =1)

    domain_df["norm_domain"] = domain_df.apply(
        lambda x: normalisation(x.Domain, list(set(domain_df.Domain))), axis=1)

    # PAM Info
    be_pam_dfs = {}
    for be, be_info in be_dict.items():
        positions = {pam: set(pos_dom_sge_df[pos_dom_sge_df[be + "_" + pam]]["aa_pos"])
                     for pam in pam_list}
        pam_dfs = {}
        for pam in pam_list:
            df = pandas.DataFrame(None, index=coding_score_piv_df.index,
                                  columns=[be+"_"+pam, "Class_" + pam])
            df = df.reset_index()
            df[be+ "_" + pam] = df.apply(
                lambda x: 1 if x[aa_position_col] in positions[pam] else 0, axis=1)
            df["Class_" + pam] = df.apply(
                lambda x: 1 if list(pos_dom_sge_df[
                                        (pos_dom_sge_df.aa_pos == x[aa_position_col]) &
                                        (pos_dom_sge_df.alt == be_info["alt"])]
                                    [func_class_col].values) == ["LOF"] else 0, axis =1)

            df = df.set_index([aa_position_col, aa_ref_col])
            pam_dfs[pam] =df

        pam_df = pandas.concat(pam_dfs.values(), axis = 1)
        pam_df[be] = pam_df.apply(lambda x: check_pam_positions(pam_df, x.name), axis=1)

        be_pam_dfs[be] = pam_df

    be_lol_piv_df = pandas.DataFrame(None, index=coding_score_piv_df.index)
    for be, be_df in be_pam_dfs.items():
        be_lol_piv_df = pandas.concat([be_lol_piv_df, be_df], axis =1)

    be_lol_piv_df["norm_cbe"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.CBE, list(set(be_lol_piv_df.CBE))), axis=1)
    be_lol_piv_df["norm_abe"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.ABE, list(set(be_lol_piv_df.ABE))), axis=1)
    be_lol_piv_df["norm_cgbe"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.CGBE, list(set(be_lol_piv_df.CGBE))), axis=1)

    coding_p_piv_df = pivot_sge(df=sge_df, index_col=[aa_position_col, aa_ref_col],
                                column_col=[aa_alt_col], value_col=nonfunctional_p_col,
                                aa_position_col=aa_position_col, coding=True)

    y_coding_perc_pivot_df = coding_p_piv_df.copy()
    x_coding_p_piv_df = coding_p_piv_df.T
    x_coding_perc_pivot_df = x_coding_p_piv_df.copy()

    # Y bar plot
    y_coding_perc_pivot_df["LOF_Percentages"] = coding_p_piv_df.apply(
        lambda x: func_class_percentages(x.name, coding_p_piv_df)["LOF"], axis=1)

    y_coding_perc_pivot_df["INT_Percentages"] = coding_p_piv_df.apply(
        lambda x: func_class_percentages(x.name, coding_p_piv_df)["INT"], axis=1)

    y_coding_perc_pivot_df["FUNC_Percentages"] = coding_p_piv_df.apply(
        lambda x: func_class_percentages(x.name, coding_p_piv_df)["FUNC"], axis=1)

    y_coding_perc_pivot_df = y_coding_perc_pivot_df[["LOF_Percentages", "INT_Percentages", "FUNC_Percentages"]]

    # X bar plot

    x_coding_perc_pivot_df["LOF_Percentages"] = x_coding_p_piv_df.apply(
        lambda x: func_class_percentages(x.name, x_coding_p_piv_df)["LOF"], axis=1)

    x_coding_perc_pivot_df["INT_Percentages"] = x_coding_p_piv_df.apply(
        lambda x: func_class_percentages(x.name, x_coding_p_piv_df)["INT"], axis=1)

    x_coding_perc_pivot_df["FUNC_Percentages"] = x_coding_p_piv_df.apply(
        lambda x: func_class_percentages(x.name, x_coding_p_piv_df)["FUNC"], axis=1)

    x_coding_perc_pivot_df = x_coding_perc_pivot_df[["LOF_Percentages", "INT_Percentages", "FUNC_Percentages"]]

    # Heatmap visual info
    pal = sns.color_palette("coolwarm_r", 3)
    color_sets = pal.as_hex()

    dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                       rgb_rgbtples((196, 56, 90)),
                                       rgb_rgbtples((99, 42, 173)),
                                       rgb_rgbtples((43, 38, 136))])

    sns.set_theme(style="ticks")

    fig = plt.figure(figsize=(18, 25))

    # Subplots
    widths = [14.2, 0.2, 0.2, 0.2, 0.2, 3]
    heights = [3, 22, 1, 1, 1, 1]

    gs = fig.add_gridspec(6, 6, wspace=0.05, hspace=0.03, width_ratios=widths,
                          height_ratios=heights)

    ax = fig.add_subplot(gs[1:, :-5])
    ax.xaxis.set_ticks_position('none')

    ax_x = fig.add_subplot(gs[0, :-5], sharex=ax)
    ax_x.xaxis.set_ticks_position('none')

    ax_y = fig.add_subplot(gs[1:, -1], sharey=ax)

    ax_y2 = fig.add_subplot(gs[1:, -5], sharey=ax)
    ax_y2.yaxis.set_ticks_position('none')
    ax_y2.xaxis.set_ticks_position('none')
    ax_y2.axes.get_xaxis().set_visible(False)

    ax_y3 = fig.add_subplot(gs[1:, -4], sharey=ax)
    ax_y3.yaxis.set_ticks_position('none')
    ax_y3.xaxis.set_ticks_position('none')
    ax_y3.axes.get_xaxis().set_visible(False)

    ax_y4 = fig.add_subplot(gs[1:, -3], sharey=ax)
    ax_y4.yaxis.set_ticks_position('none')
    ax_y4.xaxis.set_ticks_position('none')
    ax_y4.axes.get_xaxis().set_visible(False)

    ax_y5 = fig.add_subplot(gs[1:, -2], sharey=ax)
    ax_y5.yaxis.set_ticks_position('none')
    ax_y5.xaxis.set_ticks_position('none')
    ax_y5.axes.get_xaxis().set_visible(False)

    # Fill axes
    ax.pcolor(coding_score_piv_df, cmap="coolwarm_r", linewidths=.1, edgecolors="lightgrey")
    #ax_y2.pcolor(domain_df.replace(dom_value_to_int), linewidths=0, cmap="cubehelix_r")
    ax_y2.pcolor(domain_df[["norm_domain"]], linewidths=0, cmap=dom_cmap)
    #ax_y3.pcolor(be_pam_dfs["CBE"], linewidths=0, cmap="PuRd")
    ax_y3.pcolor(be_lol_piv_df[["norm_cbe"]], linewidths=0, cmap="PuRd")
    #ax_y4.pcolor(be_pam_dfs["ABE"], linewidths=0, cmap="Purples")
    ax_y4.pcolor(be_lol_piv_df[["norm_abe"]], linewidths=0, cmap="Purples")
    #ax_y5.pcolor(be_pam_dfs["CGBE"], linewidths=0, cmap="Blues")
    ax_y5.pcolor(be_lol_piv_df[["norm_cgbe"]], linewidths=0, cmap="Blues")

    y_coding_perc_pivot_df.plot.barh(stacked=True, ax=ax_y, color=color_sets, legend=False,
                                     align='edge', width=1)
    x_coding_perc_pivot_df.plot.bar(stacked=True, ax=ax_x, color=color_sets,
                                    legend=False, align='edge', width=1)

    ax.set_yticks(numpy.arange(len(coding_score_piv_df.index)) + 0.5)
    ax.set_yticklabels([str(i[0]) + "-" + i[1] for i in coding_score_piv_df.index], fontsize=4)

    ax.set_xticks(numpy.arange(len(coding_score_piv_df.columns)) + 0.5)
    ax.set_xticklabels(coding_score_piv_df.columns)

    ax.set_ylabel("Reference Amino Acid", fontsize=10)
    ax.set_xlabel("Altered Amino Acids", fontsize=10)
    plt.show()
    plt.savefig("/Users/cd7/Desktop/last/Figure_3.png")
    plt.savefig("/Users/cd7/Desktop/last/Figure_3.pdf")

    # Presentation
    """
    sns.set_theme(style="ticks")
    fig = plt.figure(figsize=(14, 8))

    # Subplots
    widths = [12, 0, 0, 0, 0, 2]
    heights = [1, 0.1, 0.1, 0.1, 0.1, 6.7]

    gs = fig.add_gridspec(6, 6, wspace=0.05, hspace=0.03, width_ratios=widths,
                          height_ratios=heights)

    ax = fig.add_subplot(gs[-1:, :-1])
    ax.xaxis.set_ticks_position('none')

    ax_x = fig.add_subplot(gs[0, :-1], sharex=ax)
    ax_x.xaxis.set_ticks_position('none')
    ax_x.axes.get_xaxis().set_visible(False)
    ax_x.yaxis.set_ticks_position('none')
    ax_x.axes.get_yaxis().set_visible(False)

    ax_x2 = fig.add_subplot(gs[1, :-1], sharex=ax)
    ax_x2.xaxis.set_ticks_position('none')
    ax_x2.axes.get_xaxis().set_visible(False)
    ax_x2.yaxis.set_ticks_position('none')
    ax_x2.axes.get_yaxis().set_visible(False)

    ax_x3 = fig.add_subplot(gs[2, :-1], sharex=ax)
    ax_x3.xaxis.set_ticks_position('none')
    ax_x3.axes.get_xaxis().set_visible(False)
    ax_x3.yaxis.set_ticks_position('none')
    ax_x3.axes.get_yaxis().set_visible(False)

    ax_x4 = fig.add_subplot(gs[3, :-1], sharex=ax)
    ax_x4.xaxis.set_ticks_position('none')
    ax_x4.axes.get_xaxis().set_visible(False)
    ax_x4.yaxis.set_ticks_position('none')
    ax_x4.axes.get_yaxis().set_visible(False)

    ax_x5 = fig.add_subplot(gs[4, :-1], sharex=ax)
    ax_x5.xaxis.set_ticks_position('none')
    ax_x5.axes.get_xaxis().set_visible(False)
    ax_x5.yaxis.set_ticks_position('none')
    ax_x5.axes.get_yaxis().set_visible(False)

    ax_y = fig.add_subplot(gs[-1:, -1], sharey=ax)

    dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                       rgb_rgbtples((235, 115, 8)),
                                       rgb_rgbtples((173, 54, 110)),
                                       rgb_rgbtples((61, 38, 166))])

    ax.pcolor(coding_score_piv_df.T, cmap="coolwarm_r", linewidths=.1, edgecolors="lightgrey")
    ax_x2.pcolor(domain_df.replace(dom_value_to_int).T, linewidths=0, cmap=dom_cmap)
    ax_x3.pcolor(be_pam_dfs["CBE"].T, linewidths=0, cmap="PuRd")
    ax_x4.pcolor(be_pam_dfs["ABE"].T, linewidths=0, cmap="Purples")
    ax_x5.pcolor(be_pam_dfs["CGBE"].T, linewidths=0, cmap="Blues")

    y_coding_perc_pivot_df.plot.bar(stacked=True, ax=ax_x, color=color_sets, legend=False,
                                    align='edge', width=1)

    ax_x.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_x.axes.get_xaxis().set_visible(False)

    x_coding_perc_pivot_df.plot.barh(stacked=True, ax=ax_y, color=color_sets, legend=False,
                                     align='edge', width=1)

    ax.set_xticks(numpy.arange(len(coding_score_piv_df.index)) + 0.5)
    ax.set_xticklabels([str(i[0]) + "-" + i[1] for i in coding_score_piv_df.index], fontsize=3, rotation=90)

    ax.set_yticks(numpy.arange(len(coding_score_piv_df.columns)) + 0.5)
    ax.set_yticklabels(coding_score_piv_df.columns)

    ax.set_ylabel("Altered Amino Acids", fontsize=10)
    ax.set_xlabel('Reference Amino Acid', fontsize=10)
    plt.show()
    plt.savefig("/Users/cd7/Desktop/last/Figure_3_presentation.png")
    plt.savefig("/Users/cd7/Desktop/last/Figure_3_presentation.pdf")
    """
    return True


def figure_lollipop(sge_df, position_col, ref_col, alt_col, aa_pos_col,
                    functional_score_col, nonfunctional_p_col, pam_list, be_dict):

    lol_piv_df = pivot_sge(df = sge_df, index_col=[position_col, ref_col],
                           column_col=[alt_col], value_col=functional_score_col,
                           aa_position_col=aa_pos_col, coding = False)

    lol_piv_df["Max"] = lol_piv_df.apply(
        lambda x: numpy.nanmax(list(lol_piv_df.loc[x.name])), axis=1)
    lol_piv_df["Min"] = lol_piv_df.apply(
        lambda x: numpy.nanmin(list(lol_piv_df.loc[x.name])), axis=1)

    lol_piv_df = lol_piv_df.reset_index()
    lol_piv_df["position"] = lol_piv_df.apply(lambda x: str(x[position_col]), axis=1)
    lol_piv_df = lol_piv_df.set_index(["position"])
    lol_piv_df = lol_piv_df[[ref_col, 'A', 'C', 'G', 'T', 'Max', 'Min']]

    # Domain

    lol_piv_df["domain"] = lol_piv_df.apply(
        lambda x: list(set(sge_df[sge_df.position == int(x.name)]["domain"]))[0],axis =1)
    non_domain_coding_positions = [str(i) for i in list(set(
        sge_df[(~pandas.isna(sge_df[aa_pos_col])) & (sge_df["domain"] == "No Domain")]["position"]))]
    splice_site_positions = [str(i) for i in list(set(
        sge_df[sge_df["consequence"].isin(["Splice region", "Canonical splice"])]["position"]))]
    lol_piv_df["Domain"] = lol_piv_df.apply(
        lambda x: 0 if x.domain == "No Domain" and x.name not in non_domain_coding_positions and x.name not in splice_site_positions else (
            1 if x.domain == "No Domain" and x.name in non_domain_coding_positions else (
                2 if x.domain == "No Domain" and x.name in splice_site_positions else (
                    3 if x.domain == "BRCT 1" else (
                        4 if x.domain == "BRCT 2" else (
                            5 if x.domain == "RING-type" else None))))), axis =1)

    # ClinVar
    """
    transitions = []
    for i in ["A", "T", "C", "G"]:
        for k in ["A", "T", "C", "G"]:
            if i != k:
                transitions.append(i + ">" + k)
    for transition in transitions:
        lol_piv_df[transition] = lol_piv_df.apply(
            lambda x: list(set(sge_df[(sge_df.position == int(x.name))&(sge_df.ref == transition.split(">")[0])&
                                 (sge_df.alt == transition.split(">")[1])]["clinvar"]))[0] if
            len(sge_df[(sge_df.position == int(x.name))&(sge_df.ref == transition.split(">")[0])&
                       (sge_df.alt == transition.split(">")[1])].index) != 0 else numpy.nan, axis = 1)

        lol_piv_df = lol_piv_df.replace(["Pathogenic", "Likely pathogenic"], 4)
        lol_piv_df = lol_piv_df.replace(["Benign", "Likely benign"], 3)
        lol_piv_df = lol_piv_df.replace("Conflicting interpretations of pathogenicity", 2)
        lol_piv_df = lol_piv_df.replace("Uncertain significance", 1)
        lol_piv_df = lol_piv_df.replace([numpy.nan, "absent"], 0)
    """
    # Base editors

    be_pam_dfs = {}
    for be, be_info in be_dict.items():
        positions = {pam: set(sge_df[sge_df[be + "_" + pam]][position_col])
                     for pam in pam_list}
        pam_dfs = {}
        for pam in pam_list:
            df = pandas.DataFrame(None, index=[lol_piv_df.index, lol_piv_df.ref],
                                  columns=[be+"_"+pam])
            df = df.reset_index()
            df[be+ "_" + pam] = df.apply(
                lambda x: 1 if int(x[position_col]) in positions[pam] else 0, axis=1)
            df = df.set_index([position_col, ref_col])
            pam_dfs[pam] =df

        pam_df = pandas.concat(pam_dfs.values(), axis = 1)
        pam_df[be] = pam_df.apply(lambda x: check_pam_positions(pam_df, x.name), axis=1)

        be_pam_dfs[be] = pam_df[[be]]

    be_lol_piv_df = lol_piv_df.reset_index()
    be_lol_piv_df = be_lol_piv_df.set_index([position_col, ref_col])
    for be, be_df in be_pam_dfs.items():
        be_lol_piv_df = pandas.concat([be_lol_piv_df, be_df], axis =1)
    be_lol_piv_df = be_lol_piv_df.reset_index()
    be_lol_piv_df = be_lol_piv_df.set_index([position_col])

    # Classification

    lol_p_piv_df = pivot_sge(df = sge_df, index_col=[position_col, ref_col],
                             column_col=[alt_col], value_col=nonfunctional_p_col,
                             aa_position_col=aa_pos_col, coding = False)

    lol_p_piv_df["Max"] = lol_p_piv_df.apply(
        lambda x: numpy.nanmax(list(lol_p_piv_df.loc[x.name])), axis=1)
    lol_p_piv_df["Min"] = lol_p_piv_df.apply(
        lambda x: numpy.nanmin(list(lol_p_piv_df.loc[x.name])), axis=1)

    lol_p_piv_df["Max_group"] = lol_p_piv_df.apply(lambda x: "LOF" if x.Max > 0.99 else (
        "INT" if 0.01 < x.Max < 0.99 else ("FUNC" if x.Max < 0.01 else None)), axis=1)

    lol_p_piv_df["Min_group"] = lol_p_piv_df.apply(lambda x: "LOF" if x.Min > 0.99 else (
        "INT" if 0.01 < x.Min < 0.99 else ("FUNC" if x.Min < 0.01 else None)), axis=1)

    lol_p_piv_df["Class"] = lol_p_piv_df.apply(
        lambda x: 1 if x["Max_group"] == x["Min_group"] and x["Max_group"] == "LOF" else (
            2 if x["Max_group"] != x["Min_group"] and x["Max_group"] == "LOF" else (
                3)), axis=1)

    lol_p_piv_df = lol_p_piv_df[["Max_group", "Min_group", "Class"]]
    lol_p_piv_df = lol_p_piv_df.reset_index()
    lol_p_piv_df[position_col] = lol_p_piv_df.apply(lambda x: str(x[position_col]), axis = 1)
    lol_p_piv_df = lol_p_piv_df.set_index([position_col, ref_col])
    lol_piv_df2 = lol_piv_df.reset_index()
    lol_piv_df2 = lol_piv_df2.set_index([position_col, ref_col])
    class_piv_df = pandas.concat([lol_piv_df2, lol_p_piv_df], axis = 1)
    class_piv_df = class_piv_df.reset_index()

    # Base editor functional scores
    lol_piv_df["C/T"] = lol_piv_df.apply(lambda x: x["T"] if x[ref_col] == "C" else numpy.nan, axis=1)
    lol_piv_df["A/G"] = lol_piv_df.apply(lambda x: x["G"] if x[ref_col] == "A" else numpy.nan, axis=1)
    lol_piv_df["G/C"] = lol_piv_df.apply(lambda x: x["C"] if x[ref_col] == "G" else numpy.nan, axis=1)

    be_lol_piv_df["norm_domain"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.Domain, list(set(be_lol_piv_df.Domain))), axis=1)
    class_piv_df["norm_class"] = class_piv_df.apply(
        lambda x: normalisation(x.Class, list(set(class_piv_df.Class))), axis=1)
    be_lol_piv_df["norm_cbe"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.CBE, list(set(be_lol_piv_df.CBE))), axis=1)
    be_lol_piv_df["norm_abe"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.ABE, list(set(be_lol_piv_df.ABE))), axis=1)
    be_lol_piv_df["norm_cgbe"] = be_lol_piv_df.apply(
        lambda x: normalisation(x.CGBE, list(set(be_lol_piv_df.CGBE))), axis=1)

    # Class colours
    cl_colours = [rgb_rgbtples((35, 139, 69)),rgb_rgbtples((214, 47, 30)),
                  rgb_rgbtples((255, 255, 255))]

    cl_nodes = [normalisation(1, list(set(class_piv_df["Class"]))),
                normalisation(2, list(set(class_piv_df["Class"]))),
                normalisation(3, list(set(class_piv_df["Class"])))]

    cl_cmap = mcolors.LinearSegmentedColormap.from_list("class_cmap", list(zip(cl_nodes, cl_colours)))

    # CBE colours
    cbe_colors = [rgb_rgbtples((255, 255, 255)), rgb_rgbtples((205, 160, 205)),
                  rgb_rgbtples((223, 33, 121))]

    cbe_nodes = [normalisation(0, list(set(be_lol_piv_df["CBE"]))),
                 normalisation(1, list(set(be_lol_piv_df["CBE"]))),
                 normalisation(2, list(set(be_lol_piv_df["CBE"])))]

    cbe_cmap = mcolors.LinearSegmentedColormap.from_list("cbe_cmap", list(zip(cbe_nodes, cbe_colors)))

    # ABE colours
    abe_colors = [rgb_rgbtples((255, 255, 255)), rgb_rgbtples((198, 199, 225)),
                  rgb_rgbtples((121, 110, 178))]

    abe_nodes = [normalisation(0, list(set(be_lol_piv_df["ABE"]))),
                 normalisation(1, list(set(be_lol_piv_df["ABE"]))),
                 normalisation(2, list(set(be_lol_piv_df["ABE"])))]

    abe_cmap = mcolors.LinearSegmentedColormap.from_list("abe_cmap", list(zip(abe_nodes, abe_colors)))

    # CGBE colours
    cgbe_colors = [rgb_rgbtples((255, 255, 255)), rgb_rgbtples((171, 208, 230)),
                   rgb_rgbtples((55, 135, 192))]

    cgbe_nodes = [normalisation(0, list(set(be_lol_piv_df["CGBE"]))),
                  normalisation(1, list(set(be_lol_piv_df["CGBE"]))),
                  normalisation(2, list(set(be_lol_piv_df["CGBE"])))]

    cgbe_cmap = mcolors.LinearSegmentedColormap.from_list("cgbe_cmap", list(zip(cgbe_nodes, cgbe_colors)))

    # Lollipop Plot

    fig = plt.figure(figsize=(40, 1.5))
    # Subplots
    widths = [40]*5
    heights = [0.3, 0.3, 0.3,0.3,0.3]
    gs = fig.add_gridspec(5, 5, wspace=0.05, hspace=0.05,
                          width_ratios=widths, height_ratios=heights)
    
    ax = fig.add_subplot(gs[-1, :])
    ax.xaxis.set_ticks_position('none')
    ax.vlines(ymin=lol_piv_df['Min'],ymax=lol_piv_df['Max'],x=lol_piv_df.index,
              color='grey', alpha=0.3)

    ax.scatter(y=lol_piv_df['Max'],x=lol_piv_df.index,
               color='lightblue', alpha=.5, label='Max', s=3)
    ax.scatter(y=lol_piv_df['Min'], x=lol_piv_df.index,
               color='red', alpha=.5, label='Min', s=3)
    
    ax.scatter(y=lol_piv_df['C/T'], x=lol_piv_df.index,
               color='deeppink', alpha=.5, label='C>T', s=4)
    ax.scatter(y=lol_piv_df['A/G'], x=lol_piv_df.index,
               color='blueviolet', alpha=.5, label='A>G', s=4)
    ax.scatter(y=lol_piv_df['G/C'], x=lol_piv_df.index,
               color='royalblue', alpha=.5, label='G>C', s=4)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(1.5)
        tick.label.set_rotation(90)
    
    dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                       rgb_rgbtples((230, 173, 18)),
                                       rgb_rgbtples((240, 101, 16)),
                                       rgb_rgbtples((196, 56, 90)),
                                       rgb_rgbtples((99, 42, 173)),
                                       rgb_rgbtples((43, 38, 136))])

    ax_d = fig.add_subplot(gs[0, :])
    ax_d.yaxis.set_ticks_position('none')
    ax_d.xaxis.set_ticks_position('none')
    ax_d.axes.get_xaxis().set_visible(False)
    ax_d.axes.get_yaxis().set_visible(False)
    ax_d.pcolor(be_lol_piv_df[["norm_domain"]].T,linewidths=0, cmap=dom_cmap)

    ax_c = fig.add_subplot(gs[1, :], sharex=ax_d)
    ax_c.yaxis.set_ticks_position('none')
    ax_c.xaxis.set_ticks_position('none')
    ax_c.axes.get_xaxis().set_visible(False)
    ax_c.axes.get_yaxis().set_visible(False)
    ax_c.pcolor(class_piv_df[["norm_class"]].T,linewidths=0, cmap=cl_cmap)

    ax_cb = fig.add_subplot(gs[2, :], sharex=ax_d)
    ax_cb.yaxis.set_ticks_position('none')
    ax_cb.xaxis.set_ticks_position('none')
    ax_cb.axes.get_xaxis().set_visible(False)
    ax_cb.axes.get_yaxis().set_visible(False)
    ax_cb.pcolor(be_lol_piv_df[["norm_cbe"]].T, linewidths=0, cmap=cbe_cmap)

    ax_ab = fig.add_subplot(gs[3, :], sharex=ax_d)
    ax_ab.yaxis.set_ticks_position('none')
    ax_ab.xaxis.set_ticks_position('none')
    ax_ab.axes.get_xaxis().set_visible(False)
    ax_ab.axes.get_yaxis().set_visible(False)
    ax_ab.pcolor(be_lol_piv_df[["norm_abe"]].T, linewidths=0, cmap=abe_cmap)

    ax_cgb = fig.add_subplot(gs[4, :], sharex=ax_d)
    ax_cgb.yaxis.set_ticks_position('none')
    ax_cgb.xaxis.set_ticks_position('none')
    ax_cgb.axes.get_xaxis().set_visible(False)
    ax_cgb.axes.get_yaxis().set_visible(False)
    ax_cgb.pcolor(be_lol_piv_df[["norm_cgbe"]].T, linewidths=0, cmap=cgbe_cmap)

    plt.savefig("/Users/cd7/Desktop/last/Figure1.png")
    plt.savefig("/Users/cd7/Desktop/last/Figure1.pdf")

    """
    fig = plt.figure(figsize=(16, 20))

    # Subplots
    widths = [20]
    widths.extend([0]*27)
    heights = [0.1,0.1,0.1,0.1,0.1,3.0,0.5]*4

    gs = fig.add_gridspec(28, 28, wspace=0.05, hspace=0.05,
                          width_ratios=widths, height_ratios=heights)
    
    
    position_ranges = list(range(0, len(lol_piv_df), int(len(lol_piv_df)/4)+1))
    position_ranges.append(len(lol_piv_df))
    p=0
    for i in range(5, 27, 7):

        ax = fig.add_subplot(gs[i, :])
        ax.xaxis.set_ticks_position('none')
        ax.vlines(ymin=lol_piv_df['Min'][position_ranges[p]:position_ranges[p+1]],
                  ymax=lol_piv_df['Max'][position_ranges[p]:position_ranges[p+1]],
                  x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]],
                  color='grey', alpha=0.3)

        ax.scatter(y=lol_piv_df['Max'][position_ranges[p]:position_ranges[p+1]],
                   x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]],
                   color='grey', alpha=.5, label='Max', s=4)
        ax.scatter(y=lol_piv_df['Min'][position_ranges[p]:position_ranges[p+1]],
                   x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]],
                   color='grey', alpha=.5, label='Min', s=4)
        
        marker_map = {0: 's', 4: "*", 2: "o", 1: "o", 3: "o"}
        for be_name, be_d in be_dict.items():
            m = [marker_map[n]
                 for n in lol_piv_df[be_d["ref"]+">"+be_d["alt"]]
                 [position_ranges[p]:position_ranges[p+1]]]
            x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]]
            y=lol_piv_df[be_d["ref"]+">"+be_d["alt"]][position_ranges[p]:position_ranges[p+1]]

            for _m, _x, _y in zip(m, x, y):
                ax.scatter(_x, _y, marker=_m, c=be_d["colour"], alpha=.5,
                           label=be_d["ref"]+">"+be_d["alt"], s=6)
        

        ax.scatter(y=lol_piv_df['C/T'][position_ranges[p]:position_ranges[p+1]],
                   x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]],
                   color='deeppink', alpha=.5, label='C>T', s=6)
        ax.scatter(y=lol_piv_df['A/G'][position_ranges[p]:position_ranges[p+1]],
                   x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]],
                   color='blueviolet', alpha=.5, label='A>G', s=6)
        ax.scatter(y=lol_piv_df['G/C'][position_ranges[p]:position_ranges[p+1]],
                   x=lol_piv_df.index[position_ranges[p]:position_ranges[p+1]],
                   color='royalblue', alpha=.5, label='G>C', s=6)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(3)
            tick.label.set_rotation(90)

        d, c, cbe, abe, cgbe = 0,0,0,0,0
        if p == 0: d, c, cbe, abe, cgbe = 0,1,2,3,4
        elif p == 1: d, c, cbe, abe, cgbe = 7, 8, 9, 10, 11
        elif p ==2: d, c, cbe, abe, cgbe = 14, 15, 16, 17, 18
        elif p ==3: d, c, cbe, abe, cgbe = 21, 22, 23, 24, 25

        if p ==0:
            dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                               rgb_rgbtples((240, 101, 16)),
                                               rgb_rgbtples((99, 42, 173))])

        elif p ==1:
            dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                               rgb_rgbtples((230, 173, 18)),
                                               rgb_rgbtples((240, 101, 16)),
                                               rgb_rgbtples((196, 56, 90)),
                                               rgb_rgbtples((99, 42, 173))])

        elif p ==2:
            dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                               rgb_rgbtples((230, 173, 18)),
                                               rgb_rgbtples((240, 101, 16)),
                                               rgb_rgbtples((196, 56, 90))])

        elif p ==3:
            dom_cmap = mcolors.ListedColormap([rgb_rgbtples((255, 255, 255)),
                                               rgb_rgbtples((230, 173, 18)),
                                               rgb_rgbtples((240, 101, 16)),
                                               rgb_rgbtples((43, 38, 136))])

        ax_d = fig.add_subplot(gs[d, :], sharex=ax)
        ax_d.yaxis.set_ticks_position('none')
        ax_d.xaxis.set_ticks_position('none')
        ax_d.axes.get_xaxis().set_visible(False)
        ax_d.axes.get_yaxis().set_visible(False)
        ax_d.pcolor(be_lol_piv_df[["norm_domain"]][position_ranges[p]:position_ranges[p+1]].T,
                    linewidths=0, cmap=dom_cmap)

        ax_c = fig.add_subplot(gs[c, :], sharex=ax)
        ax_c.yaxis.set_ticks_position('none')
        ax_c.xaxis.set_ticks_position('none')
        ax_c.axes.get_xaxis().set_visible(False)
        ax_c.axes.get_yaxis().set_visible(False)
        ax_c.pcolor(class_piv_df[["norm_class"]][position_ranges[p]:position_ranges[p+1]].T,
                    linewidths=0, cmap=cl_cmap)

        ax_cb = fig.add_subplot(gs[cbe, :], sharex=ax)
        ax_cb.yaxis.set_ticks_position('none')
        ax_cb.xaxis.set_ticks_position('none')
        ax_cb.axes.get_xaxis().set_visible(False)
        ax_cb.axes.get_yaxis().set_visible(False)
        ax_cb.pcolor(be_lol_piv_df[["norm_cbe"]][position_ranges[p]:position_ranges[p + 1]].T,
                     linewidths=0, cmap=cbe_cmap)

        ax_ab = fig.add_subplot(gs[abe, :], sharex=ax)
        ax_ab.yaxis.set_ticks_position('none')
        ax_ab.xaxis.set_ticks_position('none')
        ax_ab.axes.get_xaxis().set_visible(False)
        ax_ab.axes.get_yaxis().set_visible(False)
        ax_ab.pcolor(be_lol_piv_df[["norm_abe"]][position_ranges[p]:position_ranges[p + 1]].T,
                     linewidths=0, cmap=abe_cmap)

        ax_cgb = fig.add_subplot(gs[cgbe, :], sharex=ax)
        ax_cgb.yaxis.set_ticks_position('none')
        ax_cgb.xaxis.set_ticks_position('none')
        ax_cgb.axes.get_xaxis().set_visible(False)
        ax_cgb.axes.get_yaxis().set_visible(False)
        ax_cgb.pcolor(be_lol_piv_df[["norm_cgbe"]][position_ranges[p]:position_ranges[p + 1]].T,
                      linewidths=0, cmap=cgbe_cmap)
        p += 1
    
    plt.tight_layout()
    plt.savefig("/Users/cd7/Desktop/last/Figure2.png")
    plt.savefig("/Users/cd7/Desktop/last/Figure2.pdf")
    """
    return True
    


def main():
    be_dict={"ABE":{"ref":"A","alt":"G","colour": "purple"},
             "CBE":{"ref":"C","alt": "T","colour":"pink"},
             "CGBE":{"ref":"C","alt":"G","colour":"blue"}}
    pam_list = ["NGG", "NG"]
    uniprot_object = Uniprot(uniprotid="P38398")
    uniprot_object.extract_uniprot_info()

    sge_df = sge_file(file_path=os.getcwd() + "/venv/data/Findlay_BRCA1.csv",
                      wanted_columns=["position (hg19)", "reference", "alt", "aa_pos", "aa_ref",
                                      "aa_alt","function.score.mean", "func.class",
                                      "p.nonfunctional", "clinvar_simple", "consequence"],
                      new_columns=["position", "ref", "alt", "aa_pos", "aa_ref" ,"aa_alt",
                                   "functional score", "functional class", "nonfunctional p",
                                   "clinvar", "consequence"])

    if uniprot_seq_checking(uniprot_obj=uniprot_object, sge_df=sge_df,
                            aa_pos_column="aa_pos", aa_ref_column="aa_ref"):

        dom_sge_df = add_domain_info(sge_df=sge_df, aa_pos_column="aa_pos",
                                     uniprot_obj=uniprot_object)

        pos_dom_sge_df = dom_sge_df.copy()

        for be, be_info in be_dict.items():
            for pam in pam_list:
                f = os.getcwd() + "/venv/output/PAM_Extractor/v4/"+be+"_"+pam+"_hg19_edit_df.csv"
                pos_dom_sge_df= add_be_annotation(sge_df=pos_dom_sge_df,
                                                  ref_columns_name="ref", alt_column_name="alt",
                                                  position_column_name="position",
                                                  be_file=f, be=be, pam=pam,
                                                  ref=be_info["ref"], alt=be_info["alt"])

        figure_heatmap(sge_df= pos_dom_sge_df, aa_position_col="aa_pos", aa_ref_col="aa_ref",
                       aa_alt_col= "aa_alt", func_score_col = "functional score",
                       nonfunctional_p_col = "nonfunctional p", pam_list=pam_list,
                       be_dict=be_dict, uniprot_object = uniprot_object,
                       func_class_col="functional class")

        figure_lollipop(sge_df=pos_dom_sge_df, position_col="position", ref_col="ref", alt_col="alt",
                        aa_pos_col="aa_pos",functional_score_col="functional score",
                        nonfunctional_p_col="nonfunctional p", pam_list=pam_list, be_dict=be_dict)


###########################################################################################
# Quantification (needed to be systematical)

class_1_coding = class_piv_df[(class_piv_df.Class == 1) & (class_piv_df.Domain.isin([1,3,4,5]))]
class_1_coding = class_1_coding.reset_index()
class_1_noncoding = class_piv_df[(class_piv_df.Class == 1) & (class_piv_df.Domain.isin([0]))]
class_1_noncoding = class_1_noncoding.reset_index()
class_1_sp = class_piv_df[(class_piv_df.Class == 1) & (class_piv_df.Domain == 2)]
class_1_sp = class_1_sp.reset_index()

class_2_coding = class_piv_df[(class_piv_df.Class == 2) & (class_piv_df.Domain.isin([1,3,4,5]))]
class_2_coding = class_2_coding.reset_index()
class_2_noncoding = class_piv_df[(class_piv_df.Class == 2) & (class_piv_df.Domain.isin([0]))]
class_2_noncoding = class_2_noncoding.reset_index()
class_2_sp = class_piv_df[(class_piv_df.Class == 2) & (class_piv_df.Domain == 2)]
class_2_sp = class_2_sp.reset_index()

class_3_coding = class_piv_df[(class_piv_df.Class == 3) & (class_piv_df.Domain.isin([1,3,4,5]))]
class_3_coding = class_3_coding.reset_index()
class_3_noncoding = class_piv_df[(class_piv_df.Class == 3) & (class_piv_df.Domain.isin([0]))]
class_3_noncoding = class_3_noncoding.reset_index()
class_3_sp = class_piv_df[(class_piv_df.Class == 3) & (class_piv_df.Domain == 2)]
class_3_sp = class_3_sp.reset_index()


melt_class_1_coding = class_1_coding.melt(id_vars = ["position", "ref", "Class"],
                                          value_vars = ["A", "C", "T", "G"])
melt_class_1_noncoding = class_1_noncoding.melt(id_vars = ["position", "ref", "Class"],
                                                value_vars = ["A", "C", "T", "G"])
melt_class_1_sp = class_1_sp.melt(id_vars = ["position", "ref", "Class"],
                                  value_vars = ["A", "C", "T", "G"])

melt_class_2_coding = class_2_coding.melt(id_vars = ["position", "ref", "Class"],
                                          value_vars = ["A", "C", "T", "G"])
melt_class_2_noncoding = class_2_noncoding.melt(id_vars = ["position", "ref", "Class"],
                                                value_vars = ["A", "C", "T", "G"])
melt_class_2_sp = class_2_sp.melt(id_vars = ["position", "ref", "Class"],
                                  value_vars = ["A", "C", "T", "G"])


melt_class_3_coding = class_3_coding.melt(id_vars = ["position", "ref", "Class"],
                                          value_vars = ["A", "C", "T", "G"])
melt_class_3_noncoding = class_3_noncoding.melt(id_vars = ["position", "ref", "Class"],
                                                value_vars = ["A", "C", "T", "G"])
melt_class_3_sp = class_3_sp.melt(id_vars = ["position", "ref", "Class"],
                                  value_vars = ["A", "C", "T", "G"])

class_piv_df = class_piv_df.reset_index()
melt_all_class = class_piv_df.melt(id_vars = ["position", "ref", "Class","Domain"],
                                   value_vars = ["A", "C", "T", "G"])

"""
Check if there is a functional score difference between non coding and (coding +splice sites)
P:pvalue=2.1294825309720274e-08
"""

stats.mannwhitneyu(list(melt_all_class[melt_all_class.Domain.isin([2,3,4,5])]["value"]),
                   list(melt_all_class[melt_all_class.Domain.isin([0])]["value"]))


"""
        coding noncoding
    1
    2
    3
"""
"""
Check which class in which domain
Observed
    1   74  57  0
    2   279 41  9
    3   604 186 183
Expected
    1   87.48569435,  25.96231682,  17.55198883
    2   219.71598046,  65.20307048,  44.08094906
    3   649.79832519, 192.8346127 , 130.36706211
P: 2.922103298885345e-32
"""

class123_g, class123_p, class123_dof, class123_expctd = chi2_contingency(
    [[74, 57, 0],
     [279, 41, 9],
     [604, 186, 183]],
    lambda_="log-likelihood")

"""
Check if Class I or Class II in coding/Splice
Observed
        c   s
    2   74  57
    1   279 41
P: 5.727617754807361e-12
OR: 5.241595253790376
"""
class12_OR, class12_p = stats.fisher_exact([[279, 41],[74, 57]])

"""
Check which nucletide is prefeable for splice Class I positions 
Chi square
Observed
                        A   T   C   G
    Splice Class I      16  15  26  0
    No splice Class I   29  10  25  10
Expected
                        A   T   C   G
    Splice Class I      20  11  22  4      
    No splice Class I   25  14  29  6
P: 0.005126943506264445
"""

class_1_non_sp = class_piv_df[(class_piv_df.Class == 1) & (class_piv_df.Domain != 2)]
class_1_non_sp = class_1_non_sp.reset_index()
melt_class_1_non_sp = class_1_non_sp.melt(id_vars = ["position", "ref", "Class"],
                                          value_vars = ["A", "C", "T", "G"])


_, splice_C1_nucleotide_p, _, splice_C1_nucleotide_exp = chi2_contingency([
    [len(set(melt_class_1_sp[melt_class_1_sp.ref == "A"]["position"])),
     len(set(melt_class_1_sp[melt_class_1_sp.ref == "T"]["position"])),
     len(set(melt_class_1_sp[melt_class_1_sp.ref == "C"]["position"])),
     len(set(melt_class_1_sp[melt_class_1_sp.ref == "G"]["position"]))],
    [len(set(melt_class_1_non_sp[melt_class_1_non_sp.ref == "A"]["position"])),
     len(set(melt_class_1_non_sp[melt_class_1_non_sp.ref == "T"]["position"])),
     len(set(melt_class_1_non_sp[melt_class_1_non_sp.ref == "C"]["position"])),
     len(set(melt_class_1_non_sp[melt_class_1_non_sp.ref == "G"]["position"]))]])


stats.kruskal(list(melt_class_1_sp[melt_class_1_sp.ref == "A"]["value"]),
              list(melt_class_1_sp[melt_class_1_sp.ref == "T"]["value"]),
              list(melt_class_1_sp[melt_class_1_sp.ref == "C"]["value"]),
              list(melt_class_1_sp[melt_class_1_sp.ref == "G"]["value"]))

"""
Check which nucletide is prefeable for coding Class II positions
Chi square
Observed
                        A   T   C   G
    coding Class II     82  70  78  49
    No coding Class II  20  15  6   9
Expected
                        A   T   C   G
    coding Class II     86  72  71  49
    No coding Class II  16  13  13  9
P:0.10370250326840082
"""
class_2_non_coding = class_piv_df[(class_piv_df.Class == 2) & (class_piv_df.Domain.isin([0,2]))]
class_2_non_coding = class_2_non_coding.reset_index()
melt_class_2_non_coding = class_2_non_coding.melt(id_vars = ["position", "ref", "Class"],
                                                  value_vars = ["A", "C", "T", "G"])

_, coding_C2_nucleotide_p, _, coding_C2_nucleotide_exp = chi2_contingency([
    [len(set(melt_class_2_coding[melt_class_2_coding.ref == "A"]["position"])),
     len(set(melt_class_2_coding[melt_class_2_coding.ref == "T"]["position"])),
     len(set(melt_class_2_coding[melt_class_2_coding.ref == "C"]["position"])),
     len(set(melt_class_2_coding[melt_class_2_coding.ref == "G"]["position"]))],
    [len(set(melt_class_2_non_coding[melt_class_2_non_coding.ref == "A"]["position"])),
     len(set(melt_class_2_non_coding[melt_class_2_non_coding.ref == "T"]["position"])),
     len(set(melt_class_2_non_coding[melt_class_2_non_coding.ref == "C"]["position"])),
     len(set(melt_class_2_non_coding[melt_class_2_non_coding.ref == "G"]["position"]))]])

"""
Check if there is a difference between the median of the class II coding nucleotides
P:2.537492286176995e-05
"""
stats.kruskal(list(melt_class_2_coding[melt_class_2_coding.ref == "A"]["value"]),
              list(melt_class_2_coding[melt_class_2_coding.ref == "T"]["value"]),
              list(melt_class_2_coding[melt_class_2_coding.ref == "C"]["value"]),
              list(melt_class_2_coding[melt_class_2_coding.ref == "G"]["value"]))

"""
Posthoc --> dunn test
          1         2         3         4
1  1.000000  0.000182  0.991159  0.277919
2  0.000182  1.000000  0.000207  0.000015
3  0.991159  0.000207  1.000000  0.286816
4  0.277919  0.000015  0.286816  1.000000
T (higher value) is significantly different then others 
"""
scikit_posthocs.posthoc_dunn([list(melt_class_2_coding[melt_class_2_coding.ref == "A"]["value"]),
                              list(melt_class_2_coding[melt_class_2_coding.ref == "T"]["value"]),
                              list(melt_class_2_coding[melt_class_2_coding.ref == "C"]["value"]),
                              list(melt_class_2_coding[melt_class_2_coding.ref == "G"]["value"])])

class_piv_df = class_piv_df.set_index(["position"])
merged_df = be_lol_piv_df.copy()
merged_df["Class"] = merged_df.apply(lambda x: class_piv_df.loc[x.name]["Class"],axis =1)

# Counts and Percentages

cbe= merged_df[(merged_df.CBE != 0) & (merged_df.ABE == 0)]
abe= merged_df[(merged_df.ABE != 0) & (merged_df.CBE == 0) & (merged_df.CGBE == 0)]
cgbe= merged_df[(merged_df.CGBE != 0) & (merged_df.ABE == 0)]

cbe_ong = set(cbe[(cbe.CBE == 1)].index)
cbe_ngg = set(cbe[(cbe.CBE == 2)].index)
cbe_ngg_brct1_n = cbe[(cbe.CBE == 2) & (cbe.domain == "BRCT 1")]
cbe_ngg_brct2_n = cbe[(cbe.CBE == 2) & (cbe.domain == "BRCT 2")]
cbe_ngg_ring_n = cbe[(cbe.CBE == 2) & (cbe.domain == "RING-type")]
cbe_ong_brct1_n = cbe[(cbe.CBE == 1) & (cbe.domain == "BRCT 1")]
cbe_ong_brct2_n = cbe[(cbe.CBE == 1) & (cbe.domain == "BRCT 2")]
cbe_ong_ring_n = cbe[(cbe.CBE == 1) & (cbe.domain == "RING-type")]
cbe_ng_brct1_n = cbe[(cbe.CBE != 0) & (cbe.domain == "BRCT 1")]
cbe_ng_brct2_n = cbe[(cbe.CBE != 0) & (cbe.domain == "BRCT 2")]
cbe_ng_ring_n = cbe[(cbe.CBE != 0) & (cbe.domain == "RING-type")]

abe_ong = set(abe[(abe.ABE == 1)].index)
abe_ngg = set(abe[(abe.ABE == 2)].index)
abe_ngg_brct1_n = abe[(abe.ABE == 2) & (abe.domain == "BRCT 1")]
abe_ngg_brct2_n = abe[(abe.ABE == 2) & (abe.domain == "BRCT 2")]
abe_ngg_ring_n = abe[(abe.ABE == 2) & (abe.domain == "RING-type")]
abe_ong_brct1_n = abe[(abe.ABE == 1) & (abe.domain == "BRCT 1")]
abe_ong_brct2_n = abe[(abe.ABE == 1) & (abe.domain == "BRCT 2")]
abe_ong_ring_n = abe[(abe.ABE == 1) & (abe.domain == "RING-type")]
abe_ng_brct1_n = abe[(abe.ABE != 0) & (abe.domain == "BRCT 1")]
abe_ng_brct2_n = abe[(abe.ABE != 0) & (abe.domain == "BRCT 2")]
abe_ng_ring_n = abe[(abe.ABE != 0) & (abe.domain == "RING-type")]


cgbe_ong = set(cgbe[(cgbe.CGBE == 1)].index)
cgbe_ngg = set(cgbe[(cgbe.CGBE == 2)].index)
cgbe_ngg_brct1_n = cgbe[(cgbe.CGBE == 2) & (cgbe.domain == "BRCT 1")]
cgbe_ngg_brct2_n = cgbe[(cgbe.CGBE == 2) & (cgbe.domain == "BRCT 2")]
cgbe_ngg_ring_n = cgbe[(cgbe.CGBE == 2) & (cgbe.domain == "RING-type")]
cgbe_ong_brct1_n = cgbe[(cgbe.CGBE == 1) & (cgbe.domain == "BRCT 1")]
cgbe_ong_brct2_n = cgbe[(cgbe.CGBE == 1) & (cgbe.domain == "BRCT 2")]
cgbe_ong_ring_n = cgbe[(cgbe.CGBE == 1) & (cgbe.domain == "RING-type")]
cgbe_ng_brct1_n = cgbe[(cgbe.CGBE != 0) & (cgbe.domain == "BRCT 1")]
cgbe_ng_brct2_n = cgbe[(cgbe.CGBE != 0) & (cgbe.domain == "BRCT 2")]
cgbe_ng_ring_n = cgbe[(cgbe.CGBE != 0) & (cgbe.domain == "RING-type")]


c1_cbe_ong = set(cbe[(cbe.Class == 1) & (cbe.CBE == 1)].index)
c1_cbe_ngg = set(cbe[(cbe.Class == 1) & (cbe.CBE == 2)].index)
c1_cbe =set(cbe[(cbe.Class == 1) & (cbe.CBE != 0)].index)

c1_abe_ong = set(abe[(abe.Class == 1) & (abe.ABE == 1)].index)
c1_abe_ngg = set(abe[(abe.Class == 1) & (abe.ABE == 2)].index)
c1_abe =set(abe[(abe.Class == 1) & (abe.ABE != 0)].index)

c1_cgbe_ong = set(cgbe[(cgbe.Class == 1) & (cgbe.CGBE == 1)].index)
c1_cgbe_ngg = set(cgbe[(cgbe.Class == 1) & (cgbe.CGBE == 2)].index)
c1_cgbe =set(cgbe[(cgbe.Class == 1) & (cgbe.CGBE != 0)].index)

c2_cbe_ong = set(cbe[(cbe.Class == 2) & (cbe.CBE == 1)].index)
c2_cbe_ngg = set(cbe[(cbe.Class == 2) & (cbe.CBE == 2)].index)
c2_cbe =set(cbe[(cbe.Class == 2) & (cbe.CBE != 0)].index)

c2_abe_ong = set(abe[(abe.Class == 2) & (abe.ABE == 1)].index)
c2_abe_ngg = set(abe[(abe.Class == 2) & (abe.ABE == 2)].index)
c2_abe =set(abe[(abe.Class == 2) & (abe.ABE != 0)].index)

c2_cgbe_ong = set(cgbe[(cgbe.Class == 2) & (cgbe.CGBE == 1)].index)
c2_cgbe_ngg = set(cgbe[(cgbe.Class == 2) & (cgbe.CGBE == 2)].index)
c2_cgbe =set(cgbe[(cgbe.Class == 2) & (cgbe.CGBE != 0)].index)

c3_cbe_ong = set(cbe[(cbe.Class == 3) & (cbe.CBE == 1)].index)
c3_cbe_ngg = set(cbe[(cbe.Class == 3) & (cbe.CBE == 2)].index)
c3_cbe =set(cbe[(cbe.Class == 3) & (cbe.CBE != 0)].index)

c3_abe_ong = set(abe[(abe.Class == 3) & (abe.ABE == 1)].index)
c3_abe_ngg = set(abe[(abe.Class == 3) & (abe.ABE == 2)].index)
c3_abe =set(abe[(abe.Class == 3) & (abe.ABE != 0)].index)

c3_cgbe_ong = set(cgbe[(cgbe.Class == 3) & (cgbe.CGBE == 1)].index)
c3_cgbe_ngg = set(cgbe[(cgbe.Class == 3) & (cgbe.CGBE == 2)].index)
c3_cgbe =set(cgbe[(cgbe.Class == 3) & (cgbe.CGBE != 0)].index)


bes= merged_df[(merged_df.CBE != 0)|(merged_df.ABE != 0)|(merged_df.CGBE != 0)]
bes_pos = set(bes.index)
bes_ngg= merged_df[(merged_df.CBE == 2)|(merged_df.ABE == 2)|(merged_df.CGBE == 2)]
bes_ong= merged_df[(merged_df.CBE == 1)|(merged_df.ABE == 1)|(merged_df.CGBE == 1)]
bes_ngg_pos = set(bes_ngg.index)
others = merged_df[(merged_df.CGBE == 0)&(merged_df.CBE == 0)&(merged_df.CGBE == 0)]
others_pos = set(others.index)

bes_C1 = set(bes[bes.Class == 1].index)
bes_ngg_C1 = set(bes_ngg[bes_ngg.Class == 1].index)

bes_C2 = set(bes[bes.Class == 2].index)
bes_ngg_C2 = set(bes_ngg[bes_ngg.Class == 2].index)

bes_C3 = set(bes[bes.Class == 3].index)
bes_ngg_C3 = set(bes_ngg[bes_ngg.Class == 3].index)

# Nucleotide level clinical relevance

merged_df["clinic_relevance"] = merged_df.apply(
    lambda x: 1 if x["A>T"] == 4 or x["A>C"] == 4 or x["A>G"] == 4 or x["T>A"] == 4 or
                   x["T>C"] == 4 or x["T>G"] == 4 or x["G>T"] == 4 or x["G>C"] == 4 or
                   x["C>A"] == 4 or x["C>T"] == 4 or x["C>G"] == 4 or x["G>A"] == 4 else 0,axis =1)

pat_merged_df = merged_df[merged_df["clinic_relevance"] == 1]
pat_bes= pat_merged_df[(pat_merged_df.CBE != 0)|(pat_merged_df.ABE != 0)|(pat_merged_df.CGBE != 0)]

pat_cbe= pat_merged_df[(pat_merged_df.CBE != 0)]
pat_cbe_ngg= pat_merged_df[(pat_merged_df.CBE == 2)]
pat_abe= pat_merged_df[pat_merged_df.ABE != 0]
pat_abe_ngg= pat_merged_df[(pat_merged_df.ABE == 2)]
pat_cgbe= pat_merged_df[pat_merged_df.CGBE != 0]
pat_cgbe_ngg= pat_merged_df[(pat_merged_df.CGBE == 2)]

# Amino acid level analysis

aa_pos_dom_sge_df = pos_dom_sge_df[~pandas.isna(pos_dom_sge_df.aa_pos)]

aa_pos_dom_sge_df["ABE"] = aa_pos_dom_sge_df.apply(
    lambda x: 2 if x.ABE_NG and x.ABE_NGG else (
        1 if x.ABE_NG else 0), axis =1)

aa_pos_dom_sge_df["CBE"] = aa_pos_dom_sge_df.apply(
    lambda x: 2 if x.CBE_NG and x.CBE_NGG else (
        1 if x.CBE_NG else 0), axis =1)

aa_pos_dom_sge_df["CGBE"] = aa_pos_dom_sge_df.apply(
    lambda x: 2 if x.CGBE_NG and x.CGBE_NGG else (
        1 if x.CGBE_NG else 0), axis =1)

aa_pos_dom_sge_df = aa_pos_dom_sge_df.set_index(["aa_pos"])

aa_bes= aa_pos_dom_sge_df[(aa_pos_dom_sge_df.CBE != 0)|(aa_pos_dom_sge_df.ABE != 0)|(aa_pos_dom_sge_df.CGBE != 0)]
aa_bes_pos = set(aa_bes.index)
aa_bes_ngg= aa_pos_dom_sge_df[(aa_pos_dom_sge_df.CBE == 2)|(aa_pos_dom_sge_df.ABE == 2)|(aa_pos_dom_sge_df.CGBE == 2)]
aa_bes_ngg_pos = set(aa_bes_ngg.index)

aa_cbe= len(set(aa_bes[(aa_bes.CBE != 0)].index))
aa_cbe_ngg= len(set(aa_bes[(aa_bes.CBE == 2)].index))

aa_abe= aa_bes[(aa_bes.ABE != 0)]
aa_abe_ngg= aa_bes[(aa_bes.ABE == 2)]
aa_cgbe= aa_bes[(aa_bes.CGBE != 0)]
aa_cgbe_ngg= aa_bes[(aa_bes.CGBE == 2)]

aa_brct1 = aa_pos_dom_sge_df[aa_pos_dom_sge_df.domain == "BRCT 1"]
aa_brct1_bes= aa_brct1[(aa_brct1.CBE != 0)|(aa_brct1.ABE != 0)|(aa_brct1.CGBE != 0)]
aa_brct1_cbe = aa_brct1[aa_brct1.CBE != 0]
aa_brct1_abe = aa_brct1[aa_brct1.ABE != 0]
aa_brct1_cgbe = aa_brct1[aa_brct1.CGBE != 0]

aa_brct2 = aa_pos_dom_sge_df[aa_pos_dom_sge_df.domain == "BRCT 2"]
aa_brct2_bes= aa_brct2[(aa_brct2.CBE != 0)|(aa_brct2.ABE != 0)|(aa_brct2.CGBE != 0)]
aa_brct2_cbe = aa_brct2[aa_brct2.CBE != 0]
aa_brct2_abe = aa_brct2[aa_brct2.ABE != 0]
aa_brct2_cgbe = aa_brct2[aa_brct2.CGBE != 0]

aa_zinc = aa_pos_dom_sge_df[aa_pos_dom_sge_df.domain == "RING-type"]
aa_zinc_bes= aa_zinc[(aa_zinc.CBE != 0)|(aa_zinc.ABE != 0)|(aa_zinc.CGBE != 0)]
aa_zinc_cbe = aa_zinc[aa_zinc.CBE != 0]
aa_zinc_abe = aa_zinc[aa_zinc.ABE != 0]
aa_zinc_cgbe = aa_zinc[aa_zinc.CGBE != 0]

aa_pos_dom_sge_df["clinical_relevance"] = aa_pos_dom_sge_df.apply(
    lambda x: 1 if x.clinvar in ["Pathogenic", "Likely Pathogenic"] else 0,axis=1)

pat_aa_pos = aa_pos_dom_sge_df[aa_pos_dom_sge_df["clinical_relevance"] == 1]
pat_aa_pos_bes= pat_aa_pos[(pat_aa_pos.CBE != 0)|(pat_aa_pos.ABE != 0)|(pat_aa_pos.CGBE != 0)]

pat_aa_cbe= pat_aa_pos[(pat_aa_pos.CBE != 0)]
pat_aa_cbe_ngg= pat_aa_pos[(pat_aa_pos.CBE == 2)]
pat_aa_abe= pat_aa_pos[pat_aa_pos.ABE != 0]
pat_aa_abe_ngg= pat_aa_pos[(pat_aa_pos.ABE == 2)]
pat_aa_cgbe= pat_aa_pos[pat_aa_pos.CGBE != 0]
pat_aa_cgbe_ngg= pat_aa_pos[(pat_aa_pos.CGBE == 2)]

# heatmap run for x and y functional class percentages

lof_aa_pos = y[y.LOF_Percentages > 0.8]
lof_aa_pos_pos = [i[0] for i in lof_aa_pos.index]
aa_pos_dom_sge_df = aa_pos_dom_sge_df.set_index(["aa_pos"])
lof_aa_sge = aa_pos_dom_sge_df[aa_pos_dom_sge_df.index.isin(lof_aa_pos_pos)]
lof_aa_bes= lof_aa_sge[(lof_aa_sge.CBE != 0)|(lof_aa_sge.ABE != 0)|(lof_aa_sge.CGBE != 0)]

lof_aa_cbe= lof_aa_bes[(lof_aa_bes.CBE != 0)]
lof_aa_cbe_ngg= lof_aa_bes[(lof_aa_bes.CBE == 2)]
lof_aa_abe= lof_aa_bes[(lof_aa_bes.ABE != 0)]
lof_aa_abe_ngg= lof_aa_bes[(lof_aa_bes.ABE == 2)]
lof_aa_cgbe= lof_aa_bes[(lof_aa_bes.CGBE != 0)]
lof_aa_cgbe_ngg= lof_aa_bes[(lof_aa_bes.CGBE == 2)]


not_lof_aa_pos = y[y.FUNC_Percentages > 0.8]

###########################################################################################
# Figure 2

# A - All nucleotides and amino acid coverage

df1 = pandas.DataFrame([["CBE", "Nucleotide Level",
                         percentage(1326,56), percentage(1326,131), percentage(1326,187)],
                        ["CBE", "Amino Acid Level",
                         percentage(326,31), percentage(326,83), percentage(326,114)],
                        ["ABE", "Nucleotide Level",
                         percentage(1326,74), percentage(1326,153), percentage(1326,227)],
                        ["ABE", "Amino Acid Level",
                         percentage(326,35), percentage(326,71), percentage(326,106)],
                        ["CGBE", "Nucleotide Level",
                         percentage(1326,58), percentage(1326,137), percentage(1326,195)],
                        ["CGBE", "Amino Acid Level",
                         percentage(326,34), percentage(326,86), percentage(326,120)],
                        ["All BEs", "Nucleotide Level",
                         percentage(1326,133), percentage(1326,291), percentage(1326,424)],
                        ["All BEs", "Amino Acid Level",
                         percentage(326,56), percentage(326,132), percentage(326,188)]])

df1.columns = ["BE", "Level", "NGG", "Only NG", "NGG+NG"]
mdf1 = df1.melt(id_vars=["BE", "Level"], value_vars=["NGG", "Only NG", "NGG+NG"])
mdf1.columns = ["BE", "Level", "PAM", "Percentage"]

sns.catplot(x="BE", y = "Percentage", col="Level", data = mdf1, hue= "PAM", kind="bar",
            palette="coolwarm", legend=True, col_wrap=2)
plt.savefig("/Users/cd7/Desktop/last/Figure2A.png")

# C - Class Coverage at nucleotide level

df2 = pandas.DataFrame([["CBE", "Class-I",
                         percentage(131,10), percentage(131,17), percentage(131,27)],
                        ["CBE", "Class-II",
                         percentage(329,12), percentage(329,36), percentage(329,48)],
                        ["CBE", "Class-III",
                         percentage(866,34), percentage(866,78), percentage(866,112)],
                        ["ABE", "Class-I",
                         percentage(131,15), percentage(131,15), percentage(131,30)],
                        ["ABE", "Class-II",
                         percentage(329,13), percentage(329,43), percentage(329,56)],
                        ["ABE", "Class-III",
                         percentage(866,46), percentage(866,95), percentage(866,141)],
                        ["CGBE", "Class-I",
                         percentage(131,10), percentage(131,19), percentage(131,29)],
                        ["CGBE", "Class-II",
                         percentage(329,13), percentage(329,37), percentage(329,50)],
                        ["CGBE", "Class-III",
                         percentage(866,35), percentage(866,81), percentage(866,116)],
                        ["All BEs", "Class-I",
                         percentage(131,25), percentage(131,34), percentage(131,59)],
                        ["All BEs", "Class-II",
                         percentage(329,26), percentage(329,80), percentage(329,106)],
                        ["All BEs", "Class-III",
                         percentage(866,82), percentage(866,177), percentage(866,259)]])

df2.columns = ["BE", "Class", "NGG", "Only NG", "NGG+NG"]
mdf2 = df2.melt(id_vars=["BE", "Class"], value_vars=["NGG", "Only NG", "NGG+NG"])
mdf2.columns = ["BE", "Class", "PAM", "Percentage"]

sns.catplot(x="BE", y = "Percentage", col="Class", data = mdf2, hue= "PAM", kind="bar",
            palette="coolwarm", legend=True)
plt.savefig("/Users/cd7/Desktop/last/Figure1C.png")

# B - Domain Coverage at nucleotide level

df3 = pandas.DataFrame([["CBE", "BRCT-1", percentage(278,3), percentage(278,37), percentage(278,40)],
                        ["CBE", "BRCT-2", percentage(296,21), percentage(296,31), percentage(296,52)],
                        ["CBE", "RING", percentage(122,7), percentage(122,8), percentage(122,15)],

                        ["ABE", "BRCT-1", percentage(278,8), percentage(278,31), percentage(278,39)],
                        ["ABE", "BRCT-2", percentage(296,17), percentage(296,27), percentage(296,44)],
                        ["ABE", "RING", percentage(122,13), percentage(122,12), percentage(122,25)],

                        ["CGBE", "BRCT-1", percentage(278,3), percentage(278,39), percentage(278,42)],
                        ["CGBE", "BRCT-2", percentage(296,23), percentage(296,33), percentage(296,56)],
                        ["CGBE", "RING", percentage(122,7), percentage(122,10), percentage(122,17)],

                        ["All BEs", "BRCT-1", percentage(278,11), percentage(278,70), percentage(278,81)],
                        ["All BEs", "BRCT-2", percentage(296,40), percentage(296,61), percentage(296,101)],
                        ["All BEs", "RING", percentage(122,20), percentage(122,22), percentage(122,42)]])

df3.columns = ["BE", "Domain", "NGG", "Only NG", "NGG+NG"]
mdf3 = df3.melt(id_vars=["BE", "Domain"], value_vars=["NGG", "Only NG", "NGG+NG"])
mdf3.columns = ["BE", "Domain", "PAM", "Percentage"]

sns.catplot(x="BE", y = "Percentage", col="Domain", data = mdf3, hue= "PAM", kind="bar",
            palette="coolwarm", legend=True)
plt.savefig("/Users/cd7/Desktop/Figure1B.png")



