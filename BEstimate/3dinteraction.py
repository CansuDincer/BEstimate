#-----------------------------------------------------------------------------------------#
# 									  3D Interaction									  #
#-----------------------------------------------------------------------------------------#

# Import necessary packages
import os, pandas

#-----------------------------------------------------------------------------------------#
# YU Lab Interactome Insider

yulab = pandas.read_table("/Users/cd7/Documents/Rotations/Garnett_Group/workspace/venv/data/H_sapiens_interfaces.txt")


def extract_pis(pis):
	sites = list()
	print(pis)
	if pis != "[]":
		for site in pis.split(","):
			if site[0] == "[" and site[-1] != "]":
				s_first = site[1:]
				if len(s_first.split("-")) > 1:
					for s in list(range(int(s_first.split("-")[0]), int(s_first.split("-")[1])+1)):
						sites.append(int(s))
				else:
					sites.append(int(s_first))
			elif site[-1] == "]" and site[0] != "[":
				s_last = site[:-1]
				if len(s_last.split("-")) > 1:
					for s in list(range(int(s_last.split("-")[0]), int(s_last.split("-")[1])+1)):
						sites.append(int(s))
				else:
					sites.append(int(s_last))
			else:
				if len(site.split("-")) > 1:
					for s in list(range(int(site.split("-")[0]), int(site.split("-")[1])+1)):
						sites.append(int(s))
				else:
					sites.append(int(site))
		sites.sort()
		return sites
	else:
		return None


def collect_pis(uniprot):
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
		return True, pdb_partner_list, i3d_partner_list, eclair_partner_list
	else:
		False, None, None, None
