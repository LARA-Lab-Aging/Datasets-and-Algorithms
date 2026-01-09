import pandas as pd
import os
import numpy as np

from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.obo_parser import GODag, GOTerm
from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.semsim.termwise.wang import SsWang

#annopath= os.getcwd() + "/goa_human.gaf"
#hs_anno = GafReader("goa_human.gaf")
#teste = hs_anno.read_gaf()
#print(teste['GO:0007608'])

def main():
    data = pd.read_csv("C:/Documents/LAB R&A/GO Project/Results/level_depth_df_pval_rank_anno_v2.txt", sep = '\t')
    data = data.loc[(data["List_type"] == "gobp") & (data["adj_PVal"] < 1)]

    goids_500 = ["GO:1990166", "GO:0007616", "GO:0070301", "GO:1903599",
                    "GO:0048385", "GO:0048864", "GO:0042745", "GO:0042178",
                    "GO:0007566", "GO:0019370", "GO:0035542", "GO:0032288",
                    "GO:0060997", "GO:0007095", "GO:0042276"]
    goids_200 = ["GO:0048864", "GO:0048385","GO:0035542", "GO:1990166", "GO:1903599", "GO:0042178", "GO:0032288"]
    goids_100 = ["GO:0007616", "GO:0019370", "GO:1990166", "GO:0042178"]
    goids_50 = ["GO:0035542", "GO:0042745", "GO:0048385"]
    tools = data["Tool"].unique()
    all_goids = data["GOID"].unique()
    np.append(all_goids, goids_500)

    aux_dict = {500: goids_500, 200: goids_200, 100: goids_100, 50: goids_50}
    sizes = [500, 200, 100, 50]
## -----------------------------------------------------------------------------------
    obodag = download_go_basic_obo()
    obodag = GODag("go-basic.obo", optional_attrs={'relationship'}, load_obsolete=True)

    annopath = download_ncbi_associations()
    annopath= os.getcwd() + "/gene2go"
    annoobj = Gene2GoReader(annopath, godag=obodag, taxids = [9606])
    gene2gos = annoobj.get_id2gos()

    rels = {'part_of'}

    wang = SsWang(all_goids, obodag, rels)

    ssw = []

    for goid, size in zip(data["GOID"], data["Size"]):
        val = max_sim(goid, aux_dict[size], wang)   
        ssw.append(val)
        
    data["ssw"] = ssw

    conf_df = confusion_matrix(data, tools, sizes, aux_dict)
    print(conf_df)
    conf_df.to_csv("C:/Users/Usuario/Documents/GO Proj/Results/confusion_matrix2.txt", sep = '\t', index=False)
    
    
def max_sim(goid, list, wang):
    sim_vals = []
    for go in list:
        val = wang.get_sim(goid, go)
        sim_vals.append(val)

    max_val = max(sim_vals)
    
    return max_val

def confusion_matrix(df, tools, sizes, targ_dict):
    conf_matrix = []
    for size in sizes:
        size_df = df.loc[df["Size"] == size]
        for tool in tools:
            tmp = size_df.loc[size_df["Tool"] == tool]

            res_targets = tmp.loc[tmp["GOID"].isin(targ_dict[size])]

            if tool == "BiNGO" and size == 100:
                print(res_targets)
            #print(res_targets)
            total_targets = len(targ_dict[size])
            idd_targets = len((res_targets.loc[res_targets["adj_PVal"] <= 0.05])["GOID"].unique())
            median_rank = res_targets["Rank"].median()
            median_pval = res_targets["adj_PVal"].median()
            
            tp = tmp.loc[(tmp["adj_PVal"] <= 0.05) & (tmp["ssw"] >= 0.7)].shape[0]
            fp = tmp.loc[(tmp["adj_PVal"] <= 0.05) & (tmp["ssw"] <= 0.3)].shape[0]
            fn = tmp.loc[(tmp["adj_PVal"] >= 0.05) & (tmp["ssw"] >= 0.7)].shape[0]
            tn = tmp.loc[(tmp["adj_PVal"] >= 0.05) & (tmp["ssw"] <= 0.3)].shape[0]
            #tupla = {"Tool": tool, "size": size, "TP": tp, "FP": fp, "FN": fn, "TN": tn, "Median_rank": median_rank, "Median_pval": median_pval, "Identified OG Targets": f"{idd_targets}/{total_targets}"}
            tupla = {"Tool": tool, "Size": size, "Precision": fp/(fp+tp), "Accuracy": (tp + tn)/(tp + tn + fp + fn), "Sensitivity": tp/(tp + fn), "FPR": fp/(fp + tn), "Median_rank": median_rank, "Median_pval": median_pval, "Identified OG Targets": f"{idd_targets}/{total_targets}"}
            conf_matrix.append(tupla)

    conf_df = pd.DataFrame(conf_matrix)
    return conf_df

    
                


main()