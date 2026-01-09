import pandas as pd
import os
import numpy as np
import itertools
from datetime import datetime

from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.obo_parser import GODag, GOTerm
from goatools.semsim.termwise.wang import SsWang

def main():

    out_path = "/home/fhsoliv/Projects/GO_Proj/clustering/input"
    data = pd.read_csv("/home/fhsoliv/Projects/GO_Proj/Results/level_depth_df_pval_rank_anno_v3.txt", sep = '\t')

    df = data.loc[data["adj_PVal"] < 0.05]
    tools = np.sort(data["Tool"].unique())
    #tools = ["ENRICHR"]
    sizes = [50, 100, 200, 500]
    #sizes = [500]
    ls_types = ["hallmark", "contextual"]
    #ls_types = ["contextual"]

    obodag = GODag("go-basic.obo", optional_attrs={'relationship'}, load_obsolete=True)
    with open("/home/fhsoliv/Projects/GO_Proj/clustering/log.txt", 'w') as logfile:
        current_time = datetime.now()
        logfile.write(f"{current_time}\n")
        for tool in tools:
            for list in ls_types:
                for size in sizes:
                    print(f"{tool}: {size} - {list}")
                    logfile.write(f"{tool}: {size} - {list}\n")
                    tmp_df = df.loc[(df["Size"] == size) & (df["List_type"] == list) & (df["Tool"] == tool)].copy()
                    logfile.write(f"Sig. Ontologies: {tmp_df.shape[0]}\n")
                    input = input_builder(tmp_df, obodag)
                    logfile.write(f"Total Edges: {input.shape[0]}\n\n")
                    if tmp_df.shape[0] != 0 and input.shape[0] != 0:
                        #print(f"{tool}: {list} - {size}")
                        #print(tmp_df.shape[0])
                        #print(input.shape[0])
                        #print("-----------------------------------------------------")
                        input.to_csv("/".join([out_path, tool, f"{tool}_{list}_{size}.txt"]), sep="\t", index = False, header=False)
    
def input_builder(df, obodag):
    all_goids = df["GOID"]
    combs = list(itertools.combinations(all_goids, 2))
    #print(len(all_goids))
    #print(len(combs))
    
    rels = {'part_of'}
    wang = SsWang(all_goids, obodag, rels)

    t_score = 0
    data_rows = []
    for tup in combs:
        wang_ss = wang.get_sim(tup[0], tup[1])
        t_score+=wang_ss
        if wang_ss >= 0.5:
            row = {"GOID_A":tup[0], "GOID_B":tup[1],"Wang_SS": wang_ss}
            data_rows.append(row)
    #print(t_score/len(combs))

    input = pd.DataFrame(data_rows)

    return input

main()