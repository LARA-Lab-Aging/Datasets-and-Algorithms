from goatools.obo_parser import GODag, GOTerm
from goatools.godag.go_tasks import get_go2parents
from goatools.godag.go_tasks import get_go2children
import pandas as pd
import os

godag = GODag("go-basic.obo", load_obsolete=True)

work_dir = r"/home/fhsoliv/Projects/GO_Proj/org_res"
dirs = os.listdir(work_dir)
print(dirs)
w_dirs = [dir for dir in dirs if "random" not in dir]
print(w_dirs)

cols = ["List_type", "Size", "Tool", "GOID", "Name", "Level", "Depth"]



def df_builder():
    depth_df = pd.DataFrame()
    obsolete_df = pd.DataFrame()

    for dir in w_dirs:
        print(f"-- Now processing: {dir}")
        dir_path = "/".join([work_dir, dir])
        files = os.listdir(dir_path)
        list_type = dir.split("_")[0]
        size = dir.split("_")[1]

        for file in files:
            print(f"----- Now processing: {file}")
            tool = file.split("_")[0]
            file_path = "/".join([work_dir, dir, file])

            df = pd.read_csv(file_path, header = 0, sep = '\t')

            go_pval = get_goid(df, tool, size)
            tool_depth_df, tool_obs_df = find_depth(go_pval=go_pval, tool=tool, list_type=list_type, size=size)

            depth_df = pd.concat([depth_df, tool_depth_df])
            obsolete_df = pd.concat([obsolete_df, tool_obs_df])
        
        
    depth_df.to_csv("/home/fhsoliv/Projects/GO_Proj/Results/level_depth_df_pval_rank_anno.txt", sep = '\t', index=False)
    obsolete_df.to_csv("/home/fhsoliv/Projects/GO_Proj/Results/obsolete_count_pval.txt", sep = '\t', index=False)
    


def find_depth(go_pval, tool, list_type, size):
    data = []
    obs_count = 0
    obs_rank = [0]
    obs_pval = [0]
    idx = 0
    go_pval.columns = ["goid", "pval", "anno"]
    tam = go_pval.shape[0]

    for goid, pval, anno in zip(go_pval["goid"], go_pval["pval"], go_pval["anno"]):
        if goid != 'UNCLASSIFIED':    
            tmp = godag[goid]
            tupla = {"List_type": list_type, "Size": size,"Tool": tool, "Length": tam, "GOID": tmp.item_id, "anno": anno, "adj_PVal":pval, "Name": tmp.name, "Level": tmp.level, "Depth": tmp.depth}
            data.append(tupla)

            if tmp.is_obsolete == True:
                obs_count += 1
                obs_pval.append(pval)
                obs_rank.append(idx)

        idx += 1
    obs_tupla = [{"List_type": list_type, "Size": size, "Tool": tool, "Length": tam, "Obsolete_count": obs_count, "%": obs_count/tam, "Obsolete_ranks": obs_rank, "Obsolete_pval": obs_pval}]

    data_df = pd.DataFrame(data)
    data_df = data_df.sort_values(by=["adj_PVal"])
    data_df["Rank"] = [i for i in range(data_df.shape[0])]
    obs_df = pd.DataFrame(obs_tupla)

    return data_df, obs_df

def get_goid(res_df, tool, size):
    tool_col = {"DAVID": "GOID", "BiNGO": "GOID", "ClueGO": "ID", "clusterProfiler": "ID",
                "ENRICHR": "GOID", "goana": "GOID", "GOstats": "GOBPID", "gProfiler": "term_id", 
                "PANTHER": "GOID", "ShinyGO": "GOID", "topGO": "GO.ID", "WebGestalt": "geneSet"}
    
    pval_col = {"DAVID": "FDR", "BiNGO": "corr p-value", "ClueGO": "Term PValue Corrected with Benjamini-Hochberg", 
                "clusterProfiler": "p.adjust", "ENRICHR": "Adjusted P-value", "goana": "adj_pval", "GOstats": "FDR",
                "gProfiler": "p_value", "PANTHER": "FDR", "ShinyGO": "Enrichment FDR", "topGO": "FDR", "WebGestalt": "FDR"}
    
    anno_col = {"DAVID": "Pop Hits", "BiNGO": "n", "ClueGO": "All Associated Genes", 
                "clusterProfiler": "BgRatio", "ENRICHR": "Overlap", "goana": "N", 
                "GOstats": "Size", "gProfiler": "term_size", "PANTHER": "Ref",
                "ShinyGO": "Pathway Genes", "topGO": "Annotated", "WebGestalt": "size"}
    if tool == "BiNGO":
        if int(size) <= 100:
            res_df["tmp"] = [str(num)[0:(len(str(num))-2)] for num in res_df[anno_col[tool]]]
        else:
            res_df["tmp"] = [str(num)[0:(len(str(num))-3)] for num in res_df[anno_col[tool]]]
        cols = [tool_col[tool], pval_col[tool], "tmp"]
    elif tool == "ClueGO":
        res_df["tmp"] = [len(ls.split(", ")) for ls in res_df[anno_col[tool]]]
        cols = [tool_col[tool], pval_col[tool], "tmp"]
    elif tool == "clusterProfiler":
        res_df["tmp"] = [num.split('/')[0] for num in res_df[anno_col[tool]]]
        cols = [tool_col[tool], pval_col[tool], "tmp"]
    elif tool == "ENRICHR":
        res_df["tmp"] = [num.split('/')[1] for num in res_df[anno_col[tool]]]
        cols = [tool_col[tool], pval_col[tool], "tmp"]
    else:
        cols = [tool_col[tool], pval_col[tool], anno_col[tool]]

    #go_list = res_df[tool_col[tool]]
    new_df = res_df[cols]
    
    return new_df    

df_builder()


