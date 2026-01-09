#!/usr/bin/env python3
import os
import subprocess

input_dir = "/home/fhsoliv/Projects/GO_Proj/clustering/input"
out_dir = "/home/fhsoliv/Projects/GO_Proj/clustering/output"
#print(input_dir)

for inf in [5.0]:
    for root, dirs, files in os.walk(input_dir):
        for fname in files:
            if fname.endswith(".txt"):   
                file_path = os.path.join(root, fname)
                out_fname = "_".join([str(inf), fname])
                out_path = os.path.join(out_dir, fname.split("_")[0], out_fname)
                #print(out_path)
                print(f"Processing {file_path}")
                
                cmd = [
                    "mcl", file_path, "--abc",
                    "-I", str(inf), "-o", out_path
                ]
                #print(cmd)
                
                subprocess.run(cmd, check=True)
