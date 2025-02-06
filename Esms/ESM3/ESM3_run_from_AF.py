
import requests
from huggingface_hub import login
import os
from tools.chemo_recep_dict import dict_chemo, dict_recep
from ESM3_run import chemokine_predict_ESM3, model_list
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
import esm

token = "rAwNezEAHQducqLX5P4rw"
AlphaFold_path = "../../AlphaFold_pdbs/"
path_recep = f"../../AlphaFold_pdbs/Receptors/"
path_chemo = f"../../AlphaFold_pdbs/chemokines/"
path_predict_dir = "./from_AlphaFold_predicted"
if not os.path.exists(path_predict_dir ):os.mkdir(path_predict_dir)

#______downlaod the pdbs from AF database if you need it______
def download_AlphaFold_structures(dict_protein, type):
    for id in dict_protein.keys():
        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v4.pdb"
        file_name = f"{dict_protein[id]}.pdb"
        file_path = os.path.join(AlphaFold_path,type ,file_name)
        response = requests.get(pdb_url)
        if response.status_code == 200:
            with open(file_path, "wb") as file:
                file.write(response.content)
        else:
            print(f"Failed to download the structure file of {dict_protein[id]}.")
#download_AlphaFold_structures(dict_chemo, "Chemokines")
#download_AlphaFold_structures(dict_recep, "Receptors")
#___________________Lets predict_______________________________
for gene_id_1 in dict_recep.values():
    for gene_id_2 in dict_chemo.values():
        path_pdb1 = os.path.join(path_recep, gene_id_1 + ".pdb")
        path_pdb2 = os.path.join(path_chemo, gene_id_2 + ".pdb")
        chain1 = "A"
        chain2 = "A"
        path_predict = os.path.join(path_predict_dir, f"{gene_id_1}_{gene_id_2}")
        os.makedirs(path_predict, exist_ok=True)
        chemokine_predict_ESM3(model_list, path_pdb1, path_pdb2, chain1, chain2, path_predict)
#_____fix the pdb files________________________________
def fix_pdb_chains(input_pdb, output_pdb):
    first_chain = 'A'  
    second_chain = 'B'  
    switched_to_b = False 
    change_atom_num = False 
    output_lines = []

    with open(input_pdb, 'r') as infile:
        for i, line in enumerate(infile):
            if line.startswith(("ATOM", "HETATM")):  
                res_num = int(line[23:26].strip())  
                atom_num = int(line[7:11].strip())
                atom_type = line[77:78].strip()
                if res_num==1 and i>15 and atom_type == "N":
                    length_chainA = atom_num -1
                    change_atom_num = True
                    switched_to_b = True
                if switched_to_b and change_atom_num:
                    line = line[:21] + second_chain + line[22:] 
                    new_atom_num = atom_num - length_chainA 
                    new_atom_num_str = f"{new_atom_num:5d}"
                    line = line[:6] + new_atom_num_str + line[11:]
                else: 
                    pass
                output_lines.append(line)
    with open(output_pdb, 'w') as outfile:
        outfile.writelines(output_lines)

    print(f"Fixed PDB saved as {output_pdb}")

# Example usage:
#fix_pdb_chains("/Users/gzs260/GitHub/Chemokine_benchmark/Esms/ESM3/from_AlphaFold_predicted/CCR1_CCL1/esm3-medium-2024-08/CCR1_CCL1.pdb", "output.pdb")

path_fixed_pdbs = "./from_AF_predict_fixed_pdbs/"
if not os.path.exists(path_fixed_pdbs):os.mkdir(path_fixed_pdbs)
list_dir_CR=os.listdir(path_predict_dir)
model_medium = "esm3-medium-2024-08"
model_multi = "esm3-medium-multimer-2024-09"
def pdb_all_CR(model_n, path_predict_dir=path_predict_dir, path_fixed_pdbs=path_fixed_pdbs):
    path_model_fixed = os.path.join(path_fixed_pdbs, model_n)
    if not os.path.exists(path_model_fixed):os.mkdir(path_model_fixed)
    for dir in list_dir_CR:
        path_mess = os.path.join(path_predict_dir, dir, model_n, dir + ".pdb")
        path_fixed = os.path.join(path_model_fixed, dir + ".pdb")
        fix_pdb_chains(path_mess, path_fixed)

pdb_all_CR(model_n = model_medium)
pdb_all_CR(model_n = model_multi)