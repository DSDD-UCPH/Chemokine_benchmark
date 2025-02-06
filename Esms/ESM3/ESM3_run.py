from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
import matplotlib.pyplot as plt
import torch
from tqdm import tqdm
import matplotlib.pyplot as pl
#import py3Dmol
import torch
from esm.sdk.api import ESMProtein, GenerationConfig
from esm.sdk.forge import ESM3ForgeInferenceClient
from colorama import Fore
import os
import esm
from tools.ESM3_complex import build_protein_complex_app1, build_protein_complex_app3
from tools.chemo_recep_dict import dict_chemo, dict_recep
from tools.ESM3_util import fasta_to_dict
# you need a token from huggingface if this is your first time running esm-open.
#login()

# This is the token for using the esm
token = "rAwNezEAHQducqLX5P4rw"

#Lpath = "../../data/Ligands_uniprot_fasta/"
#Rpath = "../../data/Receptors_uniprot_fasta/"
#pdb_path = "../../pdbs/"
#path_predict = "./predicted"
#if not os.path.exists(path_predict ):os.mkdir(path_predict)

# Open-sourced ones
#model_list = ["esm3-open", "esm3_sm_open_v1", "esm3_medium-2024-03", "esm3_medium-2024-08", "esm3-medium-multimer-2024-09"]
model_list = ["esm3-medium-2024-08", "esm3-medium-multimer-2024-09"]
#___________________________CHOOSE YOUR PROTEINS HERE________________________________
#chemo_chain_dict = fasta_to_dict(Lpath, dict_chemo)
#recep_chain_dict = fasta_to_dict(Rpath, dict_recep)
#gene_id_1 = "CCR2"
#gene_id_2 = "CCL1"
#protein_chain_1 = recep_chain_dict[gene_id_1]
#protein_chain_2 = chemo_chain_dict[gene_id_2]
#chain1 = "A"
#chain2 = "A"
#__ their pdb files__ and which chains in their pdbs?
#path_pdb1 = f"../../AlphaFold_pdbs/Receptors/{gene_id_1}.pdb"
#path_pdb2 = f"../../AlphaFold_pdbs/chemokines/{gene_id_1}.pdb"

#___________________________________________________________________________________




def model_predict(model_n, complex, complex_name, model_pred_path, num_steps, temperature):
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    #device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    if model_n == "esm3-open" or model_n == "esm3_sm_open_v1":
        model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(device)
    else:
        # This need a token created on the website below to gain access.
        #model = ESM3ForgeInferenceClient(model="esm3-medium-2024-03", url="https://forge.evolutionaryscale.ai", token=token)
        model = esm.sdk.client(model_n, token= token)
    complex_predicted = model.generate(complex, GenerationConfig(track="structure", num_steps= num_steps, temperature= temperature))
    pdb_file_pred = model_pred_path + f"{complex_name}.pdb"
    print(complex_predicted)
    complex_predicted.to_pdb(pdb_file_pred)
    


def chemokine_predict_ESM3(model_list, path_pdb1, path_pdb2, chain1, chain2, path_predict, num_steps= 20, temperature= 0.70):
    for model_n in tqdm(model_list):
        os.mkdir(model_pred_path := path_predict + f"/{model_n}/") if not os.path.exists(model_pred_path := path_predict + f"/{model_n}/") else None
        protein_interaction, complex_name = build_protein_complex_app1(path_pdb1, path_pdb2, chain1, chain2)
        model_predict(model_n, protein_interaction, complex_name, model_pred_path, num_steps= num_steps, temperature= temperature)


def chemonike_predict_ESM3_from_chain(model_list, gene_id_1, gene_id_2, protein_chain_1, protein_chain_2, path_predict, num_steps= 20, temperature= 0.70):
    for model_n in tqdm(model_list):
        os.mkdir(model_pred_path := path_predict + f"/{model_n}/") if not os.path.exists(model_pred_path := path_predict + f"/{model_n}/") else None
        protein_interaction, complex_name = build_protein_complex_app3(gene_id_1, gene_id_2, protein_chain_1, protein_chain_2)
        model_predict(model_n, protein_interaction, complex_name, model_pred_path, num_steps= num_steps, temperature= temperature)



#________________This uses the structure of the monomers____________
#chemokine_predict_ESM3(model_list, path_pdb1, path_pdb2, chain1, chain2)

#________________This only uses the pure sequences_________________
#chemonike_predict_ESM3_from_chain(model_list, gene_id_1, gene_id_2, protein_chain_1, protein_chain_2)
