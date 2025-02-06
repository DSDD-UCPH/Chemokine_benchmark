import os
import numpy as np
import pypdb
import pylab as pl
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBList
import pandas as pd
import pypdb
from Bio.PDB.PDBParser import PDBParser
import colorama
from colorama import Fore, Style
from collections import defaultdict
from Bio.SeqUtils import seq3, seq1
from matplotlib import cm
from numpy.ma import masked_array
import seaborn as sns
import numpy.ma as ma
from matplotlib import cm
import matplotlib.pyplot as plt
from tqdm import tqdm

#----------download the pdb format-----------------
def get_PDBs_pdb(path,  PdbID):
    # This is the path to save the pdb files
    if not  os.path.exists(path):os.mkdir(path)
    P_problem = []
    wrong_format= []

    for P in PdbID:

        try:
            float(P)
            print(Fore.RED + f"ERROR: '{P}' :Wrong PdbID format")
            wrong_format.append(P)
            pass
        except:
            pdb_file= os.path.join(path,P+'.pdb')
            pdb_string = pypdb.get_pdb_file(P, filetype='pdb')#Source:Protien Data Bank
            print(pdb_string, file=open(os.path.join(pdb_file), 'w'))
        if P + '.pdb' not in os.listdir(path)+ wrong_format:
            P_problem.append(P)
        else:
            pass
    if not P_problem:
        print(Fore.GREEN + "All the protein structures are successfully downloaded")
        
    else:
        print(Fore.RED + f"ERROR: The proteins '{P_problem}' are not downloaded")
        return P_problem
    
#--------------download the PDB files in mmcif format-----------

def get_PDBs_cif(path, PdbID):
    if not  os.path.exists(path):os.mkdir(path)
    #dataset = pd.read_csv(path_dataset)
    #PdbID = df[pdb_id_col].unique().tolist()

    ID_problem = []
    wrong_format= []
    
    for ID in PdbID:
        print(ID)

        #try:
        #    float(ID)
        #    print(Fore.RED + f"ERROR:'{ID}' :Wrong PdbID format")
        #    wrong_format.append(ID)
        #    pass
        #except:

            
        pdb_list = PDBList() 
        pdb_filename = pdb_list.retrieve_pdb_file(ID, pdir= path, file_format="mmCif")


        if (ID.lower() + '.cif') not in (os.listdir(path) + wrong_format) and (ID.upper() + '.cif') not in (os.listdir(path) + wrong_format):
            ID_problem.append(ID)
        else:
            pass
    if not ID_problem:
        print(Fore.GREEN + "All the protein structures are successfully downloaded")
        
    else:
        print(Fore.RED + f"ERROR: The proteins '{ID_problem}' are not downloaded.\n Please rerun this function, untill the structure gets downloaded.")
        return ID_problem

#------general download function-----------
#This function doesn't care about the format!
def download_pdbs(df,path, ID_col, format):  # Add sth that if cif is not supported we can download pdbs. 

    if format.lower() == 'cif':
        problemcif = get_PDBs_cif(df, path, ID_col)
        if problemcif:
            for ID in problemcif:
                pdb_file= os.path.join(path,ID+'.pdb')
                pdb_string = pypdb.get_pdb_file(ID, filetype='pdb')
                print(pdb_string, file=open(os.path.join(pdb_file), 'w'))
                print(Fore.GREEN + f"Now '{problemcif}' is/are downloaded as well.")
        else:
            pass
            
                
         
    elif format.lower() == 'pdb':
        problempdb = get_PDBs_pdb(df, path, ID_col)
        if problempdb:
            for ID in problempdb:
                pdb_list = PDBList() 
                pdb_filename = pdb_list.retrieve_pdb_file(ID, pdir= path, file_format="mmCif")
                print(Fore.GREEN + f"Now '{problempdb}' is/are downloaded as well.")
        else:
            pass
                
        
        
    else:
        print(Fore.RED+ "Error: The format is not recognized")
        
#------------------------------these are subfunctions, so we don't work with them directly---------------------------------------Main Functions Below This Section--------------------------
def get_one_letter_code(residue):
    three_letter_code = residue.get_resname()
    one_letter_code = seq1(three_letter_code)
    return one_letter_code

def dist_cal(rdict_one, rdict_two, C_type = "CA") : #C_type = 'COM', 'CA', 'CB' // COM: Center Of Mass//
    """Returns the C-alpha distance between two residues"""
    try:
        diff_vector  = rdict_one[C_type.upper() + "_coord"] - rdict_two[C_type.upper() + "_coord"]
        dis = np.sqrt(np.sum(diff_vector * diff_vector))
    except:
        dis = None
        
    return dis

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = dist_cal(residue_one, residue_two)
    return answer


def get_residue_info_dict(res, non_AA=False): # If the non Amino Acid residues are not wanted, non_AA = False
    if not non_AA and get_one_letter_code(res) == 'X':
        return None
    else:
        res_dict = {}
        res_dict['PdbID'] = res.parent.parent.parent.id
        res_dict['Model'] = res.parent.parent.id
        res_dict['chain Id'] = res.parent.id
        res_dict['position'] = res.get_full_id()[3][1]
        if get_one_letter_code(res) == 'X':
            res_dict['residue'] = res.resname
            try:
                res_dict['COM_coord'] = res.center_of_mass()
            except AttributeError:
                res_dict['COM_coord'] = None
        
        else:
            res_dict['residue'] = get_one_letter_code(res)
            try:
                res_dict['CA_coord'] = res["CA"].coord
            except KeyError:
                res_dict['CA_coord'] = None
            try:
                res_dict['CB_coord'] = res["CB"].coord
            except KeyError:
                res_dict['CB_coord'] = None
            try:
                res_dict['COM_coord'] = res.center_of_mass()
            except AttributeError:
                res_dict['COM_coord'] = None
            
        return res_dict
    
def mask_row_and_column(matrix, row_list, col_list, masked_value=np.nan):
    masked_matrix = matrix.copy()
    for row, col in zip(row_list, col_list):
        masked_matrix[row, col] = masked_value

    return masked_matrix

# ----------The Main Functions Begin-----------

# This function makes a dataframe of all amino acid interactions in a complex
def get_DF_DIS(PdbID, PDB_dir, non_Amino_Acid=True):
    if isinstance(PdbID, list):
        pass
    else:
        PdbID = PdbID.unique().tolist()

    dataframe_list = []

    file_extensions = ['.cif', '.pdb']

    for ID in tqdm(PdbID):
        if any(os.path.isfile(os.path.join(PDB_dir, ID + ext)) for ext in file_extensions):
            file_extension = next(ext for ext in file_extensions if os.path.isfile(os.path.join(PDB_dir, ID + ext)))
            parser = MMCIFParser(QUIET=True) if file_extension == '.cif' else PDBParser(QUIET=True)
            structure = parser.get_structure(ID, os.path.join(PDB_dir, ID + file_extension))
        else:
            print(Fore.RED + f"The protein '{ID}' has not been downloaded")
            continue

        data_dict_list = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    data_dict_list.append(get_residue_info_dict(residue, non_AA=non_Amino_Acid))

        dis_list = []

        for ri_dict in data_dict_list:
            for rj_dict in data_dict_list:
                if ri_dict and rj_dict and ri_dict['Model'] == rj_dict['Model']:
                    rr_dict = {
                        'PdbID': ri_dict['PdbID'],
                        'Model': ri_dict['Model'],
                        'ChainRes1': ri_dict['chain Id'],
                        'ChainRes2': rj_dict['chain Id'],
                        'Res1': ri_dict['position'],
                        'Res2': rj_dict['position'],
                        'AA1': ri_dict['residue'],
                        'AA2': rj_dict['residue'],
                        'd_CA': dist_cal(ri_dict, rj_dict, "CA"),
                        'd_CB': dist_cal(ri_dict, rj_dict, "CB"),
                        'd_COM': dist_cal(ri_dict, rj_dict, "COM"),
                    }
                    dis_list.append(rr_dict)

        df = pd.DataFrame(dis_list)
        dataframe_list.append(df)

    print(Fore.GREEN + "The protein dataframes can be called in this format: _PdbID_df \n" +
          "Example: _6mcf_df")

    return dataframe_list

# This functions plots a contact heatmap of amino acids, the color represents the distances.
def matrix_CA(PdbID, d_f, pic_dir, model, chain1 , chain2, matrix = True, nonAA = False, plot = False, nonAA_info_table = False, annot = False):     
    #pic_dir  = PDB_dir.replace('/pdbs/', '/Distance/')
    if not  os.path.exists(pic_dir):os.mkdir(pic_dir)
    #pic_dir = pic_dir + PdbID + "/"
    #if not  os.path.exists(pic_dir):
    #    os.mkdir(pic_dir)
    row = []
    col = []
    row2 = []
    col2 = []
    non_AA = []
    data_nonAA = defaultdict(list)
    #d_f = locals()['_'+ PdbID + '_df']
    df = d_f[d_f["Model"] == model]
    df = df[(df["ChainRes1"] == chain1) & (df["ChainRes2"] == chain2)]
    num_AA1 = np.max(df.Res1) - np.min(df.Res1) + 1
    num_AA2 = np.max(df.Res2) - np.min(df.Res2) + 1
    answer = np.full((num_AA2 , num_AA1 ), np.nan)

    for res1, res2, i_, j_, k in zip(df.AA1, df.AA2, df.Res2, df.Res1, df.index):
        j = j_ - np.min(df.Res1)     # a shift to get rid of white area in the heatmap
        i = i_ - np.min(df.Res2)
        if np.isnan(df.d_CA[k]) and nonAA:
            answer[i, j] = df.d_COM[k]
            row.append(i)
            col.append(j)
            non_AA.append({'res1': res1, 'res2': res2, 'pos1': j_, 'pos2' : i_, 'Distance' : df.d_COM[k]})
     
        else:   
            answer[i, j] = df.d_CA[k]
            row2.append(i)
            col2.append(j)
    for d in non_AA:
        for key, value in d.items():
            data_nonAA[key].append(value)

    dataframe = pd.DataFrame(data_nonAA)
    
    # answer is a matrix with with non amino acid included. CA_mat is only for amino acid connections.

    
    if plot:
        fig = plt.subplots(figsize = (30, 15))
        CA_mat = mask_row_and_column(answer, row, col)
        CA = sns.heatmap(CA_mat, cbar = True, annot = annot, fmt=".2f", cbar_kws = {'label': 'C alpha distances'})
        plt.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)


        
        if nonAA:
            COM_mat = mask_row_and_column(answer, row2, col2)
            COM_cm = sns.dark_palette("#69d", reverse=False, as_cmap=True)
            COM = sns.heatmap(COM_mat, cmap = COM_cm, cbar = True, annot = annot,
                              xticklabels = np.arange(np.min(df.Res1), np.max(df.Res1) + 1),
                              yticklabels = np.arange(np.min(df.Res2), np.max(df.Res2) + 1), fmt=".2f",
                              cbar_kws = {'label': 'COM distances'})
    
        plt.xlabel("Pos_AA_1  chain:" + chain1, fontsize=16)
        plt.ylabel("Pos_AA_2  chain:" + chain2, fontsize=16)
        plt.title(f"distance matrix for protein '{PdbID}' model '{model}' between chains '{chain1}' and '{chain2}' \n", fontsize=20)
        plt.savefig(pic_dir+ PdbID +chain1+"_to_"+chain2+".png", dpi=300, bbox_inches='tight')
    
        if nonAA_info_table:
            plt.subplots_adjust(left=0.2, bottom=0.4)
            the_table = plt.table(cellText=dataframe.values,
                                  rowLabels=dataframe.index,
                                  colLabels=dataframe.columns,
                                  cellLoc = 'center', rowLoc = 'center',
                                  loc='bottom', transform=plt.gcf().transFigure,
                                  bbox = ([0.25, -0.5, 0.3, 0.8]))
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(8)
            plt.text(0.41, 0.32, "The non amino acid contacts", ha='center', va='center', transform=plt.gcf().transFigure, fontsize=16)
        
    if nonAA and matrix:
        return answer
    elif matrix:
        CA_mat = mask_row_and_column(answer, row, col)
        return CA_mat
    else:
        pass




