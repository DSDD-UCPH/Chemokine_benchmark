
from esm.sdk.api import ESMProtein, GenerationConfig
import os
import MDAnalysis as mda
from MDAnalysis.analysis import align
import py3Dmol
from pathlib import Path
from esm.sdk.api import ESMProtein, GenerationConfig
from esm.utils.structure.protein_chain import ProteinChain
from esm.sdk.forge import ESM3ForgeInferenceClient
from Bio.PDB import PDBParser
from esm.utils.structure.protein_complex import ProteinComplex


#_____________________________________________________________________
# This is the token for using the esm
token = "rAwNezEAHQducqLX5P4rw"

path_predict = "./predicted"
if not os.path.exists(path_predict):os.mkdir(path_predict)

# Which PdbID are you working on?
#PdbID1 = "1cm4"
#path_predict_pdb = f"./predicted/{PdbID1}/"
#if not os.path.exists(path_predict_pdb):os.mkdir(path_predict_pdb)


#model = ESM3ForgeInferenceClient(model="esm3-medium-2024-03", url="https://forge.evolutionaryscale.ai", token=token)
#_____________________________________________________________________


def align_pdb(PdbID1, PdbID2, path_predict_pdb):
    try:
        u1 = mda.Universe(path_predict_pdb + f"{PdbID1}.pdb")
        u2 = mda.Universe(path_predict_pdb + f"{PdbID2}.pdb")
    except Exception as e:
        raise ValueError(f"Error loading PDB files: {e}")

    print(f"Total atoms in {PdbID1}: {len(u1.atoms)}")
    print(f"Total atoms in {PdbID2}: {len(u2.atoms)}")

    try:
        ref_atoms = u1.select_atoms("protein and backbone")
        mov_atoms = u2.select_atoms("protein and backbone")
    except Exception as e:
        raise ValueError(f"Error selecting backbone atoms: {e}")

    # Fallback to all atoms if backbone atoms are missing
    if len(ref_atoms) == 0 or len(mov_atoms) == 0:
        print("No backbone atoms found, using all atoms for alignment.")
        ref_atoms = u1.select_atoms("all")
        mov_atoms = u2.select_atoms("all")

    print(f"Selected atoms in {PdbID1} for alignment: {len(ref_atoms)}")
    print(f"Selected atoms in {PdbID2} for alignment: {len(mov_atoms)}")

    # This is alignment based on RMSD
    try:
        rmsd = align.alignto(mov_atoms, ref_atoms)
        print(f"Alignment RMSD: {rmsd}")
    except Exception as e:
        raise ValueError(f"Error during alignment: {e}")

    # Save aligned files
    try:
        u1.atoms.write(path_predict_pdb + f"{PdbID1}_aligned.pdb")
        u2.atoms.write(path_predict_pdb + f"{PdbID2}_aligned.pdb")
    except Exception as e:
        raise ValueError(f"Error saving aligned PDB files: {e}")

# Function to visualize two PDBs: useful for plotting the experimental data and predicted models.
# PdbID2 is the ground truth. 
def show_two_pdbs(PdbID1, PdbID2, path_predict_pdb, do_align= False, style1="cartoon", style2="stick", color1="blue", color2="grey"):
    if do_align:
        align_pdb(PdbID1, PdbID2)
        with open(path_predict_pdb + f"{PdbID1}_aligned.pdb", "r") as pdb1_file:
            pdb1_content = pdb1_file.read()

        with open(path_predict_pdb + f"{PdbID2}_aligned.pdb", "r") as pdb2_file:
            pdb2_content = pdb2_file.read()
    else:
        with open(path_predict_pdb + f"{PdbID1}.pdb", "r") as pdb1_file:
            pdb1_content = pdb1_file.read()

        with open(path_predict_pdb + f"{PdbID2}.pdb", "r") as pdb2_file:
            pdb2_content = pdb2_file.read()

    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb1_content, "pdb")
    viewer.setStyle({style1: {"color": color1}})

    viewer.addModel(pdb2_content, "pdb")
    viewer.setStyle({"chain": "A"}, {"cartoon": {"color": "green"}})
    viewer.setStyle({"chain": "B"}, {"cartoon": {"color": "orange"}})

    
    # Set the background and zoom
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    return viewer.show()

#-------------------------building complexs from the monomers pdb----------------------------------------app1----------------------------
def build_protein_complex_app1(path_pdb1, path_pdb2, chain1, chain2):
    gene_id_1 = os.path.splitext(os.path.basename(path_pdb1))[0]
    gene_id_2 = os.path.splitext(os.path.basename(path_pdb2))[0]
    path_pdb1 = Path(path_pdb1)
    path_pdb2 = Path(path_pdb2)
    #gene_id_interaction = path_pdb1.parent.name
    #gene_id_1 = path_pdb1.parent.name
    #gene_id_2 = path_pdb2.parent.name
    protein_chain_1 = ProteinChain.from_pdb(path_pdb1, chain_id=chain1)
    protein_chain_2 = ProteinChain.from_pdb(path_pdb2, chain_id=chain2)
    protein_complex = ProteinComplex.from_chains([protein_chain_1, protein_chain_2])
    print(f"The protein complex sequence {protein_complex.sequence}")
    protein_interaction = ESMProtein.from_protein_complex(protein_complex)
    #complex_name = f"{gene_id_1}_{gene_id_2}"
    complex_name = f"{gene_id_1}_{gene_id_2}"
    return protein_interaction, complex_name
    # Predict
    #protein_interaction = model.generate(protein_interaction, GenerationConfig(track="structure", num_steps= 20, temperature= 0.0 ))
    #protein_interaction.to_pdb(path_predict_pdb + f"{gene_id_1}_{gene_id_2}.pdb")
    #visualize the complex predicted
    #name_predicted = f"{PdbID1}_app1"
    #name_real = f"{PdbID1}_real"
    #show_two_pdbs(name_predicted, name_real, do_align = do_align)


#-------------------------------------I am not sure if this function is useful :)----------------------------
def build_protein_complex_app1_pdb(protein_chain_1, protein_chain_2, model, complex_name, path_predict_pdb, do_align= False):

    protein_complex = ProteinComplex.from_chains([protein_chain_1, protein_chain_2])
    print(f"The protein complex sequence {protein_complex.sequence}")
    protein_interaction = ESMProtein.from_protein_complex(protein_complex)
    # Predict
    protein_interaction = model.generate(protein_interaction, GenerationConfig(track="structure", num_steps= 20, temperature= 0.0 ))
    protein_interaction.to_pdb(path_predict_pdb + f"{complex_name}_app1.pdb")
    #visualize the complex predicted
    name_predicted = f"{complex_name}_app1"
    name_real = f"{complex_name}_real" # this is the pdb file from the experiment
    show_two_pdbs(name_predicted, name_real, do_align = do_align)
#----------------------------------if you have only protein sequences------------------------------app3--------
def build_protein_complex_app3(gene_id_1, gene_id_2, protein_chain_1, protein_chain_2):
    complex_chain = protein_chain_1 +"|"+ protein_chain_2
    protein_interaction = ESMProtein(sequence=complex_chain)
    complex_name = gene_id_1 +"_"+ gene_id_2 + "from_seq"
    return protein_interaction, complex_name


#-------------------This is the code snippets provided by the esm GitHub for the complexs-------------------app2--------------------

def get_sample_protein_complex(PdbID1) -> ESMProtein:
    protein = ProteinComplex.from_rcsb(PdbID1)
    protein = ESMProtein.from_protein_complex(protein)
    return protein


def build_protein_complex_app2(PdbID1, path_predict_pdb , model, do_align= False):

    protein = get_sample_protein_complex(PdbID1)
    #sequence_length = len(protein.sequence)  
    folded_protein = model.generate(
        protein,
        GenerationConfig(
            track="structure", schedule="cosine", num_steps=15, temperature=0.0
        ),
    )
    assert isinstance(
        folded_protein, ESMProtein
    ), f"ESMProtein was expected but got {protein}"

    folded_protein.to_pdb(path_predict_pdb + f"{PdbID1}_app2.pdb")
    #visualize the complex predicted
    name_predicted = f"{PdbID1}_app2"
    name_real = f"{PdbID1}_real"
    show_two_pdbs(name_predicted, name_real, do_align= do_align)


#build_protein_complex_app1(PdbID1, path_predict_pdb = path_predict_pdb, do_align= False)
#build_protein_complex_app2(PdbID1, path_predict_pdb = path_predict_pdb, do_align= False)
