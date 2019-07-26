from rdkit import Chem
from rdkit.Chem import Draw

class emap():
    def __init__(self,filename,structure,chain_list,custom_residue_list,smiles_str_list):
        self.filename = filename
        self.structure=structure
        self.chains = chain_list
        self.custom_residues = custom_residue_list
        custom_residue_names=[]
        for res in custom_residue_list:
            custom_residue_names.append(res.resname)
        self.smiles_list = smiles_str_list
        self.custom_dict = dict(zip(custom_residue_names,smiles_str_list))
    def save_graph(self):
        pass
    def save_paths(self):
        pass
    def view_custom_residue(self,resname):
        mol = Chem.MolFromSmarts(self.custom_dict.get(resname))
        Draw.ShowMol(mol, kekulize=False, size=(100, 100))