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
    def save_agraph(self,user_res_list,agraph):
        self.user_res_list=user_res_list
        self.agraph = agraph
    def save_paths(self):
        pass
    def view_custom_residue(self,resname):
        #Just writes to file, would prefer if it displays to GUI
        mol = Chem.MolFromSmarts(self.custom_dict.get(resname))
        Draw.MolToFile(mol, resname+".png", kekulize=False, size=(200, 200))









