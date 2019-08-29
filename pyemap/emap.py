from rdkit import Chem
from rdkit.Chem import Draw
from .custom_residues import is_pi_bonded, dist
from .data import *
import networkx as nx
from networkx.drawing.nx_agraph import from_agraph, to_agraph
from .shortest_paths import Branch, ShortestPath
from .smiles import getSimpleSmiles
from collections import OrderedDict
from PIL import Image
import warnings
import os
import logging
from pysmiles import write_smiles, fill_valence, correct_aromatic_rings, add_explicit_hydrogens


class emap():
    '''
    Manages the data generated at all stages of eMap analysis. 

    Attributes
    ----------
    filename: str
        Name of the crystal structure file being processed by pyemap.
    structure: Bio.PDB.Structure.Structure 
        Macromolecular protein structure. Contains model which contains residues and atoms.
    eta_moieties: dict of str: Bio.PDB.Residue.Residue
        Non-protein eta moieties automatically identified at the parsing step.
    chain_list: list of str
        List of chains identified at the parsing step.
    smiles: dict of str:Bio.PDB.Residue.Residue
        List of smiles strings for non-protein eta moieties identified at the parsing step.
    residues: dict of str:Bio.PDB.Residue.Residue
        Residues included in the graph after the process step.
    ngl_strings: dict of str:str
        Formatted NGL viewer selection strings for residues included in the graph after the process step.
    user_residues: dict of str:Bio.PDB.Residue.Residue
        Custom residues specified by the user.
    init_graph: networkx.Graph 
        Graph generated after the process step.
    branches: dict of str:pyemap.Branch
        Branches found by eMap analysis
    paths: dict of str:pyemap.ShortestPath
        Paths found by emap sorted by lowest to highest score.
    paths_graph: networkx.Graph 
        Graph generated after the shortest paths step.
    '''
    def __init__(self, filename, structure, eta_moieties, chain_list):
        '''Initializes emap object.

        Parameters
        ----------
        filename: str
            Name of file
        structure: Bio.PDB.Structure.Structure object
            Macromolecular protein structure
        eta_moieties: list of Bio.PDB.Residue.Residue objects
            Customized residue objects generated for automatically detected eta moieties
        chain_list: list of str
            Chains identified by the parser
        '''
        self.filename = filename
        self.structure = structure
        self.residues = {}
        self.chains = chain_list
        self.eta_moieties = {}
        self.user_residues = {}
        self.smiles = {}
        self.paths = OrderedDict()
        self.paths_graph = []
        self.init_graph  = []
        self.ngl_strings = {}
        self.branches= OrderedDict()
        for residue in eta_moieties:
            self._add_eta_moiety(residue)
    
    def _store_initial_graph(self, graph):
        '''Stores graph representation of emap model

        Parameters
        ----------
        graph: networkx.Graph
            Graph theory representation of emap model

        Notes
        -----
        Attributes of nodes/edges store info on surface exposure, edge weights etc. 
        '''
        self.init_graph = graph

    def _store_paths(self, branches,yens=False):
        '''Stores pathways in emap object, and sets their ngl selection strings for visualization.

        Parameters
        ---------
        shortest_paths: list of pyemap.shortest_paths.Branch
            branches found by emap
        yens: boolean, optional
            True when target specified, False when only source is specified
        '''
        shortest_paths=[]
        for br in branches:
            self.branches[br.branch_id] = br
            shortest_paths+=br.paths
        shortest_paths = sorted(shortest_paths)
        for pt in shortest_paths:
            self.paths[pt.path_id] = pt
            self._visualize_pathway(pt,yens)
        
    def _store_paths_graph(self, graph):
        '''Stores graph representation of emap model with selected pathway(s) highlighted

        Parameters
        ----------
        graph: networkx.Graph 
            Graph theory representation of emap model

        Notes
        -----
        Attributes of nodes/edges store info on surface exposure, edge weights etc. 
        '''
        self.paths_graph = graph

    def _reset_process(self):
        '''Returns emap object to state it was in after parsing.
        '''
        self.residues={}
        self.user_residues={}
        self.init_graph=[]
        self.paths = OrderedDict()
        self.paths_graph = []
        self.branches = OrderedDict()
        self.ngl_strings = {}
    
    def _reset_paths(self):
        '''Returns emap object to stat it was in after the process step.
        '''
        self.paths = OrderedDict()
        self.paths_graph=[]
        self.branches = OrderedDict()

    def _get_residue_graph(self,residue):
        atoms = list(residue.get_atoms())
        arom_atoms = ['O', 'P', 'N', 'C', 'S']
        res_graph = nx.Graph()
        for i in range(len(atoms)):
            for k in range(i, len(atoms)):
                if (not i == k) and is_pi_bonded(atoms[i], atoms[k]):
                    if atoms[i].element in arom_atoms and atoms[k].element in arom_atoms:
                        res_graph.add_edge(i, k)
                        res_graph.nodes[i]['element'] = atoms[i].element
                        res_graph.nodes[k]['element'] = atoms[k].element
        return res_graph

    def _add_eta_moiety(self,residue):
        '''Gets the smiles string for an automatically identified non-protein eta moiety, 
        and adds the residue to the eta_moieties dictionary.

        Parameters
        ----------
        residue: Bio.PDB.Residue.Residue
             Customized residue object generated for automatically detected eta moiety
        '''
        res_graph = self._get_residue_graph(residue)
        #smiles_str = getSimpleSmiles(res_graph, atoms)
        fill_valence(res_graph)
        add_explicit_hydrogens(res_graph)
        correct_aromatic_rings(res_graph)
        smiles_str = write_smiles(res_graph)
        molecule = Chem.MolFromSmarts(smiles_str)
        smiles_str = Chem.MolToSmarts(molecule, True)
        residue.smiles = smiles_str
        self.smiles[residue.resname] = smiles_str
        self.eta_moieties[residue.resname]=residue

    def _visualize_pathway(self,pathway,yens):
        '''
        '''
        colors = {"F": "orange", "Y": "blue", "W": "red", "H": "green"}
        selection_strs = []
        color_list = []
        labeled_atoms = []
        label_texts = [] 
        for res in pathway.get_path_as_list()[1:-1]:
            label_texts.append(res)
            try:
                color_list.append(colors[res[0]])
                labeled_atoms.append(".CA")
            except KeyError:
                color_list.append("pink")
                labeled_atoms.append(next(self.residues[res].get_atoms()).name)
            selection_strs.append(self.residues[res].ngl_string)
        color_list[0]="yellow"
        if yens:
            color_list[-1]="turquoise"
        pathway.set_visualization(selection_strs,color_list,labeled_atoms,label_texts)

    def _add_residue(self, residue):
        '''Gets ngl string for residue, and adds the residue to the residues and ngl_strings dictionaries.
        '''
        res_graph = self._get_residue_graph(residue)
        #smiles_str = getSimpleSmiles(res_graph, atoms)
        fill_valence(res_graph,respect_hcount=False,respect_bond_order=False)
        add_explicit_hydrogens(res_graph)
        correct_aromatic_rings(res_graph)
        smiles_str = write_smiles(res_graph)
        molecule = Chem.MolFromSmarts(smiles_str)
        smiles_str = Chem.MolToSmarts(molecule, True)
        residue.smiles = smiles_str
        self.smiles[residue.node_label] = smiles_str
        residue.ngl_string = self._get_ngl_string(residue)
        self.residues[residue.node_label] = residue
        self.ngl_strings[residue.node_label] = residue.ngl_string
    
    def _get_ngl_string(self,residue):
        """Returns NGL selection string for residue

        Parameters
        ----------
        residue: Bio.PDB.Residue.Residue object
        
        Returns
        -------
        select_string: str
            NGL selection string for residue
        """
        select_string = ""
        atm_list = list(residue.get_atoms())
        first_atm = atm_list[0]
        select_string += "(" + str(first_atm.original_id[3][1]) + " and :" + str(
            first_atm.original_id[2]) + " and ." + first_atm.name + ")"
        for i in range(1, len(atm_list)):
            atm = atm_list[i]
            select_string += " or "
            select_string += "(" + str(atm.original_id[3][1]) + " and :" + str(
                atm.original_id[2]) + " and ." + atm.name + ")"
        return select_string

    def save_residue(self, resname, dest="",size=(200,200)):
        '''Saves image of residue to file

        Parameters
        ----------
        resname: str
            Name of residue (node label) to be saved to file.
        dest: str, optional
            destination to save the image
        size: (float,float), optional
            dimensions of image saved to file
        '''

        if self.residues and resname in self.residues:
            if not self.residues[resname].smiles:
                raise KeyError("Not yet implemented for standard or custom residues.")
            mol = Chem.MolFromSmarts(self.residues[resname].smiles)
        elif resname in self.eta_moieties:
            mol = Chem.MolFromSmarts(self.eta_moieties[resname].smiles)
        else:
            raise KeyError("No record of any residue by that name.")
        if dest:
            Draw.MolToFile(mol, dest, kekulize=False, size=size)
        else:
            Draw.MolToFile(mol, resname+".png", kekulize=False, size=size)
    
    def show_residue(self,resname,size=(200,200)):
        '''Opens image of chemical structure in default image viewer.

        Notes
        -----
        Uses pillow library.
        '''
        if self.residues and resname in self.residues:
            if resname not in self.smiles:
                raise KeyError("Not yet implemented for standard or custom residues.")
            mol = Chem.MolFromSmarts(self.smiles[resname])
        elif resname in self.eta_moieties:
            mol = Chem.MolFromSmarts(self.smiles[resname])
        else:
            raise KeyError("No record of any residue by that name.")
        Draw.MolToFile(mol, "tmp.png", kekulize=False, size=size)
        img = Image.open("tmp.png")
        img.show()
        os.remove("tmp.png")

    def save_init_graph(self,dest=""):
        '''Saves image of graph generated by process step to file.

        Parameters
        ----------
        dest: str
            Destination for writing to file.
        '''
        if self.init_graph:
            if dest:
                fn=dest
            else:
                fn = self.filename[:-4] + "_graph.png"
            agraph=to_agraph(self.init_graph)
            agraph.graph_attr.update(ratio=1.0, overlap="false", mode="ipsep", splines="true")
            agraph.layout(args="-Gepsilon=0.05 -Gmaxiter=50")
            agraph.draw(fn,prog='neato')
        else:
            raise RuntimeError("Nothing to draw.")

    def save_paths_graph(self,dest=""):
        '''Saves image of graph generated by pathways step to file.

        Parameters
        ----------
        dest:str
            Destination for writing to file.
        '''
        if self.paths_graph:
            if dest:
                fn=dest
            else:
                fn = self.filename[:-4] + "_graph.png"
            agraph=to_agraph(self.paths_graph)
            agraph.graph_attr.update(ratio=1.0, overlap="false", mode="ipsep", splines="true")
            agraph.layout(args="-Gepsilon=0.05 -Gmaxiter=50")
            agraph.layout(args="-n2")
            agraph.draw(fn,prog='neato')
        else:
            raise RuntimeError("Nothing to draw.")
            
    def get_surface_exposed_residues(self):
        '''Returns list of surface exposed residues.

        Returns
        -------
        surface_exposed: list of str
            List of surface exposed residues identified by pyemap
        '''
        if self.init_graph:
            surface_exposed=[]
            for n, d in self.init_graph.nodes(data=True):
                if d['shape'] == "box":
                    surface_exposed.append(n)
            return surface_exposed
        else:
            raise RuntimeError("No graph found. Please run pyemap.process(my_emap) to generate the graph.")
    
    def report(self,dest=""):
        '''Writes report of most probable pathways to console or to file.

        Parameters
        -----------
        dest: str, optional
            Destination for writing to file

        Raises
        -------
        RuntimeError
            Nothing to report
        '''
        if self.branches:
            output = ""
            for br in self.branches.values():
                output += str(br)
            if dest:
                fi = open(dest,"w")
                fi.write(output)
            else:
                return output
        else:
            raise RuntimeError("Nothing to report.")
    
    def show_init_graph(self):
        '''Opens image of graph after processing in default image viewer.

        Notes
        ------
        Uses the Pillow library.
        '''
        if self.init_graph:
            fn = "tmp.png"
            agraph=to_agraph(self.init_graph)
            agraph.graph_attr.update(ratio=1.0, overlap="ipsep", mode="ipsep", splines="true")
            agraph.layout(args="-Gepsilon=0.05 -Gmaxiter=50")
            agraph.draw(fn,prog='neato')
            img = Image.open("tmp.png")
            img.show()
            os.remove("tmp.png")
        else:
            raise RuntimeError("Nothing to draw.")
    
    def show_paths_graph(self):
        '''Opens image of paths graph in default image viewer.

        Notes
        ------
        Uses the Pillow library.
        '''
        if self.paths_graph:
            fn = "tmp.png"
            agraph=to_agraph(self.paths_graph)
            agraph.graph_attr.update(ratio=1.0, overlap="false", mode="ipsep", splines="true")
            agraph.layout(args="-Gepsilon=0.05 -Gmaxiter=50")
            agraph.draw(fn,prog='neato')
            img = Image.open("tmp.png")
            img.show()
            os.remove("tmp.png")
        else:
            raise RuntimeError("Nothing to draw.")


            
    

    
