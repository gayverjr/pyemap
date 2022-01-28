## Copyright (c) 2017-2022, James Gayvert, Ruslan Tazhigulov, Ksenia Bravaya
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from .custom_residues import is_pi_bonded
from .data import clusters, metal_ligands, res_name_to_smiles
import networkx as nx
from .utils import extract_resname
from networkx.drawing.nx_agraph import to_agraph
from collections import OrderedDict
from PIL import Image
import os
from shutil import copyfile
import tempfile
from cairosvg import svg2png
from .process_data import get_atom_list
import datetime
from pysmiles import write_smiles
from .pyemap_exceptions import *


class emap():
    '''
    Manages the data generated at all stages of PyeMap analysis.

    Attributes
    ----------
    file_path: str
        Crystal structure file being analyzed by PyeMap.
    eta_moieties: dict of str: :class:`Bio.PDB.Residue.Residue`
        Non-protein eta moieties automatically identified at the parsing step.
    chain_list: list of str
        List of chains identified at the parsing step.
    sequences: dict of str, str
        Amino acid sequence for each chain in FASTA format
    residues: dict of str: :class:`Bio.PDB.Residue.Residue`
        Residues included in the graph after the process step.
    user_residues: dict of str: :class:`Bio.PDB.Residue.Residue`
        Custom residues specified by the user.
    init_graph: :class:`networkx.Graph`
        Graph generated after the process step.
    branches: dict of int: :class:`~pyemap.Branch`
        Branches found by PyeMap analysis
    paths: dict of str: :class:`~pyemap.ShortestPath`
        Paths found by PyeMap sorted by lowest to highest score.
    paths_graph: :class:`networkx.Graph`
        Graph generated after the shortest paths step.
    '''

    def __init__(self, file_path, pdb_id, eta_moieties, chain_list, sequences):
        '''Initializes emap object.

        Parameters
        ----------
        file_path: str
            Name of file
        eta_moieties: list of :class:`Bio.PDB.Residue.Residue`
            Customized residue objects generated for automatically detected eta moieties
        chain_list: list of str
            Chains identified by the parser
        sequences: dict of str:str
            Key is chain id, value is sequence in fasta format
        '''
        self.file_path = file_path
        self.pdb_id = pdb_id
        self.residues = {}
        self.chains = chain_list
        self.sequences = sequences
        self.active_chains = {}
        self.eta_moieties = {}
        self.user_residues = {}
        self._include_residues = []
        self._process_params = {}
        self.paths = OrderedDict()
        self.paths_graph = None
        self.init_graph = None
        self.branches = OrderedDict()
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
        self.init_graph = graph.copy()

    def _store_paths(self, branches, yens=False):
        '''Stores pathways in emap object, and sets their ngl selection strings for visualization.

        Parameters
        ---------
        shortest_paths: list of pyemap.shortest_paths.Branch
            branches found by emap
        yens: boolean, optional
            True when target specified, False when only source is specified
        '''
        shortest_paths = []
        for br in branches:
            self.branches[br.branch_id] = br
            shortest_paths += br.paths
        shortest_paths = sorted(shortest_paths)
        for pt in shortest_paths:
            self.paths[pt.path_id] = pt
            self._visualize_pathway(pt, yens)

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
        self.paths_graph = graph.copy()

    def _reset_process(self):
        '''Returns emap object to state it was in after parsing.
        '''
        self.residues = {}
        self.active_chains = {}
        self.user_residues = {}
        self.init_graph = None
        self.paths = OrderedDict()
        self.paths_graph = None
        self.branches = OrderedDict()
        self._include_residues = []
        self._process_params = {}

    def _reset_paths(self):
        '''Returns emap object to state it was in after the process step.
        '''
        self.paths = OrderedDict()
        self.paths_graph = None
        self.branches = OrderedDict()

    def _get_residue_graph(self, residue):
        atoms = list(get_atom_list(residue))
        arom_atoms = ['O', 'P', 'N', 'C', 'S']
        res_graph = nx.Graph()
        for i in range(len(atoms)):
            for k in range(i, len(atoms)):
                if (not i == k) and is_pi_bonded(atoms[i], atoms[k]):
                    if atoms[i].element in arom_atoms and atoms[k].element in arom_atoms:
                        res_graph.add_edge(i, k)
                        res_graph.nodes[i]["element"] = atoms[i].element
                        res_graph.nodes[i]["coords"] = atoms[i].coord
                        res_graph.nodes[k]["element"] = atoms[k].element
                        res_graph.nodes[k]["coords"] = atoms[k].coord
        return res_graph

    def _add_eta_moiety(self, residue):
        '''Gets the smiles string for an automatically identified non-protein eta moiety,
        and adds the residue to the eta_moieties dictionary.

        Parameters
        ----------
        residue: Bio.PDB.Residue.Residue
             Customized residue object generated for automatically detected eta moiety
        '''
        if extract_resname(residue) not in clusters + list(metal_ligands.keys()) and "CUST" not in residue.resname:
            res_graph = self._get_residue_graph(residue)
            try:
                smiles_str = write_smiles(res_graph)
            except Exception as e:
                smiles_str = "unknown"
            residue.smiles = smiles_str
        elif extract_resname(residue) in metal_ligands:
            atm_name = next(residue.get_atoms()).name
            atm_name = atm_name[0] + atm_name[1].lower()
            charge = metal_ligands[extract_resname(residue)]
            smiles_str = "[" + atm_name
            if charge > 0:
                smiles_str += '+{}'.format(charge)
            elif charge < 0:
                smiles_str += '-{}'.format(charge)
            smiles_str += ']'
            residue.smiles = smiles_str
        self.eta_moieties[residue.resname] = residue

    def visualize_pathway_in_nglview(self,ptid,view):
        '''Visualize pathway in nglview widget
        
        Parameters
        ----------
        ptid: str
            Pathway ID to be visualized
        view: :class:`nglview.widget.NGLWidget`
            NGL Viewer widget 
        '''
        pt = self.paths[ptid]
        for i in range(0,len(pt.selection_strs)):
            first_atm_select = pt.selection_strs[i][:pt.selection_strs[i].index(')')+1]+"]"
            view.add_representation('ball+stick',sele=pt.selection_strs[i],color=pt.color_list[i])
            view.add_representation('label',color="black",sele=first_atm_select,labelText=pt.label_texts[i])

    def _visualize_pathway(self, pathway, yens):
        colors = {"F": "orange", "Y": "blue", "W": "red", "H": "green"}
        selection_strs = []
        color_list = []
        labeled_atoms = []
        label_texts = []
        for res in pathway.get_path_as_list()[1:-1]:
            label_texts.append(res)
            try:
                if res not in self.eta_moieties:
                    color_list.append(colors[res[0]])
                    labeled_atoms.append(".CA")
                else:
                    color_list.append("pink")
                    labeled_atoms.append(next(self.residues[res].get_atoms()).name)
            except KeyError:
                color_list.append("pink")
                labeled_atoms.append(next(self.residues[res].get_atoms()).name)
            selection_strs.append(self.residues[res].ngl_string)
        color_list[0] = "yellow"
        if yens:
            color_list[-1] = "turquoise"
        pathway.set_visualization(selection_strs, color_list, labeled_atoms, label_texts)

    def _add_residue(self, residue):
        residue.ngl_string = self._get_ngl_string(residue)
        if residue.resname in res_name_to_smiles:
            residue.smiles = res_name_to_smiles[residue.resname.upper()]
        self.residues[residue.node_label] = residue
        chain_id = residue.full_id[2]
        if chain_id in self.active_chains:
            self.active_chains[chain_id].append(residue)
        else:
            self.active_chains[chain_id] = [residue]

    def _get_ngl_string(self, residue):
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
        if hasattr(first_atm, "original_id"):
            id = first_atm.original_id
        else:
            id = first_atm.full_id
        select_string += "(" + str(id[3][1]) + " and :" + str(
            id[2]) + " and ^" + residue.id[2] + " and ." + first_atm.name + ")"
        for i in range(1, len(atm_list)):
            atm = atm_list[i]
            if hasattr(atm, "original_id"):
                id = first_atm.original_id
            else:
                id = atm.full_id
            select_string += " or "
            select_string += "(" + str(id[3][1]) + " and :" + str(
                id[2]) + " and ^" + residue.id[2] + " and ." + atm.name + ")"
        return select_string

    def residue_to_file(self, resname, dest="", size=(100, 100)):
        '''Saves image of residue to file in .svg format.

        Parameters
        ----------
        resname: str
            Name of residue (node label) to be saved to file.
        dest: str, optional
            destination to save the image
        size: (float,float), optional
            dimensions of image saved to file
        '''
        try:
            from rdkit import Chem
            from rdkit.Chem.Draw import rdMolDraw2D
        except (ModuleNotFoundError, ImportError) as e:
            raise ModuleNotFoundError(
                "RDKit is required for visualization of chemical structures. See https://www.rdkit.org/docs/Install.html"
            ) from e
        if resname in self.residues or resname in self.eta_moieties:
            if "CUST" in resname:
                raise KeyError("Not available for user defined residues.")
            elif resname[:3] in clusters:
                cluster_img_name = os.path.abspath(
                    os.path.dirname(__file__)) + '/data/clusters/' + resname[:3] + '.svg'
                if dest:
                    target_name = dest
                else:
                    target_name = resname + ".svg"
                copyfile(cluster_img_name, target_name)
            else:
                try:
                    try:
                        mol = Chem.MolFromSmiles(self.eta_moieties[resname].smiles)
                    except Exception:
                        mol = Chem.MolFromSmiles(self.residues[resname].smiles)
                    d2d = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
                    d2d.DrawMolecule(mol)
                    d2d.FinishDrawing()
                    if dest:
                        with open(dest, 'w') as f:
                            f.write(d2d.GetDrawingText())
                    else:
                        with open(resname + '.svg', 'w') as f:
                            f.write(d2d.GetDrawingText())
                except Exception as e:
                    raise PyeMapException("Could not draw residue: {}".format(resname)) from e
        else:
            raise KeyError("No record of any residue by that name.")

    def residue_to_Image(self, resname, scale=1.0):
        '''Returns PIL image of chemical structure. 

        Parameters
        -----------
        resname: str
            Name of residue
        scale: float, optional
            Output scaling factor, default dimensions are (100,100)
        Returns
        --------
        img: :class:`PIL.Image.Image`
        '''
        try:
            from rdkit import Chem
            from rdkit.Chem.Draw import rdMolDraw2D
        except (ModuleNotFoundError, ImportError) as e:
            raise ModuleNotFoundError(
                "RDKit is required for visualization of chemical structures. See https://www.rdkit.org/docs/Install.html"
            ) from e
        if resname in self.residues or resname in self.eta_moieties:
            if "CUST" in resname:
                raise KeyError("Not available for user defined residues.")
            elif resname[:3] in clusters:
                dest = tempfile.NamedTemporaryFile(suffix=".png").name
                cluster_img_name = os.path.abspath(
                    os.path.dirname(__file__)) + '/data/clusters/' + resname[:3] + '.svg'
                svg2png(url=cluster_img_name, write_to=dest, scale=scale)
                img = Image.open(dest)
                return img
            try:
                try:
                    mol = Chem.MolFromSmiles(self.eta_moieties[resname].smiles)
                except Exception:
                    mol = Chem.MolFromSmiles(self.residues[resname].smiles)
                d2d = rdMolDraw2D.MolDraw2DSVG(100, 100)
                d2d.DrawMolecule(mol)
                d2d.FinishDrawing()
                dest1 = tempfile.NamedTemporaryFile(suffix=".svg").name
                dest2 = tempfile.NamedTemporaryFile(suffix=".png").name
                with open(dest1, 'w') as f:
                    f.write(d2d.GetDrawingText())
                svg2png(url=dest1, write_to=dest2, scale=scale)
                img = Image.open(dest2)
                return img
            except Exception as e:
                raise PyeMapException("Could not draw residue: {}".format(resname)) from e
        else:
            raise KeyError("No record of any residue by that name.")

    def _graph_to_file(self, G, dest=""):
        '''Saves image of graph generated by process step to file.

        Parameters
        ----------
        dest: str
            Destination for writing to file.
        '''
        if G:
            agraph = to_agraph(G)
            agraph.graph_attr.update(ratio=1.0, overlap="rc", mode="ipsep", splines="true")
            if agraph.number_of_nodes() <= 200:
                try:
                    agraph.layout(prog='neato', args="-Gepsilon=0.01 -Gmaxiter=50")
                except Exception as e:
                    raise RuntimeError("There was a problem with Graphviz. See https://graphviz.gitlab.io/") from e
            else:
                try:
                    agraph.layout(prog='dot')
                except Exception as e:
                    raise RuntimeError("There was a problem with Graphviz. See https://graphviz.gitlab.io/") from e
            if dest:
                svg_fn = dest + '.svg'
                png_fn = dest + '.png'
                agraph.draw(svg_fn, prog='neato', args="-Gepsilon=0.01 -Gmaxiter=50")
                agraph.draw(png_fn, prog='neato', args="-Gepsilon=0.01 -Gmaxiter=50")
            else:
                png_fn = self.file_path[:-4] + "_graph.png"
                agraph.draw(png_fn, prog='neato', args="-Gepsilon=0.01 -Gmaxiter=50")
        else:
            raise RuntimeError("Nothing to draw.")

    def init_graph_to_file(self, dest=""):
        '''Saves image of graph generated by process step to file.

        Parameters
        ----------
        dest: str
            Destination for writing to file.
        '''
        return self._graph_to_file(self.init_graph, dest)

    def paths_graph_to_file(self, dest=""):
        '''Saves image of graph generated by pathways step to file.

        Parameters
        ----------
        dest:str
            Destination for writing to file.
        '''
        return self._graph_to_file(self.paths_graph, dest)

    def get_surface_exposed_residues(self):
        '''Returns list of surface exposed residues.

        Returns
        -------
        surface_exposed: list of str
            List of surface exposed residues identified by pyemap
        '''
        if self.init_graph:
            surface_exposed = []
            for n, d in self.init_graph.nodes(data=True):
                if d['shape'] == "box":
                    surface_exposed.append(n)
            return surface_exposed
        else:
            raise RuntimeError("No graph found. Please run pyemap.process(my_emap) to generate the graph.")

    def _report_header(self):
        full_str = ""
        full_str += "Generated:\n" + str(datetime.datetime.now()) + "\n"
        full_str += "Parameters:\n{}\n".format(str(self._process_params))
        full_str += "Included residues:\n"
        full_str += str(self._include_residues) + "\n"
        full_str += "Active chains:\n"
        full_str += str(list(self.active_chains.keys())) + "\n"
        full_str += "Included non protein moieties:\n"
        full_str += str(list(self.eta_moieties.keys())) + "\n"
        custom_res_atms = []
        custom_res_names = []
        for key, val in self.user_residues.items():
            custom_res_names.append(key)
            custom_res_atms.append([atm.serial_number for atm in val])
        custom_res_dict = dict(zip(custom_res_names, custom_res_atms))
        full_str += "User defined residues:\n{}\n".format(custom_res_dict)
        return full_str

    def report(self, dest=""):
        '''Returns report of most probable pathways. Writes to file if destination is specified.

        Parameters
        -----------
        dest: str, optional
            Destination for writing to file

        Returns
        --------
        output: str
            Formatted report of pathways found

        Raises
        -------
        RuntimeError
            Nothing to report
        '''
        if self.branches:
            output = self._report_header() + "\nPathways:\n"
            for br in self.branches.values():
                output += str(br) + "\n"
            if dest:
                fi = open(dest, "w")
                fi.write(output)
                return output
            else:
                return output
        else:
            raise RuntimeError("Nothing to report.")

    def _graph_to_Image(self, G):
        '''Returns PIL image of graph

        Returns
        --------
        img: :class:`PIL.Image.Image`
        '''
        if G:
            fout = tempfile.NamedTemporaryFile(suffix=".png")
            agraph = to_agraph(G)
            agraph.graph_attr.update(ratio=1.0, overlap="ipsep", mode="ipsep", splines="true")
            try:
                agraph.layout(prog='neato', args="-Gepsilon=0.01 -Gmaxiter=50")
            except Exception as e:
                raise RuntimeError("There was a problem with Graphviz. See https://graphviz.gitlab.io/") from e
            if agraph.number_of_nodes() <= 200:
                agraph.draw(fout.name, prog='neato')
            else:
                agraph.draw(fout.name, prog='dot')
            img = Image.open(fout.name)
            return img
        else:
            raise RuntimeError("Nothing to draw.")

    def init_graph_to_Image(self):
        '''Returns PIL image of initial graph

        Returns
        --------
        img: :class:`PIL.Image.Image`
        '''
        return self._graph_to_Image(self.init_graph)

    def paths_graph_to_Image(self):
        '''Returns PIL image of pathways graph
        
        Returns
        --------
        img: :class:`PIL.Image.Image`
        '''
        return self._graph_to_Image(self.paths_graph)
