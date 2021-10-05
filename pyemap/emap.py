# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
from .custom_residues import is_pi_bonded
from .data import clusters
import networkx as nx
from networkx.drawing.nx_agraph import to_agraph
from .structures import getSimpleSmarts
from collections import OrderedDict
from PIL import Image
import os
from shutil import copyfile
import tempfile
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM
from .process_data import get_atom_list
from .png2svg import convert_to_svg


class emap():
    '''
    Manages the data generated at all stages of PyeMap analysis.

    Attributes
    ----------
    filename: str
        Name of the crystal structure file being processed by PyeMap.
    structure: :class:`Bio.PDB.Structure.Structure`
        Macromolecular protein structure. Contains model which contains residues and atoms.
    eta_moieties: dict of node label(str): :class:`Bio.PDB.Residue.Residue`
        Non-protein eta moieties automatically identified at the parsing step.
    chain_list: list of str
        List of chains identified at the parsing step.
    smarts: dict of node label(str): smarts string(str)
        List of SMARTS strings for non-protein eta moieties identified at the parsing step.
    residues: dict of node label(str): :class:`Bio.PDB.Residue.Residue`
        Residues included in the graph after the process step.
    ngl_strings: dict of nodel label(str):ngl_string(str)
        Formatted NGL viewer selection strings for residues included in the graph after the process step.
    user_residues: dict of str: :class:`Bio.PDB.Residue.Residue`
        Custom residues specified by the user.
    init_graph: :class:`networkx.Graph`
        Graph generated after the process step.
    branches: dict of branch id(int): :class:`~pyemap.Branch`
        Branches found by PyeMap analysis
    paths: dict of path id(str): :class:`~pyemap.ShortestPath`
        Paths found by PyeMap sorted by lowest to highest score.
    paths_graph: :class:`networkx.Graph`
        Graph generated after the shortest paths step.
    '''

    def __init__(self, filename, structure, eta_moieties, chain_list):
        '''Initializes emap object.

        Parameters
        ----------
        filename: str
            Name of file
        structure: :class:`Bio.PDB.Structure.Structure`
            Macromolecular protein structure
        eta_moieties: list of :class:`Bio.PDB.Residue.Residue`
            Customized residue objects generated for automatically detected eta moieties
        chain_list: list of str
            Chains identified by the parser
        '''
        self.filename = filename
        self.pdb_id = filename[-8:-4].upper()
        self.structure = structure
        self.residues = {}
        self.chains = chain_list
        self.eta_moieties = {}
        self.user_residues = {}
        self.smarts = {}
        self.paths = OrderedDict()
        self.paths_graph = []
        self.init_graph = []
        self.ngl_strings = {}
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
        self.init_graph = graph

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
        self.paths_graph = graph

    def _reset_process(self):
        '''Returns emap object to state it was in after parsing.
        '''
        self.residues = {}
        self.user_residues = {}
        self.init_graph = []
        self.paths = OrderedDict()
        self.paths_graph = []
        self.branches = OrderedDict()
        self.ngl_strings = {}

    def _reset_paths(self):
        '''Returns emap object to state it was in after the process step.
        '''
        self.paths = OrderedDict()
        self.paths_graph = []
        self.branches = OrderedDict()

    def _get_residue_graph(self, residue):
        atoms = get_atom_list(residue)
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
        '''Gets the SMARTS string for an automatically identified non-protein eta moiety,
        and adds the residue to the eta_moieties dictionary.

        Parameters
        ----------
        residue: Bio.PDB.Residue.Residue
             Customized residue object generated for automatically detected eta moiety
        '''
        if not residue.resname[:3] in clusters and "CUST" not in residue.resname:
            res_graph = self._get_residue_graph(residue)
            smarts_str = getSimpleSmarts(res_graph)
            residue.smarts = smarts_str
            self.smarts[residue.resname] = smarts_str
        self.eta_moieties[residue.resname] = residue

    def _visualize_pathway(self, pathway, yens):
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
        pathway.set_visualization(
            selection_strs, color_list, labeled_atoms, label_texts)

    def _add_residue(self, residue):
        '''Gets ngl string for residue, and adds the residue to the residues and ngl_strings dictionaries.
        '''
        if not residue.resname[:3] in clusters and "CUST" not in residue.resname and residue.resname[:3] in ["PHE","HIS","TYR","TRP"]:
            res_graph = self._get_residue_graph(residue)
            smarts_str = getSimpleSmarts(res_graph)
            residue.Smarts = smarts_str
            self.smarts[residue.node_label] = smarts_str
        residue.ngl_string = self._get_ngl_string(residue)
        self.residues[residue.node_label] = residue
        self.ngl_strings[residue.node_label] = residue.ngl_string

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
        if hasattr(first_atm,"original_id") and first_atm.original_id is not None:
            id = first_atm.original_id
        else:
            id = first_atm.full_id
        select_string += "(" + str(id[3][1]) + " and :" + str(
            id[2]) + " and ^" + residue.id[2] + " and ." + first_atm.name + ")"
        for i in range(1, len(atm_list)):
            atm = atm_list[i]
            if hasattr(atm,"original_id") and atm.original_id is not None:
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
            from rdkit.Chem import Draw
        except (ModuleNotFoundError, ImportError) as e:
            raise ModuleNotFoundError(
                "RDKit is required for visualization of chemical structures. See https://www.rdkit.org/docs/Install.html"
            ) from e
        if resname in self.residues or resname in self.eta_moieties:
            if "CUST" in resname:
                raise KeyError("Not available for user defined residues.")
            elif resname[:3] in clusters:
                cluster_img_name = os.path.abspath(os.path.dirname(
                    __file__)) + '/data/clusters/' + resname[:3] + '.svg'
                if dest:
                    target_name = dest
                else:
                    target_name = resname + ".svg"
                copyfile(cluster_img_name, target_name)
            else:
                try:
                    mol = Chem.MolFromSmarts(self.smarts[resname])
                    mol.UpdatePropertyCache()
                    if dest:
                        Draw.MolToFile(mol, dest, kekulize=False, size=size)
                    else:
                        Draw.MolToFile(mol, resname + ".svg",
                                    kekulize=False, size=size)
                except Exception as e:
                    print("Warning: RDKIT couldn't draw: " + dest)
                    print("Trying to convert png to svg...")
                    try:
                        if dest[-4:]==".svg":
                            Draw.MolToFile(mol, dest[:-4]+".png", kekulize=False, size=size)
                            convert_to_svg(dest[:-4]+".png",dest)
                    except Exception as e:
                        print(e)
                        print("Warning: couldn't convert png to svg.")
                            
        else:
            raise KeyError("No record of any residue by that name.")

    def residue_to_Image(self, resname, size=(100, 100)):
        '''Returns PIL image of chemical structure

        Parameters
        -----------
        resname: str
            Name of residue
        size: (float,float), optional
            dimensions of image saved to file
        Returns
        --------
        img: :class:`PIL.Image.Image`
        '''
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
        except (ModuleNotFoundError, ImportError) as e:
            raise ModuleNotFoundError(
                "RDKit is required for visualization of chemical structures. See https://www.rdkit.org/docs/Install.html"
            ) from e
        if resname in self.residues or resname in self.eta_moieties:
            if "CUST" in resname:
                raise KeyError("Not available for user defined residues.")
            elif resname[:3] in clusters:
                cluster_img_name = os.path.abspath(os.path.dirname(
                    __file__)) + '/data/clusters/' + resname[:3] + '.svg'
                drawing = svg2rlg(cluster_img_name)
                img = renderPM.drawToPIL(drawing)
                return img
            else:
                try:
                    mol = Chem.MolFromSmarts(self.smarts[resname])
                    mol.UpdatePropertyCache()
                    img = Draw.MolToImage(mol, kekulize=False, size=size)
                    return img
                except Exception as e:
                    print("Warning: RDKIT couldn't draw: " + dest)
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
            agraph.graph_attr.update(
                ratio=1.0, overlap="ipsep", mode="ipsep", splines="true")
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
                agraph.draw(svg_fn,prog='neato',args="-Gepsilon=0.01 -Gmaxiter=50")
                agraph.draw(png_fn,prog='neato',args="-Gepsilon=0.01 -Gmaxiter=50")
            else:
                png_fn = self.filename[:-4] + "_graph.png"
                agraph.draw(png_fn,prog='neato',args="-Gepsilon=0.01 -Gmaxiter=50")
        else:
            raise RuntimeError("Nothing to draw.")


    def init_graph_to_file(self, dest=""):
        '''Saves image of graph generated by process step to file.

        Parameters
        ----------
        dest: str
            Destination for writing to file.
        '''
        return self._graph_to_file(self.init_graph,dest)

    def paths_graph_to_file(self, dest=""):
        '''Saves image of graph generated by pathways step to file.

        Parameters
        ----------
        dest:str
            Destination for writing to file.
        '''
        return self._graph_to_file(self.paths_graph,dest)

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
            raise RuntimeError(
                "No graph found. Please run pyemap.process(my_emap) to generate the graph.")

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
            output = ""
            for br in self.branches.values():
                output += str(br)
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

