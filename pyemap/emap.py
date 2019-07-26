from rdkit import Chem
from rdkit.Chem import Draw
from collections import defaultdict
from heapq import heappop, heappush
from .custom_residues import is_pi_bonded,dist
from .data import *
import networkx as nx

def buildSmiles(graph, atoms, cur, prev):
    """Two pass depth first search algorithm on chemical graph to generate smiles string.

    Parameters
    ----------
    graph: NetworkX chemical graph
    atoms: array-like
        List of atoms in structure
    cur: BioPython Atom object
    prev: BioPython Atom object
        Atoms being considered at this step of the iteration

    References
    ----------
    Varnek, A. Tutorials in Chemoinformatics; John Wiley & Sons, Inc.: Hoboken, NJ, 2017.

    """
    visited.add(cur)
    seq = ''
    seq += atoms[cur].element.lower()
    for d in closingClosures[cur]:
        seq += d
        heappush(digits, d[-1])
    for a in openingClosures[cur]:
        d = str(heappop(digits))
        seq += d
        closingClosures[a].append(d)
    branches = []
    neighbors = list(graph.neighbors(cur))
    if prev in neighbors:
        neighbors.remove(prev)
    for neighbor in neighbors:
        if neighbor not in visited:
            branches.append(buildSmiles(graph, atoms, neighbor, cur))
    for branch in branches[:-1]:
        seq += "(" + branch + ")"
    if len(branches) > 0:
        seq += branches[-1]
    return seq


def getClosures(graph, cur, prev):
    """Two pass depth first search algorithm on chemical graph to generate smiles string.

    graph: NetworkX chemical graph
    cur: BioPython Atom object
    prev: BioPython Atom object
        Atoms being considered at this step of the iteration

    References
    ----------
    Varnek, A. Tutorials in Chemoinformatics; John Wiley & Sons, Inc.: Hoboken, NJ, 2017.

    """
    ancestor.add(cur)
    visited.add(cur)
    neighbors = list(graph.neighbors(cur))
    if prev in neighbors:
        neighbors.remove(prev)
    for neighbor in neighbors:
        if neighbor in ancestor:
            openingClosures[neighbor].append(cur)
        elif neighbor not in visited:
            getClosures(graph, neighbor, cur)
    ancestor.remove(cur)


def getSimpleSmiles(graph, atoms):
    """Two pass depth first search algorithm on chemical graph to generate smiles string.

    Parameters
    ----------
    graph: NetworkX chemical graph
    atoms: array-like
        List of atoms in structure

    References
    ----------
    Varnek, A. Tutorials in Chemoinformatics; John Wiley & Sons, Inc.: Hoboken, NJ, 2017.

    """
    root = list(graph.nodes())[0]
    global visited, ancestor, openingClosures, closingClosures, digits
    visited = set()
    ancestor = set()
    openingClosures = defaultdict(list)
    getClosures(graph, root, None)
    closingClosures = defaultdict(list)
    digits = [str(x) for x in range(1, 10)]
    visited = set()
    return buildSmiles(graph, atoms, root, None)

class emap():
    def __init__(self,filename,structure,custom_residues,chain_list):
        self.filename = filename
        self.structure=structure
        self.smiles_dict={}
        self.chains = chain_list
        self.custom_residues = custom_residues
        for residue in custom_residues:
            self.add_residue(residue)
    def save_agraph(self,agraph):
        self.agraph = agraph
    def save_paths(self):
        pass
    def view_residue(self,resname):
        #Just writes to file, would prefer if it displays to GUI
        print(self.smiles_dict)
        mol = Chem.MolFromSmarts(self.smiles_dict.get(resname))
        Draw.MolToFile(mol, resname+".png", kekulize=False, size=(200, 200))
    def add_residue(self,residue):
        atoms=list(residue.get_atoms())
        #creating chemical graph structure
        arom_atoms = ['O', 'P', 'N', 'C','S']  # only these elements will be considered in our search
        res_graph = nx.Graph()
        for i in range(len(atoms)):
            for k in range(i, len(atoms)):
                if (not i == k) and is_pi_bonded(atoms[i], atoms[k]):
                    if atoms[i].element in arom_atoms and atoms[k].element in arom_atoms:
                        res_graph.add_edge(i, k)
        if len(res_graph.edges())>0:
            smiles_str = getSimpleSmiles(res_graph, atoms)
            molecule = Chem.MolFromSmarts(smiles_str)
            can_smiles_str = Chem.MolToSmarts(molecule, True)
            self.smiles_dict[residue.resname] = can_smiles_str
        #else TODO











