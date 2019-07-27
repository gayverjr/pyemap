from heapq import heappop, heappush
from collections import defaultdict


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
