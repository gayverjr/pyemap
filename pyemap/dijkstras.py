# This Python script is a part of
# eMAP: online mapping of electron transfer channels in biomolecules
# Copyright(C) 2017-2018 Ruslan Tazhigulov, James Gayvert, Melissa Wei, Ksenia Bravaya (Boston University, USA)
"""Finds shortest paths in graph given a source and optionally a target node.

Defines implementations of yen's and dijkstra's algorithms for calculating the shortest path(s) from 
the source to target/surface exposed residues. Also defines ShortestPath and Branch objects for 
organizing the pathways based on their distances and first surface exposed residue reached during
the pathway.

"""
import itertools
import os
import string
import sys
import networkx as nx
import numpy as np


class ShortestPath(object):
    """Data structure used to store shortest paths.

    Contains functions for comparison and string representation, both for the user output and
    in NGL viewer selection language for the purpose of visualization. Sorting for ShortestPath
    objects is done based on the length attribute.

    Parameters
    ----------
    path_id: array-like
        List of residues that make up the shortest path
    length: float
        Total distance from source to target

    Attributes
    ----------
    path_id: array-like
        List of residues that make up the shortest path
    length: float
        Total distance from source to target
    id: string
        Unique identifier for pathway assigned based on branch and length.

    """

    def __init__(self, path, length):
        self.path = path
        self.length = length
        self.path_id = "none"

    def __eq__(self, other):
        return self.length == other.length

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.length < other.length

    def __str__(self):
        # If only one chain is present, do not include chain id in output
        printline = self.path_id + ": " + \
            str(self.path) + " " + str('{:.2f}'.format(round(self.length, 2)))
        return printline

    def get_path_as_list(self):
        original_list = [[self.path_id], self.path, [str('{:.2f}'.format(round(self.length, 2)))]]
        merged = list(itertools.chain(*original_list))
        return merged

    def set_id(self, path_id):
        """Setter for path_id"""
        self.path_id = path_id

    def selection_str(self, customname, customnum):
        """NGL selection language representation of this pathway.

        This function returns a string that is a representation of the pathway in
        the NGL viewer selection language, which gets written to file for later access
        by the front end.

        Parameters
        ----------
        customname: array-like
            List of custom residues
        customnum: array-like
            List of atom serial numbers associated with the custom residue with corresponding
            index in customname.

        Returns
        -------
        select_arr: array-like
            Returns an array which is used to represent pathway in NGL. Each node in pathway corresponds to
            three consecutive entries in array:
            Entry 1: NGL string specifying atom by atom
            Entry 2: Color in representation
            Entry 3: Which atom to place the residue label on in representation

        """
        colors = {"F": "orange", "Y": "blue", "W": "red", "H": "green"}
        # each type of residue has a corresponding color
        select_arr = []
        # custom residues must have each atom specified individually
        for residue in self.path:
            cust = False
            for name in customname:
                if name in residue:
                    select_arr.append(customnum[customname.index(name)])
                    select_arr.append("pink")
                    atm_name_index = str.index(customnum[customname.index(name)], ".")
                    end_idx = str.index(customnum[customname.index(name)], ")", atm_name_index)
                    select_arr.append(customnum[customname.index(name)][atm_name_index:end_idx])
                    cust = True
            if not cust:
                # standard residues have only .CA specified
                res_type = residue[0]
                start = str.index(str(residue), "(")
                end = str.index(str(residue), ")")
                chain_id = residue[start + 1:end]
                out = residue[1:start]
                out += " and :" + chain_id
                select_arr.append(out)
                select_arr.append(colors.get(res_type))
                select_arr.append(".CA")
        return select_arr


class Branch(object):
    """Data structure used to group shortest paths with a common first surface exposed residue.

    Shortest paths are classified into branches based on the first surface exposed residue
    reached during the course of the pathway. For example, given a source node A, and
    surface exposed nodes B and C, the paths [A,F,B] and [A,F,B,C] will both be part of the
    "B" branch. The path [A,E,C] would be part of its own 'C" branch.

    Notes
    ----
    This class is only used if no target is specified by the user.

    Parameters
    ----------
    branch_id: str
        Unique identifier for a branch
    target: str
        Target node which a branch corresponds to

    Attributes
    ----------
    branch_id: str
        Unique identifier for a branch
    target: str
        Target node which a branch corresponds to
    paths: array-like
        List of ShortestPath objects that make up a branch

    See Also
    --------
    class ShortestPath

    """

    def __init__(self, branch_id, target):
        self.branch_id = branch_id
        self.target = target
        self.paths = []

    def add_path(self, path):
        """Adds a path to the branch and sets the path_id.

        Paths are added to a branch if self.target is the first surface exposed residue reached
        during the pathway. Each time a path is added, the paths in the branch are sorted. After sorting,
        each path is assigned a path id composed of the branch id and its location in the paths list. For
        example, the shortest path in branch 12 would be assigned the id '12a', the second shortest '12b' and
        so on.

        Parameters
        ----------
        path: ShortestPath
            A ShortestPath from source to a surface exposed residue

        """

        letters = list(string.ascii_letters)
        if len(self.paths) < len(letters):
            self.paths.append(path)
            self.paths = sorted(self.paths)
            i = 0
            while i < len(self.paths) and i < len(letters):
                self.paths[i].set_id(str(self.branch_id) + str(letters[i]))
                i += 1

    def __str__(self):
        printline = "Branch: " + str(self.target)
        for pt in self.paths:
            printline += "\n" + str(pt)
        printline += "\n"
        return printline

    def get_branch_as_list(self):
        branch_list = []
        printline = "Branch: " + str(self.target)
        branch_list.append([printline])
        for pt in self.paths:
            branch_list.append(pt.get_path_as_list())
        return branch_list


def is_parent_pathway(shortest_path, targets):
    """Returns true if ShortestPath is a parent pathway, false if not.

    A ShortestPath object is the parent of a branch if its terminal residue is the
    only surface exposed residue in the path. For example, if targets=[A,B,C] and
    the pathway is [H,I,C], then this pathway is a parent pathway. In contrast, if
    the pathway is [H,B,A], then this pathway is not a parent pathway.

    Parameters
    ----------
    shortest_path: ShortestPath
        ShortestPath object
    targets: array-like
        List of surface exposed residues

    Returns
    -------
    True if a parent pathway
    False if not a parent pathway
    """

    count = 0
    for res in shortest_path.path:
        if res in targets:
            count += 1
    return count == 1


def find_branch(pt, targets, branches):
    """Determines which branch a pathway belongs to and returns that branch.

    A ShortestPath belongs to a branch if the first surface exposed residue it reaches during the
    pathway is the target of that particular branch.

    Parameters
    ----------
    pt: ShortestPath
        ShortestPath object
    targets: array-like
        List of surface exposed residues
    branches: array-like
        List of branches

    Returns
    -------
    cur_branch: Branch
        Branch object that pt belongs to

    """
    res = pt.path[0]
    count = 0
    while res not in targets:
        count += 1
        res = pt.path[count]
    count = 0
    cur_branch = branches[0]
    while not res == cur_branch.target:
        count += 1
        cur_branch = branches[count]
    return cur_branch


def dijkstras_shortest_paths(G, start, targets):
    """Returns shortest path from source to each surface exposed residue.

    Performs Dijkstra's algorithm from the source to each surface exposed residue, finding the
    shortest path. The ShortestPath objects are organized into branches based on the first surface
    exposed residue reached during the course of the pathway. 

    Parameters
    ----------
    G: NetworkX graph object
        Undirected, weighted residue graph
    start: str
        Source node
    targets: array-like
        List of surface exposed residues

    Returns
    -------
    shortest_paths: array-like
        List of ShortestPath objects representing pathways found by eMap

    See Also
    --------
    module NetworkX.dijkstra_path
    class ShortestPath
    class Branch

    Raises
    ------
    Exception e:
        No shortest paths to surface found

    """
    shortestPaths = []
    for goal in targets:
        path = []
        try:
            path = nx.dijkstra_path(G, start, goal)
        except:
            path = []
        if not path == []:
            sum = 0
            for i in range(0, len(path) - 1):  # sum up edge weights
                sum += (G[path[i]][path[i + 1]]['weight'])
            shortestPaths.append(ShortestPath(path, sum))
    shortestPaths = sorted(shortestPaths)
    branches = []
    # find the parent pathways
    for pt in shortestPaths:
        if is_parent_pathway(pt, targets):
            path = pt.path
            for i in range(0, len(path) - 1):
                G[path[i]][path[i + 1]]['color'] = '#778899FF'
                G[path[i]][path[i + 1]]['penwidth'] = 6.0
                G[path[i]][path[i + 1]]['style'] = 'solid'
                G.node[path[i]]['penwidth'] = 6.0
                G.node[path[i + 1]]['penwidth'] = 6.0
                # make the nodes look opaque if they are connected to the source
                if len(G.node[path[i]]['fillcolor']) != 9:
                    G.node[path[i]]['fillcolor'] += 'FF'
                    G.node[path[i]]['color'] = '#708090FF'
                if len(G.node[path[i + 1]]['fillcolor']) != 9:
                    G.node[path[i + 1]]['fillcolor'] += 'FF'
                    G.node[path[i + 1]]['color'] = '#708090FF'
            br = Branch(len(branches) + 1, pt.path[-1])
            branches.append(br)
            br.add_path(pt)
    # find the sub pathways
    for pt in shortestPaths:
        if not is_parent_pathway(pt, targets):
            find_branch(pt, targets, branches).add_path(pt)
            path = pt.path
            for i in range(0, len(path) - 1):
                if G[path[i]][path[i + 1]]['color'] != '#778899FF':
                    G[path[i]][path[i + 1]]['color'] = '#7788995F'
                G[path[i]][path[i + 1]]['penwidth'] = 6.0
                G[path[i]][path[i + 1]]['style'] = 'solid'
                G.node[path[i]]['penwidth'] = 6.0
                G.node[path[i + 1]]['penwidth'] = 6.0
                # make the nodes look opaque if they are connected to the source
                if len(G.node[path[i]]['fillcolor']) != 9:
                    G.node[path[i]]['fillcolor'] += '5F'
                    G.node[path[i]]['color'] = '#7080905F'
                if len(G.node[path[i + 1]]['fillcolor']) != 9:
                    G.node[path[i + 1]]['fillcolor'] += '5F'
                    G.node[path[i + 1]]['color'] = '#7080905F'
    if len(shortestPaths) == 0:
        raise Exception("No paths to the surface found.")
    return shortestPaths


def yens_shortest_paths(G, start, target):
    """Returns top 5 shortest paths from source to target.

    Uses Yen's algorithm to calculate the 5 shortest paths from source to target, writes
    out the ShortestPath objects to file, and returns the 5 pathway IDs. In the graph, nodes and
    edges that are part of any pathways are made opaque, and the shortest path is highlighted.


    Parameters
    ----------
    G: NetworkX graph object
        Undirected, weighted residue graph
    start: str
        Source node
    target: str
        Target node

    Returns
    -------
    shortest_paths: array-like
        List of ShortestPath objects representing pathways found by eMap

    See Also
    --------
    module NetworkX.shortest_simple_paths
    class ShortestPath

    References
    ----------
    Jin Y. Yen, Finding the K Shortest Loopless Paths in a Network, Management Science,
    Vol. 17, No. 11, Theory Series (Jul., 1971), pp. 712-716.

    Raises
    ------
    Exception e:
        No shortest paths to target found

    """
    letters = list(string.ascii_letters)
    shortestPaths = []
    k = 0
    from itertools import islice
    paths = list(islice(nx.shortest_simple_paths(G, start, target), 5))
    for k in range(0, len(paths)):
        path = paths[k]
        sum = 0
        for i in range(0, len(path) - 1):  # sum up edge weights
            sum += (G[path[i]][path[i + 1]]['weight'])
        path = ShortestPath(path, sum)
        shortestPaths.append(path)
    if shortestPaths:
        shortestPaths = sorted(shortestPaths)
        for i in range(0, len(shortestPaths)):
            path = shortestPaths[i].path
            if i == 0:  # shortest path gets bolder edges
                for j in range(len(path) - 1):
                    G[path[j]][path[j + 1]]['penwidth'] = 6.0
                    G[path[j]][path[j + 1]]['style'] = 'solid'
                    G.node[path[j]]['penwidth'] = 6.0
                    G.node[path[j + 1]]['penwidth'] = 6.0
                    G[path[j]][path[j + 1]]['color'] = '#778899FF'
                    # make the nodes look opaque if they are connected to the source
                    if len(G.node[path[j]]['fillcolor']) != 9:
                        G.node[path[j]]['fillcolor'] += 'FF'
                        G.node[path[j]]['color'] = '#708090FF'
                    if len(G.node[path[j + 1]]['fillcolor']) != 9:
                        G.node[path[j + 1]]['fillcolor'] += 'FF'
                        G.node[path[j + 1]]['color'] = '#708090FF'
            else:
                for j in range(len(path) - 1):
                    G[path[j]][path[j + 1]]['penwidth'] = 6.0
                    G[path[j]][path[j + 1]]['style'] = 'solid'
                    G.node[path[j]]['penwidth'] = 6.0
                    G.node[path[j + 1]]['penwidth'] = 6.0
                    if G[path[j]][path[j + 1]]['color'] != '#778899FF':
                        G[path[j]][path[j + 1]]['color'] = '#7788997F'
                    # make the nodes look opaque if they are connected to the source
                    if len(G.node[path[j]]['fillcolor']) != 9:
                        G.node[path[j]]['fillcolor'] += '7F'
                        G.node[path[j]]['color'] = '#7080907F'
                    if len(G.node[path[j + 1]]['fillcolor']) != 9:
                        G.node[path[j + 1]]['fillcolor'] += '7F'
                        G.node[path[j + 1]]['color'] = '#7080907F'
            shortestPaths[i].set_id("1" + letters[i])
        return shortestPaths
    else:  # no paths found
        raise Exception("No paths to target found.")
