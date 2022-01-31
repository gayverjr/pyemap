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

"""Finds shortest paths in graph given a source and optionally a target node.

Defines implementations of yen's and dijkstra's algorithms for calculating the shortest path(s) from
the source to target/surface exposed residues. Also defines ShortestPath and Branch objects for
organizing the pathways based on their distances and first surface exposed residue reached during
the pathway.

"""
import itertools
import string
import networkx as nx
from functools import total_ordering
from .pyemap_exceptions import *
import numpy as np


@total_ordering
class ShortestPath(object):
    """Data structure used to store shortest paths.

    Contains functions for comparison and string representation, both for the user output and
    in NGL viewer selection language for the purpose of visualization. Sorting for ShortestPath
    objects is done based on the length attribute.

    Attributes
    ----------
    path: list of str
        List of residue names that make up the shortest path
    edges: list of float
        List of edge weights that make up the shortest path
    path_id: list of str
        List of residues that make up the shortest path
    length: float
        Total distance from source to target
    selection_strs: list of str
        NGL selection strings for visualization
    color_list: list of str
        Colors of residues in visualization
    labeled_atoms: list of str
        Atom names which are labeled in NGL visualization
    label_texts: list of str
        Labels of residues in NGL visualization
    """

    def __init__(self, path, edges, length):
        '''Initializes ShortestPath object.

        Parameters
        ----------
        path: list of str
            List of residues that make up the shortest path
        length: float
            Total distance from source to target
        '''
        self.path = path
        self.length = length
        self.edges = edges
        self.path_id = "none"
        self.selection_strs = []
        self.color_list = []
        self.labeled_atoms = []
        self.label_texts = []

    def __eq__(self, other):
        return self.length == other.length

    def __lt__(self, other):
        return self.length < other.length

    def __str__(self):
        rounded_edges = [np.round(x, 2) for x in self.edges]
        printline = self.path_id + ": " + \
            str(self.path) + " " + str('{:.2f}'.format(round(self.length, 2))) + "\n" + \
            "Edge weights: {}".format(str(rounded_edges))
        return printline

    def get_path_as_list(self):
        original_list = [[self.path_id], self.path, [str('{:.2f}'.format(round(self.length, 2)))]]
        merged = list(itertools.chain(*original_list))
        return merged

    def set_id(self, path_id):
        """Setter for path_id"""
        self.path_id = path_id

    def set_visualization(self, selection_strs, color_list, labeled_atoms, label_texts):
        '''Saves information needed for NGL visualization.'''
        self.selection_strs = selection_strs
        self.color_list = color_list
        self.labeled_atoms = labeled_atoms
        self.label_texts = label_texts


class Branch(object):
    """Data structure used to group shortest paths with a common first surface exposed residue.

    Shortest paths are classified into branches based on the first surface exposed residue
    reached during the course of the pathway. For example, given a source node A, and
    surface exposed nodes B and C, the paths [A,F,B] and [A,F,B,C] will both be part of the
    "B" branch. The path [A,E,C] would be part of its own 'C" branch.

    Attributes
    ----------
    branch_id: str
        Unique identifier for a branch
    target: str
        Target node which a branch corresponds to
    paths: list of :class:`~pyemap.ShortestPath`
        List of ShortestPath objects that make up a branch
    """

    def __init__(self, branch_id, target):
        self.branch_id = branch_id
        self.target = target
        self.paths = []

    def add_path(self, path):
        """Adds a path to the branch and sets the path_id.

        Each time a path is added, the paths in the branch are sorted. After sorting, each path is assigned a
        path id composed of the branch id and its location in the paths list. For example, the shortest path
        in branch 12 would be assigned the id '12a', the second shortest '12b' and so on.

        Parameters
        ----------
        path: :class:`~pyemap.ShortestPath`
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
        ''' String representation of Branch. First line is: "Branch: `branch_id`" and subsequent lines
        are the string representations of each ShortestPath object comprising the Branch.
        '''
        printline = "Branch: " + str(self.target)
        for pt in self.paths:
            printline += "\n" + str(pt)
        printline += "\n"
        return printline


def _is_parent_pathway(shortest_path, targets):
    """Returns true if ShortestPath is a parent pathway, false if not.

    A ShortestPath object is the parent of a branch if its terminal residue is the
    only surface exposed residue in the path. For example, if targets=[A,B,C] and
    the pathway is [H,I,C], then this pathway is a parent pathway. In contrast, if
    the pathway is [H,B,A], then this pathway is not a parent pathway.

    Parameters
    ----------
    shortest_path: ShortestPath
        ShortestPath object
    targets: list of str
        List of surface exposed residues

    Returns
    -------
    bool
        True if path is ShortestPath object is a parent pathway
    """

    count = 0
    for res in shortest_path.path:
        if res in targets:
            count += 1
    return count == 1


def _find_branch(pt, targets, branches):
    """Determines which branch a pathway belongs to and returns that branch.

    A ShortestPath belongs to a branch if the first surface exposed residue it reaches during the
    pathway is the target of that particular branch.

    Parameters
    ----------
    pt: ShortestPath
        ShortestPath object
    targets: list of str
        List of surface exposed residues
    branches: list of pyemapBranch objects
        List of branches already found

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
    G: :class:`networkx.Graph`
        Undirected, weighted residue graph
    start: str
        Source node
    targets: list of str
        List of surface exposed residues

    Returns
    -------
    branches: list of :class:`~pyemap.Branch`
        A list of Branch objects representing the groups of pathways found

    Raises
    ------
    RuntimeError:
        No shortest paths to surface found

    """
    shortestPaths = []
    for goal in targets:
        path = []
        try:
            path = nx.dijkstra_path(G, start, goal)
            weights = []
            sum = 0
            for i in range(0, len(path) - 1):  # sum up edge weights
                sum += (G[path[i]][path[i + 1]]['weight'])
                weights.append(G[path[i]][path[i + 1]]['weight'])
            shortestPaths.append(ShortestPath(path, weights, sum))
        except Exception as e:
            pass
    if len(shortestPaths) == 0:
        raise PyeMapShortestPathException("No paths to the surface from " + str(start) + " were found.")
    shortestPaths = sorted(shortestPaths)
    branches = []
    # find the parent pathways
    for pt in shortestPaths:
        if _is_parent_pathway(pt, targets):
            path = pt.path
            for i in range(0, len(path) - 1):
                G[path[i]][path[i + 1]]['color'] = '#778899FF'
                G[path[i]][path[i + 1]]['penwidth'] = 6.0
                G[path[i]][path[i + 1]]['style'] = 'solid'
                G.nodes[path[i]]['penwidth'] = 6.0
                G.nodes[path[i + 1]]['penwidth'] = 6.0
                # make the nodes look opaque if they are connected to the source
                if len(G.nodes[path[i]]['fillcolor']) != 9:
                    G.nodes[path[i]]['fillcolor'] += 'FF'
                    G.nodes[path[i]]['color'] = '#708090FF'
                if len(G.nodes[path[i + 1]]['fillcolor']) != 9:
                    G.nodes[path[i + 1]]['fillcolor'] += 'FF'
                    G.nodes[path[i + 1]]['color'] = '#708090FF'
            br = Branch(len(branches) + 1, pt.path[-1])
            branches.append(br)
            br.add_path(pt)
    # find the sub pathways
    for pt in shortestPaths:
        if not _is_parent_pathway(pt, targets):
            _find_branch(pt, targets, branches).add_path(pt)
            path = pt.path
            for i in range(0, len(path) - 1):
                if G[path[i]][path[i + 1]]['color'] != '#778899FF':
                    G[path[i]][path[i + 1]]['color'] = '#7788995F'
                G[path[i]][path[i + 1]]['penwidth'] = 6.0
                G[path[i]][path[i + 1]]['style'] = 'solid'
                G.nodes[path[i]]['penwidth'] = 6.0
                G.nodes[path[i + 1]]['penwidth'] = 6.0
                # make the nodes look opaque if they are connected to the source
                if len(G.nodes[path[i]]['fillcolor']) != 9:
                    G.nodes[path[i]]['fillcolor'] += '5F'
                    G.nodes[path[i]]['color'] = '#7080905F'
                if len(G.nodes[path[i + 1]]['fillcolor']) != 9:
                    G.nodes[path[i + 1]]['fillcolor'] += '5F'
                    G.nodes[path[i + 1]]['color'] = '#7080905F'
    return branches


def yens_shortest_paths(G, start, target, max_paths=10):
    """Returns shortest paths from source to target.

    Uses Yen's algorithm to calculate the shortest paths from source to target, writes
    out the ShortestPath objects to file, and returns the 10 pathway IDs. In the graph, nodes and
    edges that are part of any pathways are made opaque, and the shortest path is highlighted.


    Parameters
    ----------
    G: :class:`networkx.Graph` object
        Undirected, weighted residue graph
    start: str
        Source node
    target: str
        Target node
    max_paths: int, optional
        Maximum number of paths to search for

    Returns
    -------
    list of :class:`~pyemap.Branch` objects
        A list of length 1 containing a single Branch object which represents the group of pathways found.

    References
    ----------
    Jin Y. Yen, Finding the K Shortest Loopless Paths in a Network, Management Science,
    Vol. 17, No. 11, Theory Series (Jul., 1971), pp. 712-716.

    Raises
    ------
    PyeMapShortestPathException:
        No shortest paths to target found.

    """
    letters = list(string.ascii_letters)
    shortestPaths = []
    k = 0
    try:
        paths = list(itertools.islice(nx.shortest_simple_paths(G, start, target), max_paths))
    except Exception:
        raise PyeMapShortestPathException("No paths between " + str(start) + " and " + str(target) + " were found.")
    for k in range(0, len(paths)):
        path = paths[k]
        sum = 0
        weights = []
        for i in range(0, len(path) - 1):  # sum up edge weights
            sum += (G[path[i]][path[i + 1]]['weight'])
            weights.append(G[path[i]][path[i + 1]]['weight'])
        path = ShortestPath(path, weights, sum)
        shortestPaths.append(path)
    shortestPaths = sorted(shortestPaths)
    for i in range(0, len(shortestPaths)):
        path = shortestPaths[i].path
        if i == 0:  # shortest path gets bolder edges
            for j in range(len(path) - 1):
                G[path[j]][path[j + 1]]['penwidth'] = 6.0
                G[path[j]][path[j + 1]]['style'] = 'solid'
                G.nodes[path[j]]['penwidth'] = 6.0
                G.nodes[path[j + 1]]['penwidth'] = 6.0
                G[path[j]][path[j + 1]]['color'] = '#778899FF'
                # make the nodes look opaque if they are connected to the source
                if len(G.nodes[path[j]]['fillcolor']) != 9:
                    G.nodes[path[j]]['fillcolor'] += 'FF'
                    G.nodes[path[j]]['color'] = '#708090FF'
                if len(G.nodes[path[j + 1]]['fillcolor']) != 9:
                    G.nodes[path[j + 1]]['fillcolor'] += 'FF'
                    G.nodes[path[j + 1]]['color'] = '#708090FF'
        else:
            for j in range(len(path) - 1):
                G[path[j]][path[j + 1]]['penwidth'] = 6.0
                G[path[j]][path[j + 1]]['style'] = 'solid'
                G.nodes[path[j]]['penwidth'] = 6.0
                G.nodes[path[j + 1]]['penwidth'] = 6.0
                if G[path[j]][path[j + 1]]['color'] != '#778899FF':
                    G[path[j]][path[j + 1]]['color'] = '#7788997F'
                # make the nodes look opaque if they are connected to the source
                if len(G.nodes[path[j]]['fillcolor']) != 9:
                    G.nodes[path[j]]['fillcolor'] += '7F'
                    G.nodes[path[j]]['color'] = '#7080907F'
                if len(G.nodes[path[j + 1]]['fillcolor']) != 9:
                    G.nodes[path[j + 1]]['fillcolor'] += '7F'
                    G.nodes[path[j + 1]]['color'] = '#7080907F'
        shortestPaths[i].set_id("1" + letters[i])
    br = Branch(1, shortestPaths[0].path[-1])
    for pt in shortestPaths:
        br.add_path(pt)
    return [br]
