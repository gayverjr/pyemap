# This Python script is a part of
# eMAP: online mapping of electron transfer channels in biomolecules
# Copyright(C) 2017-2018 Ruslan Tazhigulov, James Gayvert, Melissa Wei, Ksenia Bravaya (Boston University, USA)
"""Finds shortest paths in graph given a source and optionally a target node.

Module used to find and visualize shortest paths. If the user specifies only a source node, dijkstra's algorithm is
used to calculate the shortest path from the source to every surface exposed residue.
If the user specifies a target as well, Yen's algorithm is used to find the 5 shortest paths from source to target.
Each of these pathways is assigned a unique ID for visualization, and a list of all of these pathways is returned.
The results of the analysis are written out to two files. "useroutput.txt" is a results
file displayed to the user on the webpage. "outputs.txt" is a list of paths in NGL Viewer's selection language, which
can be later accessed by the front end for visualization of these pathways.

Usage
-----
Called by the views module to obtain shortest paths. Returns a list of pathway ids.

"""

import itertools
import os
import string
import sys

import networkx as nx
import numpy as np
import pygraphviz as pg
from networkx.drawing.nx_agraph import from_agraph, to_agraph


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
    single_chain: boolean
        Boolean for single chain or not

    """

    def __init__(self, path, length, single_chain):
        self.path = path
        self.length = length
        self.path_id = "none"
        self.single_chain = single_chain

    def __eq__(self, other):
        return self.length == other.length

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.length < other.length

    def __str__(self):
        # If only one chain is present, do not include chain id in output
        if self.single_chain:
            for i in range(len(self.path)):
                node = self.path[i]
                if "(" in node:
                    node_str = node_str = node[:node.index(
                        "(")] + node[node.index(")") + 1:]
                    self.path[i] = node_str
        printline = self.path_id + ": " + \
            str(self.path) + " " + str('{:.2f}'.format(round(self.length, 2)))
        return printline

    def get_path_as_list(self):
        if self.single_chain:
            for i in range(len(self.path)):
                node = self.path[i]
                if "(" in node:
                    node_str = node_str = node[:node.index(
                        "(")] + node[node.index(")") + 1:]
                    self.path[i] = node_str
        original_list = [[self.path_id], self.path,
                         [str('{:.2f}'.format(round(self.length, 2)))]]
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
                    atm_name_index = str.index(
                        customnum[customname.index(name)], ".")
                    end_idx = str.index(customnum[customname.index(name)],
                                        ")", atm_name_index)
                    select_arr.append(customnum[customname.index(name)][
                        atm_name_index:end_idx])
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
    single_chain: boolean
        Boolean for single chain or not

    See Also
    --------
    class ShortestPath

    """

    def __init__(self, branch_id, target, single_chain):
        self.branch_id = branch_id
        self.target = target
        self.paths = []
        self.single_chain = single_chain

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
        if self.single_chain:
            if "(" in printline:
                printline = printline[:printline.index(
                    "(")] + printline[printline.index(")") + 1:]
        for pt in self.paths:
            printline += "\n" + str(pt)
        printline += "\n"
        return printline

    def get_branch_as_list(self):
        branch_list = []
        printline = "Branch: " + str(self.target)
        if self.single_chain:
            if "(" in printline:
                printline = printline[:printline.index(
                    "(")] + printline[printline.index(")") + 1:]
        branch_list.append([printline])
        for pt in self.paths:
            branch_list.append(pt.get_path_as_list())
        return branch_list


def yens_shortest_paths(G, start, goal, filename, single_chain):
    """Returns top 5 shortest paths from source to target.

    Uses Yen's algorithm to calculate the 5 shortest paths from source to target, writes
    out the ShortestPath objects to file, and returns the 5 pathway IDs. In the graph, nodes and
    edges that are part of any pathways are made opaque, and the shortest path is highlighted.

    Note
    ----
    This function is only called when a target is specified by the user.

    Parameters
    ----------
    G: NetworkX graph object
        Undirected, weighted residue graph
    start: str
        Source node
    goal: str
        Target node
    filename: str
        File hash for writing out to file
    single_chain: boolean
        Boolean for single chain or not

    Returns
    -------
    all_pt_ids: array-like
        List of pathway IDs corresponding to ShortestPath objects that have been written to file.
    branches_table: array-like
        List of branches and their pathway IDs corresponding to ShortestPath objects that can be displayed in a table (flask_table)

    See Also
    --------
    module NetworkX.shortest_simple_paths
    class ShortestPath

    References
    ----------
    Jin Y. Yen, Finding the K Shortest Loopless Paths in a Network, Management Science,
    Vol. 17, No. 11, Theory Series (Jul., 1971), pp. 712-716.

    """

    letters = list(string.ascii_letters)
    shortestPaths = []
    all_pt_ids = []
    done = False
    k = 0
    from itertools import islice
    paths = list(islice(nx.shortest_simple_paths(G, start, goal), 5))
    for k in range(0, len(paths)):
        path = paths[k]
        sum = 0
        for i in range(0, len(path) - 1):  # sum up edge weights
            sum += (G[path[i]][path[i + 1]]['weight'])
        path = ShortestPath(path, sum, single_chain)
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
        # write out paths to file
        useroutput = open(filename +".txt", "w")
        fi = open(filename +".txt", "r")
        strArr = fi.readlines()
        customname = []
        customnum = []
        for i in range(len(strArr)):
            tempArr = strArr[i].split(";")
            customname.append(tempArr[0])
            customnum.append(tempArr[1][:-1])
        fi.close()
        f = open(filename +".txt", "w+")
        branches_table = []
        for pt in shortestPaths:
            f.write(str(pt.path_id) + ":")
            f.write(str(pt.selection_str(customname, customnum)) + "\n")
            branches_table.append(pt.get_path_as_list())
            all_pt_ids.append(str(pt).split()[0] + ' ' + str(pt).split()[-1])
            useroutput.write(str(pt) + "\n")
        useroutput.close()
        f.close()
        fi.close()
        return all_pt_ids, branches_table
    else:  # no paths found
        useroutput = open(filename +".txt", "w")
        str_out = "No paths found."
        branches_table.append(str_out)
        useroutput.write(str_out)
        useroutput.close()
        return []


def is_parent_pathway(shortest_path, goals):
    """Returns true if ShortestPath is a parent pathway, false if not.

    A ShortestPath object is the parent of a branch if its terminal residue is the
    only surface exposed residue in the path. For example, if goals=[A,B,C] and
    the pathway is [H,I,C], then this pathway is a parent pathway. In contrast, if
    the pathway is [H,B,A], then this pathway is not a parent pathway.

    Parameters
    ----------
    shortest_path: ShortestPath
        ShortestPath object
    goals: array-like
        List of surface exposed residues

    Returns
    -------
    True if a parent pathway
    False if not a parent pathway
    """

    count = 0
    for res in shortest_path.path:
        if res in goals:
            count += 1
    return count == 1


def find_branch(pt, goals, branches):
    """Determines which branch a pathway belongs to and returns that branch.

    A ShortestPath belongs to a branch if the first surface exposed residue it reaches during the
    pathway is the target of that particular branch.

    Parameters
    ----------
    pt: ShortestPath
        ShortestPath object
    goals: array-like
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
    while res not in goals:
        count += 1
        res = pt.path[count]
    count = 0
    cur_branch = branches[0]
    while not res == cur_branch.target:
        count += 1
        cur_branch = branches[count]
    return cur_branch


def dijkstras_shortest_paths(G, start, goals, filename, single_chain):
    """Returns shortest path from source to each surface exposed residue.

    Performs Dijkstra's algorithm from the source to each surface exposed residue, finding the
    shortest path. The ShortestPath objects are organized into branches based on the first surface
    exposed residue reached during the course of the pathway. The ShortestPaths are then written out
    to file, and a list of pathway IDs returned. In the graph, nodes and edges that are part of any
    paths are made opaque.

    Note
    ----
    This function is only called when a target is not specified by the user.

    Parameters
    ----------
    G: NetworkX graph object
        Undirected, weighted residue graph
    start: str
        Source node
    goals: array-like
        List of surface exposed residues
    filename: str
        File hash for writing out to file
    single_chain: boolean
        Boolean for single chain or not

    Returns
    -------
    all_pt_ids: array-like
        List of pathway IDs corresponding to ShortestPath objects that have been written to file.
    branches_table: array-like
        List of branches and their pathway IDs corresponding to ShortestPath objects that can be displayed in a table (flask_table)

    See Also
    --------
    module NetworkX.dijkstra_path
    class ShortestPath
    class Branch

    """
    shortestPaths = []
    all_pt_ids = []
    for goal in goals:
        path = []
        try:
            path = nx.dijkstra_path(G, start, goal)
        except:
            path = []
        if not path == []:
            sum = 0
            for i in range(0, len(path) - 1):  # sum up edge weights
                sum += (G[path[i]][path[i + 1]]['weight'])
            shortestPaths.append(ShortestPath(path, sum, single_chain))
    shortestPaths = sorted(shortestPaths)
    branches = []
    # find the parent pathways
    for pt in shortestPaths:
        if is_parent_pathway(pt, goals):
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
            br = Branch(len(branches) + 1, pt.path[-1], single_chain)
            branches.append(br)
            br.add_path(pt)
    # find the sub pathways
    for pt in shortestPaths:
        if not is_parent_pathway(pt, goals):
            find_branch(pt, goals, branches).add_path(pt)
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
    '''
    # write out to file
    fi = open(filename +".txt", "r")
    strArr = fi.readlines()
    customname = []
    customnum = []
    for i in range(len(strArr)):
        tempArr = strArr[i].split(";")
        customname.append(tempArr[0])
        customnum.append(tempArr[1][:-1])
    fi.close()
    # output for NGL viewer
    f = open(filename +".txt", "w+")
    for pt in shortestPaths:
        f.write(str(pt.path_id) + ":")
        f.write(str(pt.selection_str(customname, customnum)) + "\n")
    f.close()
    '''
    # output for text file
    useroutput = open(filename + ".txt","w")

    branches_table = []
    for br in branches:
        branch_paths = str(br).strip().split('\n')
        branches_table.append(br.get_branch_as_list())
        for path in branch_paths:
            if not "Branch" in path:
                # for view pathways combobox keep only path ID and distance
                all_pt_ids.append(path.split()[0] + ' ' + path.split()[-1])
            else:
                all_pt_ids.append(path)
        useroutput.write(str(br))
    branches_table = list(itertools.chain(*branches_table))
    if len(shortestPaths) == 0:
        str_out = "No paths found."
        branches_table.append(str_out)
        useroutput.write(str_out)
        raise Exception("No paths to the surface found.")
    useroutput.close()
    return all_pt_ids, branches_table


def processName(G, name, single_chain):
    """Returns the node in G that the string name corresponds to.

    On the front end, the user selects a source/target residue. If there is only a single
    chain, that chain name is omitted in what is displayed to the user, but still included
    in G on the backend. This function ensures that the correct node is selected on
    the backend in the single chain case.

    Note
    ----
    This used to be for helping the user out for mis-spellings etc. but now that we use
    combo boxes for node selection, it's really just for the single chain case.

    Parameters
    ----------
    G: NetworkX graph
        A weighted, undirected residue graph
    name: str
        A source/target node specified by the user
    single_chain: boolean
        Boolean for single chain or not

    Raises
    ------
    Exception e:
        Invalid name (should never happen)

    """
    name = name.strip()
    for node in G.nodes():
        if single_chain:
            node_name = node[:node.index("(")] + node[node.index(")") + 1:]
        else:
            node_name = node
        if name == node_name:
            return node
    raise Exception("Invalid name")


def draw_graph(G, original_shape_start, source, filename):
    """Draws the graph with the shortest pathways highlighted, and writes them out to file for downlaod and use by
    the front end.

    Parameters
    ----------
    G: NetworkX graph object
        Residue graph
    source: str
        name of source node
    original_shape_start: str
        shape of source node
    filename: str
        file hash for writing out to file

    """
    G.node[source]['fillcolor'] = '#FFD700FF'
    G.node[source]['penwidth'] = 6.0
    G.node[source]['shape'] = original_shape_start
    # if not is not involved in pathways, make it less opaque on the graph.
    # change font color to slate gray, and change transparency of edges
    for name_node in G.nodes():
        if len(G.node[name_node]['fillcolor']) != 9:
            G.node[name_node]['fillcolor'] += '40'
            G.node[name_node]['fontcolor'] = '#708090'
    for edge in G.edges():
        name_node1, name_node2 = edge[0], edge[1]
        if G[name_node1][name_node2]['style'] == 'dashed':
            G[name_node1][name_node2]['color'] = '#7788994F'
    # draw graph
    A_new = to_agraph(G)
    try:  # no need for try catch in released version
        A_new.graph_attr.update(
            ratio=1.0, overlap="ipsep", mode="ipsep", splines="true")
        A_new.layout(args="-n2")
        #graph_filename = dir_path + \
        #    filename + "/graph" + filename + ".svg"
        downloadable = filename + ".png"
        #A_new.draw(graph_filename, prog="neato")
        A_new.draw(downloadable, prog="neato")
    except Exception as e:
        #graph_filename = dir_path + \
        #    filename + "/graph" + filename + ".svg"
        downloadable = filename + ".png"
        #A_new.draw(graph_filename, prog="neato")
        A_new.draw(downloadable, prog="neato")


def shortest_paths(emap, source, target=None, single_chain=False):
    """Main method of dijkstras module.

    Takes in input from views module, and then performs shortest path analysis
    on source and (optionally) target residues. After analysis is completed, the updated
    graph is drawn and written to file.

    Parameters
    ---------
    filename: str
        file hash for writing out to file
    source: str
        source node for analysis
    target: str
        target node for analysis. Can be [] if no target specified
    single_chain: boolean
        boolean for single chain or not

    Returns
    ------
    all_pt_ids: array-like
        List of pathway IDs corresponding to ShortestPath objects that have been written to file.


    See Also
    --------
    module Views

    Raises
    -------
    Exception e:
        Invalid name

    """
    # read in graph from file
    A = emap.agraph
    G = from_agraph(A)
    filename=emap.filename
    for u, v, d in G.edges(data=True):
        d['weight'] = np.float64(d['weight'])
    # process source and target
    source = source.strip()
    source = processName(G, source, single_chain)
    original_shape_start = G.node[source]['shape']
    G.node[source]['shape'] = 'oval'
    if target:
        target = target.strip()
        target = processName(G, target, single_chain)
        all_pt_ids, branches_table = yens_shortest_paths(
            G, source, target, filename, single_chain)
        # color target node blue
        G.node[target]['fillcolor'] = '#40e0d0FF'
        G.node[target]['penwidth'] = 6.0
    else:
        goals = []
        for n, d in G.nodes(data=True):
            if d['shape'] == "box":
                goals.append(n)
        all_pt_ids, branches_table = dijkstras_shortest_paths(
            G, source, goals, filename, single_chain)
    draw_graph(G, original_shape_start, source, filename)
    return all_pt_ids, branches_table
