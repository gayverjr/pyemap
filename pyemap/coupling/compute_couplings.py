import math
import numpy as np
import scipy as sp
from .qc_data import parent_residues
import os
from pickle import load
try:
    from tensorflow import keras
    import tensorflow as tf
    from tensorflow.keras.layers import Dense, LeakyReLU
    from tensorflow.keras.initializers import GlorotNormal
    from tensorflow.keras.models import Sequential
    n_first_layer = 500
    n_second_layer = 1000
    leaky_relu = LeakyReLU(alpha=0.01)
    initializer = GlorotNormal()
    regularizer = None
    model = Sequential()
    model.add(Dense(n_first_layer, input_dim=32, activation=leaky_relu,
            kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.add(Dense(n_second_layer, activation=leaky_relu, kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.add(Dense(n_second_layer / 2, activation='tanh', kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.add(Dense(n_second_layer / 4, activation=leaky_relu, kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.add(Dense(n_second_layer / 8, activation=leaky_relu, kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.add(Dense(n_second_layer / 16, activation=None, kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.add(Dense(1, activation=None, kernel_initializer=initializer, kernel_regularizer=regularizer))
    model.summary()
    model.load_weights(os.path.join(os.path.dirname(__file__),'full_weights'))
    feature_transform = load(open(os.path.join(os.path.dirname(__file__),'full.pkl'), 'rb'))
except Exception as e:
    print(e)
    pass
from pyscf import gto
from functools import reduce
from pandas import DataFrame as df
import os




name_to_charge = {'C':6, 'N':7, 'O':8, 'H':1}
name_to_mass = {'C': 12.0107, 'H': 1.00784, 'N': 14.0067, 'O': 15.999}

def dist(vec1, vec2):
    return np.linalg.norm(vec2 - vec1)

def equation_plane(p1, p2, p3):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    x3, y3, z3 = p3
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * x1 - b * y1 - c * z1)
    return a, b, c, d


def find_dihedral_angle_between_planes(coeffs1, coeffs2):
    a1, b1, c1 = coeffs1
    a2, b2, c2 = coeffs2
    d = (a1 * a2 + b1 * b2 + c1 * c2)
    e1 = np.sqrt(a1 * a1 + b1 * b1 + c1 * c1)
    e2 = np.sqrt(a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)
    angle = np.degrees(np.arccos(d))
    return angle

class Atom(object):    
    def __init__(self,element,name,x,y,z):
        self.element = element
        self.name = name
        self.mass = name_to_mass[self.element]
        self.charge = name_to_charge[self.element]
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.xyz = np.array([self.x, self.y, self.z])

    def __str__(self):
        return '  {:}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format(self.element, self.xyz[0], self.xyz[1], self.xyz[2])


class Molecule(object):
    def __init__(self,residue,parent):
        self.atoms = []
        for atm in parent.atms:
            res_atm = residue[atm]
            self.atoms.append(Atom(res_atm.element,atm,res_atm.coord[0],res_atm.coord[1],res_atm.coord[2]))
        self.compute_com_()
        self.parent=parent
        self.mol = parent.rotate_coords(residue)
        self.charges_p = list(self.parent.C_p_charges.values())
        self.charges = list(self.parent.C_charges.values())

    def __str__(self):
        string_output = ''
        for atom in self.atoms:
            string_output += atom.__str__()
        return string_output.rstrip()

    def compute_com_(self):
        masses = np.array([atom.mass for atom in self.atoms])[:, np.newaxis]
        all_xyz = np.array([atom.xyz for atom in self.atoms])
        sum_mass = np.sum(masses)
        sum_prod = np.sum(masses * all_xyz, axis=0)
        self.com = sum_prod / sum_mass
        self.n_atoms = len(masses)
        return self.com


    def compute_min_dist_with_mol(self, mol):
        mol1_xyz = np.array([atom.xyz for atom in mol.atoms])
        mol2_xyz = np.array([atom.xyz for atom in self.atoms])
        min_dist_matrix = np.sqrt((np.square(mol1_xyz[:, np.newaxis, :] - mol2_xyz).sum(axis=2)))
        min_dist = np.min(min_dist_matrix)
        return min_dist

    def compute_dihedral_with_mol(self, mol):
        mol1_xyz = np.array([atom.xyz for atom in self.atoms])
        mol2_xyz = np.array([atom.xyz for atom in mol.atoms])
        plane1 = equation_plane(mol1_xyz[0], mol1_xyz[1], mol1_xyz[2])
        plane2 = equation_plane(mol2_xyz[0], mol2_xyz[1], mol2_xyz[2])
        dihedral_angle = find_dihedral_angle_between_planes(plane1[:3], plane2[:3])
        return dihedral_angle

    def compute_energy_gap(self, mol):
        mol1_xyz = np.array([atom.xyz for atom in self.atoms])
        mol2_xyz = np.array([atom.xyz for atom in mol.atoms])
        charge_dist_vec1, charge_dist_vec2 = [], []
        for i in range(len(mol1_xyz)):
            for j in range(len(mol2_xyz)):
                d = dist(mol1_xyz[i], mol2_xyz[j])
                term1 = (self.charges_p[i] * mol.charges[j]) / d
                term2 = (self.charges[i] * mol.charges_p[j]) / d
                charge_dist_vec1.append(term1)
                charge_dist_vec2.append(term2)
        energy_state1 = np.array(charge_dist_vec1).sum()
        energy_state2 = np.array(charge_dist_vec2).sum()
        energy_gap = energy_state1 - energy_state2
        return np.array([energy_state1, energy_state2, energy_gap])

    def compute_distance_vec_with_mol(self, mol):
        mol1_xyz = np.array([atom.xyz for atom in mol.atoms])
        mol2_xyz = np.array([atom.xyz for atom in self.atoms])
        dist_matrix = np.sqrt(np.sum(np.square(mol1_xyz[:, np.newaxis, :] - mol2_xyz), axis=2))
        dist_vec = np.sort(dist_matrix.flatten())
        return dist_vec

    def compute_eigv_charge_dist_vec_with_mol(self, mol, n_eigv=20):
        all_atoms = mol.atoms + self.atoms
        all_atoms_xyz = np.array([atom.xyz for atom in all_atoms])
        charges1 = np.array([atom.charge for atom in all_atoms])
        cm_size = len(all_atoms)
        cm1 = np.zeros((cm_size, cm_size))
        for i in range(len(all_atoms)):
            for j in range(len(all_atoms)):
                if i == j:
                    cm1[i, i] = 0.5 * (charges1[i] ** 2.4)
                else:
                    d = dist(all_atoms_xyz[i], all_atoms_xyz[j])
                    cm1[i, j] = charges1[i] * charges1[j] / d
        eigval1, eigvec1 = sp.linalg.eigh(cm1)
        maxdim = n_eigv
        eigval1 = np.sort(np.abs(eigval1))[::-1]
        eigval1 = np.concatenate([eigval1, np.zeros(maxdim - eigval1.shape[0])])
        return eigval1

def compute_overlap(molecule1, molecule2):
    mol1 = molecule1.mol
    mol2 = molecule2.mol
    parent1 = molecule1.parent
    parent2 = molecule2.parent
    s = gto.intor_cross('int1e_ovlp', mol1, mol2)
    mo_coeff1 = parent1.rotate_mo(mol1)
    mo_coeff2 = parent2.rotate_mo(mol2)
    S_AB = reduce(np.dot, (mo_coeff1[:, parent1.nocc - 1].T, s, mo_coeff2[:, parent2.nocc - 1]))
    return S_AB

def get_feature_vec(res1,res2):
    feature_vec = []
    parent1 = parent_residues[res1.resname.upper()]
    parent2 = parent_residues[res2.resname.upper()]
    molecule1 = Molecule(res1,parent1)
    molecule2 = Molecule(res2,parent2)
    homo_overlap = compute_overlap(molecule1,molecule2)
    feature_vec = np.append(feature_vec, np.abs(homo_overlap))
    mol1_com = molecule1.compute_com_()
    mol2_com = molecule2.compute_com_()
    com_dist = dist(mol1_com, mol2_com)
    feature_vec = np.append(feature_vec, com_dist)
    dist_vec = molecule1.compute_distance_vec_with_mol(molecule2)
    mean = np.mean(dist_vec)
    median = np.median(dist_vec)
    min_dist = np.min(dist_vec)
    hopfield = 2.7/np.sqrt(molecule1.parent.natms*molecule2.parent.natms) * np.exp(-0.72 * min_dist) * 1000
    max_dist = np.max(dist_vec)
    std = np.std(dist_vec)
    feature_vec = np.append(feature_vec, [mean, median, min_dist, hopfield, max_dist, std])
    dihedral_angle = molecule1.compute_dihedral_with_mol(molecule2)
    feature_vec = np.append(feature_vec, dihedral_angle)
    energy_gap = molecule1.compute_energy_gap(molecule2)
    feature_vec = np.append(feature_vec, energy_gap)
    nuccm_eigv = molecule1.compute_eigv_charge_dist_vec_with_mol(molecule2)
    feature_vec = np.append(feature_vec, nuccm_eigv)
    return np.array(feature_vec)

def compute_couplings(G,residues):
    res_dict = {}
    for res in residues:
        res_dict[res.node_label] = res
    feature_vec = []
    edge_list = []
    for edge in G.edges():
        edge_list.append(edge)
        feature_vec.append(get_feature_vec(res_dict[edge[0]],res_dict[edge[1]]))
    couplings = model.predict(feature_transform.transform(feature_vec))
    counter = 0
    for edge in G.edges():
        if edge in edge_list:
            G.edges[edge]['coupling'] = couplings[counter]*1E-3
            counter+=1
        else:
            res1 = res_dict[edge[0]]
            res2 = res_dict[edge[1]]
            parent1 = parent_residues[res1.resname.upper()]
            parent2 = parent_residues[res2.resname.upper()]
            hopfield = 2.7/np.sqrt(parent1.natms*parent2.natms) * np.exp(-0.72 * G.edges[edge]['catm'])
            G.edges[edge]['coupling'] = hopfield
    
