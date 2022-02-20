import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from pyscf import lib, gto
from pyscf.tools.molden import load as load_molden
import scipy
import os
from pandas import DataFrame as df


def fast_build(self,geom,unit='ANG'):
    PTR_ENV_START   = 20
    mol = gto.Mole()
    mol.verbose = 0
    mol.cart = self.cart
    mol._env = self._env.copy()
    mol._atom = mol.format_atom(geom, unit=unit)
    mol._basis = mol.format_basis(self._basis)
    env = mol._env[:PTR_ENV_START]
    mol._atm, mol._bas, mol._env = \
                mol.make_env(mol._atom, mol._basis, env, self.nucmod,
                              self.nucprop)
    return mol

def rotation_matrices(mol, R):
    ANG_OF = 1
    l_max = mol._bas[:, ANG_OF].max()
    Ds = [np.ones((1, 1))]
    for l in range(1, l_max + 1):
        # All possible x,y,z combinations
        cidx = np.sort(lib.cartesian_prod([(0, 1, 2)] * l), axis=1)
        addr = 0
        affine = np.ones((1, 1))
        for i in range(l):
            nd = affine.shape[0] * 3
            affine = np.einsum('ik,jl->ijkl', affine, R).reshape(nd, nd)
            addr = addr * 3 + cidx[:, i]
        uniq_addr, rev_addr = np.unique(addr, return_inverse=True)
        ncart = (l + 1) * (l + 2) // 2
        assert ncart == uniq_addr.size
        trans = np.zeros((ncart, ncart))
        for i, k in enumerate(rev_addr):
            trans[k] += affine[i, uniq_addr]
        Ds.append(trans)
    return Ds

class ResidueData(object):

    def _rotate_mo(self,a,b):
        sup = SVDSuperimposer()
        sup.set(a, b)
        sup.run()
        rms = sup.get_rms()
        rot, tran = sup.get_rotran()
        mats = rotation_matrices(self.mol,rot)
        mat_list = []
        for i in range(0,self.mol.nbas):
            l = self.mol.bas_angular(i)
            nctr = self.mol.bas_nctr(i)
            for j in range(0,nctr):
                mat_list.append(mats[l])
        M = scipy.linalg.block_diag(*mat_list)
        new_coeff = M @ self.mo
        return new_coeff

    def rotate_mo(self, residue):
        a = np.array([atm[1] for atm in self.mol._atom if atm[0][0]!='H'])
        b = np.array([atm[1] for atm in residue._atom if atm[0][0]!='H'])
        return self._rotate_mo(a,b)

    def rotate_coords(self,residue):
        full = []
        a = []
        b = []
        for atm in self.atms:
            a.append(self.atom_map[atm])
            b.append(residue[atm].coord)
        for atm in self.atom_map.keys():
            full.append(self.atom_map[atm])
        a = np.array(a)
        b = np.array(b)
        sup = SVDSuperimposer()
        sup.set(b,a)
        sup.run()
        rms = sup.get_rms()
        rot, tran = sup.get_rotran()
        a = np.dot(full,rot) + tran
        geom_str = ''
        for i,atm in enumerate(self.atom_map.keys()):
            geom_str+= '{}{}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format(atm[0], i+1, a[i][0], a[i][1], a[i][2])
        return fast_build(self.mol,geom_str)


class Tryptophan(ResidueData):

    def __init__(self):
        self.atms = [
            'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'
        ]
        self.atom_map = {
            'CB': [2.64268823, -1.34949645, 1.25e-06],
            'CG': [1.62471728, -0.25016315, 2.5e-07],
            'CD1': [1.85536396, 1.09108068, -9e-08],
            'CD2': [0.1934185, -0.408283, 2.5e-07],
            'NE1': [0.66093584, 1.78154989, -4.6e-07],
            'CE2': [-0.37429193, 0.88051207, -3.7e-07],
            'CE3': [-0.65591677, -1.52148069, 4.9e-07],
            'CZ2': [-1.75357861, 1.08756334, -7.1e-07],
            'CZ3': [-2.0194667, -1.32349491, -3e-08],
            'CH2': [-2.56352604, -0.02752022, -6.2e-07],
            'H1': [2.15537456, -2.32333191, -9.41e-06],
            'H2': [3.28386335, -1.30837946, 0.88090513],
            'H3': [3.28387783, -1.30836746, -0.88089149],
            'H4': [2.79132717, 1.6253691, -4.9e-07],
            'H5': [0.56689533, 2.77879499, -1.35e-06],
            'H6': [-0.24891796, -2.52458898, 1.03e-06],
            'H7': [-2.1734844, 2.08508871, -1.09e-06],
            'H8': [-2.68613426, -2.17512427, -1.2e-07],
            'H9': [-3.6381292, 0.09572978, -8.2e-07]
        }
        self.natms = len(self.atom_map)
        self.C_charges = {'CB': 6.62037, 'CG': 6.1172, 'CD1': 6.02088, 'CD2': 6.03094, 'NE1': 7.52094, 'CE2': 5.8854, 'CE3': 6.21218, 'CZ2': 6.23801, 'CZ3': 6.23906, 'CH2': 6.2128}
        self.C_p_charges = {'CB': 6.65815, 'CG': 5.89543, 'CD1': 5.90018, 'CD2': 6.08476, 'NE1': 7.42766, 'CE2': 5.87708, 'CE3': 6.08411, 'CZ2': 6.1727, 'CZ3': 6.224, 'CH2': 6.09866} 
        self.h_charges = {'CB': 2.35271, 'CG': 0.0, 'CD1': 0.79461, 'CD2': 0.0, 'NE1': 0.59211, 'CE2': 0.0, 'CE3': 0.79193, 'CZ2': 0.78969, 'CZ3': 0.79071, 'CH2': 0.79046}
        self.h_p_charges = {'CB': 2.2355099999999997, 'CG': 0.0, 'CD1': 0.75357, 'CD2': 0.0, 'NE1': 0.55621, 'CE2': 0.0, 'CE3': 0.76389, 'CZ2': 0.75688, 'CZ3': 0.75471, 'CH2': 0.75648}
        self.charges = {}
        self.p_charges = {}
        for key in self.C_charges:
            self.charges[key] = self.h_charges[key] + self.C_charges[key]
            self.p_charges[key] = self.h_p_charges[key] + self.C_p_charges[key]
        self.nocc = (load_molden(os.path.join(os.path.dirname(__file__),'trp.molden'))[3] >= 2).sum()
        self.mol = load_molden(os.path.join(os.path.dirname(__file__),'trp.molden'))[0]
        self.mo =  np.load(os.path.join(os.path.dirname(__file__),'trp_mo.npy'))
        self.sigma = {'CB': 0.358141284692,
        'CG': 0.354577689820,
        'CD1': 0.355005321205,
        'CD2': 0.331414323148,
        'NE1': 0.329632525712,
        'CE2': 0.331414323148,
        'CE3': 0.354577689820,
        'CZ2': 0.354577689820,
        'CZ3': 0.355005321205,
        'CH2': 0.355005321205
        }
        self.epsilon = {'CB': 0.23430,
        'CG': 0.30543,
        'CD1': 0.29288,
        'CD2': 0.41422,
        'NE1': 0.83680,
        'CE2': 0.41422,
        'CE3': 0.30543,
        'CZ2': 0.30543,
        'CZ3': 0.29288,
        'CH2': 0.29288
        }



class Histidine(ResidueData):

    def __init__(self):
        self.atms = ['CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2']
        self.atom_map = {
            'CB': [-2.19224994, 0.05373494, 0.0],
            'H1': [-2.56147712, 0.57783875, -0.88267527],
            'H2': [-2.62539166, -0.94357792, 0.0],
            'CG': [-0.70729482, -0.04887384, 0.0],
            'ND1': [0.12115375, 1.05094163, 0.0],
            'CE1': [1.40561546, 0.59719904, 0.0],
            'H3': [2.25158141, 1.26382539, 0.0],
            'NE2': [1.45578618, -0.70546026, 0.0],
            'CD2': [0.14065963, -1.11632245, 0.0],
            'H4': [-0.12271217, -2.16038534, 0.0],
            'H5': [-0.16914902, 2.01165174, 0.0],
            'H6': [-2.56147712, 0.57783875, 0.88267527]
        }
        self.natms = len(self.atom_map)
        self.C_charges = {'CB': 6.6346, 'CG': 5.91116, 'ND1': 7.52299, 'CE1': 5.84072, 'NE2': 7.45174, 'CD2': 6.10723}
        self.C_p_charges = {'CB': 6.69157, 'CG': 5.6386, 'ND1': 7.53567, 'CE1': 5.62346, 'NE2': 7.39665, 'CD2': 5.90179}
        self.h_charges = {'CB':2.3373, 'CG': 0.0, 'ND1':0.58974, 'CE1':0.80811, 'NE2':0.0, 'CD2':0.79642}
        self.h_p_charges = {'CB':2.16663, 'CG':0.0, 'ND1':0.54635, 'CE1':0.75293, 'NE2':0.0, 'CD2':0.74636}
        self.charges = {}
        self.p_charges = {}
        for key in self.C_charges:
            self.charges[key] = self.h_charges[key] + self.C_charges[key]
            self.p_charges[key] = self.h_p_charges[key] + self.C_p_charges[key]
        self.sigma = {'CB':0.358141284692,
                     'CG':0.320723538531 ,
                     'ND1':0.329632525712,
                     'CE1':0.320723538531,
                     'NE2':0.329632525712,
                     'CD2':0.320723538531
        }
        self.epsilon = {'CB':0.23430,
                     'CG':0.20920,
                     'ND1':0.83680,
                     'CE1': 0.20920,
                     'NE2':0.83680,
                     'CD2': 0.20920
        }
        self.mo =  np.load(os.path.join(os.path.dirname(__file__),'his_mo.npy'))
        self.nocc = (load_molden(os.path.join(os.path.dirname(__file__),'his.molden'))[3] >= 2).sum()
        self.mol = load_molden(os.path.join(os.path.dirname(__file__),'his.molden'))[0]


class Phenylalanine(ResidueData):

    def __init__(self):
        self.atms = ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']
        self.atom_map = {
            'CB': [-2.44881021, -0.0, 0.00714743],
            'CG': [-0.94385642, -0.0, -0.01144184],
            'CD1': [-0.23081234, 1.19480096, -0.00838797],
            'CE1': [1.15624112, 1.19786294, 0.00352554],
            'CZ': [1.85570972, -0.0, 0.01058398],
            'CE2': [1.15624112, -1.19786294, 0.00352554],
            'CD2': [-0.23081234, -1.19480096, -0.00838797],
            'H1': [-2.8238392, -0.0, 1.03187612],
            'H2': [-2.84908222, -0.88238577, -0.48972183],
            'H3': [-0.76939713, 2.13433596, -0.01820754],
            'H4': [1.69155514, 2.13770061, 0.00388821],
            'H5': [2.93684931, -0.0, 0.01714545],
            'H6': [1.69155514, -2.13770061, 0.00388821],
            'H7': [-0.76939713, -2.13433596, -0.01820754],
            'H8': [-2.84908222, 0.88238577, -0.48972183]
        }
        self.natms = len(self.atom_map)
        self.C_charges = {'CB': 6.62583, 'CG': 6.00386, 'CD1': 6.21962, 'CE1': 6.22013, 'CZ': 6.17063, 'CE2': 6.22013, 'CD2': 6.21962}
        self.C_p_charges = {'CB': 6.67757, 'CG': 5.74335, 'CD1': 6.16757, 'CE1': 6.19323, 'CZ': 5.90275, 'CE2': 6.19323, 'CD2': 6.16757}
        self.h_charges = {'CB': 2.35271, 'CG':0.0, 'CD1':0.79392, 'CE1':0.79051, 'CZ':0.79859, 'CE2':0.79051, 'CD2':0.79392}
        self.h_p_charges = {'CB':2.19617, 'CG': 0.0, 'CD1':0.75084, 'CE1':0.74535, 'CZ':0.76618, 'CE2':0.74535, 'CD2':0.75084}
        self.charges = {}
        self.p_charges = {}
        self.sigma = { 'CB': 0.358141284692,
                        'CG':0.355005321205,
                        'CD1':0.355005321205,
                        'CE1':0.355005321205,
                        'CZ':0.355005321205,
                        'CE2':0.355005321205,
                        'CD2':0.355005321205,

        }
        self.epsilon = { 'CB': 0.23430,
                        'CG':0.29288,
                        'CD1':0.29288,
                        'CE1':0.29288,
                        'CZ':0.29288,
                        'CE2':0.29288,
                        'CD2':0.29288
        }
        for key in self.C_charges:
            self.charges[key] = self.h_charges[key] + self.C_charges[key]
            self.p_charges[key] = self.h_p_charges[key] + self.C_p_charges[key]
        self.mo =  np.load(os.path.join(os.path.dirname(__file__),'phe_mo.npy'))
        self.nocc = (load_molden(os.path.join(os.path.dirname(__file__),'phe.molden'))[3] >= 2).sum()
        self.mol = load_molden(os.path.join(os.path.dirname(__file__),'phe.molden'))[0]


class Tyrosine(ResidueData):

    def __init__(self):
        self.atms = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']
        self.atom_map = {
            'CB': [2.9355494, 0.05017324, 8.1e-07],
            'CG': [1.43013275, 0.03199901, 2e-08],
            'CD1': [0.69036761, 1.20609065, -1.7e-07],
            'CD2': [0.7305609, -1.17403155, -1.9e-07],
            'CE1': [-0.69896649, 1.18806853, -3.2e-07],
            'CE2': [-0.65188666, -1.21173495, -3.1e-07],
            'CZ': [-1.37130283, -0.02341759, -1.7e-07],
            'OH': [-2.73642182, -0.11046395, 5e-08],
            'H1': [3.31509106, 1.07069129, 5.1e-07],
            'H2': [3.33702509, -0.45550079, -0.87907029],
            'H3': [3.33702407, -0.45549983, 0.87907295],
            'H4': [1.20384621, 2.15969978, -7e-08],
            'H5': [1.2816002, -2.10703356, -2.5e-07],
            'H6': [-1.25573016, 2.1182111, -4.2e-07],
            'H7': [-1.18838502, -2.15044513, -3.6e-07],
            'H8': [-3.11095376, 0.77322914, 1.17e-06]
        }
        self.natms = len(self.atom_map)
        self.C_charges = {'CB': 6.61758, 'CG': 6.05805, 'CD1': 6.18615, 'CD2': 6.17422, 'CE1': 6.31795, 'CE2': 6.29227, 'CZ': 5.62632, 'OH': 8.68117}
        self.C_p_charges = {'CB': 6.6654, 'CG': 5.81423, 'CD1': 6.16524, 'CD2': 6.16619, 'CE1': 6.24021, 'CE2': 6.21479, 'CZ': 5.49169, 'OH': 8.53906}
        self.h_charges = {'CB': 2.3582899999999998, 'CG': 0.0, 'CD1': 0.79277, 'CD2': 0.79355, 'CE1': 0.79604, 'CE2': 0.77906, 'CZ': 0.0, 'OH': 0.52654}
        self.h_p_charges = {'CB': 2.2204699999999997, 'CG': 0.0, 'CD1': 0.7497, 'CD2': 0.75034, 'CE1': 0.75303, 'CE2': 0.73651, 'CZ': 0.0, 'OH': 0.49315}
        self.charges = {}
        self.p_charges = {}
        for key in self.C_charges:
            self.charges[key] = self.h_charges[key] + self.C_charges[key]
            self.p_charges[key] = self.h_p_charges[key] + self.C_p_charges[key]
        self.mo =  np.load(os.path.join(os.path.dirname(__file__),'tyr_mo.npy'))
        self.nocc = (load_molden(os.path.join(os.path.dirname(__file__),'tyr.molden'))[3] >= 2).sum()
        self.mol = load_molden(os.path.join(os.path.dirname(__file__),'tyr.molden'))[0]
        self.sigma = {  
            'CB': 0.358141284692,
            'CG': 0.355005321205,
            'CD1': 0.355005321205,
            'CD2': 0.355005321205,
            'CE1': 0.355005321205,
            'CE2': 0.355005321205,
            'CZ': 0.355005321205,
            'OH': 0.315378146222,
        }
        self.epsilon = {
            'CB': 0.23430,
            'CG': 0.29288,
            'CD1': 0.29288,
            'CD2': 0.29288,
            'CE1':  0.29288,
            'CE2': 0.29288,
            'CZ': 0.29288,
            'OH': 0.63639,
        }

parent_residues =  {'TRP': Tryptophan(), 'TYR':Tyrosine(), 'HIS': Histidine(), 'PHE': Phenylalanine()}