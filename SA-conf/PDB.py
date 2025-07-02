#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" MODULE TO MANIPULATE PDB FILE

"""
__author__ = "jb cheron"
__version__ = "0.2"
__date__ = "06-18-2013"

import os
import re
import sys
import time
import urllib2
import collections
from math import sqrt, pi
from random import randint

#
#
#
#
#
#
#*******************************************************************************
# CONSTANTS
#
#*******************************************************************************
#

#### PREFIX
ERR_PREF = '[ERROR]'
NFO_PREF = '[INFO]'

#### REGEX
RE_ATOM = re.compile('^ATOM')
RE_HETATM = re.compile('^HETATM')
RE_COORDLINE = re.compile('^(ATOM|HETATM)')
RE_HOH = re.compile('^HETATM.+HOH')
RE_CONECT = re.compile('^CONECT')

#### OTHERS
#3 letters code Amino Acid dictionary
AA = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'GLY': 'G',
    'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'VAL': 'V', 'GLU': 'E',
    'TYR': 'Y', 'MET': 'M', 'UNK':'X'
}
#1 letter code Amino Acid list
AA_MONO = "ACDEFGHIKLMNPQRSTVWY"

CRIS_AGENTS = [
  'GOL', 'EDO', '12P', '15P', '7PE', 'P6G', 'PE3', 'PE4', 'PE5', 'PEG', 'MPD',
  'HEZ', 'MEB', 'DMS', 'EOH', 'OXA', 'OHE', 'EOX', 'IPA', 'MOH', 'OME', 'OMB',
  'SWE', 'SUC', 'MRY', 'XYL', 'INS', 'CBU', 'BUD', 'BU3'
]
IONS = [
  '118', '1AL', '2HP', '3CO', '3MT', '3NI', '4MO', '6MO', 'ACT', 'AG', 'AL',
  'ALF', 'ATH', 'AU', 'AU3', 'AUC', 'AZI', 'BA', 'BCT', 'BEF', 'BF4', 'BO4',
  'BR', 'BS3', 'CA', 'CAC', 'CD', 'CE', 'CHT', 'CL', 'CO', 'CO3', 'CON', 'CR',
  'CS', 'CSB', 'CU', 'CU1', 'CUA', 'CUZ', 'CYN', 'DME', 'DMI', 'DSC', 'DTI',
  'EDR', 'EMC', 'ER3', 'EU', 'EU3', 'F', 'FE', 'FE2', 'FPO', 'GA', 'GD3', 'GEP',
  'HAI', 'HG', 'IN', 'IOD', 'IR', 'IR3', 'IRI', 'IUM', 'K', 'LA', 'LCO', 'LCP',
  'LI', 'LU', 'MAC', 'MG', 'MH2', 'MLI', 'MMC', 'MN', 'MN3', 'MOO', 'MOS', 'MOW',
  'NA', 'NET', 'NH4', 'NI', 'NO2', 'NO3', 'NRU', 'OAA', 'OH', 'OS', 'OS4', 'OXL',
  'PB', 'PBM', 'PD', 'PDV', 'PER', 'PI', 'PO3', 'PO4', 'PR', 'PT', 'PT4', 'PTN',
  'RB', 'RH3', 'RHD', 'RU', 'SB', 'SCN', 'SE4', 'SEK', 'SM', 'SMO', 'SO3', 'SO4',
  'SR', 'T1A', 'TB', 'TBA', 'TCN', 'TEA', 'THE', 'TL', 'TMA', 'TRA', 'UNX', 'V',
  'VN3', 'VO4', 'W', 'WO5', 'Y1', 'YB', 'YB2', 'YT3', 'ZN'
]

#
#
#
#
#
#
#*******************************************************************************
# ATOM
#
#*******************************************************************************
#

def ask_to_user(text, values):
    """ prompt: question to user, while valid answer """
    fail = True
    while fail:
        response = raw_input(text+' (choices ['+' '.join(values)+']): ')
        print '> Choice', response
        if response in values:
            fail = False
    return response

def nb_heavy_atoms(atoms):
    """ compute the heavy atom number in list """
    heavy_atoms = 0
    for atm in atoms:
        atom = Atom(atm)
        if atom.name[0] != 'H':
            heavy_atoms += 1
    return heavy_atoms

def euc_dist(line1, line2):
    """ Calculation of euclidian distance between two atom lines
    ** arguments : 2 atom lines (char)
    ** return    : distance     (float)
    """
    coords1 = [float(line1[30:38]),\
               float(line1[38:46]),\
               float(line1[46:54])]
    coords2 = [float(line2[30:38]),\
               float(line2[38:46]),\
               float(line2[46:54])]
    dist = sqrt((coords1[0] - coords2[0])**2
              + (coords1[1] - coords2[1])**2
              + (coords1[2] - coords2[2])**2)
    return dist

def angle(l, c, r):
    """ return the angle between 3 atoms """
    def num(l, c, r):
        """ numerator """
        somme = 0
        for p in xrange(3):
            somme += ((l[p]-c[p])*(r[p]-c[p]))
        return somme
    def denom(l, c, r):
        """ denominator """
        s1, s2 = 0, 0
        for p in xrange(3):
            s1 += ((l[p]-c[p])**2)
            s2 += ((r[p]-c[p])**2)
        return sqrt(s1)*sqrt(s2)
    cl = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    cc = [float(c[30:38]), float(c[38:46]), float(c[46:54])]
    cr = [float(r[30:38]), float(r[38:46]), float(r[46:54])]
    numerat = num(cl, cc, cr)
    denominat = denom(cl, cc, cr)
    return 180*acos(numerat/denominat)/pi

def fetch_pdb(pdb_id, path="./"):
    """Download PDB file from RCSB

    ** argument : PDB 4 letters code (char)
    ** return   : -
    """
    url = "http://www.rcsb.org/pdb/files/"+pdb_id+".pdb"
    content = ''
    attempts = 0
    while True and attempts < 5:
        attempts += 1
        try:
            content = urllib2.urlopen(url)
            content = content.read()
            break
        except Exception, ex:
            if re.findall("403", str(ex)):
                time.sleep(randint(5, 10)/10.0)
            else:
                break
    if content and not re.findall("requested file is not available", content):
        pdb_file = open(os.path.join(path,pdb_id+".pdb"), "w+")
        pdb_file.write(content)
        pdb_file.close()
        return os.path.join(path, pdb_id+".pdb")
    return None

def split_seq(seq, step):
    '''carriage return in sequence every n nucleotides'''
    for i in xrange(0, len(seq)+step, step+1):
        seq.insert(i, '\n')
    return ''.join(seq)[1:-1]

#
#
#
#
#
#*******************************************************************************
# ATOM
#
#*******************************************************************************
#

class Atom(object):
    """Atom Class
    """
    def __init__(self, line):
        """ Instantitate """
        self.line = line
        self.x, self.y, self.z = self.get_coords()
        self.type = str(self.line[0:6]).replace(' ', '')
        self.name = str(self.line[12:16]).replace(' ', '')
        self.num = int(self.line[6:11])
        self.resid = str(self.line[17:20]).replace(' ', '')
        self.resnum = int(self.line[22:26])
        self.code = self.get_code()
        self.chain = str(self.line[21]).replace(' ','')
        self.res_uid = self.resid+str(self.resnum)+self.code+self.chain

    def get_coords(self):
        """ retrieve coordinates """
        x = float(self.line[30:38])
        y = float(self.line[38:46])
        z = float(self.line[46:54])
        return x, y, z

    def get_code(self):
        """ alternative residue code """
        if self.line[26] != ' ':
            return str(self.line[26])
        else:
            return ''

#
#
#
#
#
#*******************************************************************************
# RESIDUE
#
#*******************************************************************************
#

class Residue(object):
    """Residue Class
    """
    def __init__(self, lines):
        """ Instantiate """
        self.lines = lines

    def test(self):
        """ test """
        self.lines += ""
        somme = 1+1
        return somme

#
#
#
#
#
#*******************************************************************************
# LIGAND
#
#*******************************************************************************
#

class Ligand(object):
    """ Ligand Class
    """
    def __init__(self, lines):
        """ Instantiate """
        self.lines = lines

    def test(self):
        """ test """
        self.lines += ""
        somme = 1+1
        return somme

#
#
#
#
#
#*******************************************************************************
# PEPTIDE
#
#*******************************************************************************
#

class Peptide(object):
    """ Peptide Class
    """

    def __init__(self, lines):
        """ Instantiate """
        self.lines = lines

    def test(self):
        """ test """
        self.lines += ""
        a = 1+1
        return a

#
#
#
#
#
#*******************************************************************************
# CHAIN
#
#*******************************************************************************
#

class Chain(object):
    """ Chain Class
    """

    def __init__(self, lines):
        """ Instantiate """
        self.lines = lines
        self.atom = []
        self.hetatm = []
        self.water = []
        self.residues = self.get_residues()
        self.get_atoms()


    def get_atoms(self):
        """ parsing to retrieve each atom type """
        for line in self.lines:
            if RE_ATOM.search(line):
                self.atom.append(line)
            elif RE_HETATM.search(line) and not RE_HOH.search(line):
                self.hetatm.append(line)
            elif RE_HOH.search(line):
                self.water.append(line)

    def nb_water(self):
        """ number O water in chain """
        return len(self.water)

    def nb_atom(self):
        """ number ATOM in chain """
        return len(self.atom)

    def nb_hetatm(self):
        """ number HETAM in chain """
        return len(self.hetatm)

    def nb_res(self):
        """ number of residues in chain """
        return len(self.residues)

    def nb_het_res(self):
        """ number of hetero residues in chain """
        return len(self.get_hetatm())

    def get_water(self):
        """ get water atoms """
        water = collections.OrderedDict()
        for hoh in self.water:
            atom = Atom(hoh)
            water[atom.resid+str(atom.resnum)] = hoh
        return water

    def get_water_list(self):
        """ get water uid list """
        return self.get_water().keys()


    def get_hetatm(self):
        """ get HETATM atoms """
        hetatm = collections.OrderedDict()
        for het in self.hetatm:
            atom = Atom(het)
            uid = atom.resid+str(atom.resnum)
            if not hetatm.has_key(uid):
                hetatm[uid] = [het]
            else:
                hetatm[uid].append(het)
        return hetatm

    def get_hetatm_list(self):
        """ get HETATM uid list """
        return self.get_hetatm().keys()

    def get_residues(self):
        """ get residues dic with corresponding atoms """
        res = collections.OrderedDict()
        for line in self.lines:
            if RE_ATOM.search(line):
                atom = Atom(line)
                if len(atom.resid) == 3:
                    key = atom.resid+str(atom.resnum)+atom.code
                    if res.has_key(key):
                        res[key].append(line)
                    else:
                        res[key] = [line]
        return res

    def get_residue_list(self):
        """ get residue uid list """
        return self.residues.keys()

    def sequence(self):
        """ get chain sequence """
        seq = []
        for res in self.residues.keys():
            seq.append(AA[res[:3]])
        return seq

    def fasta(self):
        """ get chain sequence in fasta format """
        return split_seq(self.sequence(), 80)

#
#
#
#
#
#*******************************************************************************
# PROTEIN
#
#*******************************************************************************
#

class Protein(object):
    """Protein class

    parse atomic data from PDB
    Arguments: structural data (list)
    """
    NB_ATOM_PEPTIDE = 25

    def __init__(self, pdb_id, lines):
        """ Instantiate """
        self.pdbid = pdb_id.upper()
        self.protein_atoms = []
        self.chain_list = [] # all chains
        self.prot_chain_list = [] # only protein chains
        self.lines = lines
        self.water = []
        self.__pdb_content = {}
        self.chains = self.get_chains()
        self.__extract_elem_from_chain()


    def __check_oligo_peptide(self, lines):
        """
        verify if oligopeptide is well a peptide
        ** argument : atom list
        ** return   : boolean
        """
        ca_atoms = ['HETATM'+l[6:-1] for l in lines \
              if l[12:16].replace(' ', '') in ['CA', 'C']]
        ca_atoms.sort()
        for i in xrange(0, len(ca_atoms)-2):
            if euc_dist(ca_atoms[i], ca_atoms[i+1]) > 4:
                return False
        return True


    def __extract_elem_from_chain(self):
        """ parse atom lines and order in dict """
        for ch in self.chain_list:
            chain = Chain(self.chains[ch])
            self.__pdb_content[ch] = {'atom':None, 'hetatm':None,
                                      'peptide':None, 'water':None}
            if chain.hetatm and chain.atom:
                nb_all_res = chain.nb_res() + chain.nb_het_res()
                if nb_all_res > self.NB_ATOM_PEPTIDE:
                    self.__pdb_content[ch]['atom'] = chain.atom
                    self.protein_atoms += chain.atom
                    self.__pdb_content[ch]['hetatm'] = chain.hetatm
                    self.prot_chain_list.append(ch)
                else:
                    if self.__check_oligo_peptide(chain.atom+chain.hetatm):
                        self.__pdb_content[ch]['peptide'] = chain.atom+chain.hetatm
            elif chain.atom:
                if chain.nb_res() > self.NB_ATOM_PEPTIDE:
                    self.__pdb_content[ch]['atom'] = chain.atom
                    self.protein_atoms += chain.atom
                    self.prot_chain_list.append(ch)
                else:
                    self.__pdb_content[ch]['peptide'] = chain.atom
            if chain.water:
                self.__pdb_content[ch]['water'] = chain.water
            self.water += chain.water


    def get_protein_chains(self):
        protein_chains = {}
        for ch in self.prot_chain_list:
            if self.__pdb_content[ch]['atom']:
                protein_chains[ch] = self.__pdb_content[ch]['atom']
            if self.__pdb_content[ch]['hetatm']:
                protein_chains[ch] += self.__pdb_content[ch]['hetatm']
            if self.__pdb_content[ch]['water']:
                protein_chains[ch] += self.__pdb_content[ch]['water']
        return protein_chains

    def get_peptide(self):
        """ get peptide in prot """
        peptide = {}
        for ch in self.chain_list:
            if self.__pdb_content[ch]['peptide']:
                peptide['PEP'+ch] = self.__pdb_content[ch]['peptide']
        return peptide

    def get_peptide_list(self):
        """ get peptide uid list """
        return self.get_peptide().keys()

    def get_ligand(self, min_heavy_atom=6, cris_agents=False):
        """ get ligands from pdb, def ligand """
        ligands = {}
        for ch in self.chain_list:
            if self.__pdb_content[ch]['hetatm']:
                for het in self.__pdb_content[ch]['hetatm']:
                    atom = Atom(het)
                    uid = atom.resid+'-'+str(atom.resnum)+'-'+ch
                    if not ligands.has_key(uid):
                        ligands[uid] = [het]
                    else:
                        ligands[uid].append(het)
        if min_heavy_atom:
            for uid, hetatms in ligands.items():
                if nb_heavy_atoms(hetatms) < min_heavy_atom:
                    del ligands[uid]
        if not cris_agents:
            for uid, hetatms in ligands.items():
                if uid[:3] in CRIS_AGENTS:
                    del ligands[uid]
        return ligands

    def get_ligand_list(self, min_heavy_atom=6, cris_agents=False):
        """ get ligand uid list """
        return self.get_ligand(min_heavy_atom, cris_agents).keys()

    def export_ligand(self, path="./", prefix='_lig_', pdb_code=True, peptide=True,\
                      min_heavy_atom=6, cris_agents=False):
        """ export ligand in files """
        if pdb_code:
            prefix = self.pdbid+prefix
        ligands = self.get_ligand(min_heavy_atom, cris_agents)
        paths = []
        if peptide:
            ligands = dict(ligands.items()+self.get_peptide().items())
        for uid, lines in ligands.items():
            paths.append(path+prefix+uid+'.pdb')
            ligand_file = open(path+prefix+uid+'.pdb', 'w+')
            ligand_file.write('\n'.join(lines))
            ligand_file.close()
        return paths

    def __define_pocket(self, het_lines, cut_off, water_extend):
        """ define pocket by proximity """
        pocket_atoms = []
        for atom in self.protein_atoms:
            for het in het_lines:
                if euc_dist(atom, het) <= cut_off and atom not in pocket_atoms:
                    pocket_atoms.append(atom)
        if water_extend:
            water = []
            for het in het_lines:
                for hoh in self.water:
                    if euc_dist(het, hoh) <= cut_off and hoh not in water:
                        water.append(hoh)
                        for atom in self.protein_atoms:
                            if euc_dist(atom, hoh) <= cut_off \
                            and atom not in pocket_atoms:
                                pocket_atoms.append(atom)
            pocket_atoms += water
        return pocket_atoms

    def get_pocket(self, ligand='*', peptide=True, cut_off=5,\
                   include_ligand=False, water_extend=False, \
                   min_heavy_atom=6, cris_agents=False):
        """ define pocket in prot """
        if ligand != '*' and not type(ligand) is list \
        and ligand not in self.get_ligand_list():
            ligand = ask_to_user('Invalid ligand name', self.get_ligand_list())
        ligands = self.get_ligand(min_heavy_atom, cris_agents)
        if peptide:
            ligands = dict(ligands.items()+self.get_peptide().items())
        if type(ligand) is list:
            pocket = self.__define_pocket(ligand, cut_off, water_extend)
            if include_ligand:
                pocket += ligand
            return pocket
        else:
            pockets = {}
            if ligand == '*':
                for uid, het_lines in ligands.items():
                    pockets[uid] = self.__define_pocket(het_lines, cut_off, water_extend)
                    if include_ligand:
                        pockets[uid] += het_lines
            elif ligand in ligands:
                pockets[ligand] = self.__define_pocket(ligands[ligand], cut_off, water_extend)
                if include_ligand:
                    pockets[ligand] += ligands[ligand]
        return pockets

    def export_pocket(self, ligand='*', peptide=True, cut_off=5, \
            include_ligand=False, water_extend=False, path='./', \
            prefix='_pocket_', pdb_code=True, min_heavy_atom=6, \
            cris_agents=False):
        """ export pocket in files """
        pockets = self.get_pocket(ligand, peptide, cut_off, include_ligand, \
                                  water_extend, min_heavy_atom, cris_agents)
        paths = []
        if pdb_code: prefix = self.pdbid+prefix
        for uid, lines in pockets.items():
            paths.append(path+prefix+uid+'.pdb')
            pocket_file = open(path+prefix+uid+'.pdb', 'w+')
            pocket_file.write('\n'.join(lines))
            pocket_file.close()
        return paths

    def get_chains(self):
        """ get prot chains """
        chains = {}
        for line in self.lines:
            if RE_ATOM.search(line) or RE_HETATM.search(line):
                atom = Atom(line)
                if atom.chain not in self.chain_list:
                    self.chain_list.append(atom.chain)
                    chains[atom.chain] = []
                chains[atom.chain].append(line)
        return chains

    def chain(self, ch=None):
        """ build chain class """
        if ch in self.chain_list:
            return Chain(self.chains[ch])
        else:
            ch = ask_to_user('Invalid chain, select', self.chain_list)
            return Chain(self.chains[ch])

    def export_chain(self, ch='*', path='./', prefix='_', water=True, hetatm=True):
        """ export chain in file """
        paths = []
        if ch != '*' and ch not in self.chain_list:
            ch = ask_to_user('Invalid chain, select', self.chain_list)
        if ch == '*':
            for ch in self.chain_list:
                chain = Chain(self.chains[ch])
                if not water and not hetatm:
                    atoms = chain.atom
                elif not water and hetatm:
                    atoms = chain.atom + chain.hetatm
                elif not hetatm and water:
                    atoms = chain.atom + chain.water
                else:
                    atoms = self.chains[ch]
                paths.append(os.path.join(path, self.pdbid+prefix+ch+'.pdb'))
                chain_file = open(os.path.join(path, self.pdbid+prefix+ch+'.pdb'), 'w+')
                chain_file.write('\n'.join(atoms))
                chain_file.close()
        else:
            chain = Chain(self.chains[ch])
            if not water and not hetatm:
                atoms = chain.atom
            elif not water and hetatm:
                atoms = chain.atom + chain.hetatm
            elif not hetatm and water:
                atoms = chain.atom + chain.water
            else:
                atoms = self.chains[ch]
            paths.append(os.path.join(path, self.pdbid+prefix+ch+'.pdb'))
            chain_file = open(os.path.join(path, self.pdbid+prefix+ch+'.pdb'), 'w+')
            chain_file.write('\n'.join(atoms))
            chain_file.close()
        return paths


    def protein_center(self):
        """ define center of the protein """
        x_sum, y_sum, z_sum = 0, 0, 0
        for atm in self.protein_atoms:
            atom = Atom(atm)
            x_sum += atom.x
            y_sum += atom.y
            z_sum += atom.z
        N = len(self.protein_atoms)
        return x_sum/N, y_sum/N, z_sum/N


    def renumbering(self, display=False):
        """renumber atoms and residues"""
        new_lines = []
        atom_counter = 0
        res_counter = 0
        previous_atom = None
        flag = True
        for i, line in enumerate(self.lines):
            if RE_COORDLINE.search(line):
                atom = Atom(line)
                if i < len(self.lines)-1:
                    if RE_COORDLINE.search(self.lines[i+1]):
                        next_atom = Atom(self.lines[i+1])
                        if previous_atom:
                            if atom.res_uid != previous_atom.res_uid:
                                flag = True
                            else:
                                flag = False
                        atom_counter += 1
                        if (atom.name == 'N' and next_atom.name in ['CA', 'HN']) or flag:
                            res_counter += 1
                            if display:
                                print atom.res_uid, res_counter
                    previous_atom = Atom(line)
                new_lines.append(str(line[:6]+"%5d"+line[11:22]+"%4d"+' '+line[27:])
                        % (atom_counter, res_counter))
        return new_lines

    def export_renumbering_protein(self, path='./', suffix='_renum', display=False):
        """ export renumbering protein """
        new_lines = self.renumbering(display=display)
        protein_file = open(path+self.pdbid+suffix+'.pdb', 'w+')
        protein_file.write('\n'.join(new_lines))
        protein_file.close()
        return path+self.pdbid+suffix+'.pdb'

    def centering(self):
        """ translate coordinate around center """
        new_lines = []
        x_center, y_center, z_center = self.protein_center()
        for line in self.lines:
            if RE_COORDLINE.search(line):
                atom = Atom(line)
                x_new = round(atom.x - x_center, 3)
                y_new = round(atom.y - y_center, 3)
                z_new = round(atom.z - z_center, 3)
                new_lines.append(str(line[:30]+"%8.3f%8.3f%8.3f"+line[54:])
                        % (x_new, y_new, z_new))
        return new_lines

    def export_centering_protein(self, path='./', suffix='_center', \
                                 center_point=False):
        """ export the centering coordinates """
        new_lines = self.centering()
        if center_point:
            center = "ATOM      0  C   CEN A  0        0.000  0.000   0.000   1.00"
            new_lines.insert(0, center)
        protein_file = open(path+self.pdbid+suffix+'.pdb', 'w+')
        protein_file.write('\n'.join(new_lines))
        protein_file.close()
        return path+self.pdbid+suffix+'.pdb'

    def inertia(self):
        """ find inertia moment from atoms """
        N = len(self.protein_atoms)
        max1, max2, max3 = 0, 0, 0
        inertia1, inertia2, inertia3 = [], [], []
        for i in xrange(0, N-2):
            if Atom(self.protein_atoms[i]).name == 'CA':
                for j in xrange(i, N-1):
                    if Atom(self.protein_atoms[j]).name == 'CA':
                        d = euc_dist(self.protein_atoms[i], self.protein_atoms[j])
                        if d > max1:
                            max1 = d
                            inertia1 = [self.protein_atoms[i], self.protein_atoms[j]]
        return (inertia1, inertia2, inertia3)

    def export_inertia_points(self, path='./', suffix='_inertia'):
        """export atoms which constitute inertia moemnts """
        inertia_moments = self.inertia()
        inertia_file = open(path+self.pdbid+suffix+'.pdb', 'w+')
        for inertia in inertia_moments:
            inertia_file.write('\n'.join(inertia))
            inertia_file.write('\n')
        inertia_file.close()
        return path+self.pdbid+suffix+'.pdb'


    def nb_atom(self):
        """ number of ATOM in protein """
        sum_nb_atom = 0
        for ch in self.chain_list:
            sum_nb_atom += Chain(self.chains[ch]).nb_atom()
        return sum_nb_atom

    def nb_water(self):
        """ number of O water in protein """
        return len(self.water)

    def nb_hetatm(self):
        """ number of HETATM in protein """
        sum_nb_hetatm = 0
        for ch in self.chain_list:
            sum_nb_hetatm += Chain(self.chains[ch]).nb_hetatm()
        return sum_nb_hetatm

    def nb_res(self):
        """ number of residues in protein """
        sum_nb_res = 0
        for ch in self.chain_list:
            sum_nb_res += Chain(self.chains[ch]).nb_res()
        return sum_nb_res

    def nb_het_res(self):
        """ number of hetero residues in protein """
        sum_nb_het_res = 0
        for ch in self.chain_list:
            sum_nb_het_res += Chain(self.chains[ch]).nb_het_res()
        return sum_nb_het_res

    def sequence(self, chain='*'):
        """ protein sequence """
        seq = []
        if chain == '*':
            for ch in self.chain_list:
                seq += Chain(self.chains[ch]).sequence()
        else:
            if chain in self.chain_list:
                return Chain(self.chains[chain]).sequence()
            else:
                print chain, 'not in', self.chain_list
        return seq

    def fasta(self, chain='*', concat=False, split=True):
        """ protein sequence in fasta format """
        seq = ''
        if chain == '*':
            chain_list = self.chain_list
        elif chain == 'protein':
            chain_list = self.prot_chain_list
        else:
            chain_list = list(chain)
        if concat:
            seq += '>'+self.pdbid+'\n'
        for i, ch in enumerate(chain_list):
            if concat:
                seq += ''.join(Chain(self.chains[ch]).sequence())
            else:
                if split:
                    chain_seq = Chain(self.chains[ch]).fasta()
                else:
                    chain_seq = ''.join(Chain(self.chains[ch]).sequence())
                if chain_seq:
                    if re.search("[A-Za-z0-9]{4}_[A-Za-z]", self.pdbid):
                        seq += '>'+self.pdbid+'\n'
                    else:
                        seq += '>'+self.pdbid+'_'+ch+'\n'
                    seq += chain_seq
                    if i != len(self.chain_list):
                        seq += '\n'
        if concat:
            seq += '\n'
        return seq

#
#
#
#
#
#*******************************************************************************
# PDB FILE
#
#*******************************************************************************
#

class PDB(object):
    """PDB class

    Arguments: file or pdb id (char)
    """

    def __init__(self, pdb):
        """ Instantiate """
        self.path = ''
        self.atomic_data = []
        self.conect_table = []
        self.header = []
        self.remark = []
        self.compnd = []
        self.source = []
        self.keywds = []
        self.jrnl = []
        self.title = ''
        self.nmodels = None
        if not pdb:
            sys.exit('%s Instantiate PDB object with path or pdb ID' % ERR_PREF)
        else:
            self.__load(pdb)
        self.pdbid = self.path.split('/')[-1].split('.')[0]
        self.__parse()
        self.prot = Protein(self.pdbid, self.atomic_data)


    # PRIVATE
    def __parse(self):
        """ parse pdb file """
        def format_line(line, key):
            """ replace elem and space """
            line = line.replace(key, '')
            line = re.sub('^(( )*[0-9]{1,3})?( )+', '', line)
            line = re.sub('( )+$', '', line)
            line = re.sub('( )+', ' ', line)
            return line
        begin = False
        with open(self.path) as pdbfile:
            for line in pdbfile.xreadlines():
                line = line.replace('\n', '')
                if not begin:
                    self.header.append(line)
                if RE_ATOM.search(line) or RE_HETATM.search(line):
                    begin = True
                    self.atomic_data.append(line)
                elif RE_CONECT.search(line):
                    self.atomic_data.append(line)
                else:
                    if re.search('^TITLE', line):
                        self.title += format_line(line, 'TITLE')
                    elif re.search('^KEYWDS', line):
                        for key in format_line(line, 'KEYWDS').split(' '):
                            key = re.sub(r'[\(\)\[\],;]', '', key)
                            if key and key not in self.keywds and\
                                    re.search('[A-Za-z]', key):
                                self.keywds.append(key)
                    elif re.search('^COMPND', line):
                        self.compnd.append(format_line(line, 'COMPND'))
                    elif re.search('^JRNL', line):
                        self.jrnl.append(format_line(line, 'JRNL'))
                    elif re.search('^SOURCE', line):
                        self.source.append(format_line(line, 'SOURCE'))
                    elif re.search('^COMPND', line):
                        self.compnd.append(format_line(line, 'COMPND'))
                    elif re.search('^REMARK', line):
                        self.remark.append(format_line(line, 'REMARK'))
                    elif re.search('^ENDMDL', line):
                        if not self.nmodels:
                            self.nmodels = 1
                        else:
                            self.nmodels += 1


    def __load(self, pdb):
        """ load pdb file: from path or fetch"""
        # .. if is path
        if os.path.exists(pdb) and pdb[-4:] == '.pdb':
            self.path = pdb
        elif os.path.exists(pdb+'.pdb'):
            self.path = pdb+'.pdb'
        elif len(pdb) == 4 and re.search("[A-Za-z0-9]{4}", pdb):
            # .. if is a code : try download
            dl_path = fetch_pdb(pdb)
            if dl_path:
                self.path = dl_path
                print "%s pdb %s is downloaded" % (NFO_PREF, self.path)
        # .. else : error !
        else:
            sys.exit("%s PDB file doesn't exist or incorrect path: %s" % (ERR_PREF,pdb))


    def exp_method(self):
        """ Return the experimental procedure """
        method = re.compile(r"([A-Za-z\-]+ ?)+")
        for line in self.header:
            if line[:6] == 'EXPDTA':
                expdata_line = method.search(line[6:])
                return expdata_line.group(0)[:-1]
        if self.resolution():
            return 'X-RAY DIFFRACTION'
        return None


    def resolution(self):
        """ Return X-ray resolution """
        res_value = re.compile(r'[0-9]+(\.[0-9]+)?')
        for line in self.remark:
            if re.search(r".*RESOLUTION.+([0-9]+\.?)+.+ANGSTROM", line):
                line = re.sub(r"(.*RESOLUTION|ANGSTROM.*)", "", line)
                res = res_value.search(line)
                if res:
                    return float(res.group(0))
        return None


    def date(self):
        """ Return PDB date """
        value = re.compile(r"[0-9]{1,2}\-[a-zA-Z]{3}\-[0-9]{2,4}")
        for line in self.header:
            date_line = value.search(line)
            if date_line:
                return date_line.group(0)
        return None


#
#
###### END
######
######
