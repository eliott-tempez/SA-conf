#!/usr/local/bin/python

#
# A parser for PDB files
#
# Written (2001-2003) by P. Tuffery, INSERM, France
#
# No warranty of any kind is provided
#
# Addition made by J. Martin 2005
# Addition made by L. Regad 2014
# Addition made by E. Tempez 2025

import string
import sys
import os
import copy
import math
import gzip
import types
import subprocess
import re

from FileBasics import *
from Geo3DUtils import *

#GBINPATH="/home/leslieregad/Sauvegarde/sauvegardeAmpere/Prog/ENCODE_juliette/"
#GHMMPATH="/home/leslieregad/Sauvegarde/sauvegardeAmpere/Prog/ENCODE_juliette/Bestmodels/"

#GBINPATH="/home/lregad/Research/Projects/src/srcEncodage/"
#GHMMPATH="/home/lregad/Research/Projects/src/srcEncodage/Bestmodels/"

GBINPATH="/home/eltem/Documents/Cours/M2/projet_long/srcSA-Conf3/src/"
GHMMPATH="/home/eltem/Documents/Cours/M2/projet_long/srcSA-Conf3/src/Bestmodels/"



AA1 = "ACDEFGHIKLMNPQRSTVWY"
AA3 = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS"]
AA1seq = "ACDEFGHIKLMNPQRSTVWYXXXCXXMCXXYMDXECXXY"
AA3STRICT = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

RNA3 = ["U"]
DNA3 = ["A","T","G","C"]
SOLV = ["HOH","H2O","WAT","DOD"]

BBATMS = ["N","CA","C","O"]
NCHIS  = [0,1,2,3,2,0,2,2,4,2,3,2,0,3,5,1,1,1,2,2]

CHIATMS = [ \
	[], \
	[["N","CA","CB","SG"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","OD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","OE1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[], \
	[["N","CA","CB","CG"],["CA","CB","CG","ND1"]], \
	[["N","CA","CB","CG1"],["CA","CB","CG1","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","CE"],["CG","CD","CE","NZ"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","SD"], \
	 ["CB","CG","SD","CE"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","OD1"]], \
	[], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","OE1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","NE"],["CG","CD","NE","CZ"], \
	 ["CD","NE","CZ","NH1"]], \
	[["N","CA","CB","OG"]], \
	[["N","CA","CB","OG1"]], \
	[["N","CA","CB","CG1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]]] 

AASC=[["CB"], \
      ["CB","SG"], \
      ["CB","CG","OD1","OD2"], \
      ["CB","CG","CD","OE1","OE2"], \
      ["CB","CG","CD1","CD2","CE1","CE2","CZ"], \
      [],["CB","CG","ND1","CD2","CE1","NE2"], \
      ["CB","CG1","CG2","CD1"], \
      ["CB","CG","CD","CE","NZ"], \
      ["CD","CG","CD1","CD2"], \
      ["CB","CG","SD","CE"], \
      ["CB","CG","OD1","ND2"], \
      ["CB","CG","CD"], \
      ["CB","CG","CD","OE1","NE2"], \
      ["CB","CG","CD","NE","CZ","NH1","NH2"], \
      ["CB","OG"], \
      ["CB","OG","OG1","CG2"], \
      ["CB","CG1","CG2"], \
      ["CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"], \
      ["CB","CG","CD1","CD2","CE1","CE2","CZ","OH"]]

AABB=["N","CA","C","O"]
# modif JM

def resType(aName):
	if aName == "":
		return "Error"
	if aName in AA3:
		return AA3.index(aName)
	else:
		return "Error"

def aa3Type(aName):
	if aName == "":
		return "Error"
	if aName in AA3:
		return AA3.index(aName)
	else:
		return "Error"

def aa1Type(aName):
	if aName == "":
		return "Error"
	if aName in AA1:
		return AA1.index(aName)
	else:
		return "Error"

#
# a series of AA3 separated by blanks into aa1 string
#
def SEQREStoAA1(seqres):
	seq = ""
	aList = seqres.split()
	for aRes in aList:
		# print aRes
		if aRes in AA3STRICT:
			seq += AA1[AA3STRICT.index(aRes)]
		else:
			seq += "X"
		# print seq
	return seq

#
# This will convert an alignement into a selection mask.
# s1 and s2 must be 2 strings of identical lengths.
# gaps (indels) must be represented by '-'
#
def aln2mask(s1, s2):
	res = ""
	if len(s1) != len(s2):
		return res
	for i in range(len(s1)):
		if s1[i] == '-':
			continue
		if s2[i] == '-':
			res += '-'
		else:
			res += s1[i]
	return res

# 
# Recreate python 2's apply function
#
def apply(func, args, kwargs=None):
    return func(*args) if kwargs is None else func(*args, **kwargs)

#
# any PDB line
#
class PDBLine:
	def __init__(self, aLine=""):
		self.txt = aLine

	def __getitem__(self, key):
		if isinstance(key, slice):
			return self.txt[key.start:key.stop]
		return self.txt[key]

	def __repr__(self):
		return str(self.txt)

	def __len__(self):
		return len(self.txt)

	# back to list of lines
	def flat(self):
		return self.txt

	# header de la ligne
	def header(self):
		try:
			return self.txt[0:6].split()[0]
		except:
			return ""


#
# a PDB ATOM (HETATM) line
#
class atmLine(PDBLine):
	
	def __init__(self, aLine=""):
		if isinstance(aLine, atmLine):
			## print "atmLine from atmLine"
			self.txt = aLine.txt
		elif isinstance(aLine, PDBLine):
			## print "atmLine from PDBLine"
			self.txt = aLine.txt
		elif isinstance(aLine, str):
			## print "atmLine from string"
			self.txt = aLine
		else:
			self.txt = aLine

	def atmNum(self, anum=""):
		if anum != "":
			self.txt = f"{self.txt[:6]}{anum:5d}{self.txt[11:]}"
			return anum
		try:
			anum = self.txt[6:11].split()[0]
			return anum
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return "UNK"

	def atmName(self, aname=""):
		if aname != "":
			self.txt = f"{self.txt[:12]}{aname:>4}{self.txt[16:]}"
			return aname
		try:
			rnum = self.txt[12:16].split()[0]
			return rnum
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return "UNK"

	def alt(self, acode=""):
		if acode != "":
			self.txt = f"{self.txt[:16]}{acode}{self.txt[17:]}"
			return acode
		try:
			alt=self.txt[16]
			return alt
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return " "
		
	def resName(self, rName=""):
		if rName != "":
			self.txt = f"{self.txt[:17]}{rName:>3}{self.txt[20:]}"
			return rName
		try:
			rname = self.txt[17:20].split()[0]
			return rname
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return "UNK"
		
	def chnLbl(self, lbl=""):
		if lbl != "":
			self.txt = f"{self.txt[:21]}{lbl[0]}{self.txt[22:]}"
			return lbl
		try:
			lbl = self.txt[21]
			return lbl
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return "UNK"

	def resNum(self, rnum=""):
		if rnum != "":
			self.txt = f"{self.txt[:22]}{rnum:4d}{self.txt[26:]}"
			return rnum
		try:
			rnum = self.txt[22:26].split()[0]
			return rnum
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return "UNK"

	def icode(self, thecode=""):
		if thecode != "":
			self.txt = f"{self.txt[:26]}{thecode}{self.txt[27:]}"
			return thecode
		try:
			icode = self.txt[26]
			return icode
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return " "
		
	def xyz(self):
		try:
			x=self.txt[30:38].split()[0]
			y=self.txt[38:46].split()[0]
			z=self.txt[46:54].split()[0]
			return float(x), float(y), float(z)
		except ValueError:
			print("Incorrect ATOM line format for:", self.txt)
			return 0., 0., 0.

	def crds(self):
		return self.txt[30:54]
	
	def setcrds(self, x, y, z):
		self.txt = f"{self.txt[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{self.txt[54:]}"
		return 

	def q(self, aQ=""):
		if aQ != "":
			self.txt = f"{self.txt[:54]}{aQ:7.3f}{self.txt[61:]}"
			return aQ
		try:
			occ = self.txt[54:60]
			return occ
		except ValueError:
			return "       "
		
	def r(self, aR=""):
		if aR != "":
			self.txt = f"{self.txt[:61]}{aR:7.3f}{self.txt[68:]}"
			return aR
		try:
			occ = self.txt[61:68]
			return occ
		except ValueError:
			return "       "
		
	def occ(self, aOcc=""):
		if aOcc != "":
			self.txt = f"{self.txt[:54]}{aOcc:6.2f}{self.txt[60:]}"
			return aOcc
		try:
			occ = self.txt[54:60]
			return occ
		except ValueError:
			return "      "
		
	def tfac(self):
		try:
			tfac = self.txt[60:66]
			return tfac
		except ValueError:
			return "      "
		
	def segId(self):
		try:
			segId = self.txt[72:76]
			return segId
		except ValueError:
			return "    "
		
	def ele(self):
		try:
			ele = self.txt[76:78]
			return ele
		except ValueError:
			return "  "
		
	def chrg(self):
		try:
			chrg = self.txt[78:80]
			return chrg
		except ValueError:
			return "  "


## ========================================
## A series of PDB ATOM (HETATM) lines
## Considered as a set of residues
## Tabulation of residues is achieved
##
## atmList always return atmList,
## EXCEPT for __getitem__ when requesting in 1 residue
## where it is desirable to return atmLine
##
## atom lines accessible as: x.atms
##
## With this class, we are simply manipulating text
## No semantics associated
## ========================================
class atmList(atmLine):

	# instanciate
	def __init__(self, data="", chId="", hetSkip=0, verbose=0):
		#
		# Order of parsing is important (inheritance)
		#
		
		# from PDB: just retain PDB.data field
		if isinstance(data, PDB):
			if verbose == 2:
				print("atmList from PDB")
			self.list = data.data
		# from residue: just retain residue.list field
		elif isinstance(data, residue):
			if verbose == 2:
				print("atmList from residue")
			self.list = data.atms
		# from atmList: just propagate
		elif isinstance(data, atmList):
			if verbose == 2:
				print("atmList from atmList")
			self.list = data.list
		# from atmLine: just wrap
		elif isinstance(data, atmLine):
			## We force one line as a residue
			if verbose == 2:
				print("atmList  from atmLine")
			## print data
			self.list = []
			self.list.append(data)
		# from list: suppose a list of atomic lines
		elif isinstance(data, list):
			if verbose == 2:
				print("atmList from list")
			## print len(data)
			self.list = []
			for aLine in data:
				self.list.append(atmLine(aLine))

		else:
			if verbose == 2:
				print("atmList from unknown")
			self.list  = []

	def __len__(self):
		return len(self.list)

	# return atmList
	def __add__(self, new):
		print("__add__.atmList")
		return atmList(self.list[:] + new.list[:])

	def __getitem__(self, key):
		if isinstance(key, slice):
			return atmList(self.list[key.start:key.stop])
		return self.list[key]

	# del x[i] : ready for deletions !!
	def __delitem__(self, aPos):
		# delete old series
		aDex = self.rt[aPos][0]

		for aAtm in range(self.rt[aPos][0], self.rt[aPos+1][0]):
			del self.atms[aDex]
		self.resTab(0)

	# Managing x[i] = y : ready for mutations !!
	def __setitem__(self, aPos, new):
		del self[aPos]

		# insert new
		aDex = self.rt[aPos][0]
		for aAtm in new.atms:
			self.atms.insert(aDex, aAtm)
			aDex += 1
		self.resTab(0)

	# return string for visualization
	def __repr__(self):
		res = ""
		## print "atmList repr"
		for aAtm in self.list:
			res += str(aAtm)
		return res

	# back to list of lines
	def flat(self):
		res = []
		for aAtm in self.list:
			res.append(aAtm.flat())
		return res

	# Managing x.insert(i,y) : ready for insertions !!
	def insert(self, aPos, new):
		# insert new
		aDex = self.rt[aPos][0]
		for aAtm in new.atms:
			self.list.insert(aDex,aAtm)
			aDex = aDex + 1
		self.resTab(0)



	#
	# A series of coordinates of the range
	#
	def crds(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self.list)
		res = []
		for aAtm in range(ffrom, tto):
			res.append(atmLine(self.list[aAtm]).crds())
		return res

	#
	# A series of coordinates of the range
	#
	def xyz(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self.list)
		##print "atmList xyz", ffrom , tto
		if (ffrom == 0) and (len(self.list) == 1):
			return atmLine(self.list[ffrom]).xyz()
		else:
			res = []
			for aAtm in range(ffrom, tto):
				res.append(atmLine(self.list[aAtm]).xyz())
			return res

	#
	# center of geometry of a collection of atoms
	#
	def BC(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self)
		(x,y,z) = (0.,0.,0.)
		nAtm = 0.
		for aAtm in self[ffrom:tto].atms:
			(x1,y1,z1) = aAtm.xyz()
			x = x + x1
			y = y + y1
			z = z + z1
			nAtm = nAtm + 1.
		x = x / nAtm
		y = y / nAtm
		z = z / nAtm
		return (x,y,z)

	def oneChis(self):

		resTpe = resType(self.list[0].resName())
		if resTpe == "Error":
			return

		res = [AA3[resTpe]]
		for aChi in CHIATMS[resTpe]:
			aAtm = self.theAtm(aChi[0])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[1])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[2])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[3])
			if aAtm == []:
				return res
			a = self.theAtm(aChi[0]).xyz()
			b = self.theAtm(aChi[1]).xyz()
			c = self.theAtm(aChi[2]).xyz()
			d = self.theAtm(aChi[3]).xyz()
			res.append(apply(dihedral,a+b+c+d))
		return res

	def chis(self):

		res = []
		if len(self) == 1:
			res.append(self.oneChis())
			return res
		for aRes in range(len(self)):
			res.append(self[aRes].oneChis())
		return res

	def outChis(self):
		chis = self.chis()
		for i in chis:
			print(i[0], end=' ')
			for j in i[1:]:
				print(j, end=' ')
			print("")
			
	def atmPos(self, aName):
		for aPos in range(len(self.list)):
			if self.list[aPos].atmName() == aName:
				return aPos
		return "None"
		
	def Npos(self):
		for aPos in range(len(self.list)):
			if self.list[aPos].atmName() == "N":
				# print str(self[aPos])
				return aPos
		return "None"

	def CApos(self):
		for aPos in range(len(self.list)):
			if self.list[aPos].atmName() == "CA":
				# print str(self[aPos])
				return aPos
		return "None"

	def Cpos(self):
		for aPos in range(len(self.list)):
			if self.list[aPos].atmName() == "C":
				# print str(self[aPos])
				return aPos
		return "None"

	def Opos(self):
		for aPos in range(len(self.list)):
			if self.list[aPos].atmName() == "O":
				# print str(self[aPos])
				return aPos
		return "None"

	def out(self):
		pass

	def resName(self):
		return atmLine(self.list[0]).resName()

	def theAtm(self, atmName = ""):
		for aLine in self.list:
			if atmLine(aLine).atmName() == atmName:
				return atmLine(aLine)
		return []

	def isPDB(self):
		return 1
	#
	# write PDB or PDB chain(s) to file
	#
	def write(self, outName="", label="", hetSkip=0,verbose=0):
		if outName == "":
			f = sys.stdout
		else:
			f = open(outName,"w")

		f.write(f"HEADER {label} ({len(self)} residues)\n")
		for aAtm in self.list:
			f.write(f"{aAtm}")


	def oneHMMGeo(self, aCA):
		CA1x, CA1y, CA1z = self[aCA].xyz()
		CA2x, CA2y, CA2z = self[aCA+1].xyz()
		CA3x, CA3y, CA3z = self[aCA+2].xyz()
		CA4x, CA4y, CA4z = self[aCA+3].xyz()
		d1 = distance(CA1x, CA1y, CA1z, CA3x, CA3y, CA3z)
		d2 = distance(CA1x, CA1y, CA1z, CA4x, CA4y, CA4z)
		d3 = distance(CA2x, CA2y, CA2z, CA4x, CA4y, CA4z)
		x1, y1, z1 = vecteur(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
		x2, y2, z2 = vecteur(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
		x3, y3, z3 = vecteur(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
		d4 = mixtproduct(x1, y1, z1, x2, y2, z2, x3, y3, z3)
		d5 = distance(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
		d6 = distance(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
		d7 = distance(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
		return d1,d2,d3,d4,d5,d6,d7

## ========================================
## The ONE residue class
## ========================================
class residue(atmList):
	def __init__(self, data="", verbose=0):
		if data == "":
			self.atms = []
			self.type = None
			self.name = None
		else:
			if isinstance(data, residue): # atmList instance
				self.atms = data.atms
				self.type = self.rType()
				self.name = self.rName()
			elif isinstance(data, atmList): # atmList instance
				self.atms = data
				self.type = self.rType()
				self.name = self.rName()
			elif isinstance(data, atmLine): # atmList instance
				self.atms = atmList(data)
			else:
				self.atms = atmList(data)
				self.type = self.rType()
				self.name = self.rName()

	def __len__(self):
		return len(self.atms)

	def __repr__(self):
		return self.atms.__repr__()

	# managing x[i]
	def __getitem__(self, key):
		if isinstance(key, slice):
			if len(self.atms) == 1:
				return residue(self.atms)
			if key.stop > len(self.atms):
				key.stop = len(self.atms)
			return self.atms[key.start:key.stop]
		if isinstance(key, int):
			if key > len(self.atms):
				return None
			elif key < 0:
				if key + len(self.atms) < 0:
					return None
				else:
					return self.atms[key]
			return self.atms[key]
		elif isinstance(key, str):
			for iAtm in range(len(self.atms)):
				if self.atms[iAtm].atmName() == key:
					return self.atms[iAtm]

	# back to list of lines
	def flat(self):
		return self.atms.flat()
	
	def rName(self, name="", verbose=0):
		if name == "":
			return self.atms[0].resName()
		else:
			for atm in self.atms:
				atm.resName(name)

	def rNum(self, aNum="", verbose=0):
		if aNum == "":
			return self.atms[0].resNum()
		else:
			for atm in self.atms:
				atm.resNum(aNum)

	def riCode(self, icode="", verbose=0):
		if icode == "":
			return self.atms[0].icode()
		else:
			for atm in self.atms:
				atm.icode(icode)

	def rType(self, verbose=0):
		aName = self.atms[0].resName()
		if aName in AA3:
			return "AMINO-ACID"
		elif aName in RNA3:
			return "RNA"
		elif aName in DNA3:
			return "DNA"
		elif aName in SOLV:
			return "SOLVENT"
		else:
			return "HETERO"

	def chnLbl(self, lbl="", verbose=0):
		if lbl == "":
			return self.atms[0].chnLbl()
		else:
			for atm in self.atms:
				atm.chnLbl(lbl)

	def atmPos(self, aName):
		return self.atms.atmPos(aName)

	def hasAltAtms(self, verbose=0):
		BBAltAtm = "No"
		SCAltAtm = "No"
		for iAtm in range(len(self.atms)):
			aAtm = self.atms[iAtm]
##			print aAtm
			
			alt = aAtm.alt()

			if alt != ' ':
				isAlt = 1
				if any(char.isdigit() for char in aAtm.txt[12]):
					isAlt = 0
				if aAtm.txt[12] == ' ' and aAtm.txt[13] == 'H':
					isAlt = 0
				if isAlt == 0:
					continue
				theAtmTpe = aAtm.atmName()
				if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "O":
					BBAltAtm = "Yes"
				else:
					SCAltAtm = "Yes"
		return BBAltAtm, SCAltAtm


	#
	# return a selection of atoms
	# an atmList
	#
	def select(self, awhat=[""]):

		res = atmList()
		lAltBB_ca = self.listAltBB_ca()
		if len(lAltBB_ca) > 0 :
			alt_bb = lAltBB_ca[0]
		else:
			alt_bb = " "
		for iAtm in range(len(self.atms)):
			if awhat == [""]:
				res.list.append(atmLine(self.atms[iAtm].txt))
			else:
				if awhat[0] !=  "-":
					if awhat.count(self.atms[iAtm].atmName()) > 0 and (self.atms[iAtm].alt()==alt_bb):
						res.list.append(atmLine(self.atms[iAtm].txt))
				else:
					if awhat.count(self.atms[iAtm].atmName()) == 0 and (self.atms[iAtm].alt()== alt_bb):
						res.list.append(atmLine(self.atms[iAtm].txt))
					
		return res


	def listAltBB_ca(self):
		lAltBB = []
		if self.hasAltAtms()[0]:
			for iAtm in range(len(self.atms)):
				aAtm = self.atms[iAtm]
				if aAtm.atmName()=="CA":
					lAltBB.append(aAtm[16])
		lAltBB.sort()
		return lAltBB




	def tutu(self, verbose=0):
		toto = self.atms[0].atmName()
		return toto




	def BBAtmMiss(self, verbose=0):
		missp = []
		for atms in AABB:
			if self.atms.atmPos(atms) == "None":
				missp.append(atms)
				break
		if verbose:
			print(missp)
		return missp

## ========================================
## The PDB file parser
## p.chn("A") does not work !!
## ========================================
class PDB:

	def __init__(self, fname="", chId="", model=1, hetSkip=0, verbose=0):
		self.pdb_line = PDBLine()
		self.residue = residue()
  		
		if fname != "":
			if fname == None:
				return None
			elif isinstance(fname, PDB):                 # already a PDB instance
				self.info  = fname.info
				self.id    = fname.id
				self.data  = fname.data
				self.mdls  = fname.mdls
				self.atms  = fname.atms
				self.seq   = fname.seq
				self.seq3D = fname.seq3D
				self.ss    = fname.ss
				self.s2    = fname.s2
				self.nModel = fname.nModel
				self.mdls  = fname.mdls
				self.dbref = fname.dbref
				self.chns  = fname.chns
				self.setModel(model, verbose)
				self.resTab(verbose)

			# a flat series of text lines
			elif isinstance(fname, list):    # a list of atoms
				self.parse(fname, "", chId, hetSkip, verbose)
				self.setModel(model, verbose)
				self.resTab(verbose)

			# from disk file
			elif isinstance(fname, str):  # read file from disk
				self.load(fname, chId, hetSkip, PDBDIR="/home/ionesco/pdb/data/structures/", verbose = verbose)
				self.setModel(model, verbose)
				self.resTab(verbose)

	# return residue  or PDB (chains)
	def __getitem__(self, key):
		if isinstance(key, slice):
			res = self[key.start].flat()
			for i in range(key.start+1, key.stop):
				res += self[i].flat()
			return PDB(res)
		elif isinstance(key, int):
			return self.rt[key]
		elif isinstance(key, str):
			res = []
			for i in self:
				if str.count(key, i.chnLbl()) > 0:
					res = res + i.flat()
			return PDB(res)
		
	# return number of residue 
	def __len__(self):
		return len(self.rt)

	# merge two PDB
	def __add__(self,new):

		return PDB(self.flat() + new.flat())


	def __repr__(self):

		res = ""
		for aRes in self:
			res = res + aRes.__repr__()
		return res

	# back to a string
	# an internal vital function to transit from PDB
	# to atmList etc
	def flat(self):
		res = []
		for i in self:
			res = res + i.flat()
		return res

	#
	# read PDB or PDB chain(s) from disk file
	# chainId: may be a string of several accepted Ids,
	#          or a string starting with - to indicate a list of rejected Ids.
	# hetSkip : 1 to avoid all non peptidic residuesn 0 else
	# model : the number of the model to install (from 1)
	#
	def out(self, outName="", chainId="", hetSkip=0, fmode="w", verbose=0):
		res = self.__repr__()
		if outName == "":
			f = sys.stdout
		else:
			try:
				f = open(outName,fmode)
			except:
				sys.stderr.write(f"Failed to write to {outName}\n")
				return

		f.write(f"HEADER {self.id}\n")
		f.write(f"{res}")
		f.write("TER\n")
		f.flush()
		if f != sys.stdout:
			f.close()
		else:
			sys.stdout.flush()

	def xyzout(self, outName="", chainId="", hetSkip=0, fmode="w", verbose=0):
		res = []
		for aRes in self:
			res = res + aRes.atms.crds()
			
		if outName == "":
			f = sys.stdout
		else:
			try:
				f = open(outName, fmode)
			except:
				print("Failed to write to ", outName)

		for aCrd in res:
			f.write(f"{aCrd}\n")
		if f != sys.stdout:
			f.close()

	def xyz(self, outName="", chainId="", hetSkip=0, fmode="w", verbose=0):
		res = []
		for aRes in self:
			res = res + aRes.atms.crds()
			
		return res

	def load(self, fname, chainId="", hetSkip=0, PDBDIR="/home/ionesco/pdb/data/structures/", verbose=0, model=1):

 
		try:
			if verbose:
				print("Trying: ", fname)
			allPDB = gsimpleload(fname, 0)
		except IOError:

			pdbEntry = fname[:4]
			if chainId == "":
				chainId = fname[4:]
			
			try:
			## Experimental structure
				if verbose:
					print("Trying: ",PDBDIR+"/all/pdb/pdb"+pdbEntry+".ent.Z")
				allPDB = gsimpleload(PDBDIR+"/all/pdb/pdb"+pdbEntry+".ent.Z",0)
			except IOError:
				if verbose:
					print("Failed")
			## Model structure
				try:
					if verbose:
						print("Attempting: ",PDBDIR+"/models/current/pdb/"+pdbEntry[1:3]+"/pdb"+pdbEntry+".ent.Z")
					allPDB = gsimpleload(PDBDIR+"/models/current/pdb/"+pdbEntry[1:3]+"/pdb"+pdbEntry+".ent.Z",0)
				except IOError:
					if verbose:
						print('Sorry: PDB entry ', pdbEntry, 'not found')
					raise UnboundLocalError
			
	
		# Organize series of lines
		idName = fname
		if fname.find("/") != -1:
			idName = fname[fname.rindex("/")+1:]
		self.parse(allPDB, idName[:40]+"_"+chainId, chainId, hetSkip, verbose, model)
		return None

	#
	# Flat line format to PDB format
	#
	def parse(self, allPDB, id="", chainId="", hetSkip=0, verbose=0, model=1):

		if id == "":
			id = "unkwn"
		self.info  = []
		self.id    = id
		self.data  = []   # All ATOM DATA, N MODELS
		self.mdls  = []
		self.atms  = []
		self.seq   = []
		self.seq3D = []
		self.ss    = []
		self.s2    = []
		self.nModel = 0
		self.mdls.append(0)
		#self.dbref = "" ##2020
		self.dbref = []##2020
		self.chns  = ""

		for curLine in allPDB:

			aLine = PDBLine(curLine)

			#print items
			header = aLine.header()
			if header == "ATOM" or header == "HETATM":
				aLine = atmLine(aLine)
				OK = 0
				if chainId == "":
					OK = 1
				elif chainId[0] != '-':
					if aLine.chnLbl() in chainId:
						OK = 1
				else:
					if aLine.chnLbl() not in chainId:
						OK = 1
					
				if OK:
					if hetSkip:
						if AA3.count(aLine.resName()) > 0:
							self.data.append(aLine)
						elif hetSkip == 2:
							if SOLV.count(aLine.resName()) == 0:
								self.data.append(aLine)
					else:
						self.data.append(aLine)
			elif header == "TER":
				pass
			elif header == "HEADER":
				self.info.append(curLine)
			elif header == "COMPND":
				self.info.append(curLine)
			elif header == "SOURCE":
				self.info.append(curLine)
			elif header == "REMARK":
				self.info.append(curLine)
			elif header == "SEQRES":
				self.seq.append(curLine)
			elif header == "HELIX" or header == "SHEET" or header == "TURN":
				self.s2.append(curLine)
			elif header == "SSBOND":
				self.ss.append(curLine)
			elif header == "DBREF":
				#self.dbref = aLine  ##2020
				self.dbref.append(aLine) ##2020
			elif header == "ENDMDL":
				self.mdls.append(len(self.data))
				self.nModel = self.nModel+1
			else:
				self.info.append(curLine)
				
		self.nModel = self.nModel+1
		self.mdls.append(len(self.data))
		return self


	# tabulate residues
	def resTab(self, verbose):
		"PDB.resTab"

		start   = 1
		self.rt = []
		curResNum = "-1000"
		curResName = "XXX"
		curICode  = ""
		curChn = ""
		atmFrom = 0

		if len(self.atms) == 0:
			print("Empty PDB instance")
			return
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
		
			resName = aAtm.resName()
			resNum  = aAtm.resNum()
			iCode   = aAtm.icode()
			chn     = aAtm.chnLbl()
			if resNum != curResNum or resName != curResName or iCode != curICode or chn != curChn:
				curResNum = resNum
				curResName = resName
				curICode = iCode
				curChn = chn

				if start:
					start = 0
				else:
					self.rt.append(residue(self.atms[atmFrom:iAtm]))
					atmFrom = iAtm
		self.rt.append(residue(self.atms[atmFrom:iAtm+1]))
		if verbose:
			print(" Found ", len(self.rt), ' residues')


	#
	# How many models in the PDB ?
	#
	def nModels(self):
		return self.nModel


	#
	# Install current model
	def setModel(self, model=1, verbose=0):

		if model > self.nModels():
			print("Sorry: no model number ",model," (Total of ",self.nModels,")")
			return
		self.atms = []
		if verbose:
			print("Installing model ",model," (atoms ",self.mdls[model-1]," - ",self.mdls[model],")")
		for aLine in range(self.mdls[model-1], self.mdls[model]):
			self.atms.append(atmLine(self.data[aLine]))
		self.curMdl = model
		return
	
	#
	# what chains in the PDB ?
	#
	def chnList(self):
		curChn = ""
		self.chns = ""
		for aLine in range(len(self.atms)):
			if self.atms[aLine][21] not in self.chns:
				curChn = self.atms[aLine][21]
				self.chns = self.chns + curChn
		return self.chns

	#
	# what chains in the PDB ?
	#
	def nChn(self):
		if self.chns == "":
			return len(self.chnList())
		else:
			return len(self.chns)

	#
	# is there such a chain in the PDB file ?
	#
	def hasChn(self, chnId):
		if self.chns == "":
			return str.count(self.chnList(), chnId)
		else:
			return str.count(self.chns, chnId)


	#
	# extract particular chain(s) passed in string chainId
	# the default is to return all the chains
	#
	def chn(self,chainId="", hetSkip = 0):

		res = []
		for i in self:
			if chainId == "":
				res = res + i.flat()
			elif chainId[0] != '-' and str.count(chainId, i.chnLbl()) > 0:
				res = res + i.flat()
			elif chainId[0] == '-' and str.count(chainId, i.chnLbl()) == 0:
				res = res + i.flat()
		return PDB(res, hetSkip=hetSkip)

	#
	# the molecular type of chain(s) in string chainId
	#
	def chnType(self, chainId="", verbose=0):
		if chainId == "":
			chainId = self.chnList()
		res = []
		unres = []
		for aChain in chainId:
			theChain = self.chn(aChain)
			nAA  = 0
			nRNA  = 0
			nDNA  = 0
			nHET = 0
			nH2O = 0
			for i in range(0,len(theChain)):
				resName = theChain[i].rName()
				if AA3.count(resName) > 0:
					nAA = nAA +1
				elif RNA3.count(string.split(resName)[0]) > 0:
					nRNA = nRNA +1
				elif DNA3.count(string.split(resName)[0]) > 0:
					nDNA = nDNA +1
				elif SOLV.count(string.split(resName)[0]) > 0:
					nH2O = nH2O +1
				else:
					nHET = nHET + 1
					if verbose:
						if unres.count(resName) == 0:
							unres.append(resName)
							print(unres)
							print("Unknown residue type (1)",resName)

			if verbose:
				print("nAA : ",nAA," nNA : ",nDNA + nRNA," nHET : ",nHET)
			nOTHER = nHET + nDNA + nRNA
			if nOTHER < nAA:
				res.append("Protein")
			elif nAA > nDNA + nRNA:
				res.append("Protein")
			elif nRNA + nDNA > 0:
				if nRNA > 0:
					res.append("RNA")
				else:
					res.append("DNA")
			else:
				if nH2O > nHET:
					res.append("SOLVENT")
				else:
					res.append("HETERO")
		if len(chainId) == 1:
			return res[0]
		return res

	#
	# return a selection of (sub) residues
	#
	def select(self,rwhat=[""],awhat=[""]):
		res = []
		for i in self:
			if rwhat == [""]:
				res = res + i.select(awhat).flat()
			elif rwhat[0] !=  "-":
				if rwhat.count(i.rName()) > 0:
					res = res + i.select(awhat).flat()
			else:
				if rwhat.count(i.rName()) == 0:
					res = res + i.select(awhat).flat()
		## print res
		if res == []:
			return None
		return PDB(res)

	#
	# return a selection of (sub) residues for a structure
	# the mask (if specified) is a string of length to-from
	# positions corresponding to '-' will be discarded
	#
	def mask(self,ffrom=0,tto=-1,mask=""):
		res = []
		aPos = 0
		if tto == -1:
			tto = len(self)
		if (mask != "") and (len(mask) < tto-ffrom):
			tto = ffrom + len(mask)
			
		for i in range(ffrom,tto):
			if mask == "" or ((mask != "") and (mask[aPos] != '-')):
				res = res + self[i].flat()
			aPos = aPos + 1
		if res == []:
			return None
		return PDB(res)

	# retourne le champ Header du fichier PDB
	def header(self):
		title=''
		for Line in self.info:
			if Line[:6]=='HEADER':
				items = string.split(Line)
				for aItem in items[1:]:
					if str.count(aItem,"-") == 2:
						break
					if title != '':
						title = title + " "
					title = title + aItem
		return title



	# add by Leslie REGAD
	# retourne le titre du fichier pdb (Title champs)
	def Titre_pdb(self, verbose = 0):
		titre=""
		for Line in self.info:
			#if re.search('TITLE',Line):
			if Line[0:5] == "TITLE": 
				titre = titre+Line
			titreModif = titre.replace("\n","").replace("TITLE","")
			return(titreModif)
	


	# add by Leslie REGAD
	# retourne les HETATM contenu dans le fichier
	def get_HETATM(self, verbose=0):
		hetList=[]
		occHetList=[]

		for Line in self.info:
			if Line[0:4] == "HET ":
				hetList.append(Line.split()[1])
		hetListU = list(set(hetList))
		for het in hetListU : 
			occHetList.append(hetList.count(het))
		return(hetListU, occHetList)




	# add by Leslie REGAD
	# retourne les codes UniProt pour chaque chaine
	def get_UniProtId(self, verbose = 0):
		UNPList=[]
		for Line in self.dbref :
			List1 = str(Line).split()
			if List1[5]=="UNP":
				ph = List1[2]+":"+List1[6]+":"+List1[8]+"-"+List1[9]
				UNPList.append(ph)
			return(UNPList)


	# Nature du fichier
	def compound(self):
		title=''
		for Line in self.info:
			if Line[:6]=='COMPND':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title
				

	# Provenance de la molecule
	def source(self):
		title=''
		for Line in self.info:
			if Line[:6]=='SOURCE':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title

			
	# Auteur
	def author(self):
		title=''
		for Line in self.info:
			if Line[:6]=='AUTHOR':
				print(Line)
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title
				
	# KEYWDS lines
	def  keywords(self):
		keylist = ''
		for Line in self.info:
			if Line[:6]=='KEYWDS':
				keylist=keylist+Line[10:-1]

		aPos = 0
		OK = 1
		while str.find(keylist,'\'',aPos) != -1:
			aPos = str.find(keylist,'\'',aPos)
			afunc = keylist[0:aPos]+"\\"+keylist[aPos:]
			keylist = afunc
			aPos = aPos + 1
		return keylist

	# Date de creation du fichier
	def date(self):
		date=''
		for Line in self.info:
			if Line[:6]=='HEADER':
				items = string.split(Line)
				for aItem in items:
					if str.count(aItem,"-") == 2:
						date = aItem
				break
		if date != '':
			return date

		# If no creation date, try revision date
		return self.revdate()

	# Revision date (supposes last revision is first REVDAT
	def revdate(self):
		date=''
		for Line in self.info:
			if Line[:6]=='REVDAT':
				date=string.split(Line[13:22])[0]
				break
	       ## 	print 'creation fichier',date
		return date

	# method by which crds were generated
	def expmethod(self, verbose = 0):
		for Line in self.info:
			if Line[:6]=='EXPDTA':
				if str.find(Line,'X-RAY DIFFRACTION')!=-1:
					return 'X-RAY DIFFRACTION'
				if str.find(Line,'X-RAY POWDER DIFFRACTION')!=-1:
					return 'X-RAY POWDER DIFFRACTION'
				elif str.find(Line,'NMR')!=-1:
					return 'NMR'
				elif str.find(Line,'ELECTRON DIFFRACTION')!=-1:
					return 'ELECTRON DIFFRACTION'

				elif str.find(Line,'FIBER DIFFRACTION')!=-1:
					return 'FIBER DIFFRACTION'

				elif str.find(Line,'FLUORESCENCE TRANSFER')!=-1:
					return 'FLUORESCENCE TRANSFER'

				elif str.find(Line,'NEUTRON DIFFRACTION')!=-1:
					return 'NEUTRON DIFFRACTION'


				elif str.find(Line,'THEORETICAL MODEL')!=-1:
					return 'THEORETICAL MODEL'

				elif str.find(Line,'SYNCHROTRON')!=-1:
					return 'SYNCHROTRON'

				elif str.find(Line,'ELECTRON MICROSCOPY')!=-1:
					return 'ELECTRON MICROSCOPY'

				else:
					return ''
		# Suppose if resolution set: Xray
		if self.resolution() != -1.:
			return 'X-RAY DIFFRACTION'

	# Coordinates resolution
	def resolution(self, verbose = 0):
		resol = -1.
		for Line in self.info:
			if str.find(Line,'REMARK   2 RESOLUTION')!=-1:
				posMax=str.find(Line,'ANGSTROM')-1
				posMin=str.find(Line,'RESOLUTION')+11
				if posMax!=-1:
					try:
						resol=float(Line[posMin:posMax])
					except ValueError:
						pass
		return resol


	# R Value
	def rvalue(self, verbose = 0):
		R_VALUE = "NULL"
		checkRValue = 0

		for Line in self.info:
			if str.find(Line,'REMARK   3') != -1:
				# Case where it is on the next line !!
				if R_VALUE == "NULL" and ((checkRValue == 1) or (checkRValue == 2)):
					if checkRValue == 1:
						if str.find(Line,'.') != -1:
							# print Line
							pos=str.find(Line,'.')-1
							checkRValue == 0
							try:
								R_VALUE=float(Line[pos:pos+5])
							except ValueError:
								R_VALUE   = "NULL"
					elif checkRValue == 2:
						startPos = str.find(Line,'VALUE')
						if str.find(Line,'.', startPos) != -1:
							pos=str.find(Line,'.', startPos)-1
							toPos = pos+5
							# check for cases such as: 0.20.
							if str.count(Line,'.', pos,toPos) > 1:
								toPos = str.find(Line,'.', pos+2)
							try:
								R_VALUE=float(Line[pos:toPos])
							except ValueError:
								R_VALUE   = "NULL"
					checkRValue = 0
					
				# On one line ?
				if R_VALUE == "NULL" and (str.find(Line,' R ') != -1 or str.find(Line,'R VALUE') != -1 or str.find(Line,'R-VALUE') != -1 or str.find(Line,'R-FACTOR') != -1) and str.find(Line,'TEST') == -1 and str.find(Line,'FREE') == -1 and str.find(Line,'ESTIMATE') == -1 and str.find(Line,'BIN') == -1 and  str.find(Line,'ERROR') == -1:
					startPos = str.find(Line,'R VALUE')
					if startPos == -1:
						startPos = str.find(Line,'R-VALUE')
					if startPos == -1:
						startPos = str.find(Line,'R-FACTOR')
					if startPos == -1:
						if str.find(Line,' R '):
							checkRValue = 2
					if verbose:
						print(Line[:-1])
						print(Line[startPos:-1])
					if str.find(Line,'.', startPos) != -1:
						pos=str.find(Line,'.', startPos)-1
						toPos = pos+5
						# check for cases such as: 0.20.
						if str.count(Line,'.', pos,toPos) > 1:
							toPos = str.find(Line,'.', pos+2)
						try:
							R_VALUE=float(Line[pos:toPos])
							print(R_VALUE)
						except ValueError:
							if Line[pos] == 'O':
								try:
									R_VALUE=float(Line[pos+1:toPos])
								except ValueError:
									R_VALUE   = "NULL"
							else:
								R_VALUE   = "NULL"

					else:
						checkRValue = 1

		return R_VALUE

	# Missing residues
	def missingRes(self, Chain=None):
		misRes = {}

		# Get chunk with info
		for i in range(len(self.info)):
			has_missing = False
			if self.info[i].startswith('REMARK 465   M RES C SSSEQI'):
				has_missing = True
				startPos = i + 1
				break
		while self.info[i].startswith('REMARK 465'):
			i += 1
		endPos = i

		# Get missing residues
		if has_missing:
			for i in range(startPos, endPos):
				line_split = self.info[i].strip().split()
				res3 = line_split[2]
				res1 = SEQREStoAA1(res3)
				chain = line_split[3]
				if Chain is not None and chain != Chain:
					continue
				resNum = line_split[4]
				if chain not in misRes:
					misRes[chain] = []
				misRes[chain].append([int(resNum), res1])

		return misRes			
		

	def freervalue(self):

		FREE_R_VALUE   = "NULL"
		for Line in self.info:

			if str.find(Line,'FREE R VALUE') != -1 and str.find(Line,'TEST') == -1 and str.find(Line,'ESTIMATE') == -1 and str.find(Line,'BIN') == -1 and  str.find(Line,'ERROR') == -1:
				if str.find(Line,'.') != -1:
					pos=str.find(Line,'.')-1
					try:
						FREE_R_VALUE=float(Line[pos:pos+5])
					except ValueError:
						FREE_R_VALUE   = "NULL"
		return FREE_R_VALUE

	def seqres(self,chIds='NOTSPECIFIED'):
		
		if chIds == 'NOTSPECIFIED':
			chIds = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES':
					if Line[11] not in chIds:
						chIds = chIds + Line[11]

		if len(chIds) > 1:
			rs = []
		for chId in chIds:
			aseqres = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES' and Line[11] == chId:
					aseqres = aseqres + Line[19:70]+' '

			type = self.chnType(chId)

			if type == 'Protein':
				aa1seq=SEQREStoAA1(aseqres)
			elif type == "DNA" or type == "RNA":
				curseqres = string.split(aseqres)
				aa1seq = ""
				for i in curseqres:
					aa1seq = aa1seq + i[0]
			else:
				aa1seq = aseqres

			if len(chIds) > 1:
				rs.append(aa1seq)
			else:
				rs = aa1seq
		return rs

	#
	# Does the file contain only CAs ?
	#
	def CAonly(self,verbose=0):
		res="Yes"
		for aLine in self.data:
			#print aLine[12:15]
			if str.find(aLine[12:15],"CA")==-1:
				res="No"
				if verbose:
					print('pas uniquement les CA')
				return res
				break
		if verbose:
			print('uniquement les CA')
		return res
		

	def SCatmMiss(self, verbose = 0):
		SCatmMiss=""
		status ="No"
		nSCMiss = 0
		for i in range(0,len(self)):
			resName = self[i].rName()
			if AA3STRICT.count(resName) == 0:
			## Suppose nonstandard amino-acids are OK
				continue

			aaTpe = AA3STRICT.index(resName)
			chaine = ""
			for atm in self[i].atms:
				chaine=chaine+atm.atmName()+' '
			if verbose:
				print(chaine)
			missp = 0
			for atms in AASC[aaTpe]:
				if str.find(chaine,atms)==-1:
					missp = 1
					break
			if missp:
				status ="Yes"
				nSCMiss = nSCMiss+1
				resName = self[i].rName()
				resNum = self[i].rNum()
				icode  = self[i].riCode()
				lbl  = self[i].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				SCatmMiss = SCatmMiss+Res
	
		return 	nSCMiss, SCatmMiss		

	def BBatmMiss(self, verbose = 0):
		BBatmMiss=""
		status ="No"
		nBBMiss = 0
		for i in range(0,len(self)):
			resName = self[i].rName()
			if AA3STRICT.count(resName) == 0:
			## Suppose nonstandard amino-acids are OK
				continue

			aaTpe = AA3STRICT.index(resName)
			theCheck = self[i].BBAtmMiss()
			missp = 0
			if theCheck != []:
				missp = 1
			if (missp == 1) and (theCheck[0] == "O") and (self[i].atmPos("OXT") != None):
				missp = 0
			if missp:
				status ="Yes"
				if i == 0:
					status = "Ext"
				if i == len(self) -1 and status != "Yes":
					status = "Ext"
				if i > 0 and i < len(self) -1:
					status = "Yes"
				nBBMiss = nBBMiss+1
				resName = self[i].rName()
				resNum = self[i].rNum()
				icode  = self[i].riCode()
				lbl  = self[i].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				BBatmMiss = BBatmMiss+Res
	
		return 	nBBMiss, BBatmMiss		

	def hasAltAtms(self,verbose = 0):
		BBAltAtm = "No"
		SCAltAtm = "No"
		for i in self:
			BB, SC = i.hasAltAtms()
			if BB == "Yes":
				BBAltAtm = "Yes"
			if SC == "Yes":
				SCAltAtm = "Yes"
		return BBAltAtm, SCAltAtm
	

	def altAtmsResList(self,verbose = 0):
		nBBAltAtm = 0
		nSCAltAtm = 0
		BBAltAtm  = ""
		SCAltAtm  = ""
		for i in self:
			BB, SC = i.hasAltAtms()
			if BB == "No" and SC == "No":
				continue
			resName = i.rName()
			resNum = i.rNum()
			icode  = i.riCode()
			lbl  = i.chnLbl()
			if icode == ' ':
				icode = ''
			if lbl == ' ':
				lbl = ''
			resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "

			if BB == "Yes":
				nBBAltAtm = nBBAltAtm + 1
				BBAltAtm = BBAltAtm + resLabel
			if SC == "Yes":
				nSCAltAtm = nSCAltAtm + 1
				SCAltAtm = SCAltAtm + resLabel

		return nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm


	# 
	# Check if BB peptidic geometry is correct (distance)
	# THIS WILL NOT DETECT FRAGMENTS. IF MANY, THE GAPS ARE IGNORED
	# AND DO NOT RESULT IN "Bad" RETURN.
	# This allows to scan that all the fragments are correct at once.
	#
	def geomCheck(self,verbose=0):

		aN = "None"
		aC = "None"
		Cx, Cy, Cz = 0., 0., 0.
		BBGeoOK = "Ok"
		
		for aRes in self:
			# skip heteros
			if AA3.count(aRes.rName()) == 0:
				continue
			aN = aRes.atmPos("N")
			if aN != "None":
				Nx, Ny, Nz = aRes[aN].xyz()
				theN = aRes[aN]
			if aC != "None":
				if theN.chnLbl() == theC.chnLbl():
					aDist = distance(Nx, Ny, Nz, Cx, Cy, Cz)
					if aDist > 1.50 and aDist < 3.:
						if verbose:
							print("Poor peptidic bond of ",aDist," for ", theC.resName(), theC.resNum(), theN.resName(), theN.resNum())
						if BBGeoOK == "Ok":
							BBGeoOK = "Poor"
					elif aDist > 3.:
						if verbose:
							print("Bad peptidic bond  of ",aDist," for :", theC.resName(), theC.resNum(), theN.resName(), theN.resNum())
						BBGeoOK = "Bad"
			aC  = aRes.atmPos("C")
			if aC != "None":
				Cx, Cy, Cz = aRes[aC].xyz()
				theC = aRes[aC]

		return BBGeoOK

	# 
	# Check if BB peptidic geometry is correct (distance)
	# 
	def traceCheck(self, hetSkip = 0, verbose = 0):
		theTrace = self.select(awhat=["CA"])
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Ok"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""
		tracePB = ""

		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1][0].xyz()
			except ValueError:
				if verbose:
					print("Sorry: incorrect ATOM format for:", theTrace[aRes - 1])
					return CisWarning,"No"

			try:
				x2, y2, z2 = theTrace[aRes][0].xyz()
			except ValueError:
				if verbose:
					print(" Sorry: fname incorrect ATOM format for:", theTrace[aRes])
					return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist < 3.60: # CIS peptide
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if resName != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + resLabel
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + resLabel

			if aDist > 4.20: # mauvaise geometrie
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					print("Bad Trace for ",theTrace[aRes-1])

		return traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep


	def HMMGeo(self, theId, verbose = 0):
		theTrace = atmList(self.select(awhat=["CA"]).atms)
		dst7 = []
		for aCA in range(0,len(theTrace)-3):
			d1,d2,d3,d4,d5,d6,d7 = theTrace.oneHMMGeo(aCA)
			dst7.append("%s %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %3d %s" % (atmLine(theTrace[aCA]).resNum(), d1,d2,d3,d4,d5,d6,d7, len(self)-3, theId))
		return dst7

	def HMMTrace(self, theId):
		theTrace = atmList(self.select(awhat=["CA"]).atms)
		return theTrace 
	

# JM : je rajoute la possibilite d'utiliser un autre modele
	def HMMfrgEncode(self, theId="unknwn", MODEL='none', BINPATH=GBINPATH, HMMPATH = GHMMPATH, verbose=0):

		trace, tracePb, nCISPRO, CISPRO, nCISPep, CISPep = self.traceCheck()
		if trace == "Bad":
			sys.stderr.write("HMMEncode (%s): Sorry ! incorrect alpha carbon trace. (%s)\n" % (theId, tracePb))
			return []
		dst7 = self.HMMGeo(theId)
		
		# ici choix du modle d'encodage
		if MODEL=='none':
			cmd = os.path.join(BINPATH, "HMMPred") + " -iMdl " + os.path.join(HMMPATH, "27best-2.model") + " -idst stdin -noconfmat 2> /dev/null"
		else:
			cmd = BINPATH+"/HMMPred -iMdl "+MODEL+" -idst stdin -noconfmat 2> /dev/null"
		if verbose:
			print(cmd)

		# Prepare input for the subprocess
		input_data = f"{len(dst7)}\n"
		for i in dst7:
			input_data += f"{i}\n"

		# Run the subprocess
		process = subprocess.Popen(
			cmd,
			shell=True,
			stdin=subprocess.PIPE,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
			text=True
		)

		stdout, stderr = process.communicate(input=input_data)
		if process.returncode != 0:
			sys.stderr.write(f"HMMPred error: {stderr}\n")
			return []

		# Process the output
		rs = stdout.splitlines()
		return rs


	def HMMTraceCheck(self, theId="unknwn", BINPATH=GBINPATH, HMMPATH = GHMMPATH, verbose = 0):
		chList = self.chnList()

		res = []
		for chId in chList:
			if verbose:
				print("Encoding chaing \""+chId+"\"")
			curChn = self[chId]
			if verbose:
				print(curChn[0])
				print(curChn[-1])
			nFrg, frgList = curChn.frgList()
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				curdex = curdex + 1

				curId = theId
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1]  != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep = curChn[i[0]:i[1]+1].traceCheck(verbose)
				if traceOK == "Bad":
					return traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep
		return "Ok", tracePB, nCISPRO, CISPRO, nCISPep, CISPep


# JM je rajoute la possibilite d'utiliser un autre modele
	def HMMEncode(self, theId="unknwn", MODEL='none', BINPATH=GBINPATH, HMMPATH = GHMMPATH, verbose=0):
		chList = self.chnList()
		res = []
		for chId in chList:
			if verbose:
				print("Encoding chain \""+chId+"\"")
			curChn = self[chId]
			if verbose:
				print(curChn[0])
				print(curChn[-1])
			nFrg, frgList = curChn.frgList()
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				curdex = curdex + 1

				curId = theId
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1]  != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				lrs = curChn[i[0]:i[1]+1].HMMfrgEncode(curId, MODEL, BINPATH, HMMPATH, verbose)
				res.append(lrs)
		return res

	def HMMfrgRNum(self, verbose = 0):
		rs = []
		for i in self:
			rs.append(i.rNum()+i.riCode())
		return rs

	#
	# This will format residue numbers similarly to what is achied using HMMEncode
	#
	def HMMrNum(self, theId="unknwn", verbose=0):
		chList = self.chnList()

		res = []
		for chId in chList:
			curChn = self[chId]
			nFrg, frgList = curChn.frgList()
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				curdex = curdex + 1

				curId = theId
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1]  != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				lrs = curChn[i[0]:i[1]+1].HMMfrgRNum(verbose)
				res.append(lrs)
		return res


# Function added by JM, sept 2006
	#
	# This will format CA distances similarly to what is achied using HMMEncode
	#
	def HMMDist(self, theId="unknwn", verbose=0):
		chList = self.chnList()

		res = []
		for chId in chList:
			curChn = self[chId]
			nFrg, frgList = curChn.frgList()
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				print("newfrag")
				curdex = curdex + 1

				curId = theId
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1]  != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				dst7 = self.HMMGeo(theId)
				
				lrs = curChn[i[0]:i[1]+1].HMMfrgRNum(verbose)
				res.append(lrs)
		return res




	def HMMfrgChnLbl(self, verbose = 0):
		rs = []
		for i in self:
			rs.append(i.chnLbl())
		return rs

	#
	# This will format residue numbers similarly to what is achied using HMMEncode
	#
	def HMMChnLbl(self, theId="unknwn", verbose=0):
		chList = self.chnList()

		res = []
		for chId in chList:
			curChn = self[chId]
			nFrg, frgList = curChn.frgList()
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				curdex = curdex + 1

				curId = theId
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1]  != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				lrs = curChn[i[0]:i[1]+1].HMMfrgChnLbl(verbose)
				res.append(lrs)
		return res

	#
	# This will format sequence similarly to what is achied using HMMEncode
	#
	def HMMSeq(self, theId="unknwn", maxNCDist = 1.7, maxCADist = 4.1, verbose=0):
		chList = self.chnList()
		# liste de chaines 
		res = []
		for chId in chList:
			curChn = self[chId]
			# curChn extrait PDB de la chaine courante
			nFrg, frgList = curChn.frgList()
			# nFrg : nombre de fragments
			# frgList : liste qui donne les indices de debut et de fin des fragments
			#	 exemple, pour 1bxw [[0, 21], [22, 64], [65, 68], [69, 146], [147, 157], [158, 171]]
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				curdex = curdex + 1
				curId = theId
				# curId : toujours le nom du fichier PDB (je ne sais pas pourquoi il est dans la boucle)
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1] != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				
				lrs = curChn[i[0]:i[1]+1].select(awhat=["CA"]).frgseq(maxNCDist,maxCADist)
				if lrs != []:
					header = ["> "+curId+" "+str(len(lrs[0]))]
					rs = header + lrs
				else:
					rs = lrs
				res.append(rs)
		return res

	#
	# This will format sequence similarly to what is achied using HMMEncode
	#
	def HMMxyz(self, theId="unknwn", maxNCDist = 1.7, maxCADist = 4.1, verbose=0):
		chList = self.chnList()

		res = []
		for chId in chList:
			curChn = self[chId]
			nFrg, frgList = curChn.frgList()
			if chId == " ":
				chId = ""

			curdex = 0
			for i in frgList:
				curdex = curdex + 1

				curId = theId
				if theId == "unknwn":
					curId = theId+chId
				if chId != "" and curId[-1]  != chId:
					curId = curId+chId
				if nFrg > 1:
					curId = curId+"_"+str(curdex)
				lrs = curChn[i[0]:i[1]+1].select(awhat=["CA"]).xyz()
				header = ["> "+curId+" "+str(len(lrs))]
				rs = header + lrs
				res.append(rs)	
		return res


	#
	# determine fragments based on alpha carbon inter-atomic distance alone
	# 4.10 is default threshold
	#
	def chnCAFrgList(self, chId = "", maxDist = 4.10): #x = PDB("12as")

		if chId == "" and len(self.chnList()) > 1:
			print("PDB.chnFrgList() : cannot handle several chains as \""+self.chnList()+"\"")
			return []
		res = []
		oriRes = 0
		lastRes = 0
		nFrg = 0

		for aRes in range(1,len(self)):

			try:
				aaTpe = AA3.index(self[aRes-1].rName())
			except ValueError:
			# skip non amino acid
				continue

			aC = self[aRes-1].atmPos("CA")
			if aC == "None":
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = self[aRes-1][aC].xyz()
			CchnLbl = self[aRes-1].chnLbl()
			lastRes = aRes-1

			try:
				aaTpe = AA3.index(self[aRes].rName())
			except ValueError:
				continue
			
			aN = self[aRes].atmPos("CA")
			
			if aN == "None":
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue

			Nx, Ny, Nz = self[aRes][aN].xyz()
			NchnLbl = self[aRes].chnLbl()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			lastRes = aRes
			
			if aDist > maxDist or (CchnLbl != NchnLbl):
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(lastRes)
		res.append(lRes)
		nFrg = nFrg + 1
			
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def asOneChn(self,chnId = ' '):
		for aRes in range(0,len(self)):
			self[aRes].chnLbl(chnId)
			self[aRes].rNum(aRes+1)
		return PDB(self.flat())
			

	#
	# determine fragments based on inter-atomic distance C'-N
	# 1.70 is default threshold
	#
	def chnFrgList(self, chId = "", maxDist = 1.7): #x = PDB("12as")

		if chId == "" and len(self.chnList()) > 1:
			print("PDB.chnFrgList() : cannot handle several chains as \""+self.chnList()+"\"")
			return []
		res = []
		oriRes = 0
		lastRes = 0
		nFrg = 0

		for aRes in range(1,len(self)):

			# print aRes , "/", len(self)
			try:
				aaTpe = AA3.index(self[aRes-1].rName())
			except ValueError:
			# skip non amino acid
				continue

			aC = self[aRes-1].atmPos("C")
			if aC == "None":
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = self[aRes-1][aC].xyz()
			CchnLbl = self[aRes-1].chnLbl()
			lastRes = aRes-1


			try:
				aaTpe = AA3.index(self[aRes].rName())
			except ValueError:
				continue
			
			aN = self[aRes].atmPos("N")
			
			if aN == "None":
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue

			Nx, Ny, Nz = self[aRes][aN].xyz()
			NchnLbl = self[aRes].chnLbl()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			lastRes = aRes
			
			if aDist > maxDist or (CchnLbl != NchnLbl):
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		# lRes.append(len(self) - 1)
		lRes.append(lastRes)
		res.append(lRes)
		nFrg = nFrg + 1
			
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	#
	# This will return the number of fragments
	# and their boundaries
	# This will manage different chains
	#
	def frgList(self, maxNCDist = 1.7, maxCADist = 4.1, verbose = 0): #x = PDB("12as")

		# (3, [[0, 326], [327, 655], [656, 860]])
		res = []
		oriRes = 0
		nFrg = 0

		chnIds = self.chnList()

		curDex = 0
		for chId in chnIds:
			curChn = self.chn(chId)

			if self.chnType(chId) != "Protein":
				curDex = curDex + len(curChn)
				continue

			CAonly = curChn.CAonly()
			if CAonly=="No":
				curNFrg, curFrgList = curChn.chnFrgList(maxNCDist)
			else:
				curNFrg, curFrgList = curChn.chnCAFrgList(maxCADist)

			for i in range(0,len(curFrgList)):
				curFrgList[i][0] = curFrgList[i][0] + curDex
				curFrgList[i][1] = curFrgList[i][1] + curDex
				res.append(curFrgList[i])

			nFrg = nFrg + curNFrg
			
			curDex = curDex + len(curChn)
		return nFrg, res

	def nFrg(self, maxNCDist = 1.7, maxCADist = 4.1, verbose=0):
		nFrg, frgList = self.frgList(maxNCDist,maxCADist,verbose)
		return nFrg

	#
	# This will not check for fragments
	#
	def aaseq(self, verbose=0):
		res = ""
		unres = []
		for aRes in self:
			if aRes.rName() in AA3STRICT:
				res += AA1[AA3STRICT.index(aRes.rName())]
			elif aRes.rName() in AA3:
				rName = aRes.rName()
				if verbose:
					print("Unfrequent residue type: ", end=' ')
				if rName == "MSE": # seleno MET
					res = res+"M"
				elif rName == "CSE": # seleno CYS
					res = res+"C"
				elif rName == "FGL": # Formyl GLY
					res = res+"C"
				elif rName == "CEA": # SHYDROXY-CYS
					res = res+"C"
				elif rName == "TPQ": # 2,4,5-TRIHYDROXYPHE
					res = res+"Y"
				elif rName == "CGU": # GAMMA-CARBOXY-GLU
					res = res+"E"
				elif rName == "MHO": # Hydroxy-MET
					res = res+"M"
				elif rName == "IAS": # BETA-CARBOXY ASP
					res = res+"D"
				elif rName == "TYS": # SULFONATED TYROSINE
					res = res+"Y"
				else:
					res = res+'X'
			else:
				if aRes.rName() not in unres:
					unres.append(aRes.rName())
					print("Unknown residue type (2): ", aRes.rName())
		return res


	def frgseq(self, maxNCDist = 1.7, maxCADist = 4.1, verbose=0):

		res = []
		nFrg, frgList = self.frgList(maxNCDist,maxCADist,verbose)

		for i in frgList:
			res.append( self[i[0]:i[1]+1].aaseq())
		return res

	def SGList(self):
		SGList = []
		for aRes in self:
			if aRes.rName() == "CYS":
				lSGList = []
				for aAtm in aRes.atms:
					if aAtm.atmName() == "SG":
						lSGList.append(aAtm.xyz())
				if lSGList != []:
					SGList.append(lSGList)
		return SGList
	
	def nSSIntra(self):
		nSSBond = 0
		aSGList = self.SGList()
		for aRes1 in range(0,len(aSGList)):
			for aSG1 in range(0,len(aSGList[aRes1])):
				for aRes2 in range(aRes1+1,len(aSGList)):
					for aSG2 in range(0,len(aSGList[aRes2])):
						if apply(distance, aSGList[aRes1][aSG1]+aSGList[aRes2][aSG2]) < 2.35:
							nSSBond = nSSBond + 1
							break


		return nSSBond


	def isHalfCys(self, aRes):	
		if self[aRes].rName() != "CYS":
			return 0,0,0
		x = 0.
		y = 0.
		z = 0.
		isSet = 0
		for aAtm in range(0,len(self[aRes].atms)):
			if self[aRes].atms[aAtm].atmName() == "SG":
				x,y,z = self[aRes].atms[aAtm].xyz()
				isSet = 1
		if isSet == 0:
			return 0,0,0
		for aPos in range(0,len(self)):
			if self[aPos].rName() != "CYS":
				continue
			if aPos == aRes:
				continue
			for aAtm in range(0,len(self[aPos].atms)):
				if self[aPos].atms[aAtm].atmName() == "SG":
					x1,y1,z1 = self[aPos].atms[aAtm].xyz()
					if distance(x,y,z,x1,y1,z1) < 2.35:
						return 1, aPos, distance(x,y,z,x1,y1,z1)
		return 0,0,0


	def getAllRes(self, Chain=None):
		res = {}
		solved_dict = {}
		bfact_dict = {}
     
		# Extract List of solved residues
		if Chain is not None:
			chains = Chain
		else:
			chains = self.chnList()
		for chain in chains:
			solved_dict[chain] = []
			chain_data = self[chain]
   
			for residue in chain_data:
				res_name = residue.rName()
				res_num = residue.rNum()
				solved_dict[chain].append([int(res_num), SEQREStoAA1(res_name)])

		# Extract dict of missing residues
		missing_dict = self.missingRes(Chain)
  
		# Merge the dictionaries
		for key in set(solved_dict.keys()).union(missing_dict.keys()):
			solved_list = solved_dict.get(key, [])
			missing_list = missing_dict.get(key, [])
			solved_entries = [entry + ['solved'] for entry in solved_list]
			missing_entries = [entry + ['missing'] for entry in missing_list]

			# Combine the lists and sort by the first element (index)
			res[key] = sorted(solved_entries + missing_entries, key=lambda x: x[0])

		return dict(sorted(res.items()))
        

	def get_bfact(self, Chain=None):
    # Get the B-factors for the CA of each residue in the chains
		bfact_dict = {}
		if Chain is not None:
			chains = Chain
		else:
			chains = self.chnList()

		for chain in chains:
			bfact_dict[chain] = {}
			chain_data = self[chain]
   
			for residue in chain_data:
				res_num = int(residue.rNum())
				for atom in residue:
					if atom.atmName() == "CA":
						tfact = float(atom.tfac().strip())
						bfact_dict[chain][res_num] = tfact

		return bfact_dict





## ========================================
## Protein specific tools
## y = protein(x.chn("A"))
## ========================================
	
class protein(PDB):

	def __init__(self, data, chId = "", model = 1, hetSkip = 0, verbose = 0):
		if data != "":
			if isinstance(data,PDB):
				self.atms = data.atms
				self.rt   = data.rt
				self.nFrg    = 0
				self.resTypes(verbose)
				self.frgs = []
				self.chns  = data.chns

			elif isinstance(data,types.ListType):
				PDB.__init__(self,data, chId, model, hetSkip, verbose)
				self.resTab(verbose)
				self.resTypes(verbose)
				self.nFrg    = 0
				self.frgs = []
				self.chns  = ""
				self.id    = ""
				self.dbref = ""
			elif isinstance(data,atmList):
				self.atms = []
				for aLine in data.atms:
					self.atms.append(aLine)
				PDB.resTab(self,verbose)
				self.resTypes(verbose)
				self.nFrg    = 0
				self.frgs = []
				self.chns  = ""
				self.id    = ""
				self.dbref = ""
			elif isinstance(data,types.StringType):
				self.atms  = []
				self.info  = []
				self.seq   = []
				self.seq3D = []
				self.ss    = []
				self.s2    = []
				self.id    = ""
				self.dbref = ""
				self.chns  = ""
				self.nFrg    = 0
				self.frgs = []
				self.nModel = 0
				PDB.__init__(self,data, chId, model, hetSkip, verbose)
				self.resTab(verbose)
				self.resTypes(verbose)
	def resTypes(self, verbose = 0):
		self.tpe = []
		unres = []
		for aRes in range(0,len(self.rt) -1):
			aAtm = self.rt[aRes][0]
			aLine = self.atms[aAtm]
			if AA3.count(atmLine(aLine).resName()) != 0:
				idex = AA3.index(atmLine(aLine).resName())
				self.tpe.append(idex)
			else:
				if unres.count(atmLine(aLine).resName()) == 0:
					print("Unknown residue type (3): ",atmLine(aLine).resName())
					unres.append(atmLine(aLine).resName())
				self.tpe.append(-1)
				

	def frgList(self):
		res = []
		theBB = self.BB()
		oriRes = 0
		nFrg = 0
		for aRes in range(1,len(theBB)):
			if self.tpe[aRes-1] == -1 and atmList(theBB[aRes-1].atms).theAtm("C") == []:
				continue
			if atmList(theBB[aRes-1].atms).theAtm("C") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = atmList(theBB[aRes-1].atms).theAtm("C").xyz()
			if self.tpe[aRes] == -1 and atmList(theBB[aRes].atms).theAtm("N") == []:
				continue
			if atmList(theBB[aRes].atms).theAtm("N") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx,Ny,Nz = atmList(theBB[aRes].atms).theAtm("N").xyz()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			if aDist > 1.7:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(len(theBB) - 1)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def nFrgs(self):
		res = []
		theBB = self.BB()
		oriRes = 0
		nFrg = 0
		for aRes in range(1,len(theBB)):
			if self.tpe[aRes-1] == -1 and atmList(theBB[aRes-1].atms).theAtm("C") == []:
				continue
			if atmList(theBB[aRes-1].atms).theAtm("C") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = atmList(theBB[aRes-1].atms).theAtm("C").xyz()
			if self.tpe[aRes] == -1 and atmList(theBB[aRes].atms).theAtm("N") == []:
				continue
			if atmList(theBB[aRes].atms).theAtm("N") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx,Ny,Nz = atmList(theBB[aRes].atms).theAtm("N").xyz()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			if aDist > 1.7:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(len(theBB) - 1)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		return nFrg

	def trace(self,fname = "", chId = "", hetSkip = 0, altSel = " "):
		self.resTab(0)
		res = []
		for aRes in range(0,len(self.rt) -1):
			for aAtm in range(self.rt[aRes][0], self.rt[aRes+1][0]):
				aLine = self.atms[aAtm]
				if (hetSkip == 2) and (AA3STRICT.count(atmLine(aLine).resName()) == 0):
					break
				if hetSkip and (AA3.count(atmLine(aLine).resName()) == 0):
					break
				if atmLine(aLine).atmName() == "CA":
					res.append(aLine)
					break
		return atmList(res)


	def outSeq(self, label, hetSkip = 0, verbose = 0):
		theTrace = self.trace("","",hetSkip, verbose)
		seq = protein(theTrace).aaseq()
		print("> ",label,len(seq))
		while len(seq) > 0:
			print(seq[:80])
			seq = seq[80:]

	def outRawSeq(self, hetSkip = 0, verbose = 0):
		theTrace = self.trace("","",hetSkip, verbose)
		seq = protein(theTrace).aaseq()
		print(seq)

	def aaseq(self, verbose = 0):
		res = ""
		unres = []
		for aRes in self:
			if AA3STRICT.count(aRes[0].resName()):
				res = res + AA1[AA3STRICT.index(aRes[0].resName())]
			elif AA3.count(aRes[0].resName()):
				if verbose:
					print("Unfrequent residue type: ",aRes[0].resName())
				if aRes[0].resName() == "MSE": # seleno MET
					res = res+"M"
				elif aRes[0].resName() == "CSE": # seleno CYS
					res = res+"C"
				elif aRes[0].resName() == "FGL": # Formyl GLY
					res = res+"C"
				elif aRes[0].resName() == "CEA": # SHYDROXY-CYS
					res = res+"C"
				elif aRes[0].resName() == "TPQ": # 2,4,5-TRIHYDROXYPHE
					res = res+"Y"
				elif aRes[0].resName() == "CGU": # GAMMA-CARBOXY-GLU
					res = res+"E"
				elif aRes[0].resName() == "MHO": # Hydroxy-MET
					res = res+"M"
				elif aRes[0].resName() == "IAS": # BETA-CARBOXY ASP
					res = res+"D"
				elif aRes[0].resName() == "TYS": # SULFONATED TYROSINE
					res = res+"Y"
				else:
					res = res+'X'
			else:
				if unres.count(aRes[0].resName()) == 0:
					unres.append(aRes[0].resName())
					print("Unknown residue type (2): ",aRes[0].resName())
					print(unres)
		return res
	
	def frg(self,whatFrg, frgs = []):
		if frgs == [] and self.frgs == []:
			self.nFrg, self.frgs = self.frgList()
		return protein(self[self.frgs[whatFrg][0]:self.frgs[whatFrg][1]+1])
	
	def hasAltAtms(self,verbose):
		BBAltAtm = "No"
		SCAltAtm = "No"
		for aLine in self.atms:
			if aLine[16] != ' ':
				isAlt = 1
				if any(char.isdigit() for char in aLine[12]):
					isAlt = 0
				if aLine[12] == ' ' and aLine[13] == 'H':
					isAlt = 0
				if isAlt == 0:
					continue
				theAtmTpe = string.split(aLine[12:15])[0]
				if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "C":
					BBAltAtm = "Yes"
				else:
					SCAltAtm = "Yes"
		return BBAltAtm, SCAltAtm
	
	def altAtmsResList(self,verbose):
		nBBAltAtm = 0
		nSCAltAtm = 0
		BBAltAtm  = ""
		SCAltAtm  = ""
		for aPos in range(0,len(self)):
			curRes = self[aPos]
			for aLine in curRes.atms:
				if aLine[16] != ' ':
					isAlt = 1
					if any(char.isdigit() for char in aLine[12]):
						isAlt = 0
					if aLine[12] == ' ' and aLine[13] == 'H':
						isAlt = 0
					if isAlt == 0:
						continue
					theAtmTpe = string.split(aLine[12:15])[0]
					res    = aLine.resName()
					resNum = aLine.resNum()
					icode  = aLine.icode()
					lbl    = aLine.chnLbl()
					if icode == ' ':
						icode = ''
					if lbl == ' ':
						lbl = ''
					Res=res+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
					if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "C":
						nBBAltAtm = nBBAltAtm + 1
						BBAltAtm = BBAltAtm + Res
						break
					else:
						nSCAltAtm = nSCAltAtm + 1
						SCAltAtm = SCAltAtm + Res
						break

		return nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm

	# Check if all BB atoms are present
	def hasAllBBAtms(self,verbose):
		CAWarning = 0
		CWarning  = 0
		OWarning  = 0
		NWarning  = 0
		residuNameManquant=[]
		cp=0
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			if aRes.Npos() == "None":
				if aPos == 0:
					NWarning  = 1
				elif aPos == len(self) - 1:
					if NWarning < 1:
						NWarning  = 1
				else:
					NWarning  = 2
					cpt=1
					residuNameManquant.append(aPos)
					
			if aRes.CApos() == "None":
				if aPos == 0:
					CAWarning  = 1
				elif aPos == len(self) - 1:
					if CAWarning < 1:
						CAWarning  = 1
				else:
					CAWarning  = 2
					if cp==0:
						cp=1
						residuNameManquant.append(aPos)
			if aRes.Cpos() == "None":
				if aPos == 0:
					CWarning  = 1
				elif aPos == len(self) - 1:
					if CWarning < 1:
						CWarning  = 1
				else:
					CWarning  = 2
					if cp==0:
						cp=1
						residuNameManquant.append(aPos)
			if aRes.Opos() == "None":
				if aPos == 0:
					OWarning  = 1
				elif aPos == len(self) - 1:
					if OWarning < 1:
						OWarning  = 1
				else:
					OWarning  = 2
					if cp==0:
						cp=1
						residuNameManquant.append(aPos)
			
			cp=0

			
		if OWarning == 2 or NWarning == 2 or CAWarning == 2 or CWarning == 2:
			BBAtmMiss = "Yes"
		elif OWarning == 1 or NWarning == 1 or CAWarning == 1 or CWarning == 1:
			BBAtmMiss = "Ext"
		else:
			BBAtmMiss = "No"

		return BBAtmMiss,residuNameManquant


	# Check if BB peptidic geometry is correct (distance)
	def geomCheck(self,verbose):

		aN = "None"
		aC = "None"
		Cx, Cy, Cz = 0., 0., 0.
		BBGeoOK = "Ok"
		
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			aN = aRes.Npos()
			if aN != "None":
				Nx, Ny, Nz = aRes[aN].xyz()
				if aC != "None":
					aDist = distance(Nx, Ny, Nz, Cx, Cy, Cz)
					if aDist > 1.50 and aDist < 3.:
						if verbose:
							print("Poor peptidic bond of ",aDist," for ", aRes.resName(), aRes.resNum(), aRes.resName(), aRes.resNum())
						if BBGeoOK == "Ok":
							BBGeoOK = "Poor"
					elif aDist > 3.:
						if verbose:
							print("Bad peptidic bond  of ",aDist," for :", aRes.resName(), aRes.resNum(), aRes.resName(), aRes.resNum())
							BBGeoOK = "Bad"
			aC  = aRes.Cpos()
			if aC != "None":
				Cx, Cy, Cz = aRes[aC].xyz()

		return BBGeoOK

	# Check if BB peptidic geometry is correct (distance)
	def traceCheck(self, hetSkip =0, verbose=0):
		theTrace = self.trace("","",hetSkip, verbose)
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Yes"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""

		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1].xyz()
			except ValueError:
				if verbose:
					print(" Sorry: fname incorrect ATOM format for:", theTrace[aRes - 1])
					return CisWarning,"No"

			try:
				x2, y2, z2 = theTrace[aRes].xyz()
			except ValueError:
				if verbose:
					print(" Sorry: fname incorrect ATOM format for:", theTrace[aRes])
					return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist < 3.60: # CIS peptide
				res    = atmLine(self[aRes].atms[0]).resName()
				resNum = atmLine(self[aRes].atms[0]).resNum()
				#print i, resNum
				icode  = atmLine(self[aRes].atms[0]).icode()
				lbl  = atmLine(self[aRes].atms[0]).chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=res+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if theTrace[aRes][17:20] != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + Res
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + Res

			if aDist > 4.10: # mauvaise geometrie
				traceOK = "No"
				if verbose:
					print("Bad Trace for ",theTrace[aRes-1])

		return CisWarning, hasCisPRO, hasCisPEP, traceOK, nCISPRO, CISPRO, nCISPep, CISPep


	def BBAngles(self,aRes = -1000):
		res = []
		if aRes == -1000:
			rFrom = 0
			rTo = len(self)
		else:
			rFrom = aRes
			rTo = aRes+1
		for aPos in range(rFrom,rTo):
			phi = -1000.
			psi = -1000.
			ome = -1000.

			if aPos > 0:
				OK = 1
				aAtm = self[aPos-1].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos-1].theAtm("C").xyz()
					b = self[aPos].theAtm("N").xyz()
					c = self[aPos].theAtm("CA").xyz()
					d = self[aPos].theAtm("C").xyz()
					phi = apply(dihedral,a+b+c+d)
			if aPos < len(self) - 1:
				OK = 1
				aAtm = self[aPos].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("N")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos].theAtm("N").xyz()
					b = self[aPos].theAtm("CA").xyz()
					c = self[aPos].theAtm("C").xyz()
					d = self[aPos+1].theAtm("N").xyz()
					psi = apply(dihedral,a+b+c+d)
			if aPos < len(self) - 1:
				OK = 1
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("CA")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos].theAtm("CA").xyz()
					b = self[aPos].theAtm("C").xyz()
					c = self[aPos+1].theAtm("N").xyz()
					d = self[aPos+1].theAtm("CA").xyz()
					ome = apply(dihedral,a+b+c+d)
			res.append([phi,psi,ome])
		return res


	def SGList(self):
		SGList = []
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			if aRes[0].resName() == "CYS":
				lSGList = []
				for aAtm in range(0,len(aRes.atms)):
					if atmLine(aRes.atms[aAtm]).atmName() == "SG":
						lSGList.append(atmLine(aRes.atms[aAtm]).xyz())
				if lSGList != []:
					SGList.append(lSGList)
		return SGList
		
	def nSSIntra(self):
		
		nSSBond = 0
		aSGList = self.SGList()
		for aRes1 in range(0,len(aSGList)):
			for aSG1 in range(0,len(aSGList[aRes1])):
				for aRes2 in range(aRes1+1,len(aSGList)):
					for aSG2 in range(0,len(aSGList[aRes2])):
						if apply(distance, aSGList[aRes1][aSG1]+aSGList[aRes2][aSG2]) < 2.35:
							nSSBond = nSSBond + 1
							break


		return nSSBond


	def BB(self):
		res = []
		for aLine in self.atms:
			theName = aLine.atmName()
			if BBATMS.count(theName) > 0:
				res.append(aLine)
		return  protein(res)
			
	def SC(self):
		res = []
		for aLine in self.atms:
			theName = atmLine(aLine).atmName()
			if BBATMS.count(theName) == 0:
				res.append(aLine)
		return  PDB(res)			

	def HMMGeo(self, theId):
		print("appel fonction 2 ")
		theTrace = self.trace()
		self.dst7 = []
		for aCA in range(0,len(theTrace)-3):
			d1,d2,d3,d4,d5,d6,d7 = theTrace.oneHMMGeo(aCA)
			self.dst7.append("%s %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %3d %s" % (atmLine(theTrace[aCA]).resNum(), d1,d2,d3,d4,d5,d6,d7, len(self)-3, theId))
		return self.dst7


