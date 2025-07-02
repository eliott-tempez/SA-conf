#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import string

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# To modify
#-------------------------------------------------------------------------------

folder_out = 'saconf_out'


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

PDB_SEQ_FILE = 'pdb_list_encode.id'
POS_ALIGN_FILE = 'position_alignment.fasta2'
CSV_OUT = 'Position_description.csv'
PML_OUT = 'script_pymol.pml'


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------



def extractAAPos(file_param, vect_pos) :
   listPosAA = []
   entete = file_param[0].split(";")
   idx_neqAA = entete.index("neqAA")
   for ligne in file_param[1:]:
      liste = ligne.split()[0].split(";")
      neqAA = liste[idx_neqAA]
      if neqAA != "-":
          if float(neqAA) > 1 :
              pos = int(liste[0])
              real_pos = vect_pos[pos-1]
              if real_pos != '-':
                  listPosAA.append(real_pos)
   return(listPosAA)


def extractSLPos(file_param, vect_pos) :
   listPosSL = []
   entete = file_param[0].split(";")
   idx_neqSL = entete.index("neqSL")
   for ligne in file_param[1:]:
      liste = ligne.split()[0].split(";")
      neqSL = liste[idx_neqSL]
      if neqSL != "-":
          if float(neqSL) >= 1.5 :
              pos = int(liste[0])
              real_pos = vect_pos[pos-1]
              if real_pos != '-':
                  listPosSL.append(real_pos)
   return(listPosSL)


def extractSLPosCons(file_param, vect_pos) :
   listPosSLc = []
   entete = file_param[0].split(";")
   idx_neqSL = entete.index("neqSL")
   for ligne in file_param[1:]:
      liste = ligne.split()[0].split(";")
      neqSL = float(liste[idx_neqSL])
      if neqSL == 1 :
          pos = int(liste[0])
          real_pos = vect_pos[pos-1]
          if real_pos != '-':
              listPosSLc.append(real_pos)
   return(listPosSLc)



def extractSLPosWeakVar(file_param, vect_pos) :
   listPosSLwv = []
   entete = file_param[0].split(";")
   idx_neqSL = entete.index("neqSL")
   for ligne in file_param[1:]:
      liste = ligne.split()[0].split(";")
      neqSL = float(liste[idx_neqSL])
      if ( (neqSL >= 1) and (neqSL < 1.5)) :
          pos = int(liste[0])
          real_pos = vect_pos[pos-1]
          if real_pos != '-':
              listPosSLwv.append(real_pos)
   return(listPosSLwv)



def def_important_pos(file_param, vect_pos) :
   listPos = [[],[],[],[]]
   entete = file_param[0].split(";")
   idx_neqSL = entete.index("neqSL")
   idx_SS = entete.index("SS")
   for ligne in file_param[1:]:
      liste = ligne.split()[0].split(";")
      pos = int(liste[0])
      real_pos = vect_pos[pos-1]
      if real_pos != '-':
          neqSL = float(liste[idx_neqSL])
          ss = liste[idx_SS] 
          if neqSL == 1 :
              listPos[0].append(real_pos)
          if neqSL >= 1.5:
              if ( (ss=="h") or (ss=="s") or (ss=="l") ):
                  listPos[1].append(real_pos)
              if ( (ss=="l-h-s") or (ss=="h-l") or (ss=="s-l") or (ss=="h-s") ):
                  listPos[2].append(real_pos)
   return(listPos)






def gen_pml(args):
    file1 = open(os.path.join(folder_out, CSV_OUT),"r")
    file_param = file1.readlines()
    file1.close()

    file2 = open(os.path.join(folder_out, PDB_SEQ_FILE),"r")
    prot_ref = file2.readlines()[0].split()[0]
    file2.close()

    prot = prot_ref.split("_")[0]  ##LR: add
    ch = prot_ref.split("_")[1].lower()

    path_pdb = os.path.join(folder_out, "PDB")

    file3 = open(os.path.join(folder_out, POS_ALIGN_FILE),"r")
    vect_pos = file3.readlines()[1].split()[0].split('.')
    file3.close()

    posAA =  extractAAPos(file_param,  vect_pos)
    #posSL =  extractSLPos(file_param,  vect_pos)
    #posSLc =  extractSLPosCons(file_param,  vect_pos)
    #posSLwv =  extractSLPosWeakVar(file_param,  vect_pos)
    
    posSLc =  def_important_pos(file_param, vect_pos)[0]
    posSLv_noSS =  def_important_pos(file_param, vect_pos)[1]
    posSLv_withSS =  def_important_pos(file_param, vect_pos)[2]

    file_out = open(os.path.join(folder_out, PML_OUT),"w")
    #file_out.write("load "+os.path.join('PDB', prot_ref+".pdb")+", " + prot_ref +" \n")
    if os.path.isfile(os.path.join(path_pdb,prot+".pdb")):
        prot_ref = prot
        file_out.write("load "+os.path.join(path_pdb, prot+".pdb")+", " + prot_ref +" \n")          ##LR: cp+modif
        file_out.write("color grey, "+ prot_ref + "\n")
        file_out.write("select protch, chain "+ ch +", " +prot_ref +"\n")
        file_out.write("color cyan, protch" + "\n")
        file_out.write("hide line\n")
        file_out.write("show cartoon, "+ prot_ref +"\n")
        file_out.write("show surface, protch" +"\n")
        
    else:
        file_out.write("load "+os.path.join(path_pdb, prot_ref+".pdb")+", " + prot_ref +" \n")          ##LR: cp+modif
        file_out.write("color cyan, "+ prot_ref + "\n")
        file_out.write("show cartoon, "+ prot_ref +"\n")
        file_out.write("show surface, "+ prot_ref +"\n")

    file_out.write("remove resname hoh\n")
    file_out.write("set transparency, 0.3\n")


    #visualization of ligand
    #file_out.write("select ligand, not polymer, mut \n")     ##LR: comment
    file_out.write("select ligand, het \n")                   ##LR: copy and modif
    file_out.write("show sticks, ligand\n")
    file_out.write('util.cbay("ligand")\n')
    
    #visualization of the mutated aligned positions
    if len(posAA)!= 0:    ## LR: add
        file_out.write("select mut, chain "+ ch + " and resid " + str.join(posAA,"+") +" \n")
        file_out.write("show sticks, mut\n")                  ##LR: modify lines by sticks
        #file_out.write("color red, mut\n")                   ##LR: comment

    #visualization of the structurally conserved aligned positions
    if len(posSLc)!= 0:    ## LR: add
        file_out.write("select consPos, chain "+ ch + " and  resid " + str.join(posSLc,"+") + " \n")
        file_out.write("color green, consPos\n")

    #visualization of the structurally variable aligned positions without SS changes
    if len(posSLv_noSS)!= 0:    ## LR: add
        file_out.write("select varPos_NoSS, chain "+ ch + " and  resid " + str.join(posSLv_noSS,"+") + " \n")
        file_out.write("color slate, varPos_NoSS\n")

    #visualization of the structurally variable aligned positions with SS changes
    if len(posSLv_withSS)!= 0:    ## LR: add
        file_out.write("select varPos_withSS, chain "+ ch + " and  resid " + str.join(posSLv_withSS,"+") + " \n")
        file_out.write("color blue, varPos_withSS\n")


    #visualization of the weak structurally variable aligned positions
    #if len(posSLwv)!= 0:    ## LR: add
    #    file_out.write("select WeakvarPos, chain "+ ch + " and  resid " + str.join(posSLwv,"+") + " \n")
    #    file_out.write("color slate, WeakvarPos\n")


    #visualization of the strongly structurally variable aligned positions
    #if len(posSL)!= 0:    ## LR: add
    #    file_out.write("select varPos, chain "+ ch + " and  resid " + str.join(posSL,"+") + " \n")
    #    file_out.write("color blue, varPos\n")



    file_out.write("bg_colour white \n")
    #file_out.write("set line_width, 4 \n")
    file_out.write("set surface_quality, 2 \n")

    file_out.write("set cartoon_fancy_helices, 1 \n")
    file_out.write("set cartoon_fancy_sheets, 1 \n")
    file_out.write("set stick_radius, 0.2 \n")    ##LR: add
    #file_out.write("center \n")
    #file_out.write("ray 1200,1200 \n")
    #file_out.write("png "+ prot_ref+".png, dpi=600 \n")
    #file_out.write("quit \n")
    file_out.close()
    #os.chdir(args.output)
    #os.system('pymol script_pymol.pml &')
    #os.chdir(PATH)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    gen_pml(folder_out)
