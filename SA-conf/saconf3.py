#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" SA-conf
    a tool to analyse to compare fastly protein sequences and structures

"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SA-conf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# By Leslie Regad
# MTi, Molecule Therapeutique in silico
# Computational approaches applied to pharmacological profiling
#
# 25/06/2014
#
#
#
# Additions made by E. Tempez 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import os
import re
import sys
import time
import string
import argparse
import platform
import shutil
import urllib.request
import gzip
from shutil import which as find_executable

import PDB
import PDB6

from Bio import SeqIO
from Bio import ExPASy
from Bio import SwissProt
from shutil import copy2, rmtree
from shutil import which as find_executable
import subprocess
import xml.etree.ElementTree as ET


#pathsrcDesc ="/home/leslieregad/Dropbox/Projet/uroTest/sa-conf-beta/"
#sys.path.append(pathsrcDesc)
#import PDB6
#import PDB



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Arguments
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='SA-conf')
'''
parser.add_argument(
    "--install",
    help="Install missing dependencies",
    action='store_true'
)
'''
parser.add_argument(
    "-i", "--IDfile", required=True,
    help="Text file containing all pdb id",
    type=str
)
parser.add_argument(
    "-m", "--method", default="clustalw",
    help="Multiple alignment method (clustalw [default], tcoffee)",
    type=str
)
parser.add_argument(
    "-p", "--pdbpath",
    help="Path of the directory containing the pdb (default: download)",
    type=str
)
parser.add_argument(
    "-a", "--alignFile",
    help="Path to multiple alignment file",
    type=str
)
#parser.add_argument(
#    "-py", "--pymol",
#    help="Generate pymol (.pml) visualisation script",
#    action='store_true'
#)
parser.add_argument(
    "-o", "--output", default="./saconf_out",
    help="Path to output",
    type=str
)

parser.add_argument(
    "-u", "--uniprot",
    help="Text file containing the different submitted UniProt id",
    type=str
)



args = parser.parse_args()


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
PATH = os.path.abspath(os.path.dirname(sys.argv[0]))
GBINPATH = os.path.join(PATH,"src")
GHMMPATH = os.path.join(PATH, "src", "Bestmodels")
METHODS = ['clustalw', 'tcoffee']
FASTA_FILE = 'seq.aa'
AA_ALIGN_SWSS = 'AA_alignment_withUniProt.fasta'
AA_ALIGN_SWSS_TMP = 'AA_alignment_withUniProt_tmp.fasta'
AA_ALIGN_SWSS_TMP2 = 'AA_alignment.clustal'
AA_ALIGN_FILE = 'AA_alignment.fasta'
AA_ALIGN_FILE2 = 'AA_alignment.fasta2'
SL_ALIGN_FILE = 'SL_alignment.fasta2'
POS_ALIGN_FILE = 'position_alignment.fasta2'
PDB_ID_FILE = 'Listpdb_extSeq.id'
PDB_SEQ_FILE = 'pdb_list_encode.id'
POSDAT_FILE = 'list_seqPOS.dat'
AADAT_FILE = 'list_seqAA.dat'
CSV_OUT = 'Position_description.csv'
PML_OUT = 'script_pymol.pml'
CORRES_SWSS_ALI = 'corresp_Swss_align.csv'
COMPO_DATASET = "dataset_composition.csv"



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def install_dependencies():
    """ Install missing packages """
    print('Install missing packages: python-dev python-biopython clustalw')
    os.system('sudo apt-get install python-dev')
    os.system('sudo apt-get install clustalw')
    os.system('sudo easy_install -f http://biopython.org/DIST/ biopython')
    sys.exit("Installed!")



def check_args(args):
    """ Check input arguments """
    error = ""
    if not os.path.isfile(args.IDfile):
        error += str("[ERROR] IDfile %s\n --File doesn't exist\n" % args.IDfile)
    if args.method and args.method not in METHODS:
        error += str("[ERROR] method %s\n --method doesn't exist\n" % args.IDfile)
    if args.pdbpath and not os.path.isdir(args.pdbpath):
        error += str("[ERROR] IDfile %s\n --File doesn't exist\n" % args.IDfile)
    if args.alignFile and not os.path.isfile(args.alignFile):
        error += str("[ERROR] IDfile %s\n --File doesn't exist\n" % args.IDfile)
    if error:
        print("#", "-"*78)
        print("#", "-"*78)
        print("#", " "*34, "ERROR")
        print("#", "-"*78)
        sys.exit(error)



def header(args):
    """ Print user info """
    print("")
    print("#", "-"*78)
    print("#", "-"*78)
    print("#", " "*35, "SA-conf")
    print("#     ", end=' ')
    print("a tool to analyse to compare fastly protein sequences and structures")
    print("#", "-"*78)
    print("#")
#    print("#", os.getlogin(), "on", platform.platform())
    print("#")
    print("#", time.strftime("%d/%m/%Y"))
    print("#", time.strftime("%H:%M:%S"))
    print("#")
    print("# Arguments")
    print("# ---------")
    print("#", "IDfile:", args.IDfile)
    print("#", "method:", args.method)
    print("#", "pdbpath:", args.pdbpath)
    print("#", "alignFile:", args.alignFile)
    print("#", "output:", args.output)
    print("#", "uniprot:", args.uniprot)
    print("#")
    print("#", "-"*78)
    print("#", "-"*78)
    print("")



def getname(path):
    return os.path.basename(path).split('.')[0]



def get_pdb(pdb_id, pdb_chain, pdb_path, output, fileout):
    if pdb_chain:
        extract_chain = True
    else:
        extract_chain = False
    exist = False
    out = None
    if pdb_path:
        if pdb_chain:
            src = os.path.join(pdb_path, pdb_id+"_"+pdb_chain+".pdb")
            if os.path.isfile(src):
                copy2(src, os.path.join(output, pdb_id+"_"+pdb_chain+".pdb"))
                out = os.path.join(output, pdb_id+"_"+pdb_chain+".pdb")
                toto = get_pdb_info(out)
                fileout.write(";".join([pdb_id] + get_pdb_info(out)) + "\n")##06/03/15-LR : extract pdb information
                extract_chain = False
                exist = True
        if not pdb_chain or not exist:
            src = os.path.join(pdb_path, pdb_id+".pdb")
            if os.path.isfile(src):
                copy2(src, os.path.join(output, pdb_id+".pdb"))
                out = os.path.join(output, pdb_id+".pdb")
                fileout.write(";".join([pdb_id]+ get_pdb_info(out))+"\n")##06/03/15-LR : extract pdb information
    if not pdb_path or not out:
        out = PDB.fetch_pdb(pdb_id, output)
        if out != None :
            fileout.write(";".join([pdb_id]+ get_pdb_info(out))+"\n")##06/03/15-LR : extract pdb information
    if out and extract_chain:
        if pdb_chain:
            pdb_obj = PDB.PDB(out)
            #fileout.write(str.join([pdb_id]+ get_pdb_info(out),";")+"\n")##06/03/15-LR : extract pdb information
            out_ch = pdb_obj.prot.export_chain(
                ch=pdb_chain, path=output
                #,water=False, hetatm=False,
            )[0]
            #os.remove(out)            ##LR : comment this line to conserve the pdb
            out = out_ch
            
    return out




def get_pdb_info(pdb):
    pdb_obj = PDB6.PDB(pdb)
    #determination de la method experimentale
    method = pdb_obj.expmethod()
    if method == None:
        method = "none"
    #pour les structures cristallo determine la resolution
    if method == "X-RAY DIFFRACTION":
        reso = str(pdb_obj.resolution())
    else:
        reso = "NA"
    #pour les structures rmn extrait le nombre de modèle
    if method == "NMR": 
        nbrModel =  str(pdb_obj.nModels()-1)
    else:
        nbrModel = "NA"
    info = [method, reso,nbrModel]
    ###determination of the number of chain
    #allch = pdb_obj.chnList().split()
    #listch = list([0])
    listch= pdb_obj.chnList().split()
    ###determination of the AA sequence of each protein chain
    if len(listch) != 0:
        nbrChain = len(list(listch[0]))
        listch_taille = []
        info.append(str(nbrChain))
        for ch in list(listch[0]):
            AAseq = pdb_obj.chn(ch).aaseq() 
            listch_taille.append(ch+":"+str(len(AAseq)))
        info.append("/".join(listch_taille))
    else:
        info.append("0")
        info.append("NA")
    ###determine si il y a des HETATM
    hetListU, occHetList = pdb_obj.get_HETATM()
    if len(hetListU) == 0:
         info.append("NA")
    else:
         phr_list = []
         for het_i in range(len(hetListU)):
             phr_list.append(hetListU[het_i]+":"+str(occHetList[het_i]))
         info.append("/".join(phr_list))
    ###determine le code uniprot
    UNPList = pdb_obj.get_UniProtId()
    if len(UNPList) == 0:
        info.append("NA")
    else:
         info.append("/".join(UNPList))
        
    
    return(info)



def extract_fasta(pdb, output):
    fasta_file = open(os.path.join(output, FASTA_FILE), "a")
    pdb_obj = PDB.PDB(pdb)
    fasta = pdb_obj.prot.fasta(split=False)
    fasta_file.write(fasta)
    fasta_file.close()
    return(fasta)





#def get_pdb_list(pdb_file, pdb_path, output):
#    """ get pdb list and dl """
#    fasta_file = open(os.path.join(output, FASTA_FILE),"w")
#    pdb_list = []
#    uniprot_list = []
#    write_list_uniprot = []    #add by LR
#    write_list_pdb = []    #add by LR
#    output = os.path.join(output, "PDB")
#    os.mkdir(output)
#    with open(pdb_file) as pdbs:
#        for ligne in pdbs.xreadlines():   #add by LR
#            liste = ligne.split()         #add by LR
#            print liste
#            pdb = liste[0]                #add by LR
#            pdb = re.sub('[ \n]', '', pdb)
#            pdb = pdb.upper()
#            pdb_split = pdb.split('_')
#            if len(pdb_split) == 1:
#                pdb_id, pdb_chain = pdb, None
#            elif len(pdb_split) == 2:
#                pdb_id, pdb_chain = pdb_split
#            if len(pdb_id) == 4 and re.search("[A-Za-z0-9]{4}", pdb_id):
#                out = get_pdb(pdb_id, pdb_chain, pdb_path, output)
#                if out:
#                    pdb_list.append(out)
#                    print "[DL]", out
#                    fasta = extract_fasta(out, args)
#                    write_list_pdb.append(fasta)
#                else:
#                    print "[WARNING]", pdb_id, "not downloaded"
#            else:
#                print "[WARNING]", pdb_id, "incorrect code"
#          
#            if len(liste) > 1 :           #add by LR
#                uniprot = liste[1]      #add by LR
#                if re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",uniprot) : 
#                    if not uniprot in uniprot_list:
#                        uniprot_list.append(uniprot)
#                        phrase = extract_fasta_UniProt(uniprot)
#                        if phrase != "" : 
#                            write_list_uniprot.append(phrase)

#    if len(write_list_uniprot) > 0:
#        for ph in write_list_uniprot :
#            fasta_file.write(ph) 
#    else:
#        print "[WARNING] no uniProt id was precised"
#    if len(write_list_pdb) > 0:
#        for ph_pdb in write_list_pdb :
#            fasta_file.write(ph_pdb) 
#    return pdb_list





def get_pdb_list(pdb_file, pdb_path, output):
    """ get pdb list and dl """

    ###preparation du fichier de sortie dataset_composition.csv
    fileout = open(os.path.join(output, COMPO_DATASET), "a")   ##06/03/15 : add LR
    entete = ["pdb_id", "experimental method", "resolution", "number of models", "nbr of chain", "chain:length", "HETATM:Occ", "chn:UniProtId:pos"] ##06/03/15-LR :
    fileout.write(";".join(entete) + "\n")   ##06/03/15-LR : extract pdb information
    
    pdb_list = []
    output2 = os.path.join(output, "PDB")
    os.mkdir(output2)
    with open(pdb_file) as pdbs:
         for pdb in pdbs.readlines():
#        for ligne in pdbs.xreadlines():   #add by LR
#            liste = ligne.split()         #add by LR
#            pdb = liste[0]                #add by LR
            pdb = re.sub('[ \n]', '', pdb)
            pdb = pdb.upper()
            pdb_split = pdb.split('_')
            if len(pdb_split) == 1:
                pdb_id, pdb_chain = pdb, None
            elif len(pdb_split) == 2:
                pdb_id, pdb_chain = pdb_split
            if len(pdb_id) == 4 and re.search("[A-Za-z0-9]{4}", pdb_id):
                out = get_pdb(pdb_id, pdb_chain, pdb_path, output2, fileout)
                if out:
                    pdb_list.append(out)
                    print("[DL]", out)
                    extract_fasta(out, output)
                else:
                    print("[WARNING]", pdb_id, "not downloaded")
            else:
                print("[WARNING]", pdb_id, "incorrect code")
    fileout.close()
    return pdb_list



def extract_pdb(args):
    ###preparation du fichier de sortie seq.aa
    if args.uniprot:
        prepa_seq_uniprot(args.uniprot, args.output)

    pdb_list = get_pdb_list(args.IDfile, args.pdbpath, args.output)
    return pdb_list




def prepa_seq_uniprot(file_uniprot,output):
    uniprot_list = []
    write_list_uniprot = []    #add by LR
    write_list_pdb = []    #add by LR
    fasta_file = open(os.path.join(output, FASTA_FILE),"w")
    with open(file_uniprot,"r") as uniprotList :
        for uniprotId in uniprotList.xreadlines():
            uniprot = uniprotId.split()[0]
            if re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",uniprot) : 
                if not uniprot in uniprot_list:
                    uniprot_list.append(uniprot)
                    phrase = extract_fasta_UniProt(uniprot)
                    print(phrase)
                    if phrase != "" : 
                        write_list_uniprot.append(phrase)
    if len(write_list_uniprot) > 0:
        for ph in write_list_uniprot :
            fasta_file.write(ph) 
    fasta_file.close()






#add by Leslie
def extract_fasta_UniProt(uniprot):
    handle = ExPASy.get_sprot_raw(uniprot)
    record = SwissProt.read(handle)
    fasta = record.sequence
    if len(fasta) != 0 : 
        phrase = "> "+uniprot+"\n"+fasta+"\n"
        #print phrase
    else:
        print("[WARNING] no uniprot id is given")
        phrase = ""          
    return phrase





def convert_fasta(fasta_ali, fasta_ali_out, args):
    """
    filein : file that contains fasta sequences (80 characters by lines)
    fileout : output file that contains fasta sequence in only one line
    """
    file_id = open(os.path.join(args.output, PDB_ID_FILE),"w")
    filefasta2 = open(fasta_ali_out,'w')

    with open(fasta_ali) as fastas:
        lines = fastas.readlines()

    for i, line in enumerate(lines):
        if line[0] == ">":
            file_id.write(line[1:])
            if i == 0:
                filefasta2.write(line)
            else:
                filefasta2.write(current_seq+"\n")
                filefasta2.write(line)
            current_seq = ""
        else:
            current_seq += line.replace("\n","")
        if i+1 == len(lines):
            filefasta2.write(current_seq)

    file_id.close()
    filefasta2.close()



#def align_aa_seq(args):
#    fasta_file = os.path.join(args.output, FASTA_FILE)
#    output = os.path.join(args.output, "AA_alignment."+args.method)
#    output2 = os.path.join(args.output, "AA_alignment2."+args.method)

#    if args.method == "clustalw":
#        #align = ClustalwCommandline(os.path.join(PATH,"src/clustalw2"), infile=fasta_file, outfile = output)
#        align = ClustalwCommandline("clustalw", infile=fasta_file, outfile = output)
#        stdout, stderr = align()
#    elif args.method == "tcoffee":
#        align = TCoffeeCommandline("t_coffee", infile=fasta_file, outfile = output)
#        stdout, stderr = align()
#    #os.remove("seq.dnd")

    #supprime les prot UniProt dans alignement ###add by LR
#    listeUniProt,liste_Seq_uniProt = supp_Unirot_Seq(output, output2) ###add by LR

#    #convert to fasta
#    fasta_out = os.path.join(args.output, AA_ALIGN_FILE)
#    sequences = SeqIO.parse(output, "clustal")
#    count = SeqIO.write(sequences, fasta_out, "fasta")
#    #os.remove(output+".html")
#    #one line
#    convert_fasta(fasta_out,args)
#    pos_uniProt_conserv = supp_gap_align(args)  ###add LR

#    #creation du fichier qui contient la correspondance des codes UniProt et aligne
#    corresp_pos_ali_Swss(args, listeUniProt, liste_Seq_uniProt, pos_uniProt_conserv)
    


def align_aa_seq(args):
    fasta_file = os.path.join(args.output, FASTA_FILE)
    output = os.path.join(args.output, "AA_alignment."+args.method)
    output2 = os.path.join(args.output, "AA_alignment2."+args.method)
    output3 = os.path.join(args.output, AA_ALIGN_SWSS)
    output4 = os.path.join(args.output, AA_ALIGN_SWSS_TMP )  
    output5 = os.path.join(args.output, AA_ALIGN_SWSS_TMP2) 

    if args.method == "clustalw":
        #align = ClustalwCommandline(os.path.join(PATH,"src/clustalw2"), infile=fasta_file, outfile = output)
        cmd = ["clustalw", "-INFILE=" + fasta_file, "-OUTFILE=" + output]
        result = subprocess.run(cmd, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
        if result.returncode != 0:
            print("Error in ClustalW:", result.stderr)
    elif args.method == "tcoffee":
        cmd = ["t_coffee", "-infile=" + fasta_file, "-outfile=" + output]
        result = subprocess.run(cmd, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
        if result.returncode != 0:
            print("Error in tCoffee:", result.stderr)
    #os.remove("seq.dnd")

    sequences = SeqIO.parse(output, "clustal")
    count2 = SeqIO.write(sequences, output5 , "clustal")
    count = SeqIO.write(sequences, output3 , "fasta")
    convert_fasta(output3,output4,args)
    os.rename(output4,output3)


    motif_align(args, output, output2)




def motif_align(args, output, output2):
    """
    output : clustalw alignment avec seq uniprot
    output2 : clustalw alignement ss seq uniprot
    """
    #supprime les prot UniProt dans alignement ###add by LR
    listeUniProt,liste_Seq_uniProt = supp_Unirot_Seq(output, output2) ###add by LR

    #convert to fasta
    fasta_ali = os.path.join(args.output, AA_ALIGN_FILE)
    fasta_ali_out = os.path.join(args.output, AA_ALIGN_FILE2)   

    sequences = SeqIO.parse(output, "clustal")

    count = SeqIO.write(sequences, fasta_ali, "fasta")
    #os.remove(output+".html")
    #one line
    convert_fasta(fasta_ali, fasta_ali_out, args)
    pos_uniProt_conserv = supp_gap_align(args)  ###add LR

    #creation du fichier qui contient la correspondance des codes UniProt et aligne
    corresp_pos_ali_Swss(args, listeUniProt, liste_Seq_uniProt, pos_uniProt_conserv)
    
    #output.close()
    sequences.close()


###add by Leslie
def supp_Unirot_Seq(output, output2):
    """
    output : clustalw alignment avec seq uniprot
    output2 : clustalw alignement ss seq uniprot
    """
    fileout = open(output2,"w")
  
    file1 = open(output,"r")
    filein = file1.readlines()
    file1.close()

    listeUniProt = []
    liste_Num = []

    for ligne in filein: 
        liste = ligne.split()
        if len(liste) != 2:
            fileout.write(ligne)
        else:
          if ligne == "".join([" "]*76)+'\n':
            fileout.write(ligne)
          else:
            id_prot = ligne.split()[0]
            if not re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",id_prot):
                fileout.write(ligne)  
            else:
                seq = ligne.split()[1]
                if id_prot in listeUniProt:
                    idx = listeUniProt.index(id_prot)
                    liste_Num[idx] = liste_Num[idx]+seq
                else:
                    listeUniProt.append(id_prot)
                    liste_Num.append(seq)                 
    fileout.close()

    os.rename(output2, output)

    return(listeUniProt,liste_Num)


    

def supp_gap_align(args):
    
    file_nm1 = os.path.join(args.output, AA_ALIGN_FILE2)

    filein = open(file_nm1,'r')
    filefasta2 = filein.readlines()
    filein.close()

    #supprime les noms des prot
    list_seq = []
    list_prot = []
    for i in filefasta2:
        if not re.search(">",i):
            list_seq.append(i)
        else:
            list_prot.append(i)

    pos_supp = []
    pos_uniProt_conserv = []

    taille_seq = len(list_seq[0])
    motif_ref = "".join(["-"]*len(list_seq))
    for pos in range((taille_seq-1)):
        tmp = ""
        for seq in list_seq:
            tmp =  tmp + seq[pos] 
        if tmp == motif_ref:          
            pos_supp.append(pos)
        else:
            pos_uniProt_conserv.append(pos)   

    #reecris le fichier

    file_nm2 = os.path.join(args.output, "tmp.tmp")
    file_tmp = open(file_nm2,"w")
    
    for i in range(len(list_seq)):
        file_tmp.write(list_prot[i])
        seq = list(list_seq[i])
        newseq = ""
        for j in pos_uniProt_conserv  : 
            newseq = newseq+seq[j]
        file_tmp.write(newseq+"\n")

    file_tmp.close()

    os.rename(file_nm2,file_nm1) 

    return(pos_uniProt_conserv)


           
    
def corresp_pos_ali_Swss(args, listeUniProt, liste_Seq_uniProt, pos_uniProt_conserv):
    file_out = open(os.path.join(args.output, CORRES_SWSS_ALI),"w")

    file_out.write(",".join(["Aligned"]+listeUniProt)+"\n")

    list_Num_Swss = []
    #creation des numero de positions
    for seq in liste_Seq_uniProt:
        list_Num_Swss.append(range(1,(len(seq)+1)))

    for i in range(len(pos_uniProt_conserv)):
        pos = pos_uniProt_conserv[i]
        tmp = [str(i+1)]
        for seq_num in list_Num_Swss : 
            tmp.append( str(seq_num[pos]))
        file_out.write(",".join(tmp)+"\n")

    file_out.close()







def encode_file_to_seq(pdb_code, encode_path, args):
    file_aa = open(os.path.join(args.output, "list_seqAA.dat"), "a")
    file_sl = open(os.path.join(args.output, "list_seqLS.dat"), "a")
    file_pos = open(os.path.join(args.output, "list_seqPOS.dat"), "a")

    encode_file = open(encode_path, "r")
    encode_data = encode_file.readlines()
    encode_file.close()

    if len(encode_data) != 0:
        i = 0
        line = encode_data[i]
        while i <= (len(encode_data)-2):
            if line[0] == ">":
                list_aa = []
                list_sl = []
                list_pos = []
                entete = line
                i += 1
                line = encode_data[i]
                while (line[0] != ">") and (i <= (len(encode_data)-2)):
                    if i <= (len(encode_data)-2):
                        liste = line.split()
                        list_aa.append(liste[4])
                        list_sl.append(liste[5])
                        list_pos.append(liste[0])
                        i += 1
                        line = encode_data[i]
                    if i == (len(encode_data)-1):
                        liste = line.split()
                        list_aa.append(liste[4])
                        list_sl.append(liste[5])
                        list_pos.append(liste[0])
                if len(list_sl) != 0:
                    file_aa.write(entete)
                    file_sl.write(entete)
                    file_pos.write(entete)
                    file_aa.write("".join(list_aa)+"\n")
                    file_sl.write("".join(list_sl)+"\n")
                    file_pos.write(".".join(list_pos)+"\n")
    file_aa.close()
    file_sl.close()
    file_pos.close()



#def HMMSAencode(args, pdb_list):
#    output = os.path.join(args.output, "HMM-SA")
#    os.mkdir(output)
#    for pdb_path in pdb_list:
#        pdb_file = os.path.basename(pdb_path)
#        pdb_code = pdb_file.split(".")[0]
#        print "[Encode]", pdb_code,
#        pdb = PDB6.PDB(pdb_path, hetSkip=1)
#        geo = pdb.HMMGeo(pdb_path)
#        try:
#            encode = pdb.HMMEncode(theId=pdb_path, BINPATH=GBINPATH, HMMPATH=GHMMPATH)
#        except:
#            print("pbl dans l'encodage")
#        seq = pdb.HMMSeq(pdb_path)
#        ca = pdb.HMMxyz(pdb_path)
#        #print ca
#        index = pdb.HMMrNum(pdb_path)
#        encode_nfo = open(os.path.join(output, pdb_code+".info"), "w")
#        frag = 0
#        while frag < len(ca):
#            temp_ca = ca[frag] # coordonnees du fragment
#            temp_index = index[frag] # indices du fragment
#            temp_fasta = ''
#            j = 1
#            while j < len(seq[frag]):
#                temp_fasta = temp_fasta + seq[frag][j]
#                j += 1
#            temp_encodage = '--'
#            try:
#                temp_encodage = temp_encodage+encode[frag][1]
#            except:
#                print 'empty case'
#            temp_encodage = temp_encodage + '-'
#            if len(temp_fasta) < 4 or (len(temp_fasta) >= 4 and len(temp_fasta) == len(temp_encodage)):
#                i = 0
#                encode_nfo.write("%s\n" % temp_ca[i])
#                i = 1
#                while i < len(ca[frag]):
#                    try:
#                        encode_nfo.write("%5s " % temp_index[i-1])
#                        # decalage de -1 car index n'intere pas le nom du fragment
#                    except:
#                        print 'empty case in all info'
#                    encode_nfo.write(
#                        "%s %s %s %s %s\n"
#                        % (temp_ca[i][0:8], temp_ca[i][8:16], temp_ca[i][16:], temp_fasta[i-1], temp_encodage[i-1])
#                    )
#                    i += 1
#            else:
#                print 'length error'
#            frag += 1
#        encode_nfo.close()
#        encode_file_to_seq(pdb_code, os.path.join(output, pdb_code+".info"), args)


def HMMSAencode(args, pdb_list):    
    output = os.path.join(args.output, "HMM-SA")
    os.mkdir(output)
    for pdb_path in pdb_list:
        pdb_file = os.path.basename(pdb_path)
        pdb_code = pdb_file.split(".")[0]
        print("[Encode]", pdb_code)
        try:
            pdb = PDB6.PDB(pdb_path, hetSkip=1)
            geo = pdb.HMMGeo(pdb_path)
            encode = pdb.HMMEncode(theId=pdb_path, BINPATH=GBINPATH, HMMPATH=GHMMPATH)
            seq = pdb.HMMSeq(pdb_path)
            ca = pdb.HMMxyz(pdb_path)
            index = pdb.HMMrNum(pdb_path)
            encode_nfo = open(os.path.join(output, pdb_code+".info"), "w")
            frag = 0
            while frag < len(ca):
                temp_ca = ca[frag] # coordonnees du fragment
                temp_index = index[frag] # indices du fragment
                temp_fasta = ''
                j = 1
                while j < len(seq[frag]):
                    temp_fasta = temp_fasta + seq[frag][j]    
                    j += 1
                temp_encodage = '--'
                try:
                    temp_encodage = temp_encodage+encode[frag][1]
                except IndexError:
                    print('empty case')
                temp_encodage = temp_encodage + '-'
                if len(temp_fasta) < 4 or (len(temp_fasta) >= 4 and len(temp_fasta) == len(temp_encodage)):
                    i = 0
                    encode_nfo.write("%s\n" % temp_ca[i])
                    i = 1
                    while i < len(ca[frag]):
                        try:
                            encode_nfo.write("%5s " % temp_index[i-1])
                            # decalage de -1 car index n'intere pas le nom du fragment
                        except IndexError:
                            print('empty case in all info')
                        encode_nfo.write(
                            f"{temp_ca[i][0:8]} {temp_ca[i][8:16]} {temp_ca[i][16:]} {temp_fasta[i-1]} {temp_encodage[i-1]}\n"
                        )
                        i += 1
                else:
                    print('length error')
                frag += 1
            encode_nfo.close()
            encode_file_to_seq(pdb_code, os.path.join(output, pdb_code+".info"), args)
        except Exception as e:
            print(f"erreur dans l'encodage: {e}")


##LR: comment this function- this function must be removed

#def ExtractNumRes(pdbObj):
#  ###creation d'un vecteur qui contient les positions des residus
#  vectRes=[]
#  for res in pdbObj:
#     numres = res[0][23:27].split()[0]
#     vectRes.append(numres)
#  return(vectRes)


##LR: add a new function that replaces ExtractNumRes
def ExtractNumRes(ProtNm, pdbFile_nm):
    prot = ProtNm.split("_")[0]
    ch = ProtNm.split("_")[1]
    pdb_obj = PDB.PDB(pdbFile_nm)
    ch_obj = pdb_obj.prot.chain(ch)
    listRes = ch_obj.get_residue_list()
    listRes2 = []
    for res in listRes:
        listRes2.append(res[3:])
    return(listRes2)



def HmmfileParse(HMMfile,ch):
   vectpos= []
   vectLS = []
   for li in HMMfile:
     if re.search(">",li):
       pn = li.split()[1]
       protnm = pn.split('.pdb')[0]
       chnm = pn.split('.pdb')[1].split("_")[0]
     if not re.search(">",li):
       if chnm == ch:
          liste = li.split()
          vectpos.append(liste[0])
          vectLS.append(liste[5])
   return(vectpos,vectLS)



##LR: the following functions must be removed

#def ExtractLS(idx,vectRes , vectpos, vectLS):
#         #determine la position de ce residus dans la sequence alignee
#         posseq = vectRes[idx]
#         #recherche la position correspondante dans HMM seq et la lettre correspondante
#         if posseq in vectpos:
#           idxHMM = vectpos.index(posseq)
#           LS = vectLS[idxHMM]
#         else:
#           LS="-"
#         return LS





#def ExtractPOS(idx,vectRes , vectpos):
#    """
#    seqaa = sequence en AA
#    vectRes = liste avec le numero des positions du fichier pdb
#    vectpos = liste avec le numero des positions  du fichier hmm
#    vectLS = liste avec les lettres structuruales
#    """

#    #determine la position de ce residus dans la sequence alignee
#    posseq = vectRes[idx]
#    #recherche la position correspondante dans HMM seq et la lettre correspondante
#    if posseq in vectpos:
#        idxHMM = vectpos.index(posseq)
#        pos = vectpos[idxHMM]
#    else:
#        pos="-"
#    return pos





#def convertAALS(seqaa,vectRes , vectpos, vectLS):
#    """
#    seqaa = sequence en AA
#    vectRes = liste avec le numero des positions du fichier pdb
#    vectpos = liste avec le numero des positions  du fichier hmm
#    vectLS = liste avec les lettres structuruales

#    seqaa = seq
#    vectRes = NumRes
#    vectpos = PosVect
#    vectLS = LSVect
#    """

#    newseq=""
#    idx = -1
#    for ele in seqaa:
#      if ele =="-":
#         newseq = newseq+"-"
#         #idx += 1
#      if ele != "-":
#         idx += 1
#         BS = ExtractLS(idx,vectRes , vectpos, vectLS)
#         newseq = newseq+BS
#    return newseq



#def convertAAPOS(seqaa,vectRes , vectpos):
#    """
#    seqaa = sequence en AA
#    vectRes = liste avec le numero des positions du fichier pdb
#    vectpos = liste avec le numero des positions  du fichier hmm

#    seqaa = seq
#    vectRes = NumRes
#    vectpos = PosVect

#    """
#    newseq=[]
#    idx = -1
#    for ele in seqaa:
#      if ele =="-":
#        newseq.append("-")
#      if ele != "-":
#        idx += 1
#        BS = ExtractPOS(idx,vectRes , vectpos)
#        newseq.append(BS)
#    newseq2 = str.join(newseq,".")
#    return newseq2




##LR: I added this new function

def convertAAPOSpdb(seqaa,vectRes):
    """
    seqaa = sequence en AA
    vectRes = liste avec le numero des positions du fichier pdb
    vectpos = liste avec le numero des positions  du fichier hmm

    seqaa = seq
    vectRes = NumRes
    vectpos = PosVect

    """
    newseq=[]
    idx = 0
    for ele in seqaa:
      if ele =="-":
        newseq.append("-")
      if ele != "-":
        newseq.append(vectRes[idx])
        idx = idx +1
    return newseq


##LR: I added this new function

def convertAA_POS_SL(seqaa,NumRes,PosVect, LSVect):
    """
    seqaa = sequence en AA
    NumRes = liste avec le numero des positions extrait du fichier pdb
    PosVect = liste avec le numero des positions  du fichier hmm
    LSVect  =  liste des SL extraites du fichier hmm
    """
    newseq = convertAAPOSpdb(seqaa,NumRes)
    alignSL_vect = []
    alignPos_vect = []
    for num in newseq:
        if num != "-" :
            if (num in PosVect) :
                idx = PosVect.index(num)
                alignPos_vect.append(PosVect[idx])
                alignSL_vect.append(LSVect[idx])
            else :
                alignPos_vect.append("-")
                alignSL_vect.append("-")
        else :
            alignPos_vect.append("-")
            alignSL_vect.append("-")


    alignSL = "".join(alignSL_vect)
    alignPos = ".".join(alignPos_vect)
    return(alignPos, alignSL)



def align_sl_seq(args):
    #if args.alignFile:
    #    copy2(args.alignFile, os.path.join(args.output, AA_ALIGN_FILE2))
    file_nm = os.path.join(args.output, AA_ALIGN_FILE2)
    fileout_nm = os.path.join(args.output, SL_ALIGN_FILE)
    fileout_nm2 = os.path.join(args.output, POS_ALIGN_FILE)
    pdb_dir = os.path.join(args.output, "PDB")
    dir_hmmsa = os.path.join(args.output, "HMM-SA")

    file1 = open(file_nm,"r")
    aliFile = file1.readlines()
    file1.close()

    fileres=open(fileout_nm,"w")
    fileres2=open(fileout_nm+"_wo_header.aln","w")
    fileres3=open(fileout_nm2,"w")
    fileres4=open(file_nm+"_wo_header.aln","w")
    pdb_id_out = open(os.path.join(args.output, PDB_SEQ_FILE), "w")
    for i in range(len(aliFile)):
        ligne = aliFile[i][:-1]
        if re.search(">",ligne):
          ProtNm = re.sub("[>\n]", "", ligne)
          if not re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",ProtNm) : 
            pdb_id_out.write(ProtNm+"\n")
            prot = ProtNm.split("_")[0]
            ch = ProtNm.split("_")[1]
            if os.path.isfile(os.path.join(pdb_dir, ProtNm+".pdb")):
                pdbFile_nm = os.path.join(pdb_dir, ProtNm+".pdb")
            elif os.path.isfile(os.path.join(pdb_dir, prot+".pdb")):
                pdbFile_nm = os.path.join(pdb_dir, prot+".pdb")
            else:
                print("[ERROR] path not found")
            #print pdbFile_nm, ch
            if os.path.isfile(pdbFile_nm):
                ###transfo le pdbFile en objet PDB
                pdb_obj = PDB6.PDB(pdbFile_nm)
                ###extraction du pdb de la chain
                pdb_objch = pdb_obj.chn(ch)
                ###extraction du vecteur des numero de res
                #NumRes = ExtractNumRes(pdb_objch)
                NumRes = ExtractNumRes(ProtNm, pdbFile_nm)
                ###ouverture du fichier HMM27
                #determine la nom de la prot ss la chaine
                if os.path.isfile(os.path.join(dir_hmmsa, prot+".info")):
                    fileHMM_nm = os.path.join(dir_hmmsa, prot+".info")
                elif os.path.isfile(os.path.join(dir_hmmsa, ProtNm+".info")):
                    fileHMM_nm = os.path.join(dir_hmmsa, ProtNm+".info")
                if os.path.isfile(fileHMM_nm) :

                    file2 = open(fileHMM_nm,"r")
                    fileHMM = file2.readlines()
                    file2.close()

                    ###extraction de la sequence de LS et de pos
                    PosVect, LSVect = HmmfileParse(fileHMM,ch)
                    fileres.write(ligne+"\n")
                    fileres3.write(ligne+"\n")
                    ligne = aliFile[i+1]
                    #Parcours maintenant l'alignement
                    seq = ligne.split()[0]
                    seqPOS, seqLS = convertAA_POS_SL(seq,NumRes,PosVect, LSVect)
                    fileres.write(seqLS+"\n")
                    fileres2.write(seqLS+"\n")
                    fileres3.write(seqPOS+"\n")
                    fileres4.write(seq+"\n")
    fileres.close()
    fileres2.close()
    fileres3.close()
    fileres4.close()
    pdb_id_out.close()



#def extract_sequences(alignment_file):
#    pdbs = {}
#    with open(alignment_file, 'r') as fastas:
#        fasta_lines = fastas.readlines()
#    for i, line in enumerate(fasta_lines):
#        if line[0] == ">":
#            if i != 0:
#                seq = align.replace("-", "")
#                pdbs[seq_name] = {'align':align,'seq':seq}
#            seq_name = re.sub("[>\n]", "", line)
#            align = ''
#        else:
#            align += line.replace('\n', '')
#        if i+1 == len(fasta_lines):
#            seq = align.replace("-", "")
#            pdbs[seq_name] = {'align':align,'seq':seq}
#    return pdbs




def extract_sequences(alignment_file):
    pdbs = {}
    with open(alignment_file, 'r') as fastas:
        fasta_lines = fastas.readlines()
    for i, line in enumerate(fasta_lines):
        if line[0] == ">":
            if i != 0:
                seq = align.replace("-", "")
                if not re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",seq_name):
                    pdbs[seq_name] = {'align':align,'seq':seq}
            seq_name = re.sub("[>\n]", "", line)
            align = ''
        else:
            align += line.replace('\n', '')
        if i+1 == len(fasta_lines):
            seq = align.replace("-", "")
            if not re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",seq_name):
                pdbs[seq_name] = {'align':align,'seq':seq}
    return pdbs





def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]



def adiff(a, b):
    print('[ERROR] PDB list and alignement file length mismatch')
    diff1 = diff(a, b)
    if diff1:
        print('Fasta file', end=' ')
        print(' '.join(diff1), end=' ')
    diff2 = diff(b, a)
    if diff2:
        print('PDB list', end=' ')
        print(' '.join(diff2))
    sys.exit()



def check_size(sizes):
    maj_size = 0
    if len(set(sizes.values())) > 1:
        for i in set(sizes.values()):
            if sizes.values().count(i) > maj_size:
                maj_size = i
        print('[ERROR] Alignement length mismatch')
        print('[BAD LEN] usual', maj_size)
        for i in sizes.keys():
            if sizes[i] != maj_size:
                print('[BAD LEN]', i, sizes[i])
        sys.exit()



def check_seq(seq1, seq2):
    print('[SEQ MISMATCH] Fasta seq')
    print(PDB.split_seq(list(seq1), 80))
    print('[SEQ MISMATCH] PDB seq  ')
    print(PDB.split_seq(list(seq2), 80))
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            print('[SEQ MISMATCH] PDB:', seq1[i], i+1, seq2[i], ':Fasta')
    sys.exit()



#def check_alignement_input(pdb_list, args):

#    pdbs = extract_sequences(args.alignFile)
#    if len(pdbs) != len(pdb_list):
#        pdb_list_name = [getname(pdb) for pdb in pdb_list]
#        adiff(pdbs.keys(), pdb_list_name)
#    pdb_db = os.path.join(args.output, "PDB")
#    sizes = {}
#    for pdb_id in pdbs:
#        print '[CHECK] aligment input', pdb_id,
#        sizes[pdb_id] = len(pdbs[pdb_id]['align'])
#        pdb, ch = pdb_id.split("_")
#        if os.path.isfile(os.path.join(pdb_db, pdb_id+".pdb")):
#            error = False
#            pdb_path = os.path.join(pdb_db, pdb_id+".pdb")
#        elif os.path.isfile(os.path.join(pdb_db, pdb+".pdb")):
#            error = False
#            pdb_path = os.path.join(pdb_db, pdb+".pdb")
#        else:
#            error = True
#        pdb_obj = PDB.PDB(pdb_path)
#        seq = ''.join(pdb_obj.prot.sequence(chain=ch))
#        if seq == pdbs[pdb_id]['seq']:
#            print 'sequence OK.',
#        else:
#            print '[ERROR] sequence mismatch', pdb_id
#            check_seq(seq, pdbs[pdb_id]['seq'])
#        print 'Alignment length', len(pdbs[pdb_id]['align'])
#    check_size(sizes)
#    copy2(args.alignFile, os.path.join(args.output, "AA_alignment."+args.method))

#    output = os.path.join(args.output, "AA_alignment."+args.method)   ###add by LR
#    output2 = os.path.join(args.output, "AA_alignment2."+args.method)   ###add by LR
#    motif_align(args, output, output2) 




def check_alignement_input(pdb_list, args):

    pdbs = extract_sequences(args.alignFile)
    if len(pdbs) != len(pdb_list):
        pdb_list_name = [getname(pdb) for pdb in pdb_list]
        adiff(pdbs.keys(), pdb_list_name)
    pdb_db = os.path.join(args.output, "PDB")
    sizes = {}
    for pdb_id in pdbs:
        print('[CHECK] aligment input', pdb_id, end=' ')
        sizes[pdb_id] = len(pdbs[pdb_id]['align'])
        pdb, ch = pdb_id.split("_")
        if os.path.isfile(os.path.join(pdb_db, pdb_id+".pdb")):
            error = False
            pdb_path = os.path.join(pdb_db, pdb_id+".pdb")
        elif os.path.isfile(os.path.join(pdb_db, pdb+".pdb")):
            error = False
            pdb_path = os.path.join(pdb_db, pdb+".pdb")
        else:
            error = True
        pdb_obj = PDB.PDB(pdb_path)
        seq = ''.join(pdb_obj.prot.sequence(chain=ch))
        if seq == pdbs[pdb_id]['seq']:
            print('sequence OK.', end=' ')
        else:
            print('[ERROR] sequence mismatch', pdb_id)
            check_seq(seq, pdbs[pdb_id]['seq'])
        print('Alignment length', len(pdbs[pdb_id]['align']))
    check_size(sizes)

    file_input = os.path.join(args.output, args.alignFile)

    cmd = "cp "+args.alignFile+" "+file_input
    os.system(cmd)

    #copy2(args.alignFile, file_input)

    sequences = SeqIO.parse(file_input, "fasta")

    filein = os.path.join(args.output, "AA_alignment.clustalw")
    count = SeqIO.write(sequences, filein , "clustal")

    output2 = os.path.join(args.output, "AA_alignment2.clustalw")

    motif_align(args, filein, output2) 

    #filein.close()
    sequences.close()




def graph(args):
    os.mkdir(os.path.join(args.output, 'graph'))
    param_file = open(os.path.join(args.output, "params.R"), "w")
    param_file.write("pathdata = '"+os.path.join(args.output, "PDB")+"'\n")
    param_file.write("pathsaconf = '"+args.output+"'\n")
    param_file.write("pathFig = '"+os.path.join(args.output, "graph")+"'\n")
    param_file.write("file1 = '"+os.path.join(args.output, PDB_SEQ_FILE)+"'")
    param_file.close()
    #run Rscript
    rscript = find_executable('Rscript')
    if rscript:
        os.system(rscript+" "+os.path.join(PATH,"src/script.R")+" "+\
                  os.path.join(args.output, "params.R"))
    else:
        sys.exit('[ERROR] Rscript not found')
    #run R
    #os.system("R CMD BATCH "+os.path.join(PATH,"src/script.R"))


##LR : I modified the three following functions


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


def fetch_xml(pdb_code, args):
    """Fetch the validation XML file from the PDB"""
    out_dir = os.path.join(args.output, "XML")
    os.makedirs(out_dir, exist_ok=True)
    dir_structure = pdb_code.lower()[1:3]
    url = f"https://files.rcsb.org/pub/pdb/validation_reports/{dir_structure}/{pdb_code.lower()}/{pdb_code.lower()}_validation.xml.gz"
    file_path_gz = os.path.join(out_dir, f"{pdb_code}_validation.xml.gz")
    file_path = file_path_gz[:-3]
    if not os.path.isfile(file_path):
        try:
            urllib.request.urlretrieve(url, file_path_gz)
            with gzip.open(file_path_gz, 'rb') as gz_file:
                with open(file_path, 'wb') as out_file:
                    shutil.copyfileobj(gz_file, out_file)
            os.remove(file_path_gz)            
        except Exception as e:
            print(f"[ERROR] Downloading XML file : {e}")
            

def get_rsrz(xml_file, Chain=None):
    """Extract RSRZ for the alpha carbons from the validation XML file"""
    rsrz_dict = {}
    try:
        # Parse the XML file
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Iterate through all <ModelledSubgroup> elements
        for subgroup in root.findall(".//ModelledSubgroup"):
            model = subgroup.attrib.get("model")
            resnum = subgroup.attrib.get("resnum")
            chain = subgroup.attrib.get("chain")
            rsrz = subgroup.attrib.get("rsrz")
            # Ensure we have the right chain
            if Chain is not None:
                if chain != Chain:
                    continue
            # Ensure it's model="1" and contains the required attributes
            if model == "1" and resnum and chain and rsrz:
                resnum = int(resnum)
                rsrz = float(rsrz)
                if chain not in rsrz_dict:
                    rsrz_dict[chain] = {}
                rsrz_dict[chain][resnum] = rsrz
        return rsrz_dict
    except FileNotFoundError:
        print(f"[ERROR] File not found: {xml_file}")
        return {}
    except ET.ParseError:
        print(f"[ERROR] Failed to parse XML file: {xml_file}")
        return {}


def extractMissingRSRZ(args, pdb_list):
    """Create a file that contains info for each residue (missing/RSRZ/b-factor)"""
    fname_indiv = os.path.join(args.output, "res_info.csv")
    with open (fname_indiv, "w") as file_indiv:
        file_indiv.write("pdb,chain,resnum,resname,missing,rsrz,bfact\n")
        for pdb_path in pdb_list:
            pdb_file = os.path.basename(pdb_path)
            pdb_code = pdb_file.split(".")[0]
            pdb_code_list = pdb_code.split("_")
            if len(pdb_code_list) > 1:
                pdb_path = re.sub(f"{pdb_code}.pdb", f"{pdb_code_list[0]}.pdb", pdb_path)
                chain_code = pdb_code_list[1]
                pdb = PDB6.PDB(pdb_path, hetSkip=1)
            else:
                pdb = PDB6.PDB(pdb_path, hetSkip=1)
                chain_code = None
            # Download and read validation XML file
            xml_file = os.path.join(args.output, "XML", f"{pdb_code_list[0]}_validation.xml")
            fetch_xml(pdb_code_list[0], args)
            # Read RSRZ
            rsrz_dict = get_rsrz(xml_file, chain_code)
            # Read b-factor
            bfact_dict = pdb.get_bfact(chain_code)
            # Read all other info
            res_dict = pdb.getAllRes(chain_code)
            # Write to file
            for chain in res_dict:
                for res in res_dict[chain]:
                    resnum = res[0]
                    resname = res[1]
                    missing = int(res[2] == "missing")
                    rsrz = rsrz_dict.get(chain, {}).get(resnum, "")
                    bfact = bfact_dict.get(chain, {}).get(resnum, "")
                    file_indiv.write(f"{pdb_code_list[0]},{chain},{resnum},{resname},{missing},{rsrz},{bfact}\n")


def gen_pml(args):
    file1 = open(os.path.join(args.output, CSV_OUT),"r")
    file_param = file1.readlines()
    file1.close()

    file2 = open(os.path.join(args.output, PDB_SEQ_FILE),"r")
    prot_ref = file2.readlines()[0].split()[0]
    file2.close()

    prot = prot_ref.split("_")[0]  ##LR: add
    ch = prot_ref.split("_")[1].lower()

    path_pdb = os.path.join(args.output, "PDB")

    file3 = open(os.path.join(args.output, POS_ALIGN_FILE),"r")
    vect_pos = file3.readlines()[1].split()[0].split('.')
    file3.close()

    posAA =  extractAAPos(file_param,  vect_pos)
    #posSL =  extractSLPos(file_param,  vect_pos)
    #posSLc =  extractSLPosCons(file_param,  vect_pos)
    #posSLwv =  extractSLPosWeakVar(file_param,  vect_pos)
    
    posSLc =  def_important_pos(file_param, vect_pos)[0]
    posSLv_noSS =  def_important_pos(file_param, vect_pos)[1]
    posSLv_withSS =  def_important_pos(file_param, vect_pos)[2]

    file_out = open(os.path.join(args.output, PML_OUT),"w")
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
        file_out.write("select mut, chain "+ ch + " and resid " + "+".join(posAA) +" \n")
        file_out.write("show sticks, mut\n")                  ##LR: modify lines by sticks
        #file_out.write("color red, mut\n")                   ##LR: comment

    #visualization of the structurally conserved aligned positions
    if len(posSLc)!= 0:    ## LR: add
        file_out.write("select consPos, chain "+ ch + " and  resid " + "+".join(posSLc) + " \n")
        file_out.write("color magenta, consPos\n")

    #visualization of the structurally variable aligned positions without SS changes
    if len(posSLv_noSS)!= 0:    ## LR: add
        file_out.write("select varPos_NoSS, chain "+ ch + " and  resid " + "+".join(posSLv_noSS) + " \n")
        file_out.write("color slate, varPos_NoSS\n")

    #visualization of the structurally variable aligned positions with SS changes
    if len(posSLv_withSS)!= 0:    ## LR: add
        file_out.write("select varPos_withSS, chain "+ ch + " and  resid " + "+".join(posSLv_withSS) + " \n")
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


def clean(args):
    rmtree(os.path.join(args.output,"HMM-SA"))
    os.remove(os.path.join(args.output, POSDAT_FILE))
    os.remove(os.path.join(args.output, AADAT_FILE))
    if os.path.isfile(os.path.join(args.output, AA_ALIGN_FILE)):
        os.remove(os.path.join(args.output, AA_ALIGN_FILE))
    if args.alignFile:
        if os.path.isfile(os.path.join(args.output, args.alignFile)):
            os.remove(os.path.join(args.output, args.alignFile))
    if os.path.isfile(os.path.join(args.output, AA_ALIGN_SWSS_TMP2)):
        os.remove(os.path.join(args.output, AA_ALIGN_SWSS_TMP2))        
    if os.path.isfile(os.path.join(args.output, FASTA_FILE.split('.')[0]+'.dnd')):
        os.remove(os.path.join(args.output, FASTA_FILE.split('.')[0]+'.dnd'))
    if os.path.isfile(os.path.join(args.output, "AA_alignment."+args.method)):
        os.remove(os.path.join(args.output, "AA_alignment."+args.method))
    if os.path.isfile(os.path.join(args.output, "list_seqLS.dat")):
        os.remove(os.path.join(args.output, "list_seqLS.dat"))
    if os.path.isfile(os.path.join(args.output, AA_ALIGN_SWSS)):
        os.remove(os.path.join(args.output, AA_ALIGN_SWSS))
    if os.path.isfile(os.path.join(args.output, FASTA_FILE)):
        os.remove(os.path.join(args.output, FASTA_FILE))
    if os.path.isfile(os.path.join(args.output, "Listpdb_extSeq.id")):
        os.remove(os.path.join(args.output, "Listpdb_extSeq.id"))        
    if os.path.isfile("seq.dnd"):
            os.remove("seq.dnd")
            
            
def moveProgR(args):
    path_Rfiles = os.path.join(args.output, 'R_files_tmp')
    print(path_Rfiles)
    os.mkdir(path_Rfiles)
    if os.path.isfile(os.path.join(args.output, AA_ALIGN_FILE2)+"_wo_header.aln"):
        os.rename(os.path.join(args.output, AA_ALIGN_FILE2)+"_wo_header.aln", 
                  os.path.join(path_Rfiles, AA_ALIGN_FILE2)+"_wo_header.aln")

    if os.path.isfile(os.path.join(args.output, SL_ALIGN_FILE)+"_wo_header.aln"):
        os.rename(os.path.join(args.output, SL_ALIGN_FILE)+"_wo_header.aln",
                  os.path.join(path_Rfiles, SL_ALIGN_FILE)+"_wo_header.aln")
        
    if os.path.isfile(os.path.join(args.output, "params.R")):
        os.rename(os.path.join(args.output, "params.R"),os.path.join(path_Rfiles, "params.R"))
        
    if os.path.isfile(os.path.join(args.output, CORRES_SWSS_ALI)):
        os.rename(os.path.join(args.output, CORRES_SWSS_ALI), os.path.join(path_Rfiles, CORRES_SWSS_ALI))
        
    if os.path.isfile(os.path.join(args.output, POS_ALIGN_FILE)):
        os.rename(os.path.join(args.output, POS_ALIGN_FILE), os.path.join(path_Rfiles, POS_ALIGN_FILE))
        
    if os.path.isfile(os.path.join(args.output, "res_info.csv")):
        os.rename(os.path.join(args.output, "res_info.csv"), os.path.join(path_Rfiles, "res_info.csv"))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
if __name__ == "__main__":

    #if args.install:
    #    install_dependencies()
    check_args(args)
    header(args)

    start_time = time.time()

    #- START
    if os.path.isdir(args.output):
        fail = True
        while fail:
            response = input(f"[Warning] {args.output} exists. Erase folder? (y/n) ")
            if response in ['y', 'n']:
                fail = False
        if response == 'n':
            sys.exit()
        else:
            rmtree(args.output)
    os.mkdir(args.output)

    print("- GET PDB")
    pdb_list = extract_pdb(args)

    print("- HMM-SA ENCODING")
    HMMSAencode(args, pdb_list)
    

    if not args.alignFile:
        print("- ALIGN")
        align_aa_seq(args)
    else:
        check_alignement_input(pdb_list, args)
    extractMissingRSRZ(args, pdb_list)


    print("- SL ALIGNMENT")
    align_sl_seq(args)

    print("- GRAPH")
    graph(args)

    print("- Generate pymol script")
    gen_pml(args)

    #clean(args)

    moveProgR(args)
    
    total_time = time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))
    print("\n-- Time:", str(total_time))
    print('')
    print('*'*80)
    print('*'*80)
    print('')
    
    
    
    #source("~/Research/Projects/src/srcSA-Conf/saconf/src/USER_script.R")

    

