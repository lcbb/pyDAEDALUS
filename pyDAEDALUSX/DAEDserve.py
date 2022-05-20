from SimpleXMLRPCServer import SimpleXMLRPCServer
import xmlrpclib
import shutil
import sys
from math import floor
from os import listdir, path, makedirs
from Automated_Design.ply_to_input import ply_to_input
from Automated_Design.DX_cage_design import DX_cage_design
from Automated_Design.gen_PDB import pdbgen

# class daedalusRPC(object):

def calc(pName, helicalForm, helicalTurns, plyfile, seqfile):
    if not path.exists(pName):
        makedirs(pName)
    hMult=int(helicalTurns)
    hF=str(helicalForm)
    pFile=str(plyfile)
    sFile=str(seqfile)
    if hF=='Aform': 
        # staple crossover asymmetry, 11 nt/helical turn
        minEdgeLen=hMult*11
        hForm=True
        twist=1
    elif hF=='Bform': 
        # no asymmetry in staple or scaffold crossovers, 10.5 nt/helical turn
        minEdgeLen=floor(hMult*10.5)
        hForm=False
        twist=1
    elif hF=='Hybrid':
        # no asymmetry in staple or scaffold crossovers, 10.5 nt/helical turn
        minEdgeLen=hMult*11
        hForm=True
        twist=2
    elif hF=="Twisted":
        # scaffold crossover asymmetry in opposite direction, 11 nt/helical turn
        minEdgeLen=hMult*11
        hForm=True
        twist=3
    coordinates, edges, faces, edge_length_vec, file_name, \
        staple_name, singleXOs = ply_to_input(
            str(pFile), str(pName), minEdgeLen, hForm)
    if (sFile=='M13.txt'):
        scaf_seq = []
        scaf_name = []
    else:
        scaf_name = str(pName)
        fSeq=open(sFile,'r')
        scaf_seq=''
        for lines in fSeq:
            scaf_seq=scaf_seq+lines.strip()
        scaf_seq = scaf_seq.upper() # Force scaffold sequence to be uppercase
    full_file_name = DX_cage_design(  # noqa: F841
        coordinates, edges, faces, edge_length_vec, file_name,
        staple_name, singleXOs, scaf_seq, scaf_name, hForm, str(pName),
        twist, print_to_console=False)
    pdbout = pdbgen(full_file_name, hForm, str(pName))
    return "Finished!"
#
server = SimpleXMLRPCServer(("localhost", 4242))
server.register_multicall_functions()
server.register_function(calc, "calc")
server.serve_forever()
