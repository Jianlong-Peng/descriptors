#!/usr/bin/env python
'''
#=============================================================================
#     FileName: runMopac.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-01-03 16:41:09
#   LastChange: 2014-02-18 14:33:18
#      History:
#=============================================================================
'''
import sys
import os
import glob
import pybel
import openbabel as ob

global mopac7
if sys.platform == "win32":
    mopac7 = "mopac7.exe"
else:
    mopac7 = "mopac7"

def main(argv=sys.argv):
    if len(argv) < 2:
        print "\n  Usage: runMopac7.py mols[...]\n"
        print "  Attention:"
        print "  1. molecule title is needed for each input molecule, and"
        print "     it will be used to name the output file!"
        print "  2. if the input molecule is charged or radical, you may need"
        print "     to add DOUBLET to the keyword list passed to Mopac"
        print "  3. .mopout file will be generated in the same dir as input!"
        print "  4. by default, the program 'mopac7' can be found in PATH. Otherwise,"
        print "     you need to specify the global variable 'mopac7' to tell the script"
        print "     where to find mopac7!"
        print "  5. the input file name can contain simple shell-style wildcards!"
        print ""
        sys.exit(1)

    for candidate in argv[1:]:
        names = glob.glob(candidate)
        for name in names:
            run(name)

def genMopFile(mol):
    conv = ob.OBConversion()
    conv.SetOutFormat("mop")
    if os.sep in mol.title:
        mop_name = mol.title[mol.title.rfind(os.sep)+1:] + ".mop"
    else:
        mop_name = mol.title + ".mop"
    if not mol.OBMol.HasHydrogensAdded():
        mol.addh()
    if not mol.OBMol.Has3D():
        mol.make3D()
    key = "ESP MMOK VECTORS BONDS PI PRECISE ENPART AM1 GNORM=0.01 EF XYZ MULLIK +\nCharge=%d"%mol.charge
    conv.AddOption("k",conv.OUTOPTIONS,key)
    status = conv.WriteFile(mol.OBMol, mop_name)
    if not status:
        print "Error: failed to write mop file for "%mol.title
        return ""
    else:
        return mop_name

def runMopac(mop_name):
    lines = open(mop_name,"r").readlines()
    fi_tmp = open("FOR005","w")
    fi_tmp.writelines(lines)
    fi_tmp.close()
    basename = mop_name[:mop_name.rfind(".")]
    mopout_name = basename+".mopout"
    cmd = "%s > %s"%(mopac7,mopout_name)
    status = os.system(cmd)
    elim_list=["FOR005","FOR009","FOR010","FOR011","FOR012","SHUTDOWN"]
    for i in elim_list:
        if os.path.exists(i):
            os.unlink(i)
    if status:
        print "Error: failed to run Mopac7 for molecule \"%s\""%basename
        return False
    lines = open(mopout_name,"r").readlines()
    if "DONE" not in lines[-1]:
        print "Error: failed to run Mopac7 for molecule \"%s\""%basename
        return False
    else:
        os.unlink(mop_name)
        print "\"%s\" has been generated!"%mopout_name
        return True


def run(name):
    dirname = os.path.dirname(name)
    basename = os.path.basename(name)
    _format = basename[basename.rfind(".")+1:]
    if _format not in ("mol","mol2","sdf"):
        print "Error: invalid format %s, only \"mol, mol2, and sdf\" supported"%_format
        return
    prev_path = os.getcwd()
    if dirname != "":
        os.chdir(dirname)
    error_list = []
    print "to optimize molecule(s) from",name
    for mol in pybel.readfile(_format,basename):
        print "> molecule titled \"%s\""%mol.title
        mop_name = genMopFile(mol)
        if mop_name == "":
            error_list.append(mol.title)
            continue
        if not runMopac(mop_name):
            error_list.append(mol.title)
    if dirname != "":
        os.chdir(prev_path)
    
    if len(error_list) != 0:
        print "\nAttention: failed to run Mopac for the following molecules:"
        for item in error_list:
            print item
    else:
        print "\nAll finished successfully"

main()

