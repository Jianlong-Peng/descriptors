#!/usr/bin/env python
'''
#=============================================================================
#     FileName: fingerprints.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-15 19:10:30
#   LastChange: 2014-01-24 14:29:10
#      History:
#=============================================================================
'''

import sys
from getopt import getopt
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import Chem

class ReadMol2:
    def __init__(self,infile):
        _format = infile[infile.rfind(".")+1:]
        if _format != "mol2":
            raise TypeError("invalid format <%s>"%_format)
        self.inf = open(infile,"r")
        self.line = self.inf.readline()
        while self.line!="" and self.line.strip()!="@<TRIPOS>MOLECULE":
            self.line = self.inf.readline()

    def __iter__(self):
        while self.line != "":
            content = self.line
            self.line = self.inf.readline()
            name = self.line.strip()
            while self.line!="" and self.line.strip()!="@<TRIPOS>MOLECULE":
                content += self.line
                self.line = self.inf.readline()
            mol = AllChem.MolFromMol2Block(content)
            if self.line == "":
                self.inf.close()
            yield (name,mol)

class ReadSdf:
    def __init__(self,infile):
        _format = infile[infile.rfind(".")+1:]
        if _format != "sdf":
            raise TypeError("invalid format <%s>"%_format)
        self.inf = open(infile,"r")

    def __iter__(self):
        line = self.inf.readline()
        while line != "":
            content = ""
            name = line.strip()
            while line!="" and line.strip()!="M  END":
                content += line
                line = self.inf.readline()
            content += line
            mol = AllChem.MolFromMolBlock(content)
            while line.strip() != "$$$$":
                line = self.inf.readline()
            line = self.inf.readline()
            if line == "":
                self.inf.close()
            yield (name,mol)

def exit_with_help():
    print "\nUsage"
    print "  fingerprints.py [options] mols[,...]\n"
    print "[options]"
    print "  --type [ecfp/fcfp/maccs]"
    print "  --out file: where to save results"
    print "  --diameter int: optional <default: 4>"
    print "    specify the diamter of the atom environments considered!"
    print "    only used when type is either 'ecfp' or 'fcfp'"
    print "  --nbits int: optional <default: 1024>"
    print "    only used when type is either 'ecfp' or 'fcfp'"
    print "  --compress: optional <default False>"
    print "    if passed, there are no spaces between each bit!\n"
    print "Attention"
    print "  1. More than one molecule file can be passed!"
    print "  2. Supported formats: mol, sdf, mol2!\n"
    sys.exit(1)

def calculate_and_write_fp(title,mol,outf,_type,n,nbits,compress):
    if _type == 1: #ecfp
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,n,nBits=nbits)
        bit_string = fp.ToBitString()
    elif _type == 2: #fcfp
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,n,nBits=nbits,useFeatures=True)
        bit_string = fp.ToBitString()
    elif _type == 3: #MACCS
        mol = Chem.AddHs(mol)
        fp = MACCSkeys.GenMACCSKeys(mol)
        bit_string = fp.ToBitString()[1:]
    else:
        pass
    outf.write(title)
    if compress:
        outf.write("\t"+bit_string)
    else:
        for bit in bit_string:
            outf.write("\t"+bit)
    outf.write("\n")

def main(argv=sys.argv):
    if len(argv) < 6:
        exit_with_help()

    options,args = getopt(argv[1:],"",["type=","out=","diameter=","nbits=","compress"])
    if len(args) == 0:
        print >>sys.stderr, "Error: at least one molecule file should be passed!"
        sys.exit(1)
    type_key = {"ecfp":1,"fcfp":2,"maccs":3}
    _type = 0
    out_file = ""
    n = 2
    nbits = 1024
    compress = False
    for opt,value in options:
        if opt == "--type":
            try:
                _type = type_key[value]
            except KeyError:
                print >>sys.stderr, "Error: type should be one of 'ecfp,fcfp,maccs', but %s given"%value
                sys.exit(1)
        elif opt == "--out":
            out_file = value
        elif opt == "--diameter":
            n = int(value)
            assert n%2==0
            n /= 2
        elif opt == "--nbits":
            nbits = int(value)
        elif opt == "--compress":
            compress = True
        else:
            print >>sys.stderr, "Error: invalid options <%s>"%opt
    assert _type!=0 and out_file!=""
    #output title
    outf = open(out_file,"w")
    outf.write("Name")
    if _type == 3:
        nbits = 166
    for i in xrange(nbits):
        outf.write("\tpro_%d"%(i+1))
    outf.write("\n")
    #read and calculate fingerprints
    num_success = 0
    fails = []
    num_fail = 0
    for name in args:
        _format = name[name.rfind(".")+1:]
        if _format == "mol":
            mol = AllChem.MolFromMolFile(name)
            if mol == None:
                fails.append(name)
                print "Error: failed to read from <%s>"%name
            else:
                num_success += 1
                calculate_and_write_fp(name,mol,outf,_type,n,nbits,compress)
        elif _format == "sdf":
            i = 0
            for title,mol in ReadSdf(name):
                i += 1
                if title == "":
                    title = "Mol_%d"%i
                if mol == None:
                    fails.append(title+" from " +name)
                    print "Error: failed to read <%s> from <%s>"%(title,name)
                else:
                    num_success += 1
                    calculate_and_write_fp(title,mol,outf,_type,n,nbits,compress)
        elif _format == "mol2":
            i = 0
            for title,mol in ReadMol2(name):
                i += 1
                if title == "":
                    title = "Mol_%d"%i
                if mol == None:
                    fails.append(title+" from "+name)
                    print "Error: failed to read <%s> from <%s>"%(title,name)
                else:
                    num_success += 1
                    calculate_and_write_fp(title,mol,outf,_type,n,nbits,compress)
        else:
            fails.append(name+": unsupport format")
            print "Error: unsupport format <%s> of file <%s>"%(_format,name)
    #END of "for name in args"
    outf.close()
    print "\n  Of molecules from %d file(s), %d(success), %d(failure)"%(len(args),num_success,len(fails))
    for fail_info in fails:
        print "   ",fail_info

main()

