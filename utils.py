import os, sys
from os.path import basename
import glob

def get_grafts(dir):
    grafts = []
    if not os.path.isdir(dir):
        with open(dir) as tempf:
            for line in tempf:
                grafts.append(line.strip())
    else:
        if dir[-1] != '/':
            dir = dir+'/'
        for f in glob.glob(dir+"*idealized.pdb"):
            grafts.append(f)
        for f in glob.glob(dir+"*/*idealized.pdb"):
            grafts.append(f)
        for f in glob.glob(dir+"*/*/*idealized.pdb"):
            grafts.append(f)
    return grafts

def write_lines(file,lines):
    if len(lines) == 0:
        return
    with open(file, "w") as tempf:
        for l in lines:
            tempf.write(l)

def get_chain_pdbs(pdb,output_dir,prefix="",onlychain=""):
    prevresn = ""
    prevchain = ""
    pos = 0
    outpdbf = output_dir+prefix+os.path.basename(pdb).split(".pdb")[0]+".chain"
    outpdbs = []
    out = []
    with open(pdb) as tempf:
        for line in tempf: 
            if line[0:4] == "ATOM":
                chain = line[21:22]
                if chain != prevchain:
                    if len(out):
                        out.append("TER")
                        if len(onlychain) > 0:
                            if prevchain.strip() == onlychain:
                                write_lines(outpdbf+prevchain.strip()+".pdb",out)
                                outpdbs.append(outpdbf+prevchain.strip()+".pdb")
                        else:
                            write_lines(outpdbf+prevchain.strip()+".pdb",out)
                            outpdbs.append(outpdbf+prevchain.strip()+".pdb")
                    out = []
                    pos = 0
                    prevchain = chain
                resn =  line[23:26]
                if resn != prevresn:
                    pos = pos + 1
                    prevresn = resn
                out.append(line[0:22]+"%+4.4s"%(pos)+line[26:])

    if len(out):
        out.append("TER")
        if len(onlychain) > 0:
            if prevchain.strip() == onlychain:
                write_lines(outpdbf+prevchain.strip()+".pdb",out)
                outpdbs.append(outpdbf+prevchain.strip()+".pdb")
        else:
            write_lines(outpdbf+prevchain.strip()+".pdb",out)
            outpdbs.append(outpdbf+prevchain.strip()+".pdb")
    return outpdbs

