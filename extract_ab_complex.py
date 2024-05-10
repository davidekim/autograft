import os, sys
import pyrosetta
from pyrosetta import *

contact_dist = 30
remove_ab_dimer_chain = True

complex_pdb = ""
if len(sys.argv) == 3:
    abchain = sys.argv[1]
    complex_pdb = sys.argv[2]
elif len(sys.argv) == 4 and sys.argv[3].isnumeric():
    abchain = sys.argv[1]
    complex_pdb = sys.argv[2]
    contact_dist = float(sys.argv[3])
else:
    print(f'USAGE: {sys.argv[0]} <ab chain> <pdb> [contact distance default: {contact_dist}]')
    print()
    sys.exit()


name = os.path.basename(complex_pdb).split(".pdb")[0]

nf = open(name+f'_tmp.pdb',"w")
with open(complex_pdb) as tmpf:
    for line in tmpf:
        if line[0:6] == "ANISOU":
            continue
        nf.write(line)
nf.close()

pyrosetta.init('-mute all -detect_disulf false -out::file::renumber_pdb -in::ignore_waters -corrections::beta_nov16')

print(f'Using contact dist: {contact_dist}')
p = pyrosetta.pose_from_file(name+f'_tmp.pdb')
print(f'input number of chains: {p.num_chains()}')
abcid = 0
for cid in range(1, p.num_chains()+1):
    skip = False
    for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
        if not j.is_protein():
           skip = True
           break
    if skip:
        continue

    if pyrosetta.rosetta.core.pose.get_chain_from_chain_id( cid, p) == abchain:
        abcid = cid

tmpchainstr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz"
for i in range(1,20):
    tmpchainstr = tmpchainstr + tmpchainstr

abrosettachain = tmpchainstr[abcid-1]


p.dump_pdb(name+f'_rosetta_tmp.pdb')

p = pyrosetta.pose_from_file(name+f'_rosetta_tmp.pdb')

chain_isprotein = {}
for cid in range(1, p.num_chains()+1):
    skip = False
    for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
        if not j.is_protein():
           chain_isprotein[cid] = False
           skip = True
           break
    if skip:
        continue
    chain_isprotein[cid] = True
    if pyrosetta.rosetta.core.pose.get_chain_from_chain_id( cid, p) == abrosettachain:
        abcid = cid
contacts = {}


chainseq = p.chain_sequence(abcid)

for cid in range(1, p.num_chains()+1):
    if cid == abcid:
        continue
    contacts[cid] = []

resnums_to_keep = []
chains_to_keep = [abrosettachain]
chainids_to_keep = [abcid]
for i in pyrosetta.rosetta.core.pose.get_chain_residues(p, abcid):
    resnums_to_keep.append(p.pdb_info().number(i.seqpos()))
    for cid in range(1, p.num_chains()+1):
        if cid == abcid:
            continue
        atomid = 1
        if chain_isprotein[cid]:
            atomid = 2
        for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
            dist = (i.atom(2).xyz() - j.atom(atomid).xyz()).norm()
            if dist < contact_dist:
                contacts[cid].append(j.seqpos())

contacts = dict(sorted(contacts.items(), key=lambda item: len(item[1]), reverse=True))

cnt = 0
for cid in contacts:
    cnt = cnt + 1
    # assume the chain with most contacts is the dimer chain of the antibody
    if remove_ab_dimer_chain and (len(contacts[cid]) == 0 or cnt == 1):
        continue
    print(f'{pyrosetta.rosetta.core.pose.get_chain_from_chain_id(cid, p)} contacts: {len(contacts[cid])}')
    for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
        resnums_to_keep.append(p.pdb_info().number(j.seqpos()))
    chains_to_keep.append(pyrosetta.rosetta.core.pose.get_chain_from_chain_id(cid, p))
    chainids_to_keep.append(cid)

nf = open(name+f'_{abrosettachain}_tmp.pdb',"w")
prevchain = ""
with open(name+f'_rosetta_tmp.pdb') as tmpf:
    for line in tmpf:
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            resn = int(line[22:26].strip())
            if resn in resnums_to_keep:
                chain = line[21:22].strip()
                if chain in chains_to_keep:
                    nf.write(line)
                    if len(prevchain) and prevchain !=  chain:
                        nf.write("TER\n")
                        prevchain = chain
        else:
            nf.write(line)
    nf.write("TER\n")
nf.close()

p = pyrosetta.pose_from_file(name+f'_{abrosettachain}_tmp.pdb')
newabrosettachain = ""
for cid in range(1, p.num_chains()+1):
    if chainseq == p.chain_sequence(cid):
        newabrosettachain = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(cid, p)
        break
print(f'output number of chains: {p.num_chains()}')
print(f'New chain id: {newabrosettachain}')
os.rename(name+f'_{abrosettachain}_tmp.pdb', name+f'_{abchain}_to_{newabrosettachain}_rosetta.pdb')
os.remove(name+f'_tmp.pdb')
os.remove(name+f'_rosetta_tmp.pdb')





