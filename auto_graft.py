import os, sys
import numpy as np
from utils import *
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import glob
import re
import json
import time
import tempfile

#
# David Kim 2022
#

temp_runid = next(tempfile._get_candidate_names())
verbose = True
slurm_queue = "cpu"
slurmtimelimit = "03:00:00"
CONF = {}

checkpoint = True
checkpoints = {}
def add_checkpoint(checkpointfile,checkpointid):
    with open(checkpointfile, "a+") as f:
        f.write(checkpointid+"\n")

aa_three_to_one = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

def init_conf():
    global CONF

    # filter criteria
    if 'graft_chain' not in CONF:
        CONF['graft_chain'] = "A" # the chain of the interacting loop to graft

    if 'MAX_RMSD' not in CONF:
        CONF['MAX_RMSD'] = 1.5 # Maximum rmsd for strand pair superposition

    if 'MAX_SCAF_REPLACEMENT_LEN' not in CONF:
        CONF['MAX_SCAF_REPLACEMENT_LEN'] = 10 # maximum number of scaffold residues that can be replaced by a graft 

    if 'MAX_STRAND_PAIRS' not in CONF:
        CONF['MAX_STRAND_PAIRS'] = 3 # maximum strand pairs (differing by a residue shift) to superimpose, starting from 
                                     # closest pair to target

    if 'MAX_STRAND_PAIR_CONTACT_DISTANCE' not in CONF:
        CONF['MAX_STRAND_PAIR_CONTACT_DISTANCE'] = 15 # distance threshold to target for finding potential strand pairs

    if 'MAX_GRAFT_POINT_CONTACTS' not in CONF:
        CONF['MAX_GRAFT_POINT_CONTACTS'] = 5 # maximum contacts (<MAX_CONTACT_DISTANCE) strand pair residues have with target 
                                 # for choosing strand pairs, higher values get pairs closer to the target

    if 'MAX_CONTACT_DISTANCE' not in CONF:
        CONF['MAX_CONTACT_DISTANCE'] = 8 # higher values get pairs further from the target

    if 'MIN_GRAFT_CONTACTS' not in CONF:
        CONF['MIN_GRAFT_CONTACTS'] = 10

    if 'SHIFT_STRAND_PAIRS_BY' not in CONF:
        CONF['SHIFT_STRAND_PAIRS_BY'] = 0 # shift the strand pairs closer to the end by this number of residues (to 
                                          # get grafts closer or further from the target)

    if 'MAX_GRAFT_LEN' not in CONF:
        CONF['MAX_GRAFT_LEN'] = 20 


    if 'GRAFT_STRAND_LEN' not in CONF:
        CONF['GRAFT_STRAND_LEN'] = 3

init_conf()

def report_failure(message):
    print('\033[91m'+f'FAILED {message}'+'\x1b[0m')

# get chain id of graft chain
def get_graft_chain_id(pose):
    for cid in range(1, pose.num_chains()+1):
        if pyrosetta.rosetta.core.pose.get_chain_from_chain_id(cid, pose) == CONF['graft_chain']:
            return cid

# get all candidate scaffold strand pairs
def get_scaffold_strand_pairs(pose):
    pairs = []
    pairset = pyrosetta.rosetta.core.scoring.dssp.StrandPairingSet(pose)
    for i in range(1,pairset.size()+1):
        pairings = pairset.strand_pairing(i)
        if not pairings.antiparallel():
            continue  
        pairinglist = pyrosetta.rosetta.utility.vector1_core_scoring_dssp_Pairing()
        pairings.get_beta_pairs(pairinglist)
        for j in range(len(pairinglist), 0, -1):
            pos1 = pairinglist[j].Pos1()
            pos2 = pairinglist[j].Pos2()
            if pos2-pos1+1 > CONF['MAX_SCAF_REPLACEMENT_LEN']:
                continue
            pairs.append([pos1,pos2, pairinglist[j].Pleating()])
    return pairs

# get strand pairs and pleating in the interacting hairpin loop for grafting
def get_strand_pairs_for_grafting(pose, chain_id):
    pairset = pyrosetta.rosetta.core.scoring.dssp.StrandPairingSet(pose)
    pairsetindex_to_graft = 0
    strandpair_contacts = {}
    contacts = {}
    distances = {}
    for j in pyrosetta.rosetta.core.pose.get_chain_residues(pose, chain_id):
        contacts[j.seqpos()] = []
    for cid in range(1, pose.num_chains()+1):
        for i in pyrosetta.rosetta.core.pose.get_chain_residues(pose, cid):
            if cid == chain_id:
                continue
            atomidi = 2
            if not i.is_protein():
                atomidi = 1
            for j in pyrosetta.rosetta.core.pose.get_chain_residues(pose, chain_id):
                atomidj = 2
                if not j.is_protein():
                    atomidj = 1
                dist = (i.atom(atomidi).xyz() - j.atom(atomidj).xyz()).norm()
                if dist <= CONF['MAX_STRAND_PAIR_CONTACT_DISTANCE']:
                    contacts[j.seqpos()].append(i.seqpos())
                    distances[(i.seqpos(),j.seqpos())] = dist

    for i in range(1,pairset.size()+1):
        pairings = pairset.strand_pairing(i)
        if not pairings.antiparallel():
            continue 
        pairinglist = pyrosetta.rosetta.utility.vector1_core_scoring_dssp_Pairing()
        pairings.get_beta_pairs(pairinglist)
        if pose.residue(pairinglist[1].Pos1()).chain() == chain_id and pose.residue(pairinglist[1].Pos2()).chain() == chain_id:
            totalcontacts = 0
            for j in range(1,len(pairinglist)+1):
                pos1 = pairinglist[j].Pos1()
                pos2 = pairinglist[j].Pos2()
                if pos2-pos1+1 > CONF['MAX_GRAFT_LEN']:
                    continue                
                if pose.residue(pos1).chain() == chain_id and pose.residue(pos2).chain() == chain_id:
                    totalcontacts = totalcontacts+len(contacts[pos1])+len(contacts[pos2])
            if totalcontacts > 0:
                strandpair_contacts[i] = totalcontacts

    # get strand pairs with loop graft contacts
    strandpairs = []
    strandpairs_graftcontacts = {}
    for pairsetindex_to_graft in strandpair_contacts:
        if strandpair_contacts[pairsetindex_to_graft] == 0:
            continue
        pairings = pairset.strand_pairing(pairsetindex_to_graft)
        pairinglist = pyrosetta.rosetta.utility.vector1_core_scoring_dssp_Pairing()
        pairings.get_beta_pairs(pairinglist)
        for j in range(len(pairinglist),0,-1):
            pos1 = pairinglist[j].Pos1()
            pos2 = pairinglist[j].Pos2()
            if pos2-pos1+1 > CONF['MAX_GRAFT_LEN']:
                continue
            # get total graft contacts to target
            totalgraftcontacts = 0
            for k in range(pos1,pos2+1):
                for targetpos in contacts[k]:
                    if distances[(targetpos,k)] < CONF['MAX_CONTACT_DISTANCE']:
                        totalgraftcontacts = totalgraftcontacts + 1
            # get total graft point (residue pair) contacts to target
            totalgraftpointcontacts = 0
            for targetpos in contacts[pos1]:
                if distances[(targetpos,pos1)] < CONF['MAX_CONTACT_DISTANCE']:
                    totalgraftpointcontacts = totalgraftpointcontacts + 1
            for targetpos in contacts[pos2]:
                if distances[(targetpos,pos2)] < CONF['MAX_CONTACT_DISTANCE']:
                    totalgraftpointcontacts = totalgraftpointcontacts + 1
            if totalgraftcontacts >= CONF['MIN_GRAFT_CONTACTS'] and totalgraftpointcontacts <= CONF['MAX_GRAFT_POINT_CONTACTS']:
                # strands in strand pair to superimpose
                if verbose:
                    print(f'Considering strand pair: {pos1} {pos2} pleating: {pairinglist[j].Pleating()} graftpointcontacts: {totalgraftpointcontacts} graftcontacts: {totalgraftcontacts} pairings: {pairings}')
                strandpairs_graftcontacts[len(strandpairs)] = totalgraftcontacts
                strandpairs.append([pos1, pos2, pairinglist[j].Pleating()])

    # sort by number of graft contacts
    strandpairs_graftcontacts = dict(sorted(strandpairs_graftcontacts.items(), key=lambda item: item[1], reverse=True))

    # get beta strand pairs for superpositions
    strand_pairs_for_grafting = []
    for i in strandpairs_graftcontacts:
        if len(strand_pairs_for_grafting) >= CONF['MAX_STRAND_PAIRS']:
            break
        if verbose:
            print(f'Using: {i} graftcontacts: {strandpairs_graftcontacts[i]} pos1,pos2,pleating: {strandpairs[i]}')
        strand_pairs_for_grafting.append(strandpairs[i])
    return strand_pairs_for_grafting

# get resnums of strands to superimpose from strand pair
graft_strand_half_len = int(round((CONF['GRAFT_STRAND_LEN']-1)/2))
#if verbose:
#    print(f' graft half len: {graft_strand_half_len}')
def get_superposition_resnums_from_pair(pair):
    superposition_resnums = []
    starts = [pair[0]-graft_strand_half_len, pair[1]+graft_strand_half_len-CONF['GRAFT_STRAND_LEN']+1]
    for i in starts:
        for j in range(i, i+CONF['GRAFT_STRAND_LEN']):
            superposition_resnums.append(j)
    return superposition_resnums


# idealize
def idealize_graft(pose, scaffold, chain_id, graft_seq, idealized_outputfile):
    if verbose:
        print(f'Creating: {idealized_outputfile}')
    try:
        # idealize
        idealize = rosetta.protocols.idealize.IdealizeMover()
        to_idealize = rosetta.utility.vector1_unsigned_long()
        scorefxn=create_score_function('empty')
        scorefxn.set_weight(rosetta.core.scoring.cart_bonded, 1.0)
        scorefxn.score(pose)
        emap = pose.energies()
        bfactors = []
        mmap = MoveMap()
        mmap.set_bb(False)
        mmap.set_chi(False)
        mmap.set_jump(False)
        for r in range(1,pose.total_residue()+1):
            bfactors.append(pose.pdb_info().bfactor(r,1))
            if pose.residue(r).chain() != chain_id:
                continue
            cart = emap.residue_total_energy(r)
            if cart > 20:
                to_idealize.append(r)
                mmap.set_bb(r, True)
                if verbose:
                    print(f'idealize {r} {cart}')
        if len(to_idealize) > 0:
            idealize.set_pos_list(to_idealize)
            idealize.apply(pose)

            # cart-minimize
            scorefxn_min=create_score_function('ref2015_cart')
            min_mover = rosetta.protocols.minimization_packing.MinMover(mmap, scorefxn_min, 'lbfgs_armijo_nonmonotone', 0.00001, True)
            min_mover.max_iter(1000)
            min_mover.cartesian(True)
            if verbose:
                print("minimize...")
            min_mover.apply(pose)

        # add HOTSPOT remarks (grafted interacting loop)
        for resn in range(graft_point_resnums[0],graft_point_resnums[1]+1):
            pose.pdb_info().add_reslabel(resn,"HOTSPOT")

        pose.dump_pdb(idealized_outputfile+".tmp")
        # can't figure out how to set bfactors in pose, pdb_info().bfactor is not working so just update the pdb file
        # lets abuse REMARK
        nf = open(idealized_outputfile,"w")
        nf.write(f'REMARK  99'+"\n")
        nf.write(f'REMARK  99 scaffold: {scaffold}'+"\n")

        # do not modify this since post processing scripts depend on this format and info
        nf.write(f'REMARK  99 CHAIN: {CONF["graft_chain"]} GRAFT: {graft_point_resnums[0]}-{graft_point_resnums[1]} GRAFT_SEQ: {graft_seq}'+"\n")

        nf.write(f'REMARK  99 {CONF}'+"\n")
        nf.write(f'REMARK  99 rmsd: {rmsd:.2f}'+"\n")
        nf.write(f'REMARK  99 {scaffold_strand_pair} {graft_points}'+"\n")
        with open(idealized_outputfile+".tmp") as tmpf:
            for line in tmpf:
                if line[0:4] == "ATOM":
                    resn = int(line[22:26].strip())

                    # update to original chains
                    origchain = CONF['graft_chain']
                    if resn in chainmap:
                        origchain = chainmap[resn]
                    line = line[0:20]+"%+2.2s"%(origchain)+line[22:]

                    # update b-factor for visualization
                    if bfactors[resn-1] > 0:
                        line = line[0:60]+"100.00 "+line[67:]
                nf.write(line)
        nf.close()
        os.remove(idealized_outputfile+".tmp")
        return True
    except:
        report_failure("Idealize failed!")
        return False

def pdb_resnum(pose, rosettares):
    return pose.pdb_info().number(rosettares)

def pdb_name( pdb ):
    return os.path.basename(pdb).split(".pdb")[0]

### MAIN ###

slurm_scaff_per_job = 0
slurm_scaff_set = 0

if len(sys.argv) >= 2:
    f = open(sys.argv[1])
    CONF = json.load(f)
    init_conf()
else:
    print(f'USAGE: {sys.argv[0]} <json_input_file> [scaffolds per sbatch {slurm_queue} queue array job (min 100)]')
    print()
    print("json params:")
    print(CONF)
    print()
    sys.exit()

if len(sys.argv) == 3 and sys.argv[2].isnumeric():
    slurm_scaff_per_job = int(sys.argv[2])
elif len(sys.argv) == 4 and sys.argv[2].isnumeric() and sys.argv[3].isnumeric():
    slurm_scaff_per_job = int(sys.argv[2])
    slurm_scaff_set = int(sys.argv[3])

## required
if 'target_complex_pdb' not in CONF or not os.path.exists(CONF['target_complex_pdb']):
    print("Error! missing 'target_complex_pdb' in json input\n")
    sys.exit()
if 'scaffolds_paths' not in CONF:
    print("Error! missing 'scaffolds_paths' in json input\n")
    sys.exit()
if 'output_path' not in CONF:
    print("Error! missing 'output_path' in json input\n")
    sys.exit()

if CONF['output_path'][-1] != "/":
    CONF['output_path'] = CONF['output_path']+"/"
if not os.path.isdir(CONF['output_path']):
    os.mkdir(CONF['output_path'])
print(f'Config: {CONF}')
print(sys.argv)

# make json input filename the job id
jobid = os.path.basename(sys.argv[1]).split(".json")[0]
scaffolds_file = jobid+"_scaffolds.txt"
scaffolds = []
if os.path.exists(scaffolds_file):
    scaffoldsdict = {}
    with open(scaffolds_file) as tempf:
        for line in tempf:
            scaff = line.strip()
            if scaff not in scaffoldsdict:
                scaffolds.append(scaff)
                scaffoldsdict[scaff] = True
else:
    # create scaffolds list file
    sfl = open(scaffolds_file, "w")
    scaffoldsdict = {}
    for scaffolds_path in CONF['scaffolds_paths']:
        if not os.path.isdir(scaffolds_path):
            with open(scaffolds_path) as tempf:
                for line in tempf:
                    scaff = line.strip()
                    if scaff not in scaffoldsdict:
                        scaffolds.append(scaff)
                        sfl.write(scaff+"\n")
                        scaffoldsdict[scaff] = True
        else:
            for scaff in glob.glob(scaffolds_path+"/*.pdb"):
                if scaff not in scaffoldsdict:                
                    scaffolds.append(scaff)
                    sfl.write(scaff+"\n")
                    scaffoldsdict[scaff] = True
            for scaff in glob.glob(scaffolds_path+"/*/*.pdb"):
                if scaff not in scaffoldsdict:
                    scaffolds.append(scaff)
                    sfl.write(scaff+"\n")
                    scaffoldsdict[scaff] = True
            for scaff in glob.glob(scaffolds_path+"/*/*/*.pdb"):
                if scaff not in scaffoldsdict:
                    scaffolds.append(scaff)
                    sfl.write(scaff+"\n")
                    scaffoldsdict[scaff] = True
    sfl.close()
if verbose:
    print(f'Scaffolds to check: {len(scaffolds)}')
if slurm_scaff_per_job > 0:
    if slurm_scaff_per_job < 100:
        print("Error! scaffolds per sbatch minimum is 100")
        sys.exit()
    if slurm_scaff_per_job > len(scaffolds):
        print(f'Error! scaffolds per sbatch {slurm_scaff_per_job} is greater than the number of scaffolds {len(scaffolds)}')
        sys.exit()

if slurm_scaff_per_job > 0 and slurm_scaff_set > 0:
    # run a subset of the scaffolds (for a slurm array job)
    xindex = (slurm_scaff_per_job * slurm_scaff_set) - slurm_scaff_per_job
    yindex = slurm_scaff_per_job * slurm_scaff_set
    if yindex > len(scaffolds):
        scaffolds = scaffolds[xindex:]
    else:
        scaffolds = scaffolds[xindex:yindex]
elif slurm_scaff_per_job > 0:
    # create slurm array tasks and submit
    tasksf = open(jobid+"_tasks.txt", "w")
    cnt = 1
    while slurm_scaff_per_job*cnt < len(scaffolds)+slurm_scaff_per_job:
        tasksf.write(f'python {sys.argv[0]} {sys.argv[1]} {slurm_scaff_per_job} {cnt}'+"\n")
        if slurm_scaff_per_job*cnt >= len(scaffolds):
            break 
        cnt = cnt + 1
    tasksf.close()
    # create slurm sbatch file
    sbf = open(jobid+"_sbatch.sh", "w")
    sbf.write("#!/bin/bash\n#SBATCH -p "+slurm_queue+"\n#SBATCH --mem=8g\n")
    sbf.write(f'#SBATCH -t {slurmtimelimit}'+"\n\n")
    sbf.write("source activate pyrosetta\n")
    sbf.write('CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" '+jobid+'_tasks.txt)'+"\n")
    sbf.write('echo "${CMD}" | bash'+"\n")
    sbf.close()
    # submit
    os.system('sbatch -a 1-$(cat '+jobid+'_tasks.txt|wc -l) '+jobid+'_sbatch.sh')
    sys.exit(0)

pyrosetta.init('-mute all -detect_disulf false -out::file::renumber_pdb -in::ignore_unrecognized_res -in::ignore_waters -corrections::beta_nov16')

complex_name = pdb_name(CONF['target_complex_pdb'])
target_complex_pose = pyrosetta.pose_from_file(CONF['target_complex_pdb'])
chainid = get_graft_chain_id(target_complex_pose)
strandpairs_for_grafting = get_strand_pairs_for_grafting(target_complex_pose, chainid)
if verbose:
    print(f'strandpairs_for_grafting: {strandpairs_for_grafting}')

totalgrafts = 0
totalidealgrafts = 0

# iterate strand pairs
for graft_points in strandpairs_for_grafting:
    if CONF['SHIFT_STRAND_PAIRS_BY'] != 0:
        graft_points[0] = graft_points[0] + CONF['SHIFT_STRAND_PAIRS_BY']
        graft_points[1] = graft_points[1] - CONF['SHIFT_STRAND_PAIRS_BY']
    graft_points_str = f'{graft_points[0]}_{graft_points[1]}'
    graft_name = f'{complex_name}_{graft_points_str}'

    interacting_hairpin_loop_superposition_resnums = get_superposition_resnums_from_pair(graft_points)

    graft_points_pdb = [ pdb_resnum(target_complex_pose, graft_points[0]), pdb_resnum(target_complex_pose, graft_points[1]) ]

    # iterate scaffolds
    for s in scaffolds:
        scaffold_name = pdb_name(s)
        graftname = f'{complex_name}_{graft_points_str}'+'_'+os.path.basename(os.path.dirname(os.path.abspath(s)))
        output_dir = f'{CONF["output_path"]}grafts_{graftname}_rmsd{CONF["MAX_RMSD"]}/'
        checkpoint_file = f'{output_dir}checkpoint_{slurm_scaff_set}'
        if checkpoint:
            if os.path.exists(checkpoint_file):
                with open(checkpoint_file) as cpf:
                    for l in cpf:
                        checkpoints[l.strip()] = True
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        single_chain_scaffolds = get_chain_pdbs(s,output_dir,temp_runid+'_')
        for scaffold in single_chain_scaffolds:
           single_chain_scaffold_name = pdb_name(scaffold)
           scaffold_pose = pyrosetta.pose_from_file(scaffold)
           for scaffold_strand_pair in get_scaffold_strand_pairs(scaffold_pose):
                checkpointrunid = single_chain_scaffold_name[len(temp_runid+'_'):]+f'_{graft_points[0]}_{graft_points[1]}_{scaffold_strand_pair[0]}_{scaffold_strand_pair[1]}'
                if checkpoint:
                    if checkpointrunid in checkpoints:
                        if verbose:
                            print(f'checkpoint {checkpointrunid} exists.')
                        continue

                # make sure strand pleating at graft point is the same between scaffold and interacting loop hairpin
                if graft_points[2] != scaffold_strand_pair[2]:
                    if verbose:
                        report_failure(f'{scaffold_name} {scaffold_strand_pair} pleating mismatch to {graft_points}')
                    if checkpoint:
                        add_checkpoint(checkpoint_file,checkpointrunid)
                    continue

                # Superimpose strand pairs and get rmsd
                scaffold_superposition_resnums = get_superposition_resnums_from_pair(scaffold_strand_pair)
                ca_map = map_core_id_AtomID_core_id_AtomID()
                for i in range(0, len(interacting_hairpin_loop_superposition_resnums)):
                    ca_map[AtomID(scaffold_pose.residue(scaffold_superposition_resnums[i]).atom_index("CA"), scaffold_superposition_resnums[i])] = \
                            AtomID(target_complex_pose.residue(interacting_hairpin_loop_superposition_resnums[i]).atom_index("CA"), interacting_hairpin_loop_superposition_resnums[i])
                rmsd = pyrosetta.rosetta.core.scoring.superimpose_pose(scaffold_pose, target_complex_pose, ca_map)
                if rmsd > CONF['MAX_RMSD']:
                    if verbose:
                        report_failure(f'{scaffold_name} {scaffold_strand_pair} rmsd too high {rmsd} > max allowed ({CONF["MAX_RMSD"]}) to {graft_points}')
                    if checkpoint:
                        add_checkpoint(checkpoint_file,checkpointrunid)
                    continue

                # Create grafted pdb (to do: use pose manipulations instead)
                #   The superposition strand pair residues of the interacting hairpin loop are not included in the graft,
                #   the superposition strand pair residues of the scaffold are kept.
                #
                #   The start and end residues of the grafted interacting loop are labeled with bfactor = 100.0.
                #   The grafted residues are pdb_info labeled as HOTSPOT residues.

                tmp_scaffold = f'{single_chain_scaffold_name}_{graft_points[0]}_{graft_points[1]}_{scaffold_strand_pair[0]}_{scaffold_strand_pair[1]}.pdb' 
                scaffold_pose.dump_pdb(tmp_scaffold)
                tmp_scaffold_lines = []
                with open(tmp_scaffold) as tmps:
                    for l in tmps:
                        tmp_scaffold_lines.append(l) 
                os.remove(tmp_scaffold)

                grafted_name_prefix = scaffold_name+'_'+f'{pdb_resnum(scaffold_pose, scaffold_strand_pair[0])}_{pdb_resnum(scaffold_pose, scaffold_strand_pair[1])}'+'_'+ \
                    graft_name+'_'+f'{pdb_resnum(target_complex_pose,graft_points[0])}_{pdb_resnum(target_complex_pose,graft_points[1])}'+"_grafted"
                grafted_name = grafted_name_prefix+".pdb"

                graftedf = open(output_dir+grafted_name, "w")
                graftedf.write(f'REMARK  99'+"\n")
                graftedf.write(f'REMARK  99 {CONF}'+"\n")
                graftedf.write(f'REMARK  99 rmsd: {rmsd:.2f}'+"\n")
                graft_seq = ""
                n_added = False
                c_added = False
                with open(CONF['target_complex_pdb']) as tempf:
                    for line in tempf:
                        if line[0:6] == "HETATM" and CONF['graft_chain'] == line[21:22].strip():
                            continue # skip HETATM in graft chain
                        elif line[0:4] == "ATOM":
                            line = line[0:60]+"  0.00 "+line[67:] # zero b-factor
                            if line[21:22].strip() == CONF['graft_chain']: # loop chain
                                resn = int(line[22:26].strip())
                                if resn > graft_points_pdb[0] and resn < graft_points_pdb[1]:
                                    # N-term of scaffold before graft
                                    if resn == graft_points_pdb[0]+1:
                                        # update b-factor for visualization
                                        line = line[0:60]+"100.00 "+line[67:]
                                        if not n_added:
                                            for sline in tmp_scaffold_lines:
                                                if sline[0:4] == "ATOM":
                                                    sresn = int(sline[22:26].strip())
                                                    # add N term up to graft point of scaffold original coords (including scaffold graft point)
                                                    if sresn <= scaffold_strand_pair[0]:
                                                        sline = sline[0:60]+"  0.00 "+sline[67:] # zero b-factor
                                                        sline = sline[0:20]+"%+2.2s"%(CONF['graft_chain'])+sline[22:]
                                                        graftedf.write(sline)
                                            n_added = True
                                    if resn == graft_points_pdb[1]-1:
                                        # update b-factor for visualization
                                        line = line[0:60]+"100.00 "+line[67:]
                                    # add loop graft!
                                    line = line[0:22]+"%+4.4s"%(scaffold_strand_pair[0]+resn-graft_points_pdb[0])+line[26:]
                                    # save loop graft sequence
                                    if line[13:15].strip() == "CA":
                                        graft_seq = graft_seq + aa_three_to_one[line[17:20]]
                                    graftedf.write(line)
                                
                                # C-term of scaffold after graft
                                elif resn == graft_points_pdb[1]:
                                    if not c_added:
                                        for sline in tmp_scaffold_lines:
                                            if sline[0:4] == "ATOM":
                                                sresn = int(sline[22:26].strip())
                                                # add C term starting at graft point of scaffold original coords (including scaffold graft point)
                                                if sresn >= scaffold_strand_pair[1]:
                                                    sline = sline[0:60]+"  0.00 "+sline[67:] # zero b-factor
                                                    sline = sline[0:20]+"%+2.2s"%(CONF['graft_chain'])+sline[22:]
                                                    sline = sline[0:22]+"%+4.4s"%(graft_points_pdb[1]-graft_points_pdb[0]+sresn-(scaffold_strand_pair[1]-scaffold_strand_pair[0]+1)+1)+sline[26:]
                                                    graftedf.write(sline)
                                        c_added = True
                                        graftedf.write("TER\n")
                            else:
                                 graftedf.write(line)
                        else:
                            graftedf.write(line)
                graftedf.close()

                p = pyrosetta.pose_from_file(output_dir+grafted_name)

                os.remove(output_dir+grafted_name)

                # get start and stop resnums of grafted residues
                graft_point_resnums = []
                for r in range(1,p.total_residue()+1):
                    if p.pdb_info().bfactor(r,1) > 0:
                        graft_point_resnums.append(r)
                if len(graft_point_resnums) != 2:
                    report_failure(f'cannot get graft point resnums from {grafted_name}')
                    if checkpoint:
                        add_checkpoint(checkpoint_file,checkpointrunid)
                    totalgrafts = totalgrafts + 1
                    continue

                idealizedfile = output_dir+grafted_name_prefix+f'_{CONF["graft_chain"]}_{graft_point_resnums[0]}_{graft_point_resnums[1]}_idealized.pdb'
                if os.path.exists(idealizedfile):
                    if verbose:
                        print(f'Continue since {idealizedfile} already exists.')
                    totalidealgrafts = totalidealgrafts + 1
                    totalgrafts = totalgrafts + 1
                    continue

                # save chain mapping and check for backbone clash
                chainmap = {}
                clashing = False
                for cid in range(1, p.num_chains()+1):
                    for i in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
                        chainmap[i.seqpos()] = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(cid, p)
                    if clashing:
                        break
                    # check for intra chain clashing
                    if cid == chainid:
                        for i in pyrosetta.rosetta.core.pose.get_chain_residues(p, chainid):
                            if clashing:
                                break
                            iseqpos = i.seqpos()+2
                            for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, chainid):
                                if iseqpos < j.seqpos():
                                    if (i.atom("CA").xyz() - j.atom("CA").xyz()).norm() < 3.4: # 3.5:
                                        clashing = True
                                    elif (i.atom("O").xyz() - j.atom("O").xyz()).norm() < 2.6: #2.7:
                                        clashing = True
                                    elif (i.atom("C").xyz() - j.atom("C").xyz()).norm() < 3.6: #3.7:
                                        clashing = True
                                    elif (i.atom("N").xyz() - j.atom("N").xyz()).norm() < 3.3: #3.4:
                                        clashing = True
                                    if clashing:
                                        break
                        continue
                    # check for inter chain clashing
                    for i in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
                        if clashing:
                            break
                        if i.is_protein():
                            for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, chainid):
                                if (i.atom("CA").xyz() - j.atom("CA").xyz()).norm() < 3.4: # 3.5:
                                    clashing = True
                                elif (i.atom("O").xyz() - j.atom("O").xyz()).norm() < 2.6: #2.7:
                                    clashing = True
                                elif (i.atom("C").xyz() - j.atom("C").xyz()).norm() < 3.6: #3.7:
                                    clashing = True
                                elif (i.atom("N").xyz() - j.atom("N").xyz()).norm() < 3.3: #3.4:
                                    clashing = True
                                if clashing:
                                    break
                        else:
                            for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, chainid):
                                if clashing:
                                    break
                                for k in range(1, i.natoms()+1):
                                    if (i.atom(k).xyz() - j.atom("CA").xyz()).norm() < 3.4: # 3.5:  what distances should we use?
                                        clashing = True
                                    elif (i.atom(k).xyz() - j.atom("O").xyz()).norm() < 2.6: #2.7:
                                        clashing = True
                                    elif (i.atom(k).xyz() - j.atom("C").xyz()).norm() < 3.6: #3.7:
                                        clashing = True
                                    elif (i.atom(k).xyz() - j.atom("N").xyz()).norm() < 3.3: #3.4:
                                        clashing = True
                                    if clashing:
                                        break
                if clashing:
                    if verbose:
                        report_failure(f'{scaffold_name} {scaffold_strand_pair} clashing in graft to {graft_points}')
                    if checkpoint:
                        add_checkpoint(checkpoint_file,checkpointrunid)
                elif idealize_graft(p, s, chainid, graft_seq, idealizedfile):
                    totalidealgrafts = totalidealgrafts + 1
                    if verbose:
                        pair1_O_dist = (target_complex_pose.residue(graft_points[0]).atom("O").xyz() - scaffold_pose.residue(scaffold_strand_pair[0]).atom("O").xyz()).norm()              
                        pair2_O_dist = (target_complex_pose.residue(graft_points[1]).atom("O").xyz() - scaffold_pose.residue(scaffold_strand_pair[1]).atom("O").xyz()).norm()              
                        pair1_CA_dist = (target_complex_pose.residue(graft_points[0]).atom("CA").xyz() - scaffold_pose.residue(scaffold_strand_pair[0]).atom("CA").xyz()).norm()
                        pair2_CA_dist = (target_complex_pose.residue(graft_points[1]).atom("CA").xyz() - scaffold_pose.residue(scaffold_strand_pair[1]).atom("CA").xyz()).norm()
                        print(f'Graft rmsd: {rmsd:.2f} pair1_O_dist: {pair1_O_dist:.2f} pair2_O_dist: {pair2_O_dist:.2f} pair1_CA_dist: {pair1_CA_dist:.2f} pair2_CA_dist: {pair2_CA_dist:.2f}') 

                totalgrafts = totalgrafts + 1

        for f in single_chain_scaffolds:
            if os.path.exists(f):
                os.remove(f)
        
print(f'Total grafts: {totalgrafts}')
print(f'Total ideal grafts: {totalidealgrafts}')


