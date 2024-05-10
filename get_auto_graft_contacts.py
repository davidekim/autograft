import os, sys
import pyrosetta
from utils import *
import glob
import re
import time
import subprocess
from pyrosetta import *

contact_dist = 10
grafts_per_job = 500
slurm_queue = "cpu"
slurm_time_limit = "04:00:00"

jobindex = 0
graftdir = "./"
if len(sys.argv) == 2:
    graftdir = sys.argv[1]
elif len(sys.argv) == 3 and sys.argv[2].isnumeric():
    graftdir = sys.argv[1]
    jobindex = int(sys.argv[2])
else:
    print(f'USAGE: {sys.argv[0]} <grafts dir or list>')
    print()
    print(f'Submits {grafts_per_job} grafts per job to slurm,')
    print(f'waits for jobs to complete (time limit {slurm_time_limit}),')
    print(f'and then ranks the grafts by contacts.')
    print()
    print(f'Job id is the grafts directory name.')
    print(f'<job id>.grafts and <job id>.contacts is output.')
    print()
    sys.exit()


# make graft dir name the job id
if os.path.isdir(graftdir):
    if graftdir[-1] != '/':
        graftdir = graftdir + '/'
    jobid = os.path.dirname(graftdir).split('/')[-1]
else:
    jobid = os.path.basename(graftdir)
 
graftsfile = jobid+".grafts"
grafts = []
if os.path.exists(graftsfile):
    with open(graftsfile) as tempf:
        for line in tempf:
            grafts.append(line.strip())
else:
    grafts = get_grafts(graftdir)
    if not os.path.exists(graftsfile):
        gf = open(graftsfile, "w")
        for graft in grafts:
            gf.write(graft+"\n")
        gf.close()

total_grafts = len(grafts)
print(f'total grafts in {graftdir}: {total_grafts}')

def is_done():
    donecnt = 0
    for f in glob.glob(jobid + '_*_tmp.contacts'):
        with open(f) as tempf:
            for line in tempf:
                donecnt = donecnt + 1
    if donecnt >= total_grafts:
        return True
    else:
        return False

def get_results():
    graft_contacts = {}
    graft_avgdegree = {}
    total_graft_contacts = 0
    for f in glob.glob(jobid + '_*_tmp.contacts'):
        print(f'Reading {f}')
        with open(f) as tempf:
            for line in tempf:
                cols = line.strip().split(" ")
                graft_contacts[cols[0]] = int(cols[1])
                graft_avgdegree[cols[0]] = float(cols[2])
                total_graft_contacts = total_graft_contacts + 1
    if total_graft_contacts != total_grafts:
        print("Error! incomplete")
        exit(1)
    graft_contacts = dict(sorted(graft_contacts.items(), key=lambda item: item[1], reverse=True))
    outputfile = jobid + ".contacts"
    f = open(outputfile, "w")
    for graft in graft_contacts:
        f.write(f'{graft} {graft_contacts[graft]} {graft_avgdegree[graft]}'+"\n")
    f.close()
    for f in glob.glob(jobid + '_*_tmp.contacts'):
        os.remove(f)
    print("Done!")
    exit(0)

if jobindex > 0:
    jobid = jobid + '_' + str(jobindex) + "_tmp"
    # run a subset of the grafts (for a slurm array job)
    xindex = (grafts_per_job * jobindex) - grafts_per_job
    yindex = grafts_per_job * jobindex
    if yindex > len(grafts):
        grafts = grafts[xindex:]
    else:
        grafts = grafts[xindex:yindex]
else:
    # setup array job and submit
    # create slurm array tasks and submit
    tasksf = open(jobid+"_contact_tasks.txt", "w")
    cnt = 1
    while grafts_per_job*cnt < len(grafts)+grafts_per_job:
        tasksf.write(f'python {sys.argv[0]} {sys.argv[1]} {cnt}'+"\n")
        if grafts_per_job*cnt >= len(grafts):
            break
        cnt = cnt + 1
    tasksf.close()

    if len(glob.glob(jobid + '_*_tmp.contacts')) == cnt:
        get_results() # exits

    # create slurm sbatch file
    sbf = open(jobid+"_contact_sbatch.sh", "w")
    sbf.write("#!/bin/bash\n#SBATCH -p "+slurm_queue+"\n#SBATCH --mem=8g\n#SBATCH -t "+slurm_time_limit+"\n")
    sbf.write("source activate /home/nrbennet/miniconda3/envs/ampere\n")
    sbf.write('CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" '+jobid+'_contact_tasks.txt)'+"\n")
    sbf.write('echo "${CMD}" | bash'+"\n")
    sbf.close()
    # submit
    sbatchout = subprocess.check_output('sbatch -a 1-$(cat '+jobid+'_contact_tasks.txt|wc -l) '+jobid+'_contact_sbatch.sh', shell=True)
    match = re.search("batch\s+job\s+(\d+)", sbatchout.decode("utf-8"))
    batchid = 0
    if match is not None:
        batchid = match.group(1)
    if batchid == 0:
        print("Error! cannot get batch id")
        exit(1)
    print(f'Submitted sbatch id: {batchid}')
    time.sleep(60)
    while not is_done():
        time.sleep(60)
        done = 1
        try:
            squeue = subprocess.check_output('squeue --job '+batchid, shell=True)
            for l in squeue.decode("utf-8").split("\n"):
                if l.strip().startswith(+batchid+"_"):
                    done = 0
                    break
        except:
            done = 1
        if done:
            break

    time.sleep(15)

    if not is_done():
        print("Error! could not finish")
        exit(1)
    get_results() # exits


pyrosetta.init('-mute all -detect_disulf false -out::file::renumber_pdb -in::ignore_waters -corrections::beta_nov16')

print(f'Ranking {len(grafts)} grafts.')
graft_contacts = {}
graft_avgdegree = {}
count = 0
for graft in grafts:
    if count and count % 100 == 0:
        print(count)
    graft_chain = 'A'
    graft_start = 0
    graft_end = 0
    with open(graft) as tmpf:
        for line in tmpf:
            if line[0:6] == "REMARK":
                match = re.search("CHAIN:\s+(\S+)\s+GRAFT:\s+(\d+)\-(\d+)", line)
                if match is not None:
                    graft_chain = match.group(1)
                    graft_start = int(match.group(2))
                    graft_end = int(match.group(3))

    print(f'graft_chain: {graft_chain} start: {graft_start} end: {graft_end}')

    assert(graft_end > 0), f'Cannot parse graft info from {graft}'

    p = pyrosetta.pose_from_file(graft)

    graft_chainid = 0
    chain_isprotein = {}
    for cid in range(1, p.num_chains()+1):
        if pyrosetta.rosetta.core.pose.get_chain_from_chain_id( cid, p) == graft_chain:
            graft_chainid = cid
            break

    # count contacts
    contacts = []
    rescontacts_sum = 0
    interface_rescnt = 0
    for i in pyrosetta.rosetta.core.pose.get_chain_residues(p, graft_chainid):
        rescontacts = 0
        atomidi = 2
        if not i.is_protein():
            atomidi = 1
        for cid in range(1, p.num_chains()+1):
            if cid == graft_chainid:
                continue
            for j in pyrosetta.rosetta.core.pose.get_chain_residues(p, cid):
                if not j.is_protein():
                    for k in range(1, j.natoms()+1):
                        if (j.atom(k).xyz() - i.atom(atomidi).xyz()).norm() < contact_dist:
                            contacts.append((i.seqpos,j.seqpos()))
                            rescontacts = rescontacts + 1
                            break
                else:
                    if (j.atom(2).xyz() - i.atom(atomidi).xyz()).norm() < contact_dist:
                        contacts.append((i.seqpos,j.seqpos()))
                        rescontacts = rescontacts + 1

        if rescontacts:
            interface_rescnt = interface_rescnt + 1
            rescontacts_sum = rescontacts_sum + rescontacts
        
    graft_avgdegree[graft] = rescontacts_sum/interface_rescnt 
    graft_contacts[graft] = contacts
    count = count + 1

graft_contacts = dict(sorted(graft_contacts.items(), key=lambda item: len(item[1]), reverse=True))

outputfile = jobid + ".contacts"
f = open(outputfile, "w")
f.write("graft contacts avgdegree\n")
for graft in graft_contacts:
    f.write(f'{graft} {len(graft_contacts[graft])} {graft_avgdegree[graft]}'+"\n")
f.close()




