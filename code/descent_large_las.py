#!/usr/bin/env sage

from sage.all import *
from cado.scripts import descent
from helpers import timeprint, TakenLineMissing
import multiprocessing
import argparse
import json
import json_custom
import re
import itertools
import os
import sys
import random
import glob
from helpers import silent_remove, masked_target, construct_S, CadoExplainRenumberFile, generate_or_load_target
from helpers import truncate_S, sanity_check_descent_outfile, handle_very_large_special_q, extract_taken_relations
from contextlib import redirect_stdout,redirect_stderr
from candy import print_command_line,warning_message,major_message,error_message
from cado_sage import CadoPolyFile
from cado_nfs_binaries import CadoNFS,CadoNFSBinaries
import time
import descent_ecm_utils
import functools
import subprocess
from pathlib import Path
from helpers import slurm_wait
from timing import *


def get_las_descent_cmd_str(params, todofile):

    las_descent_lpb0 = params.parameters.get('LAS_DESCENT_UNTIL_LPB0', 'LPB0')
    las_descent_lpb1 = params.parameters.get('LAS_DESCENT_UNTIL_LPB1', 'LPB1')

    commandlist = [
        params.dirs['CADO_BUILD_DIR'] + 'sieve/las_descent',
        '--recursive-descent',
        '--allow-largesq',
        '--never-discard',
        '--adjust-strategy 2',
        '--fb1', params.files['CAPPED_FBGZ'],
        '-poly', params.files['POLYFILE'],
        '--descent-max-increase-A', str(params.parameters.get('descent_max_increase_A', 2)),
        '--descent-max-increase-lpb', str(params.parameters.get('descent_max_increase_lpb', 6)),
        '--descent-hint-table', params.files['HINTFILE'],
        '--I', str(params.parameters['I_sieving']),
        '--lim0', str(min(params.BOUNDR, 2**31)),
        '--lim1', str(min(params.BOUNDA, 2**31)),
        '--lpb0', str(las_descent_lpb0),
        '--mfb0', str(params.parameters['desc.mfb0']),
        '--lpb1', str(las_descent_lpb1),
        '--mfb1', str(params.parameters['desc.mfb1']),
        '--bkthresh1', str(params.parameters.get('desc.bkthresh1', min(2**31, params.BOUNDR, params.BOUNDA))),
        "-bkmult", params.parameters.get('desc.bkmult', "1s:1.1"),
        '-t', params.parameters['desc.thr'],
        '--B 16',
        '--todo', todofile
    ]
    command = ' '.join(commandlist)
    return command


# Yes, it's annoying to have a second slurmit.
# However we want to save these files in the desc/ directory, not the slurm/ one.
def descent_slurmit(params, command, job_jobfile_name, job_name, job_outfile_name, job_errfile_name):
    with open(job_jobfile_name,"w") as f:
        f.writelines("#!/bin/bash\n")
        f.writelines("#SBATCH -n 1\n")
        f.writelines(f"#SBATCH --job-name {job_name}\n")
        f.writelines("#SBATCH --requeue\n")
        f.writelines(f"#SBATCH --exclude={params.slurm_exclude}\n")
        f.writelines("#SBATCH --exclusive\n")
        partition = params.slurm_job_partition
        f.writelines(f"#SBATCH --partition={partition}\n")
        f.writelines(f"#SBATCH --output={job_outfile_name}\n")
        f.writelines(f"#SBATCH -e {job_errfile_name}\n")
        f.writelines("\n")
        f.writelines("export DOT_SAGE=/tmp/$USER.sage/\n")
        f.writelines("set -e\n")
        f.writelines("\n")
        f.writelines('srun time -p ' + command)
        f.writelines("\n")

    result = subprocess.run(["sbatch", job_jobfile_name],stdout=subprocess.PIPE)
    timeprint("Submitted slurm job:")
    major_message(command)

    pattern = rb'Submitted batch job (\d+)'
    match = re.search(pattern,result.stdout)
    if match:
        jobnum = int(match.group(1))
        timeprint("job number is",jobnum)
    return jobnum, job_errfile_name


@timing
def start_descent_large_las(params, specified_descent_init_file):
    with open(specified_descent_init_file, "r") as fp:
        init_dict = json.load(fp)

    seed_ten = str(init_dict['seed'])[:10]
    working_dir = params.dirs['DESC'] + 'seed' + seed_ten
    Path(working_dir).mkdir(exist_ok=True)
    working_dir += "/"

    numjobs = params.parameters['desc.large_las.numjobs']
    all_todo_filename = init_dict['todofilename']
    all_todo_qs = []

    las_descent_lpb0 = params.parameters.get('LAS_DESCENT_UNTIL_LPB0', 'LPB0')
    las_descent_lpb1 = params.parameters.get('LAS_DESCENT_UNTIL_LPB1', 'LPB1')

    with open(all_todo_filename, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            lineinfo = line.split()
            if lineinfo[0] == '0':
                qsize = Integer(lineinfo[1]).nbits()
                if qsize > las_descent_lpb0:
                    all_todo_qs.append((0, int(lineinfo[1])))
                else:
                    # This happens if the q is between LPB0 and LAS_DESCENT_UNTIL_LPB0.
                    # In this case we want to handle the q later, in the bottom,
                    # rather than during las_descent.
                    continue
            elif lineinfo[0] == '1':
                qsize = Integer(lineinfo[1]).nbits()
                if qsize > las_descent_lpb1:
                    all_todo_qs.append((1, int(lineinfo[1]), int(lineinfo[2])))
                else:
                    continue

    todo_jobs = dict()

    numjobs = len(all_todo_qs)  # makes things simpler

    if len(all_todo_qs) <= numjobs:
        # Easy: 1 q for each job (and todofile)
        for i in range(len(all_todo_qs)):
            todo_jobs[i] = [ all_todo_qs[i] ]

    if len(all_todo_qs) > numjobs:
        qs_per_job = ceil(len(all_todo_qs) / numjobs)

        for j in range(numjobs):
            todo_jobs[j] = all_todo_qs[ j*qs_per_job : min(len(all_todo_qs), (j+1)*qs_per_job) ]

    slurmlist = []

    for jobi in todo_jobs.keys():
        job_todofile_name = working_dir + 'startlas.job.' + str(jobi) + '.todo'
        job_outfile_name = working_dir + 'startlas.job.' + str(jobi) + '.out'
        job_errfile_name = working_dir + 'startlas.job.' + str(jobi) + '.err'
        job_jobfile_name = working_dir + 'startlas.job.' + str(jobi) + '.job'
        job_name = 'llas-' + str(jobi)

        with open(job_todofile_name, 'w') as tdf:
            for todoitem in todo_jobs[jobi]:
                line = " ".join([str(x) for x in todoitem])
                tdf.write(line)
                tdf.write("\n")

        command = get_las_descent_cmd_str(params, job_todofile_name)
        jobnum, errfile = descent_slurmit(params, command, job_jobfile_name, job_name, job_outfile_name, job_errfile_name)
        slurmlist.append((jobnum, errfile))

    finished_processes, cputime_slurm = slurm_wait(slurmlist)
    overall_cputime.add(cputime_slurm)


@timing
def rebalance_descent_large_las(params, seed_ten, counter, ONLY_RETURN_TODOS=False):
    # seed_ten is the 10-digit seed (string).
    # The counter is 1 or greater. 1 indicates that we look at "start" files,
    # and will write to "1" files. 2 indicates that we look at "1" files, and will
    # write to "2" files. And so on.

    # Assume slurm jobs have been killed, and we have a bunch of
    # partial *.out files, each with an associated todo file.
    # It's possible some of the files are already done.

    numjobs = params.parameters['desc.large_las.numjobs']
    working_dir = params.dirs['DESC'] + 'seed' + seed_ten + "/"

    counter = int(counter)

    if counter == 1:
        prior_outfile_regex = working_dir + 'startlas.*.out'
    else:
        prior_outfile_regex = working_dir + f'rebalance.ctr{counter-1}.job*.out'

    prior_outfiles = glob.glob(prior_outfile_regex)
    latest_output = [
        ( xx[ : len(xx)-3 ] + 'todo', xx ) for xx in prior_outfiles
    ]

    next_outfile_regex = working_dir + f'rebalance.ctr{counter}.job*.out'
    assert len(glob.glob(next_outfile_regex)) == 0

    global_todolist = []

    for (jobtodo, jobout) in latest_output:
        remaining_todolist = []
        timeprint("Looking at output: " + str(jobout))
        num_taken = 0

        with open(jobtodo, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                lineinfo = line.split()
                if lineinfo[0] == '0':
                    remaining_todolist.append((0, int(lineinfo[1])))
                elif lineinfo[0] == '1':
                    remaining_todolist.append((1, int(lineinfo[1]), int(lineinfo[2])))

        with open(jobout, 'r') as f:
            for line in f.readlines():
                line = line.strip()

                if m := re.match(r".* pushing side-(\d+) q=(\d+); rho=(\d+) .* to todo list .*", line):
                    side, sq, rho = (int(c) for c in m.groups())
                    if side == 1:
                        remaining_todolist.append((side, sq, rho))
                    elif side == 0:
                        remaining_todolist.append((side, sq))

                if m := re.match(r"# Taking decision on .* side-(\d+) q=(\d+); rho=(\d+)", line):
                    # To be stored if next line matches if condition below
                    side, sq, rho = (int(c) for c in m.groups())

                if (m := re.match("^Taken: (.*)", line)):
                    if side == 1:
                        remaining_todolist.remove((side, sq, rho))
                    elif side == 0:
                        remaining_todolist.remove((side, sq))
                    num_taken += 1

        timeprint("Taken relations: " + str(num_taken))
        timeprint("Remaining todolist: " + str(len(remaining_todolist)))

        global_todolist += remaining_todolist

    timeprint("In total, the size of our remaining todolist is: " + str(len(global_todolist)))
    timeprint("We will redistribute it over " + str(numjobs) + " slurm jobs.")

    if ONLY_RETURN_TODOS:
        return global_todolist

    #global_todolist2 = [ qq for qq in global_todolist if Integer(qq[1]).nbits() > 32 ]
    #global_todolist = global_todolist2
    #timeprint("Removed all 32-bit or smaller qs.")
    #timeprint("Now, the size of our remaining todolist is: " + str(len(global_todolist)))

    job_allocation = dict()

    if len(global_todolist) <= numjobs:
        for i in range(len(global_todolist)):
            job_allocation[i] = [ global_todolist[i] ]
    else:
        # TODO: could be smarter, based on bitsize.
        qs_per_job = ceil(1.0 * len(global_todolist) / numjobs)
        for j in range(numjobs):
            job_allocation[j] = global_todolist[j*qs_per_job : min(len(global_todolist), (j+1)*qs_per_job) ]

    slurmlist = []

    for jobi in job_allocation.keys():
        job_todofile_name = working_dir + 'rebalance.ctr' + str(counter) + '.job' + str(jobi) + '.todo'
        job_outfile_name = working_dir + 'rebalance.ctr' + str(counter) + '.job' + str(jobi) + '.out'
        job_errfile_name = working_dir + 'rebalance.ctr' + str(counter) + '.job' + str(jobi) + '.err'
        job_jobfile_name = working_dir + 'rebalance.ctr' + str(counter) + '.job' + str(jobi) + '.job'
        job_name = 'rebal-' + str(counter) + '-' + str(jobi)

        with open(job_todofile_name, 'w') as tdf:
            for todoitem in job_allocation[jobi]:
                line = " ".join([str(x) for x in todoitem])
                tdf.write(line)
                tdf.write("\n")

        command = get_las_descent_cmd_str(params, job_todofile_name)
        jobnum, errfile = descent_slurmit(params, command, job_jobfile_name, job_name, job_outfile_name, job_errfile_name)
        slurmlist.append((jobnum, errfile))

    finished_processes, cputime_slurm = slurm_wait(slurmlist)
    overall_cputime.add(cputime_slurm)


@timing
def todofile_descent_large_las(params, given_todofile_glob):
    # We look at todos in the given todofile.
    # We'll write to todofile + 'rels.out'. It could be appended somewhere after this call.
    # The purpose of this is really not to launch a ton of jobs;
    # it's to launch a few with tweaked parameters, for example, or on the beefy machines only.

    slurmlist = []

    for given_todofile in glob.glob(given_todofile_glob):

        job_todofile_name = given_todofile
        job_outfile_name = given_todofile + ".rels.out"
        job_errfile_name = given_todofile + ".err"
        job_jobfile_name = given_todofile + ".job"
        job_name = 'desc-todo'

        command = get_las_descent_cmd_str(params, job_todofile_name)
        jobnum, errfile = descent_slurmit(
            params, command, job_jobfile_name, job_name, job_outfile_name, job_errfile_name
        )
        slurmlist.append((jobnum, errfile))

    finished_processes, cputime_slurm = slurm_wait(slurmlist)
    overall_cputime.add(cputime_slurm)


@timing
def finish_descent_large_las(params, specified_descent_init_file):
    # We write all the Taken lines to one file. We need to look at:
    # 1. The large q relations
    # 2. The "start" las_descent relations
    # 3. The "rebalance" las_descent relations, from 1 to the last counter

    with open(specified_descent_init_file, "r") as fp:
        init_dict = json.load(fp)

    seed_ten = str(init_dict['seed'])[:10]
    working_dir = params.dirs['DESC'] + 'seed' + seed_ten + "/"

    largeq_rels_file = init_dict['largeq_rels_file']
    total_file = init_dict['DRELS_FILE']

    writing_file = open(total_file, "w")

    num_largeq_taken = 0

    with open(largeq_rels_file, "r") as f:
        for line in f.readlines():
            writing_file.write(line)
            num_largeq_taken += 0.5

    timeprint(f"Wrote {num_largeq_taken} relations from {largeq_rels_file} to {total_file}")

    start_files = glob.glob(working_dir + 'startlas.*.out')

    for ff in start_files:
        num_ff = 0
        with open(ff, "r") as f:
            for line in f.readlines():
                writing_file.write(line)
                if "Taken:" in line:
                    num_ff += 1
        timeprint(f"Wrote {num_ff} relations from {ff} to {total_file}")

    rebalance_files = glob.glob(working_dir + 'rebalance.ctr*.job*.out')

    for ff in rebalance_files:
        num_ff = 0
        with open(ff, "r") as f:
            for line in f.readlines():
                writing_file.write(line)
                if "Taken:" in line:
                    num_ff += 1
        timeprint(f"Wrote {num_ff} relations from {ff} to {total_file}")

    writing_file.close()


@timing
def descent_rock_bottom(params, specified_descent_init_file, strategy, prev_desc_file=None):
    # We assume all descent relations are in one file (e.g. desc.total.rels from calling ..._finish).
    # If no file is given we assume desc/seed{seed_ten}/desc.total.rels.
    # We extract all non-taken q's that are above our LPBs.
    # For each we launch a -P (embedding polynomial) or a -C (composites) job. Only one needs to work!
    # Finally, we look at all those found relations.
    # We APPEND them to the same prev_desc_file.
    # Now, there may still be some outstanding q's. We write these to a file so they
    # can be looked at.
    # It should be possible to update parameters and run this function again, to hopefully get those
    # remaining q's.

    assert strategy in ['P','C','C1','C2']

    with open(specified_descent_init_file, "r") as fp:
        init_dict = json.load(fp)

    seed_ten = str(init_dict['seed'])[:10]
    working_dir = params.dirs['DESC'] + 'seed' + seed_ten + "/rbjobs/"

    # file where we will append the new relations
    total_file = init_dict['DRELS_FILE']

    # file where we read the existing relations
    if prev_desc_file is None:
        prev_desc_file = total_file

    desired_LPB0 = params.parameters.get('current_desc_rb_lpb0', 'LPB0')
    desired_LPB1 = params.parameters.get('current_desc_rb_lpb1', 'LPB1')

    global_todolist = extract_outstanding_qs(
        params, desired_LPB0, desired_LPB1, prev_desc_file, ONLY_RETURN_TODOS=True
    )

    # If we're calling this function, these better be defined
    las_descent_lpb0 = params.parameters.get('LAS_DESCENT_UNTIL_LPB0', 'LPB0')
    las_descent_lpb1 = params.parameters.get('LAS_DESCENT_UNTIL_LPB1', 'LPB1')

    # Also need to look for q's that were ignored in our original todo
    # list, due to being in between LPB and LAS_DESCENT_UNTIL_LPB
    original_todolist = init_dict['todofilename']
    with open(original_todolist, "r") as og_todos:
        for line in og_todos.readlines():
            line = line.strip()
            lineinfo = line.split()
            if lineinfo[0] == '0':
                qsize = Integer(lineinfo[1]).nbits()
                if qsize <= las_descent_lpb0:
                    global_todolist.add((0, int(lineinfo[1])))
                    timeprint(f"Adding a side-0 q from the original todolist ({qsize})")
                else:
                    continue
            elif lineinfo[0] == '1':
                qsize = Integer(lineinfo[1]).nbits()
                if qsize <= las_descent_lpb1:
                    global_todolist.add((1, int(lineinfo[1]), int(lineinfo[2])))
                    timeprint(f"Adding a side-1 q from the original todolist ({qsize})")
                else:
                    continue

    # Exactly one job per q, we are not combining qs into one job.
    slurmlist = []
    specialq_to_file = dict()

    shuffle_list = list(global_todolist)
    random.shuffle(shuffle_list)        # sometimes stuff is killed early,
                                        # depending on progress. it'd be nice
                                        # to be representative

    for specialq in shuffle_list:
        side = specialq[0]
        q = specialq[1]

        # TEMPORARY
        size_q = Integer(q).nbits()
        #if side == 0 and size_q > 36:
        #    continue
        #if side == 1 and size_q > 38:
        #    continue

        if side == 0:
            rho = None
            job_outfile = working_dir + f'bottomq.{side}.{q}.{strategy}.job.out'
            job_errfile = working_dir + f'bottomq.{side}.{q}.{strategy}.job.err'
            job_jobfile = working_dir + f'bottomq.{side}.{q}.{strategy}.job.slurm'
            jobname = f'bot-{side}-{q}'
        elif side == 1:
            rho = specialq[2]
            job_outfile = working_dir + f'bottomq.{side}.{q}.{rho}.{strategy}.job.out'
            job_errfile = working_dir + f'bottomq.{side}.{q}.{rho}.{strategy}.job.err'
            job_jobfile = working_dir + f'bottomq.{side}.{q}.{rho}.{strategy}.job.slurm'
            jobname = f'bot-{side}-{q}-{rho}'

        specialq_to_file[specialq] = job_outfile + '.rb'

        commandlist = [
            params.files['SAGE'],
            "descent_rb_helper.py",
            "--params", params.files['PARAMS'],
            "--outfile", job_outfile + '.rb',
            "--side", str(side),
            "--q", str(q),
            "--seed", str(seed_ten),
            "--strategy", strategy,
            "--overwrite-lpb0", str(desired_LPB0),
            "--overwrite-lpb1", str(desired_LPB1)
        ]

        if side == 1:
            commandlist.append("--rho")
            commandlist.append(str(rho))

        command_str = ' '.join(commandlist)
        jobnum, errfile = descent_slurmit(params, command_str, job_jobfile, jobname, job_outfile, job_errfile)
        slurmlist.append((jobnum, errfile))

    finished_processes, cputime_slurm = slurm_wait(slurmlist)
    overall_cputime.add(cputime_slurm)

    appending_file = open(total_file, "a")
    num_written = 0

    for specialq in global_todolist:
        side = specialq[0]
        q = specialq[1]

        if side == 0:
            rho = None
        elif side == 1:
            rho = specialq[2]

        if specialq not in specialq_to_file.keys():
            continue

        expected_file = specialq_to_file[specialq]

        if not os.path.exists(expected_file):
            print(f"Note: I don't see the expected file: {expected_file}")
            print("That may be an outstanding q to try again.")
            continue
        else:
            appending_file.write("#\n")
            with open(expected_file, "r") as infile:
                for line in infile.readlines():
                    appending_file.write(line)

            num_written += 1

    appending_file.close()
    timeprint(f"Wrote {num_written} relations to file: {appending_file}")
    timeprint(f"That was out of a total of {len(global_todolist)} qs.")
    timeprint("If there are some left over you can run this command again with different parameters.")


@timing
def extract_outstanding_qs(params, desired_LPB0, desired_LPB1, prev_desc_file, ONLY_RETURN_TODOS=False):
    # This is to be run when we have a file, e.g. desc.total.rels,
    # that is smooth up to a bigger LPB than we ultimately want.
    # We can't generate our todo list from the old descent files,
    # but we can mimic construct_S and make a list of the remaining qs.

    # First look for outstanding qs that were put on the todo list.
    # Probably there aren't any but we may as well check.

    og_poly = CadoPolyFile(params.files['POLYFILE']); og_poly.read()
    g = og_poly.f[0]
    f = og_poly.f[1]
    missed_qs = []

    try:
        sanity_check_descent_outfile(g, prev_desc_file)
    except TakenLineMissing as ex:
        for side,q,rho in ex.missed:
            missed_qs.append((side,q,rho))

    timeprint(f"After a sanity check, the number of missing qs is: {len(missed_qs)}")
    timeprint("This is coming only from the las_descent logs.")

    taken = list(extract_taken_relations(prev_desc_file))

    # At this point, the logged special qs are accounted for.
    # We are now only concerned with large qs that show up in the (righthand)
    # factorizations, that were previously considered done.

    side0_covered_special_qs = set()
    side1_covered_special_qs = set()
    side0_rough_qs = set()
    side1_rough_qs = set()

    # TEMPORARY EDIT:
    # We consider the extension unlinked ideals to be "rough".
    # For 768 we don't have algebraic unlinked ideals to worry about, but we
    # could do a similar thing if needed in the future.

    ext_unlinked_ideals = set()
    ext_unlinked_files = [
        # none for now
    ]
    for filename in ext_unlinked_files:
        with open(filename, 'r') as f:
            for line in f:
                lineinfo = line.strip().split()
                assert lineinfo[0] == '0'       # 1 sided sieving
                p = Integer(lineinfo[1])
                r = Integer(lineinfo[2])
                ext_unlinked_ideals.add( (1,p,r) )

    #with open(params.files['EXT_UNLINKED_TODOS'], "r") as extfile:
    #    for line in extfile.readlines():
            # it's a todo file
    #        lineinfo = line.strip().split(" ")
    #        assert lineinfo[0] == '0'
    #        q = Integer(lineinfo[1])
    #        r = Integer(lineinfo[2])
    #        ext_unlinked_ideals.add((1,q,r))

    for rel, (side, sq, rho) in taken:
        a = rel.a
        b = rel.b

        if side == 0:
            side0_covered_special_qs.add((side,sq))
        elif side == 1:
            side1_covered_special_qs.add((side,sq,rho))

        if side == 0:
            fac0 = [ p for p in rel.sides[0] if p != sq ]
            fac1 = rel.sides[1]
        elif side == 1:
            fac0 = rel.sides[0]
            fac1 = [ p for p in rel.sides[1] if p != sq ]

        for p in fac0:
            p = Integer(p)
            if p.nbits() > desired_LPB0:
                r = ZZ(g.roots(GF(p))[0][0])
                side0_rough_qs.add((0,p))

        for p in fac1:
            p = Integer(p)

            if gcd(b,p) != 1:
                r = p
            else:
                r = (a * inverse_mod(b,p)) % p

            if p.nbits() > desired_LPB1:
                side1_rough_qs.add((1,p,r))

            if (1,p,r) in ext_unlinked_ideals:
                side1_rough_qs.add((1,p,r))

    final_side0_rough = side0_rough_qs - side0_covered_special_qs
    final_side1_rough = side1_rough_qs - side1_covered_special_qs

    side0_stats = dict()
    side1_stats = dict()

    for (side,q) in final_side0_rough:
        nb = Integer(q).nbits()
        if nb not in side0_stats.keys():
            side0_stats[nb] = 1
        else:
            side0_stats[nb] += 1

    for (side,q,rho) in final_side1_rough:
        nb = Integer(q).nbits()
        if nb not in side1_stats.keys():
            side1_stats[nb] = 1
        else:
            side1_stats[nb] += 1

    timeprint("Stats for side 0:\n" + str(side0_stats))
    timeprint("Stats for side 1:\n" + str(side1_stats))

    if ONLY_RETURN_TODOS:
        return final_side0_rough | final_side1_rough

    new_todofile = prev_desc_file + f".lpbs.{desired_LPB0}.{desired_LPB1}.stilltodo"
    num0 = 0
    num1 = 0

    with open(new_todofile, "w") as outfile:

        for (side,q) in final_side0_rough:
            outfile.write(f"0 {q}\n")
            num0 += 1

        for (side,q,rho) in final_side1_rough:
            outfile.write(f"1 {q} {rho}\n")
            num1 += 1

    timeprint(f"Wrote {num0} rational and {num1} algebraic outstanding qs to {new_todofile}")
