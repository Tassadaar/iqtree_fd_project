#!/usr/bin/env python

"""
Copyright 2024 Jason Shao & Joran Martijn.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.
"""

# standard modules
import io, os, sys, re
import argparse
import subprocess
import itertools
import logging
import shutil
import time
import concurrent.futures
import psutil

# stuff you need to install
from ete3 import Tree, TreeNode #, TreeStyle
from ete3.parser.newick import NewickError
from Bio import AlignIO
from dask_jobqueue import SGECluster
from dask.distributed import Client
from dask import delayed

# for type annotation
from typing import List, Set, Tuple

# argument parser
parser = argparse.ArgumentParser(description="Wrapper for FunDi fix")

# input
parser.add_argument( "-te", "--tree",
    required=True,
    help="Tree file in newick format, must be rooted"
)
parser.add_argument( "-s", "--alignment",
    required=True,
    help="Alignment in a supported format"
)
parser.add_argument( "-d", "--definition",
    required=True,
    help="Definition file that splits the tree and alignment by FunDi branch"
)
parser.add_argument( "-i", "--increment",
    required=False, default=0.05, type=float,
    help="Partition top branches of subtrees in sections of X subs/site when stitching subtrees together in alternate full trees, default is 0.05"
)
parser.add_argument( "-m", "--model",
    required=False, default="LG+C10+G",
    help="Model to be used with iqtree. Syntax: RateMatrix+MixtureModel+RateHeterogeneity. Default is LG+C10+G"
)
parser.add_argument( "-mdef", "--nexus",
    required=False, default=None,
    help="Custom model file in NEXUS format"
)
parser.add_argument( "-strat", "--strategy",
    required=False, type=str, default='subtree_branch_lengths',
    help="Are we using branch lengths from the full tree ('fulltree_branch_lengths') or the subtrees (subtree_branch_lengths)?"
)
parser.add_argument( "--seed",
    required=False, default=654321, type=int,
    help="Seed given to all IQTREE runs in this wrapper. Defaults to 654321"
)
parser.add_argument( "--iqtree_path",
    required=False, default='iqtree', type=str,
    help="Path to iqtree binary you would like to use. Defaults simply to 'iqtree'"
)

# performance
parser.add_argument( "-nt", "--cores",
    required=False, type=int, default=2,
    help="Number of CPU cores to use"
)

# running trees in parallel on a computer cluster
parser.add_argument( "-q", "--jobqueue",
    required=False, type=str, default=None,
    help="If running on a computer cluster, the jobqueue to submit to. Can specify multiple queues delimited by comma"
)
parser.add_argument( "-j", "--jobs",
    required=False, type=int, default=None,
    help="If running on a computer cluster, the maximum number of jobs to distribute the workload to. Choose a number that prevents memory overload"
)

# running trees in parallel on a single node/workstation
parser.add_argument( "-par", "--parallelization",
    required=False, action="store_true",
    help="If running on a single workstation, choose to run the funDi analyses in parallel"
)
parser.add_argument( "-mem", "--memory",
    required=False, type=int, default=8,
    help="Amount of memory in gigabytes to use"
)

# miscellaneous
parser.add_argument( "-k", "--keep",
    required=False, action="store_true",
    help="Keeps files associated with trees which are not the best"
)

# output
parser.add_argument( "-o", "--outdir",
    required=True,
    type=str,
    help="Name of output directory"
)

arguments = parser.parse_args()


# create out directory
try:
    os.makedirs(arguments.outdir, exist_ok=True)
except PermissionError:
    print(f"Permission denied: Unable to create '{arguments.outdir}'")

# configure logger
# strip directory path from filename
base = os.path.basename(arguments.alignment)
logging.basicConfig(
        filename=f'{arguments.outdir}/{base}.fundi_wrapper.log',
        filemode='w',
        level=logging.DEBUG,
        # level=logging.INFO,
        format='%(asctime)s -- %(levelname)s: %(message)s'
)


# log how fundi_wrapper was called
logging.info('fundi_wrapper.py was called with the following parameters:\n')
for msg in (
        f'--iqtree_path     {arguments.iqtree_path}\t'
        # f'--iqtree_path     iqtree2\t'
        f'-te/--tree        {arguments.tree}\t'
        f'-s/--alignment    {arguments.alignment}\t'
        f'-d/--definition   {arguments.definition}\t'
        f'-i/--increment    {arguments.increment}\t'
        f'-m/--model        {arguments.model}\t'
        f'-mdef/--nexus     {arguments.nexus}\t'
        f'-strat/--strategy {arguments.strategy}\t'
        f'--seed            {arguments.seed}\t'
        f'-o/--outdir       {arguments.outdir}\t'
        f'-nt/--cores       {arguments.cores}\t'
        f'-k/--keep         {arguments.keep}\n'
).split('\t'):
    logging.info(msg)
# logging.info('\n')

if arguments.jobqueue:
    jobid  = os.environ.get('JOB_ID')
    taskid = os.environ.get('SGE_TASK_ID')
    for msg in (
        'fundi_wrapper.py was invoked to be run parallel on a computer cluster\t'
        f'-q/--jobqueue   {arguments.jobqueue}\t'
        f'-j/--jobs       {arguments.jobs}\t'
        f'JOBID:          {jobid}\t'
        f'TASKID:         {taskid}\n'
    ).split('\t'):
        logging.info(msg)

else:
    logging.info(
        'fundi_wrapper.py was invoked to be run sequentially, either locally or on a computer cluster\n'
    )
# logging.info('\n')


def main(args):

    # --- PROCESS COMMAND LINE ARGUMENTS ---
    logging.info( 'Loading and validating input tree, alignment and definition file' )

    # load tree
    try:
        master_tree = Tree(args.tree)
    except NewickError:
        logging.critical(f'{args.tree} can not be parsed by ETE3, please reformat your tree')
        sys.exit()

    # validate and load alignment and definition file
    alignment = validate_alignment(master_tree, args.alignment)
    defined_groups: List[Set[str]] = validate_def_file(master_tree, args.definition)

    # parse model. Syntax: Rate Matrix+Mixture Model+Rate Heterogeneity
    mixture_model = False
    models: List[str] = args.model.split("+")
    if any(re.match(r'C[1-6]0', m) for m in models):
        mixture_model = True
        logging.info( 'Mixture model detected. After subtree calculation, new mixture models will be calculated and a NEXUS file will be generated' )

    # --- SPLIT TREE AND ALIGNMENT BASED ON DEFINITION FILE ---
    logging.info( 'Splitting input tree and alignment into two subtrees and subalignments\n' )


    # split tree into two subtrees
    a_tree, b_tree = get_ref_subtrees(master_tree, defined_groups)

    # write subtrees into newick files
    a_tree.write(format=1, outfile=f"{args.outdir}/subtree_a.newick")
    b_tree.write(format=1, outfile=f"{args.outdir}/subtree_b.newick")

    # split alignment into two sub-alignments
    a_alignment, b_alignment = write_alignment_partitions(alignment, a_tree, b_tree)
    AlignIO.write(a_alignment, f"{args.outdir}/subtree_a.aln", "fasta")
    AlignIO.write(b_alignment, f"{args.outdir}/subtree_b.aln", "fasta")

    # --- RE-OPTIMIZE ALPHA AND MIXTURE WEIGHTS ON FIXED SUBTREES ---
    logging.info( 'Re-optimizing alpha, mixture weights and branch lengths per subtree' )

    # collect iqtree tasks
    iqtree_tasks = []
    for subtree in ["subtree_a", "subtree_b"]:
        iqtree_command = [
            # "iqtree2",
            args.iqtree_path,
            "-seed", str(args.seed),
            "-s", f"{args.outdir}/{subtree}.aln",
            "-te", f"{args.outdir}/{subtree}.newick",
            "-m", "+".join(models),
            "--prefix", f'{args.outdir}/{subtree}',
            "-nt", str(args.cores),
            "-mem", "200G",
            "-prec", "10",
            "--quiet",
            "--keep-ident",
        ]
        # sometimes the user wants to use a custom model defined in a NEXUS file
        if args.nexus:
            iqtree_command = iqtree_command + ["--mdef", args.nexus]
        if mixture_model:
            iqtree_command = iqtree_command + ["-mwopt"]
        logging.info( ' '.join(iqtree_command) )
        iqtree_tasks.append(iqtree_command)

    # if we've provided a jobqueue argument,
    # the user wants to run the fundi_wrapper
    # on a computer cluster. Hence, use dask
    if args.jobqueue:
        # logging.info( iqtree_tasks )

        # determine network interface
        net_ifs = list( psutil.net_if_addrs().keys() )
        net_ifs.remove('lo')
        scheduler_interface = net_ifs[0]

        # Configure the SGE cluster
        # the cluster object contains a scheduler!
        # see cluster.scheduler
        cluster = SGECluster(
            # worker resources
            processes=1,
            cores=args.cores,
            memory="10GB",
            # job resources
            # job_extra=['-V', f'-pe threaded {str(args.cores)}', '-o /dev/null'],
            job_extra=['-V', f'-pe threaded {str(args.cores)}' ],
            walltime="99:99:99",
            #resource_spec='h_vmem=70G',
            queue=args.jobqueue,
            # interface='eno1',
            scheduler_options={'interface':scheduler_interface},
            death_timeout=120,
            log_directory='dask_logs',
        )
        logging.debug(f'Cluster Name: {cluster.name}')
        logging.debug(f'Cluster Scheduler Address: {cluster.scheduler_address}')
        logging.debug(f'Cluster Scheduler Interface: {scheduler_interface}')
        # Create a Dask client
        # and attach it to the cluster/scheduler
        client = Client(cluster)
        # Scale the cluster to the desired number of workers
        cluster.scale(jobs=2)
        logging.info('SGECluster launched and scaled to two workers')

        # collect tasks to be submitted
        delayed_tasks = [ 
            # delayed(subprocess.run)(task, stderr=subprocess.DEVNULL)
            delayed(subprocess.run)(task, check=True)
            for task in iqtree_tasks
        ]

        client.wait_for_workers(n_workers=2)

        # Compute the delayed tasks in parallel
        futures = client.compute(delayed_tasks)
        logging.debug('Statement after client.compute()')
        # Gather results
        client.gather(futures)
        logging.debug('Statement after client.gather()')
        # Close the two workers
        cluster.scale(jobs=0)

    # otherwise simply run
    # one after the other
    else:
        for task in iqtree_tasks:
            subtree = task[6]
            logging.info(f'Running iqtree for {subtree}')
            subprocess.run(task, stderr=subprocess.DEVNULL)
            logging.info(f'Completed running {subtree}')

    logging.info('Done with the two subtrees\n')

    # check if the new trees generated by iqtree have the same topology
    a_tree = conform_iqtree_tree(f"{args.outdir}/subtree_a.treefile", a_tree.get_children()[0].get_leaf_names())
    b_tree = conform_iqtree_tree(f"{args.outdir}/subtree_b.treefile", b_tree.get_children()[0].get_leaf_names())

    # re-estimate alpha and mixture weights by taking a
    # weighted average of alpha and mixture weights of the first two subtrees.
    # weighted using (n subtree_taxa / n all_taxa)
    a_weight, b_weight  = calculate_weights(a_tree, b_tree)
    avg_alpha: float    = calculate_weighted_average_alpha(f"{args.outdir}/subtree_a.iqtree", f"{args.outdir}/subtree_b.iqtree", a_weight, b_weight)
    if mixture_model:
        avg_mixture_weights = calculate_weighted_average_mixture_weights(f"{args.outdir}/subtree_a.iqtree", f"{args.outdir}/subtree_b.iqtree", a_weight, b_weight)
    #: dict[str,float]

    logging.info(f'Re-optimized alpha: {str(avg_alpha)}')

    # only log mixture weights and create
    # NEXUS file only if we are actually using a mixture model
    nexus_file = None
    if mixture_model:

        logging.info('Re-optimized class weights:')
        logging.info('Class No.\tRe-optimized weight:')
        for class_no, new_weight in avg_mixture_weights.items():
            logging.info(f'{class_no}\t{new_weight}')
            # class_name = mixture_model_freqs[class_no][0]
            # weight_statements.append( f"fundi_{class_name}:1:{str(new_weight)}" )


        # --- CREATE NEXUS MODEL FILE FROM RE-OPTIMIZED ALPHA AND MIXTURE WEIGHTS ---
        logging.info( 'Creating "re-optimized-model.nexus" using re-optimized alpha and mixture weights' )

        # if we are using a custom model,
        # defined in a custom NEXUS file,
        # provided on the command line 
        if args.nexus:
            file_handle = open(args.nexus, "r")
            key_phrase = "begin models;"
        # if we are using a standard LG+Cx0 model
        else:
            file_handle = io.StringIO(MODEL_MIXTURE_DB)
            key_phrase = f"CAT-{models[1]}"

        # combine class frequencies from custom model or LG+Cx0 model with 
        # re-optimized alpha and mixture weights to make a new nexus file, called 're-optimized-model.nexus'
        # nexus_address, models[1] = create_custom_nexus_file(avg_mixture_weights, file_handle, models[1], key_phrase)
        create_custom_nexus_file(avg_mixture_weights, file_handle, models[1], key_phrase, args.outdir)
        nexus_file = f'{args.outdir}/re-optimized-model.nexus'

        # update model name that will refer to the model
        # described in 're-optimized-model.nexus'
        models[1] = f"fundi_{models[1]}"


    # if we use branch lengths from the subtrees,
    # we need to stitch trees together
    if args.strategy == 'subtree_branch_lengths':

        # --- GENERATE LIST OF STITCHED-TOGETHER TREES ---
        logging.info( 'Stitching re-optimized subtrees together into many new full trees' )

        # calculate length of top branches of re-estimated subtrees
        a_branch: float = a_tree.get_children()[0].dist + a_tree.get_children()[1].dist
        b_branch: float = b_tree.get_children()[0].dist + b_tree.get_children()[1].dist
        logging.info(f'Length of top branch subtree a: {a_branch}')
        logging.info(f'Length of top branch subtree b: {b_branch}')

        # generate steps (0.01, 0.02, 0.03 ... branch_length)
        a_steps = float_range(start=args.increment, stop=a_branch, step=args.increment)
        b_steps = float_range(start=args.increment, stop=b_branch, step=args.increment)

        trees = []
        # get alpha and beta from a cartesian product of proportions
        for alpha, beta in itertools.product(a_steps, b_steps):
            # set new branch lengths for a
            a_tree.get_children()[0].dist = alpha
            a_tree.get_children()[1].dist = a_branch - alpha
            # set new branch lengths for b
            b_tree.get_children()[0].dist = beta
            b_tree.get_children()[1].dist = b_branch - beta
            # reconstruct master tree
            new_master_tree = TreeNode(dist=0.1)
            new_master_tree.add_child(a_tree.copy("deepcopy"))
            new_master_tree.add_child(b_tree.copy("deepcopy"))
            assert new_master_tree.robinson_foulds(master_tree)[0] == 0, "The new master tree is not the same!"
            trees.append(new_master_tree)
        # print number of trees generated
        logging.info(f'{len(trees)} trees were generated')


        # --- RUN IQTREE FUNDI WITH THE RE-OPTIMIZED MODEL ON ALL STITCHED TREES --- 
        logging.info( 'Running IQTREE with --fundi on all newly created stitched full trees with the re-optimized alpha, mixture weights and subtrees' )

        # second iqtree (funDi) execution

        ## parallelization (on a single node/workstation) turned out not be more cost effective,
        ## but we keep it here because its a nice example of how to implement parallilzation
        if args.parallelization:
            run_iqtrees_par(
                trees=trees,
                alignment_address=args.alignment,
                avg_alpha=avg_alpha,
                model="+".join(models),
                nexus_file=f'{args.outdir}/re-optimized-model.nexus',
                all_cores=args.cores,
                memory=args.memory,
                leaves=a_tree.get_leaf_names(),
                iqtree_path=args.iqtree_path,
            )

        # if we want to submit trees to a computer cluster
        elif args.jobqueue:
            #nexus_file = f'{args.outdir}/re-optimized-model.nexus' if mixture_model else None
            run_iqtrees(
                cluster=cluster,
                client=client,
                trees=trees,
                alignment_address=args.alignment,
                model="+".join(models),
                avg_alpha=avg_alpha,
                nexus_file=nexus_file,
                seed=args.seed,
                outdir=args.outdir,
                cores=args.cores,
                leaves=a_tree.get_leaf_names(),
                iqtree_path=args.iqtree_path,
            )

        ## sequential mode
        else:
            run_iqtrees(
                trees=trees,
                alignment_address=args.alignment,
                model="+".join(models),
                avg_alpha=avg_alpha,
                nexus_file=f'{args.outdir}/re-optimized-model.nexus',
                seed=args.seed,
                cores=args.cores,
                leaves=a_tree.get_leaf_names(),
                outdir=args.outdir,
                # iqtree_path='iqtree2',
                iqtree_path=args.iqtree_path,
            )

        # --- GENERATE SUMMARY. WHICH STITCHED TREE HAS THE HIGHEST LIKELIHOOD? ---
        logging.info( 'Generating summary file' )

        # To get the number of trees generated, we take number of proportions to the power of 2
        # summary_name = args.alignment.replace('.aln','') + '.summary.txt'
        # summary_name = args.alignment.replace('.fasta','') + '.summary.txt'
        summary_name = os.path.basename(args.alignment) + '.FunDiWrapper.Summary.txt'
        generate_summary(
            summary_filename=summary_name, 
            tree_count=len(trees), 
            model="+".join(models), 
            increment=args.increment, 
            taxa_groups=defined_groups, 
            keep=args.keep, 
            outdir=args.outdir
        )
        # generate summary also creates '___.best_tree.log file'

        # report key results to log
        # best_log_file = args.alignment.replace('.aln','') + '.best_tree.log'
        # best_log_file = args.alignment.replace('.fasta','') + '.best_tree.log'
        best_log_file = os.path.basename(args.alignment) + '.FunDiWrapper.best_tree.log'
        logging.info(f'Best Tree Log File: {best_log_file}')
        alpha, rho, fundi_brlen, fundi_ll, total_tree_len = parse_iqtree_log(f'{args.outdir}/{best_log_file}')

        logging.info('\n')
        for result in (
                f'Model                {"+".join(models)}\t'
                f'Strategy             SubTree\t'
                f'Seed                 {args.seed}\t'
                f'Alpha                {alpha:.3f}\t'
                f'Rho                  {rho:.3f}\t'
                f'FundiBranchLength    {fundi_brlen}\t'
                f'FundiLogLikelihood   {fundi_ll}\t'
                f'TotalTreeLength      {total_tree_len}\t'
        ).split('\t'):
            logging.info(result)

        # logging.info(f'Model\t{"+".join(models)}'
        # logging.info('Model Strategy Seed Alpha Rho FunDiBranchLength FunDiLogLikelihood')
        # logging.info(f'{"+".join(models)} SubTree,{args.seed} {alpha:.3f} {rho:.3f} {fundi_brlen} {fundi_ll}')

    # if however we want to estimate branch lengths
    # during the fundi, we don't need to stitch trees together
    elif args.strategy == 'fulltree_branch_lengths':

        logging.debug('Starting FullTree strategy')

        iqtree_tasks = []
        # we simply just insert the
        # re-optimized alpha and mixture weights
        # from the subtree analysis
        iqtree_command = [
            # "iqtree2",
            args.iqtree_path,
            "-seed", str(args.seed),
            "-s", args.alignment,
            "-te", args.tree,
            "-m", "+".join(models),
            # "-mdef", f'{args.outdir}/re-optimized-model.nexus',
            "-a", str(avg_alpha),
            "--fundi", f"{','.join(defined_groups[0])},estimate",
            "--fundi-init-rho", str(0.5),
            "--fundi-epsilon", '1e-6',
            "--prefix", f'{args.outdir}/fundi_full_tree',
            "-nt", str(args.cores),
            "-prec", "10",
            "--quiet",
        ]
        if mixture_model:
            iqtree_command = iqtree_command + ["-mdef", nexus_file]
        logging.info( ' '.join(iqtree_command) )
        iqtree_tasks.append(iqtree_command)
        logging.debug(iqtree_tasks)

        cluster.scale(jobs=2)
        # collect tasks to be submitted
        delayed_tasks = [
            # delayed(subprocess.run)(task, stderr=subprocess.DEVNULL)
            delayed(subprocess.run)(task)
            for task in iqtree_tasks
        ]

        # Compute the delayed tasks in parallel
        logging.debug('About to compute')
        futures = client.compute(delayed_tasks)
        # futures = client.compute(delayed_task)
        # Gather results
        logging.debug('About to gather')
        client.gather(futures)
        client.shutdown()

        # report key results to log
        best_log_file = 'fundi_full_tree.log'
        alpha, rho, fundi_brlen, fundi_ll, total_tree_len = parse_iqtree_log(f'{args.outdir}/{best_log_file}')

        for result in (
                f'Model                {"+".join(models)}\t'
                f'Strategy             FullTree\t'
                f'Seed                 {args.seed}\t'
                f'Alpha                {alpha:.3f}\t'
                f'Rho                  {rho:.3f}\t'
                f'FundiBranchLength    {fundi_brlen}\t'
                f'FundiLogLikelihood   {fundi_ll}\t'
                f'TotalTreeLength      {total_tree_len}\t'
        ).split('\t'):
            logging.info(result)

        # alpha, rho, fundi_brlen, fundi_ll = parse_iqtree_log(f'{args.outdir}/{best_log_file}')
        # logging.info('Model Strategy Seed Alpha Rho FunDiBranchLength FunDiLogLikelihood')
        # logging.info(f'{"+".join(models)} FullTree {args.seed} {alpha:.3f} {rho:.3f} {fundi_brlen} {fundi_ll}')


def float_range(start: float, stop: float, step: float):
    while start < stop:
        yield start
        start += step


def validate_alignment(tree, alignment_file):

    # alignment formats supported by AlignIO
    ALIGNMENT_FORMATS = [
        "fasta", "phylip-relaxed", "phylip", "clustal", "emboss", "fasta-m10", "ig", 
        "maf", "mauve", "msf", "nexus", "phylip-sequential", "stockholm"
    ]

    # Does the alignment have the exact same set of taxa as the tree?
    def detect_alignment_format():
        for alignment_format in ALIGNMENT_FORMATS:
            try:
                alignment_obj = AlignIO.read(alignment_file, alignment_format)
            except:
                continue
            if isinstance(alignment_obj, AlignIO.MultipleSeqAlignment):
                logging.info(f"Your alignment is of the format {alignment_format}")
                return alignment_obj
        return None

    # parsing fasta file
    alignment = detect_alignment_format()

    if not alignment:
        raise ValueError(f"Unrecognized alignment format! Supported formats are {ALIGNMENT_FORMATS}.")

    seq_ids = set(record.id for record in alignment)
    seen_taxa = set()

    for taxon in tree.get_leaf_names():

        if taxon not in seq_ids:
            raise ValueError(f"Tree taxon {taxon} does not exist in the alignment file")

        if taxon in seen_taxa:
            raise ValueError(f"Tree taxon {taxon} is duplicated in the alignment file")

        seen_taxa.add(taxon)

    return alignment


def validate_def_file(tree, def_file) -> List[Set[str]]:
    # store definition file in memory
    # as a list of two sets of taxa
    # this ensures uniqueness, assuming that order does not matter with def files

    with open(def_file, "r") as file:
        taxa_groups: List[Set[str]] = []
        # taxa_groups = []

        for line in file:
            # clean line by removing trailing whitespace, comma or newline characters
            clean_line: str = line.rstrip(" ,\n")

            # ignore empty sets
            if clean_line == "":
                continue

            # parse clean line delimited by comma or whitespace
            taxon_group: set[str] = set(re.split(r"[,\s]+", clean_line))
            taxa_groups.append(taxon_group)

    # check 1: definition file must not be blank
    if not taxa_groups:
        raise ValueError("Definition file is blank!")

    if len(taxa_groups) > 2:
        raise ValueError("Definition file contains more than 2 taxa groups!")

    leaves: set[str] = set(tree.get_leaf_names())  # assuming leaves are correct as reference

    # auto-complete the list if only one group of taxa is provided
    if len(taxa_groups) == 1:
        taxa_groups.append(leaves - taxa_groups[0])

    # check 2: the two groups of taxa must be non-overlapping, i.e. is there a taxon in both groups?
    if taxa_groups[0] & taxa_groups[1]:
        raise ValueError(f"Defined groups have overlapping items: {taxa_groups[0] & taxa_groups[1]}.")

    # check 3: leaves in the tree must all be in the definition file
    all_taxa = taxa_groups[0] | taxa_groups[1]
    unique_leaves = leaves - all_taxa

    if unique_leaves:
        raise ValueError(f"Some tree leaves do not exist in the definition file: {unique_leaves}.")

    # check 4: are there any taxa in the definition file that do not exist in the tree file?
    unique_taxa = all_taxa - leaves

    if unique_taxa:
        raise ValueError(f"Some taxa do not exist in the tree: {unique_taxa}.")

    # check 5: either definition file partition (taxa group) must have more than 1 taxon
    if len(taxa_groups[0]) == 1 or len(taxa_groups[1]) == 1:
        raise ValueError("Definition file has at least one partition with only 1 taxon! You can't make a subtree from an subalignment with 1 taxon!")


    return taxa_groups


def get_ref_subtrees(master_tree, leaf_groups: List[Set[str]]):
    # this is good practice, but it currently breaks the code
    # master_tree_copy = master_tree.copy("deepcopy")

    # we declare these as None, so that if
    # in the code below they are not replaced with actual objects,
    # we can catch the error later
    a_tree = None
    b_tree = None

    # probe for non-root subtree
    for leaf_group in leaf_groups:
        common_ancestor = master_tree.get_common_ancestor(*leaf_group)

        if not common_ancestor.is_root():
            master_tree.set_outgroup(common_ancestor)
            a_tree, b_tree = master_tree.get_children()

    return a_tree, b_tree


def write_alignment_partitions(alignment, a_tree, b_tree):
    a_leaves = a_tree.get_leaf_names()
    b_leaves = b_tree.get_leaf_names()
    a_alignment = AlignIO.MultipleSeqAlignment([])
    b_alignment = AlignIO.MultipleSeqAlignment([])

    for record in alignment:

        if record.id in a_leaves:
            a_alignment.append(record)

        if record.id in b_leaves:
            b_alignment.append(record)

    return a_alignment, b_alignment


def conform_iqtree_tree(iqtree_treefile, ref_outgroup_leaves):
    # extracting subtree
    iqtree_tree = Tree(iqtree_treefile)

    # if outgroup is a single taxon
    if len(ref_outgroup_leaves) == 1:
        iqtree_tree.set_outgroup(ref_outgroup_leaves[0])

        return iqtree_tree

    # if outgroup is multiple taxa
    outgroup_lca = iqtree_tree.get_common_ancestor(*ref_outgroup_leaves)

    # if the lca of the ref tree outgroup taxa is
    # the root of the new iqtree tree, 
    # the ref outgroup taxa are not monophyletic
    if outgroup_lca.is_root():

        for leaf in iqtree_tree.get_leaf_names():

            # we do an initial reroot here,
            # to ensure that the ref outgroup taxa are monophyletic
            if leaf not in ref_outgroup_leaves:
                iqtree_tree.set_outgroup(leaf)
                break

    # now that they are monophyletic,
    # we can reroot with the reference outgroup taxa
    new_outgroup_lca = iqtree_tree.get_common_ancestor(*ref_outgroup_leaves)
    iqtree_tree.set_outgroup(new_outgroup_lca)

    return iqtree_tree


# tuple[float, float] syntax only workes from 3.9 onwards
# def calculate_weights(a_tree, b_tree) -> tuple[float, float]:
def calculate_weights(a_tree, b_tree):
    a_taxa_count: int = len(a_tree.get_leaf_names())
    b_taxa_count: int = len(b_tree.get_leaf_names())
    total_taxa_count: int = a_taxa_count + b_taxa_count
    return a_taxa_count / total_taxa_count, b_taxa_count / total_taxa_count


def calculate_weighted_average_mixture_weights(
    a_iqtree_file,
    b_iqtree_file,
    a_weight: float,
    b_weight: float,
):
# ) -> dict[str, float]:

    def get_mixture_weights(iqtree_file):
        mixture_weights = {}
        if not os.path.isfile(iqtree_file):
            raise FileNotFoundError(f"'{iqtree_file}' does not exist!")
        with open(iqtree_file, "r") as file:
            for line in file:
                # NOTE: this will fail if the model provided does not exist
                if "No  Component      Rate    Weight   Parameters" in line:
                    next_line = next(file)
                    while next_line != "\n":
                        words = next_line.split()
                        class_no, class_weight = words[0], words[3]
                        mixture_weights[class_no] = float(class_weight)
                        next_line = next(file)
                    break
        return mixture_weights

    # get weights from iqtree log file, returns empty set if not found
    # a_mixture_weights: dict[str, float] = get_mixture_weights(a_iqtree_file)
    # b_mixture_weights: dict[str, float] = get_mixture_weights(b_iqtree_file)
    a_mixture_weights = get_mixture_weights(a_iqtree_file)
    b_mixture_weights = get_mixture_weights(b_iqtree_file)

    # usefulness of these safety nets questionable,
    # but doesn't hurt
    if not a_mixture_weights:
        raise ValueError(f"Cannot extract mixture weights, check if '{a_iqtree_file}' is formatted correctly!")
    if not b_mixture_weights:
        raise ValueError(f"Cannot extract mixture weights, check if '{a_iqtree_file}' is formatted correctly!")

    # mixture_weight = weight of the mixture model class
    # weight = weight used to calculate the weighted average
    avg_mixture_weights = {}

    # a_mixture_weights and b_mixture_weights have the same keys (class numbers)
    for class_no, a_mixture_weight in a_mixture_weights.items():
        weighted_avg_weight = a_mixture_weight * a_weight + b_mixture_weights[class_no] * b_weight
        avg_mixture_weights[class_no] = weighted_avg_weight

    return avg_mixture_weights


def calculate_weighted_average_alpha(a_iqtree_file, b_iqtree_file, a_weight, b_weight) -> float:

    # get alpha from iqtree log file, returns None if not found
    def get_alpha(iqtree_file):
        alpha = None

        if not os.path.isfile(iqtree_file):
            raise FileNotFoundError(f"'{iqtree_file}' does not exist!")

        with open(iqtree_file, "r") as file:
            for line in file:
                if "Gamma shape alpha:" in line:
                    words = line.split()
                    alpha = float(words[3])
                    return alpha
        return alpha

    a_alpha: float = get_alpha(a_iqtree_file)
    b_alpha: float = get_alpha(b_iqtree_file)

    if not a_alpha:
        raise ValueError(f"Cannot extract weights, check if '{a_iqtree_file}' is formatted correctly!")
    if not b_alpha:
        raise ValueError(f"Cannot extract weights, check if '{a_iqtree_file}' is formatted correctly!")

    return a_alpha * a_weight + b_alpha * b_weight


def create_custom_nexus_file(
    # new_weights: dict[str, float],
    new_weights,
    nexus_fh,
    mixture_model: str,
    key_phrase: str,
    outdir: str,
) -> None:

    # container for stock class frequencies
    # mixture_model_freqs: dict[str, tuple[str, list[str]]] = {}
    mixture_model_freqs = {}

    i: int = 0
    # load user specified NEXUS file
    # we need its classes' frequencies to build our 're-optimized-nexus.nexus'
    with nexus_fh as nexus:
        section: bool = False
        for line in nexus:
            # start parsing only after
            # 'begin models;' or desired 'CAT-Cx0' is encountered
            if line.startswith(key_phrase):
                section = True
                continue
            # skip empty lines after begin models
            if section and line.startswith("frequency"):
                # link class no. to class name and aa frequencies
                i += 1
                _, class_name, _, *frequencies = line.split()
                mixture_model_freqs[str(i)] = (class_name, frequencies)
                continue
            if section and (line.startswith("model") or 'FMIX' in line):
                break


    # link up newly optimized mixture weights to
    # corresponding class names (C60pi1, C60pi2, etc)
    weight_statements = []
    for class_no, new_weight in new_weights.items():
        class_name = mixture_model_freqs[class_no][0]
        weight_statements.append( f"fundi_{class_name}:1:{str(new_weight)}" )

    # define the model line
    weight_line = f"model fundi_{mixture_model} = FMIX{{{','.join(weight_statements)}}};"

    # write model with new, re-optimized mixture weights to
    # 're-optimized-model.nexus' file
    with open(f"{outdir}/re-optimized-model.nexus", "w") as nex_file:
        nex_file.write("#nexus\nbegin models;\n")
        for class_name, frequencies in mixture_model_freqs.values():
            class_line = f"frequency fundi_{class_name} = {' '.join(frequencies)}\n"
            nex_file.write(class_line)
        nex_file.write(weight_line + "\n")
        nex_file.write("end;")


def run_iqtrees(
        trees,
        alignment_address,
        model: str,
        avg_alpha: float,
        nexus_file,
        seed: int,
        cores: int,
        leaves,
        outdir: str,
        iqtree_path: str,
        cluster=None,
        client=None,
):
    # collect iqtree tasks
    iqtree_tasks = []
    for index, tree in enumerate(trees, start=1):
        formatted_index = f"{index:04d}"
        # render image of stitched-together-tree
        # tree.render(f"subtree_{i}.png") not supported on perun
        tree.write(format=1, outfile=f"{outdir}/tree_{formatted_index}.tree")
        iqtree_cmd = [
            # "iqtree2",
            iqtree_path,
            "-seed", str(seed),
            "-s", alignment_address,
            "-te", f"{outdir}/tree_{formatted_index}.tree",
            "-m", model,
            "-a", str(avg_alpha),
            "-blfix",
            "--fundi-init-rho", str(0.5),
            "--fundi-init-branch", str(0.01),
            "--fundi-epsilon", '1e-6',
            "-nt", str(cores),
            "-mem", "40G",
            "--prefix", f"{outdir}/tree_{formatted_index}",
            "-prec", "10",
            "--quiet",
            "--keep-ident",
            "--fundi", f"{','.join(leaves)},estimate",
            "-redo",
        ]
        if nexus_file is not None:
            iqtree_cmd = iqtree_cmd + ["--mdef", nexus_file]

        iqtree_tasks.append(iqtree_cmd)
        logging.info( ' '.join(iqtree_cmd) )
        # logging.info(iqtree_cmd)

    # submit iqtrees to cluster
    if cluster and client:
        # upscale dask to 10 workers
        # each of which will run one or more trees
        cluster.scale(jobs=arguments.jobs)
        delayed_tasks = [
            delayed(subprocess.run)(task, stderr=subprocess.DEVNULL)
            for task in iqtree_tasks
        ]
        # delayed_tasks.append(delayed(subprocess.run)(iqtree_cmd, stderr=subprocess.DEVNULL))
        # Compute the delayed tasks in parallel
        futures = client.compute(delayed_tasks)
        # Gather results
        client.gather(futures)
        # Optionally, wait for all workers to fully stop
        client.close()
        # Stop all Dask workers
        cluster.scale(jobs=0)  # This should stop all workers
        cluster.close()

    # run iqtrees sequenctially
    else:
        for task in iqtree_tasks:
            current_tree = task[4]
            logging.info(f"Running iqtree funDi on {current_tree} out of {len(trees)} total trees")
            subprocess.run(task, stderr=subprocess.DEVNULL)
            logging.info(f"Completed running {current_tree}")


def run_iqtrees_par(
        iqtree_path,
        trees,
        alignment_address,
        avg_alpha,
        model,
        nexus_file,
        all_cores,
        memory,
        leaves
):
    def get_memory_requirement(log_file):

        if not os.path.isfile(log_file):
            raise FileNotFoundError(f"'{log_file}' does not exist!")

        with open(log_file, "r") as file:

            for line in file:

                if "NOTE: " in line:
                    return int(line.split(" ")[1])

    # generate first iqtree command
    trees[0].write(format=1, outfile="tree_001.tree")
    first_iqtree_command = [
        # "iqtree2",
        iqtree_path,
        "-s", alignment_address,
        "-te", "tree_001.tree",
        "-m", f"{model}{{{avg_alpha}}}",
        "--mdef", nexus_file,
        "-nt", f"{all_cores}",
        "--prefix", "tree_001",
        "-prec", "10",
        "-blfix",
        "--fundi", f"{','.join(leaves)},estimate",
        "--fundi-init-rho", str(0.5),
        "--fundi-init-branch", str(0.01),
        "--fundi-epsilon", '1e-6',
        "--quiet"
        # "-redo",
    ]

    # get the memory requirement for a reconstructed tree, convert to gigabytes
    print(f"Running iqtree funDi for Tree 01 out of {len(trees)} and determining memory requirement.\n")
    subprocess.run(first_iqtree_command, stderr=subprocess.DEVNULL)
    print("Completed running Tree 01 and memory determination.\n")
    mem_req = get_memory_requirement("tree_001.log") * 0.001

    # determine the max numbers of concurrent processes and their cores
    if all_cores < memory / mem_req:
        max_workers = all_cores
        cores = 1
    else:
        max_workers = int(memory / mem_req)
        cores = int(all_cores / max_workers)

    # generate the rest of the iqtree commands
    iqtree_commands = {}
    for index, tree in enumerate(trees[1:], start=2):
        formatted_index = f"{index:02d}"
        # render image of stitched-together-tree
        # tree.render(f"tree_{i}.png") not supported on perun
        tree.write(format=1, outfile=f"tree_{formatted_index}.tree")
        iqtree_command = [
            # "iqtree2",
            iqtree_path,
            "-s", alignment_address,
            "-te", f"tree_{formatted_index}.tree",
            "-m", f"{model}{{{avg_alpha}}}",
            "--mdef", nexus_file,
            "-nt", f"{cores}",
            "--prefix", f"tree_{formatted_index}",
            "-prec", "10",
            "-blfix",
            "--fundi", f"{','.join(leaves)},estimate",
            "--fundi-init-rho", str(0.5),
            "--fundi-init-branch", str(0.01),
            "--fundi-epsilon", '1e-6',
            "-redo",
            "--quiet"
        ]

        iqtree_commands[formatted_index] = iqtree_command

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []

        # running iqtree funDi here by submitting subprocess to executor
        for index, iqtree_command in iqtree_commands.items():
            print(f"Running iqtree funDi for Tree {index} out of {len(trees)} in parallel...")
            futures.append(executor.submit(subprocess.run, iqtree_command, stderr=subprocess.DEVNULL))

        print()  # divider

        # obtaining confirmation of completion
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            tree_index = result.args[4][5:7]
            print(f"Completed running Tree {tree_index}.")

    # delimiter
    print()


def parse_iqtree_log(logfile: str) -> Tuple[float]:
    with open(logfile, "r") as iqtree_log:
        for line in iqtree_log:
            if "Gamma shape alpha" in line:
                words = line.split()
                alpha = float(words[3])
                continue
            if "Best FunDi parameter rho:" in line:
                words = line.split()
                rho_value = float(words[4])
                continue
            if "Best FunDi central branch length:" in line:
                words = line.split()
                central_branch_length = float(words[5])
                continue
            if "FunDi log-likelihood:" in line:
                words = line.split()
                log_likelihood = float(words[2])
                continue
            if "Total tree length:" in line:
                words = line.split()
                total_tree_len = float(words[3])
                break
    return alpha, rho_value, central_branch_length, log_likelihood, total_tree_len


def generate_summary(
        summary_filename,
        tree_count,
        model,
        increment,
        taxa_groups,
        keep,
        outdir: str
) -> None:

    # gather for each generated stitched tree,
    # the rho, central branch length and log likelihood
    tree_properties = []
    for index in range(1, tree_count + 1):
        formatted_index = f"{index:04d}"
        iqtree_log = f'{outdir}/tree_{formatted_index}.log'
        alpha, rho_value, central_branch_length, log_likelihood, total_tree_len = parse_iqtree_log(iqtree_log)
        # with open(f"{outdir}/tree_{formatted_index}.log", "r") as iqtree_log:
        # this is a problem with large datasets discovered and debugged on perun
        if not log_likelihood:
            raise NameError("Insufficient memory allocation!")
        tree_properties.append((formatted_index, alpha, log_likelihood, rho_value, central_branch_length, total_tree_len))

    # sort the trees based on largest funDi log-likelihood
    sorted_tree_properties = sorted(tree_properties, key=lambda attr_tuple: attr_tuple[2], reverse=True)

    # we're overwriting the initial .png image,
    # but only for the best tree
    #
    # style the tree
    # tree_style = TreeStyle()
    # tree_style.show_leaf_name = True
    # tree_style.show_branch_length = True
    # tree_style.branch_vertical_margin = 10
    #
    # # give each internal node an explicit name
    # for node in best_tree.traverse():
    #
    #     if not node.is_leaf():
    #         node.name = "node"
    # best_tree.render(file_name=f"tree_{best_tree_index}.png", tree_style=tree_style, units="px", w=800, h=1000)

    # reroot best tree and overwrite its treefile ?
    best_index: str = sorted_tree_properties[0][0]
    best_tree = conform_iqtree_tree(f"{outdir}/tree_{best_index}.treefile", taxa_groups[0])
    best_tree.write(format=1, outfile=f"{outdir}/tree_{best_index}.treefile")

    # create and write to some_tree.Summary.txt
    with open(f'{outdir}/{summary_filename}', "w") as summary_file:

        summary_file.write(f"Trees generated: {tree_count}\n")
        summary_file.write(f"Model used: {model}\n")
        summary_file.write(f"Increment used: {increment}\n")

        best_alpha: float            = sorted_tree_properties[0][1]
        best_fundi_ll: float         = sorted_tree_properties[0][2]
        best_rho: float              = sorted_tree_properties[0][3]
        best_fundi_branch_len: float = sorted_tree_properties[0][4]
        best_total_tree_len: float   = sorted_tree_properties[0][5]

        summary_file.write(
            f"Best tree: Tree {best_index}\n"
            f"Best gamma shape alpha: {best_alpha:.3f}\n"
            f"Best funDi log-likelihood: {best_fundi_ll:.3f}\n"
            f"Best rho: {best_rho:.3f}\n"
            f"Best Central branch length: {best_fundi_branch_len:.3f}\n"
            f"Best Total Tree Length: {best_total_tree_len:.3f}\n"
        )

        summary_file.write("\nRooted best tree:\n")

        # print tree with branch lengths
        summary_file.write(f"{best_tree.get_ascii(attributes=['name', 'dist'], show_internal=True)}\n\n")
        # summary_file.write(f"See \"tree_{best_tree_index}.png\" for a tree illustration.\n\n")
        summary_file.write("The following list ranks the remaining trees based on best log-likelihood:\n")

        for index in range(1, len(sorted_tree_properties)):
            summary_file.write(f"Tree {sorted_tree_properties[index][0]}:\t"
                               f"FunDi log-likelihood: {sorted_tree_properties[index][2]};\t"
                               f"Rho: {sorted_tree_properties[index][3]};\t"
                               f"Central Branch Length: {sorted_tree_properties[index][4]}\t"
                               f"Total Tree Length: {sorted_tree_properties[index][5]}\n")


    alignment_index = summary_filename.replace('.Summary.txt','')
    shutil.copy(f"{outdir}/tree_{best_index}.treefile", f"{outdir}/{alignment_index}.best_tree.treefile")
    shutil.copy(f"{outdir}/tree_{best_index}.iqtree",   f"{outdir}/{alignment_index}.best_tree.iqtree")
    shutil.copy(f"{outdir}/tree_{best_index}.log",      f"{outdir}/{alignment_index}.best_tree.log")
    shutil.copy(f"{outdir}/tree_{best_index}.ckp.gz",   f"{outdir}/{alignment_index}.best_tree.ckp.gz")

    if not keep:
        # shutil.rmtree(outdir)
        for filename in os.listdir(outdir):
            if filename.startswith('tree_'):
                file_path = os.path.join(outdir, filename)
                os.remove(file_path)



MODEL_MIXTURE_DB = """
#nexus

begin models;

CAT-C10 profile mixture model of Le, Gascuel & Lartillot (2008)
frequency C10pi1 = 0.4082573125 0.0081783015 0.0096285438 0.0069870889 0.0349388179 0.0075279735 0.0097846653 0.1221613215 0.0039151830 0.0125784287 0.0158338663 0.0059670150 0.0081313216 0.0061604332 0.0394155867 0.1682450664 0.0658132542 0.0018751587 0.0041579747 0.0604426865;
frequency C10pi2 = 0.1027763487 0.0418664491 0.0213272051 0.0155943616 0.0149663448 0.0440685478 0.0419667447 0.0138805792 0.0158864807 0.1066076641 0.1131944125 0.0436343681 0.0437800327 0.0180729309 0.0223250701 0.0529608087 0.1081741005 0.0045147205 0.0137373857 0.1606654446;
frequency C10pi3 = 0.0351766018 0.0019678632 0.0016591476 0.0006768741 0.0078706538 0.0016559557 0.0019686768 0.0022420602 0.0012878339 0.3515819591 0.1278183107 0.0018856550 0.0242631753 0.0126221329 0.0029771559 0.0049998099 0.0255378034 0.0011907778 0.0037539283 0.3888636245;
frequency C10pi4 = 0.0408513927 0.0269887074 0.2185648186 0.2333814790 0.0037602852 0.0380451418 0.0901238869 0.1158332065 0.0373197176 0.0025523644 0.0052164616 0.0485017266 0.0022571778 0.0025108218 0.0108333610 0.0804527209 0.0302879995 0.0010815260 0.0069890931 0.0044481118;
frequency C10pi5 = 0.0185492661 0.0062362395 0.0024895723 0.0009775062 0.0070416514 0.0083539447 0.0024891617 0.0028952913 0.0040103982 0.1632422345 0.4443079409 0.0043570878 0.1202815687 0.0733329781 0.0048827648 0.0051642443 0.0131806647 0.0068759784 0.0144734420 0.0968580644;
frequency C10pi6 = 0.1106750119 0.0352190043 0.0405186210 0.1636437899 0.0014834855 0.0877962201 0.2638456592 0.0325228293 0.0163803600 0.0068334902 0.0140679579 0.0677158208 0.0048988133 0.0023256777 0.0298982139 0.0562887953 0.0426922497 0.0010338979 0.0040522304 0.0181078719;
frequency C10pi7 = 0.0522657662 0.0668294648 0.0714836849 0.0297745257 0.0143324928 0.0736540298 0.0388386669 0.0228101108 0.1551638111 0.0187406149 0.0653779932 0.0439469345 0.0207189121 0.0624033021 0.0145475497 0.0549017631 0.0370140058 0.0193756900 0.1110694548 0.0267512268;
frequency C10pi8 = 0.0116587342 0.0050990142 0.0064011054 0.0021742457 0.0105340743 0.0040203734 0.0024251112 0.0034709143 0.0366787049 0.0187185330 0.0676489746 0.0026694717 0.0143534813 0.3650985596 0.0031159927 0.0094848536 0.0073713920 0.0509564551 0.3574858593 0.0206341497;
frequency C10pi9 = 0.0627195947 0.2038782162 0.0428629162 0.0236193294 0.0052662886 0.1098111767 0.0686284994 0.0256174957 0.0332612124 0.0128968249 0.0305627740 0.2270839355 0.0124036991 0.0039181841 0.0140440613 0.0483152469 0.0463378087 0.0025143473 0.0065521118 0.0197062770;
frequency C10pi10 = 0.1145518598 0.0324008908 0.0750614981 0.0416192189 0.0098549497 0.0339624663 0.0364907910 0.0503817581 0.0165233329 0.0092949460 0.0139153707 0.0423026886 0.0082240805 0.0046605982 0.0379221548 0.2610647896 0.1845829279 0.0017548981 0.0058538316 0.0195769483;

model C10 = POISSON+G+FMIX{C10pi1:1:0.1191344178,C10pi2:1:0.0874372456,C10pi3:1:0.1037105070,C10pi4:1:0.0922584809,C10pi5:1:0.1070492801,C10pi6:1:0.1329945166,C10pi7:1:0.0538028458,C10pi8:1:0.0691986212,C10pi9:1:0.1319937434,C10pi10:1:0.1024203429};
model C10Opt = POISSON+G+FMIX{C10pi1,C10pi2,C10pi3,C10pi4,C10pi5,C10pi6,C10pi7,C10pi8,C10pi9,C10pi10};


CAT-C20 profile mixture model of Le, Gascuel & Lartillot (2008)
frequency C20pi1 = 0.0862412505 0.0171943793 0.0791293376 0.0329908619 0.0130504558 0.0169046938 0.0184526503 0.0366905299 0.0108013340 0.0097907148 0.0112826424 0.0220195221 0.0087821483 0.0044155335 0.0189273201 0.3178152357 0.2711700523 0.0015317305 0.0048342853 0.0179753220 ;
frequency C20pi2 = 0.2035582865 0.0050980810 0.0077052407 0.0031656079 0.0348667285 0.0064044073 0.0070859400 0.0195235515 0.0024392035 0.1152573291 0.0789777393 0.0042380850 0.0309187017 0.0112429356 0.0164189221 0.0496777139 0.1118946615 0.0017762569 0.0048448213 0.2849057867 ;
frequency C20pi3 = 0.0211547413 0.0014946177 0.0012755030 0.0005492865 0.0048188557 0.0012328812 0.0014539632 0.0011430874 0.0011346394 0.3928460626 0.1250644210 0.0013579946 0.0209788805 0.0128251737 0.0020247248 0.0026240726 0.0171914121 0.0011591071 0.0036027969 0.3860677787 ;
frequency C20pi4 = 0.0376903543 0.2885196153 0.0365411474 0.0109469400 0.0064073829 0.0893564381 0.0358365464 0.0191106776 0.0329513951 0.0101711878 0.0237495504 0.2897626974 0.0096528870 0.0036349802 0.0105337370 0.0356313768 0.0355926500 0.0027925238 0.0066557222 0.0144621902 ;
frequency C20pi5 = 0.0084597802 0.0053589922 0.0072525884 0.0024487852 0.0084909000 0.0042781483 0.0025055486 0.0024277107 0.0433214027 0.0097713028 0.0380507037 0.0026741007 0.0080724771 0.3420463838 0.0021418673 0.0080418935 0.0055322116 0.0494840193 0.4375001561 0.0121410277 ;
frequency C20pi6 = 0.1759898886 0.0290429175 0.0332845569 0.1301263816 0.0017558693 0.0707183953 0.2182166681 0.0409535143 0.0130708195 0.0085622087 0.0159530702 0.0542946169 0.0054045759 0.0025276980 0.0371020404 0.0793480500 0.0540083424 0.0010592104 0.0036259116 0.0249552645 ;
frequency C20pi7 = 0.1634397322 0.0195541184 0.0438701833 0.0374272612 0.0088659891 0.0137554758 0.0220611924 0.5296717726 0.0090006141 0.0017569353 0.0061156267 0.0167117975 0.0029390787 0.0030641349 0.0126457766 0.0829342776 0.0142835614 0.0028640685 0.0032398299 0.0057985736 ;
frequency C20pi8 = 0.0917468761 0.0265853306 0.0290699087 0.0133818895 0.0284015012 0.0255084506 0.0196875685 0.0249898794 0.0449766405 0.0583555688 0.1155009222 0.0164915955 0.0395994595 0.0998479096 0.0209916159 0.0736482742 0.0661518462 0.0246463919 0.0972327226 0.0831856483 ;
frequency C20pi9 = 0.0646700714 0.0988015996 0.0228907308 0.0168733856 0.0077117603 0.0996414875 0.0544977962 0.0148893975 0.0313851988 0.0505983315 0.1844282999 0.0907931290 0.0774839960 0.0219148172 0.0105004469 0.0321196170 0.0411766062 0.0084303030 0.0206106035 0.0505824221 ;
frequency C20pi10 = 0.0135993865 0.0043408375 0.0018469375 0.0007951703 0.0100090240 0.0046420778 0.0018011758 0.0026794645 0.0072401918 0.0814026713 0.3661422246 0.0025158135 0.0734965132 0.2640965246 0.0038994134 0.0043668760 0.0075248451 0.0261564898 0.0660970801 0.0573472826 ;
frequency C20pi11 = 0.1478036236 0.0842845089 0.0726630217 0.0534743238 0.0048825808 0.0757166156 0.0727246460 0.0907725939 0.0262288856 0.0035781075 0.0126777221 0.1051660098 0.0059621792 0.0029903868 0.0156558198 0.1459903343 0.0634877444 0.0015928454 0.0050760739 0.0092719768 ;
frequency C20pi12 = 0.0186377412 0.0042055165 0.0019865236 0.0008329696 0.0054968852 0.0065890091 0.0020248504 0.0021713483 0.0023665991 0.2020809776 0.4370381920 0.0029120653 0.1241860384 0.0385383157 0.0040672279 0.0046177381 0.0149904396 0.0026871667 0.0056324117 0.1189379840 ;
frequency C20pi13 = 0.0477624336 0.0505742667 0.0209574273 0.0141349161 0.0075791708 0.0429296799 0.0462688073 0.0052327914 0.0165351815 0.1741496627 0.1121253570 0.0577575020 0.0330288046 0.0130691347 0.0124374733 0.0264988925 0.0951754678 0.0031660482 0.0112465746 0.2093704079 ;
frequency C20pi14 = 0.4164189845 0.0056100821 0.0091701381 0.0045131748 0.0406937949 0.0061320495 0.0063229801 0.0946185184 0.0031057404 0.0076443223 0.0099885414 0.0038941773 0.0069323155 0.0048438356 0.0187840756 0.2360774301 0.0746274607 0.0012172579 0.0034825786 0.0459225422 ;
frequency C20pi15 = 0.0402295888 0.0735203003 0.1036647193 0.0365523994 0.0124782975 0.0826558132 0.0372197283 0.0233618081 0.2108307125 0.0093478727 0.0360561493 0.0482410586 0.0100289536 0.0459094917 0.0098503973 0.0533383445 0.0310209005 0.0140076639 0.1064377821 0.0152480184 ;
frequency C20pi16 = 0.0323453034 0.0236282995 0.2520448083 0.2431495959 0.0035976296 0.0330831153 0.0710274499 0.1016074562 0.0366225082 0.0031410809 0.0051980542 0.0470129351 0.0024028744 0.0024429276 0.0094837826 0.0848355278 0.0359083275 0.0008730928 0.0067247672 0.0048704638 ;
frequency C20pi17 = 0.1476256642 0.0334506604 0.0211972524 0.0403051550 0.0032327194 0.0371554480 0.0576893391 0.0330850942 0.0146392559 0.0108267008 0.0256200793 0.0451350877 0.0058651400 0.0047177179 0.3473710507 0.0892065279 0.0485899446 0.0016358749 0.0044177191 0.0282335685 ;
frequency C20pi18 = 0.1031448143 0.0717747663 0.0435172139 0.0386401502 0.0061762467 0.0786603123 0.0923369140 0.0202338419 0.0246761899 0.0376904275 0.0376283678 0.0921698920 0.0161883318 0.0067666433 0.0128302120 0.0951450188 0.1378566702 0.0022144738 0.0083041573 0.0740453560 ;
frequency C20pi19 = 0.0837542823 0.0899383244 0.0518811417 0.0804870571 0.0020735078 0.1456497470 0.1947759184 0.0229030361 0.0268458796 0.0074079756 0.0190249576 0.1459287407 0.0067395241 0.0023063393 0.0085616014 0.0455739585 0.0451080843 0.0010771349 0.0049325333 0.0150302559 ;
frequency C20pi20 = 0.0578735570 0.0138313604 0.0491421636 0.2946738942 0.0011130839 0.0598250358 0.3402102668 0.0293911435 0.0139817004 0.0030525663 0.0062611922 0.0363365043 0.0027295976 0.0017034884 0.0156106390 0.0358044639 0.0249941878 0.0008664342 0.0038312977 0.0087674229 ;

[ C20 with fixed weights ]
model C20 = POISSON+G+FMIX{C20pi1:1:0.0559910600,C20pi2:1:0.0514824870,C20pi3:1:0.0812922124,C20pi4:1:0.0721976867,C20pi5:1:0.0556718858,C20pi6:1:0.0331003080,C20pi7:1:0.0589501763,C20pi8:1:0.0263756889,C20pi9:1:0.0307584220,C20pi10:1:0.0376701125,C20pi11:1:0.0303058290,C20pi12:1:0.0808775576,C20pi13:1:0.0263349134,C20pi14:1:0.0579101455,C20pi15:1:0.0371248064,C20pi16:1:0.0586867766,C20pi17:1:0.0561479138,C20pi18:1:0.0349810886,C20pi19:1:0.0544937394,C20pi20:1:0.0596471901};
[ C20 to weights to be optimized ]
model C20Opt = POISSON+G+FMIX{C20pi1,C20pi2,C20pi3,C20pi4,C20pi5,C20pi6,C20pi7,C20pi8,C20pi9,C20pi10,C20pi11,C20pi12,C20pi13,C20pi14,C20pi15,C20pi16,C20pi17,C20pi18,C20pi19,C20pi20};

model C20Test = POISSON+G+FMIX{C20pi1:1:0.089485,C20pi2:1:0.021281,C20pi3:1:0.119676,C20pi4:1:0.080933,C20pi5:1:0.064054,C20pi6:1:0.021848,C20pi7:1:0.063392,C20pi8:1:0.003629,C20pi9:1:0.007174,C20pi10:1:0.006256,C20pi11:1:0.023424,C20pi12:1:0.086825,C20pi13:1:0.038495,C20pi14:1:0.090028,C20pi15:1:0.020025,C20pi16:1:0.043484,C20pi17:1:0.076864,C20pi18:1:0.031347,C20pi19:1:0.047749,C20pi20:1:0.064031};

CAT-C30 profile mixture model of Le, Gascuel & Lartillot (2008)
frequency C30pi1 = 0.1100453954 0.0171294861 0.0640338464 0.1595411459 0.0019047235 0.0310187088 0.1098958823 0.0684301540 0.0137950707 0.0026283074 0.0073396531 0.0358553674 0.0024706414 0.0016629473 0.1669356820 0.1381790473 0.0568342547 0.0004661120 0.0035970152 0.0082365591;
frequency C30pi2 = 0.0874125465 0.0806320385 0.0382152368 0.0326119879 0.0049826376 0.0798168854 0.0951700809 0.0144042708 0.0210626652 0.0399884450 0.0301585074 0.1147200015 0.0126488911 0.0048996596 0.0137397028 0.0873769666 0.1558616621 0.0015122843 0.0053974463 0.0793880836;
frequency C30pi3 = 0.0225477414 0.0014900535 0.0013034594 0.0005959279 0.0050018158 0.0011436556 0.0015030529 0.0011570953 0.0009374322 0.3944689167 0.0889573138 0.0013600872 0.0189102669 0.0089216031 0.0018312028 0.0028336408 0.0189813395 0.0006693746 0.0023303726 0.4250556480;
frequency C30pi4 = 0.0602158209 0.0136833299 0.0414987935 0.2900084105 0.0009525462 0.0621611083 0.3610869026 0.0281925621 0.0130500799 0.0030516237 0.0060401889 0.0352704692 0.0027460635 0.0014625624 0.0127175499 0.0318109377 0.0225279521 0.0007948027 0.0034024563 0.0093258397;
frequency C30pi5 = 0.0101223637 0.0028344920 0.0012928910 0.0006379191 0.0085989355 0.0035028551 0.0011249625 0.0024085229 0.0047753376 0.0701153131 0.4135913903 0.0016748492 0.0744862631 0.2785384406 0.0040466582 0.0037087155 0.0052379329 0.0200222636 0.0523938808 0.0408860135;
frequency C30pi6 = 0.1335831781 0.0284789590 0.0213891629 0.1125775537 0.0010514541 0.0565844323 0.2099572968 0.0207551870 0.0121330488 0.0073526522 0.0133278240 0.0771772013 0.0030571689 0.0016793592 0.1890195131 0.0484054108 0.0373318180 0.0009266995 0.0026946425 0.0225174379;
frequency C30pi7 = 0.0408277374 0.0124491768 0.0080464869 0.0030634898 0.0153918410 0.0102922098 0.0066010880 0.0058113137 0.0245211764 0.1487514547 0.1637802160 0.0075923232 0.0385527359 0.1575049888 0.0058352224 0.0151578617 0.0332220362 0.0264937109 0.1213342989 0.1547706314;
frequency C30pi8 = 0.2469059247 0.0106278945 0.0168929681 0.0027418266 0.1039406309 0.0103988197 0.0054944756 0.0373263209 0.0085752319 0.0292403793 0.0535091180 0.0056123053 0.0302246485 0.0251775640 0.0078098946 0.1642352274 0.1239889705 0.0053155877 0.0163953993 0.0955868125;
frequency C30pi9 = 0.0549428629 0.1305426495 0.0202957532 0.0092915274 0.0099280995 0.0906036344 0.0417085054 0.0105563869 0.0363512470 0.0569584863 0.1681833183 0.1152521806 0.0592328363 0.0243860149 0.0083055411 0.0283778833 0.0412594019 0.0096355359 0.0249780472 0.0592100878;
frequency C30pi10 = 0.0462773303 0.0362984274 0.0412365193 0.0182504174 0.0172727117 0.0348990852 0.0224266258 0.0160971397 0.1357852215 0.0164966886 0.0598936127 0.0239396241 0.0164507129 0.1336320854 0.0117413009 0.0454156401 0.0304387749 0.0330338410 0.2350163763 0.0253978649;
frequency C30pi11 = 0.0474379955 0.0410179935 0.0222453982 0.0112116958 0.0082332447 0.0374051414 0.0388100853 0.0055998598 0.0149156570 0.1832173840 0.1100691114 0.0467850545 0.0356443791 0.0116643783 0.0100244663 0.0317171100 0.1114352326 0.0026685586 0.0099660086 0.2199312452;
frequency C30pi12 = 0.0213607696 0.0069976154 0.0039878996 0.0012941246 0.0061024858 0.0139566033 0.0036297282 0.0030017014 0.0038425894 0.1309465785 0.4566988203 0.0054567760 0.1947837355 0.0371808169 0.0040747282 0.0076991487 0.0198018718 0.0034086391 0.0064545692 0.0693207986;
frequency C30pi13 = 0.0919632044 0.0160004872 0.0764682386 0.0306717360 0.0117031014 0.0160060006 0.0171907654 0.0370684649 0.0100792697 0.0093123713 0.0097240970 0.0205385908 0.0075767282 0.0041589440 0.0179686194 0.3254471625 0.2744377258 0.0013887442 0.0044739725 0.0178217761;
frequency C30pi14 = 0.4649246103 0.0043013249 0.0075304815 0.0050731691 0.0233328752 0.0043571322 0.0057994247 0.1495242047 0.0023298425 0.0043361190 0.0055995530 0.0028525398 0.0039313170 0.0025588185 0.0186467246 0.2150194771 0.0477030158 0.0009038096 0.0020087184 0.0292668421;
frequency C30pi15 = 0.2051329382 0.0439661329 0.0339418395 0.1070980865 0.0020915940 0.0822742346 0.1989733497 0.0487574293 0.0127143076 0.0058124693 0.0133471767 0.0667787412 0.0043783406 0.0018235059 0.0110997761 0.0873961609 0.0519781961 0.0007361603 0.0023821404 0.0193174204;
frequency C30pi16 = 0.0263689890 0.0133613622 0.2727158135 0.3117715371 0.0039462429 0.0218978778 0.0694354212 0.0799842408 0.0309615130 0.0027521242 0.0038579661 0.0288630708 0.0018363656 0.0023351927 0.0062457560 0.0798729385 0.0324143174 0.0007229656 0.0063857732 0.0042705326;
frequency C30pi17 = 0.1526502637 0.0332784464 0.0168229991 0.0237392180 0.0040215287 0.0341733672 0.0377949108 0.0306214335 0.0141929803 0.0123317972 0.0290062362 0.0375543022 0.0064473224 0.0058584416 0.3864504800 0.0880336410 0.0489543188 0.0018252558 0.0048877798 0.0313552773;
frequency C30pi18 = 0.0080247558 0.0017408595 0.0006327403 0.0003385965 0.0023412143 0.0015507896 0.0007818945 0.0005403825 0.0010026402 0.3177056649 0.3737894172 0.0012598254 0.0488212345 0.0311968471 0.0020687549 0.0012095129 0.0065696791 0.0016309208 0.0043343553 0.1944599147;
frequency C30pi19 = 0.0599950319 0.1000540567 0.1334918892 0.0889730776 0.0016884984 0.0864856169 0.0962700957 0.0588796388 0.0327277145 0.0021467269 0.0070876372 0.1825860579 0.0033979446 0.0011800742 0.0141408084 0.0779002375 0.0448817374 0.0006249028 0.0032641120 0.0042241415;
frequency C30pi20 = 0.0393520657 0.0838170642 0.1425481600 0.0431197671 0.0099071945 0.1019786610 0.0394639510 0.0282866471 0.2095718357 0.0076101442 0.0258339558 0.0596434088 0.0084586675 0.0188680789 0.0096840517 0.0624998643 0.0347087967 0.0054645779 0.0564145251 0.0127685828;
frequency C30pi21 = 0.0072715487 0.0140998918 0.0019756795 0.0027603830 0.0067852535 0.0043339290 0.0025069369 0.0080834718 0.0113217919 0.0056609640 0.0394199644 0.0017735096 0.0079866080 0.1271475634 0.0041098092 0.0052244365 0.0043022271 0.6273570153 0.1084563767 0.0094226397;
frequency C30pi22 = 0.0907070068 0.0290062335 0.0860677696 0.0745872716 0.0063699858 0.0259377035 0.0386802115 0.4750046194 0.0168090013 0.0014721054 0.0055149849 0.0343855535 0.0024692074 0.0028859215 0.0112150781 0.0731110371 0.0153705714 0.0022914775 0.0041860660 0.0039281943;
frequency C30pi23 = 0.0055291882 0.0024626303 0.0046086594 0.0011413426 0.0072105915 0.0022692184 0.0009683043 0.0016070950 0.0325831191 0.0082918400 0.0353677882 0.0013849437 0.0074486804 0.3744093753 0.0013374573 0.0057402692 0.0037279636 0.0330334445 0.4609978298 0.0098802591;
frequency C30pi24 = 0.2443263138 0.0045386562 0.0062422652 0.0031590902 0.0273880205 0.0053593950 0.0076715636 0.0196089609 0.0020189401 0.1017435067 0.0468424225 0.0045492259 0.0201286022 0.0060619450 0.0185219126 0.0497753825 0.1170795523 0.0009577255 0.0035333687 0.3104931504;
frequency C30pi25 = 0.0863111274 0.0984811895 0.0313963115 0.0600902926 0.0024419845 0.1672351286 0.2036096150 0.0175221435 0.0245245046 0.0105994220 0.0271209781 0.1485789590 0.0095824358 0.0029393105 0.0068276769 0.0347800318 0.0408210979 0.0014001253 0.0055105388 0.0202271268;
frequency C30pi26 = 0.0643926114 0.0369048739 0.1031213278 0.1628208462 0.0023165895 0.0752534859 0.1762701353 0.0297139006 0.0303503732 0.0088163033 0.0148016812 0.0727140107 0.0056748403 0.0043066715 0.0099270322 0.0926433867 0.0833129915 0.0011237109 0.0093801464 0.0161550816;
frequency C30pi27 = 0.1736682858 0.0943628709 0.0520404980 0.0285984935 0.0083596568 0.0722446698 0.0483894060 0.0781901497 0.0266134684 0.0068641911 0.0219499324 0.0964011794 0.0112303313 0.0058273974 0.0169661076 0.1547802460 0.0751701930 0.0028774511 0.0082130397 0.0172524320;
frequency C30pi28 = 0.0347856579 0.3075984538 0.0314157384 0.0092355245 0.0062754891 0.0861073155 0.0323568406 0.0170288127 0.0306438905 0.0091932292 0.0224428556 0.3020845818 0.0093720833 0.0034303536 0.0104447169 0.0326882932 0.0328713449 0.0025244855 0.0064171317 0.0130832013;
frequency C30pi29 = 0.1087737102 0.0051781020 0.0032679768 0.0015823203 0.0247877480 0.0057932006 0.0041769888 0.0134703172 0.0024765788 0.1643462917 0.2337152707 0.0027000391 0.0539213396 0.0316523420 0.0154886946 0.0188187787 0.0474912345 0.0037656478 0.0073106362 0.2512827825;
frequency C30pi30 = 0.1101008748 0.0324324597 0.0435098681 0.0579268520 0.0072699765 0.0615196630 0.0828181488 0.0314463068 0.0308557019 0.0530865813 0.1096787834 0.0293860426 0.0458728977 0.0269153699 0.0296430687 0.0715887866 0.0685882454 0.0062324120 0.0257237601 0.0754042006;

model C30 = POISSON+G+FMIX{C30pi1:1:0.0095783264,C30pi2:1:0.0248476365,C30pi3:1:0.0636309366,C30pi4:1:0.0537939225,C30pi5:1:0.0295885587,C30pi6:1:0.0117587936,C30pi7:1:0.0132013428,C30pi8:1:0.0236868805,C30pi9:1:0.0261687659,C30pi10:1:0.0239821974,C30pi11:1:0.0257100906,C30pi12:1:0.0465072425,C30pi13:1:0.0546794546,C30pi14:1:0.0536085131,C30pi15:1:0.0270622670,C30pi16:1:0.0403913593,C30pi17:1:0.0474212700,C30pi18:1:0.0458816478,C30pi19:1:0.0214036510,C30pi20:1:0.0290385981,C30pi21:1:0.0123391793,C30pi22:1:0.0569350229,C30pi23:1:0.0419687568,C30pi24:1:0.0339027062,C30pi25:1:0.0388777376,C30pi26:1:0.0196343766,C30pi27:1:0.0233086174,C30pi28:1:0.0622722654,C30pi29:1:0.0184803385,C30pi30:1:0.0203395454};


CAT-C40 profile mixture model of Le, Gascuel & Lartillot (2008)
frequency C40pi1 = 0.0660259814 0.0231861755 0.1599815873 0.1054473175 0.0056586745 0.0273928499 0.0440360794 0.0711238664 0.0168194755 0.0039088727 0.0055316013 0.0366689617 0.0037412416 0.0013104807 0.0176359169 0.2497687201 0.1507079582 0.0006723214 0.0038290224 0.0065528958;
frequency C40pi2 = 0.0232377444 0.0122683027 0.2759650991 0.3532087982 0.0037987468 0.0197339134 0.0739378219 0.0576668030 0.0315866952 0.0031092806 0.0038711609 0.0259363304 0.0017355634 0.0024032103 0.0063116881 0.0657067704 0.0270483653 0.0007602894 0.0069602476 0.0047531689;
frequency C40pi3 = 0.0166486809 0.0012594763 0.0012622242 0.0005651446 0.0036665719 0.0010669784 0.0013356251 0.0008894749 0.0008231853 0.4129367561 0.0884689295 0.0011904105 0.0186054583 0.0082775676 0.0014029981 0.0021339439 0.0162167380 0.0006082049 0.0019553200 0.4206863114;
frequency C40pi4 = 0.2394741986 0.0072901253 0.0120536943 0.0044741726 0.0283811727 0.0086558850 0.0105529632 0.0135109628 0.0038929844 0.0765957115 0.0358494908 0.0071093014 0.0199496319 0.0055991131 0.0114265585 0.0847798773 0.1797284519 0.0009838000 0.0042240671 0.2454678377;
frequency C40pi5 = 0.1194613086 0.0233255669 0.0294552140 0.0134272792 0.0150526644 0.0301537796 0.0192173037 0.0337675998 0.0214746045 0.0579001821 0.1446308373 0.0147261337 0.0561242940 0.0550467421 0.0631355418 0.0925266727 0.0831230185 0.0131636136 0.0331118002 0.0811758434;
frequency C40pi6 = 0.0567043710 0.0117359330 0.0364734454 0.2955500969 0.0008924801 0.0609516515 0.3795154126 0.0230469606 0.0118360971 0.0031182036 0.0060137466 0.0314205689 0.0028584065 0.0012972333 0.0124745819 0.0300334889 0.0227051137 0.0007738758 0.0031343761 0.0094639563;
frequency C40pi7 = 0.0179027412 0.0040967133 0.0035697688 0.0008870412 0.0160760340 0.0045395474 0.0023182113 0.0039829808 0.0127292680 0.0404650518 0.1676143477 0.0027994718 0.0424172255 0.3344862590 0.0020115128 0.0075841581 0.0068227293 0.0518381385 0.2452542553 0.0326045442;
frequency C40pi8 = 0.2712170094 0.0056480837 0.0141045260 0.0021017036 0.2003830179 0.0048264059 0.0023229984 0.0502501222 0.0053727960 0.0150684657 0.0330003443 0.0020646283 0.0154811217 0.0202990358 0.0045351023 0.1764198412 0.0839578061 0.0046265242 0.0141271048 0.0741933626;
frequency C40pi9 = 0.0894736584 0.1040026384 0.0190192153 0.0272183085 0.0045538316 0.1168091917 0.1275076663 0.0115685734 0.0215746293 0.0469424171 0.0512035100 0.1382047308 0.0147656854 0.0056590176 0.0095546504 0.0383953611 0.0836652641 0.0017079427 0.0062181292 0.0819555787;
frequency C40pi10 = 0.0495441385 0.0375345822 0.0315863530 0.0143641284 0.0182505609 0.0316504100 0.0215379122 0.0140199913 0.1108543799 0.0247065801 0.0700287927 0.0258142032 0.0188271760 0.1418048822 0.0112101202 0.0456094427 0.0361427973 0.0371985427 0.2223972375 0.0369177689;
frequency C40pi11 = 0.1704314254 0.0415784004 0.0271109259 0.1098556600 0.0009747331 0.0917299929 0.2536458944 0.0249846466 0.0101389736 0.0058749399 0.0116526350 0.0903324267 0.0036512738 0.0013321301 0.0293613681 0.0561765645 0.0479045729 0.0006696817 0.0022637316 0.0203300232;
frequency C40pi12 = 0.0162725399 0.0054826071 0.0021876158 0.0010182101 0.0050614097 0.0104414465 0.0025141347 0.0021935389 0.0029914328 0.1328173512 0.4904441779 0.0040120394 0.1929931280 0.0376245580 0.0034333187 0.0040122105 0.0127074428 0.0032107554 0.0058100621 0.0647720205;
frequency C40pi13 = 0.0823765743 0.0734226431 0.0598389731 0.0311745159 0.0065694304 0.0686451074 0.0675530778 0.0178961594 0.0251143622 0.0291161743 0.0287904106 0.0982301674 0.0168022878 0.0064717899 0.0114044922 0.1302995288 0.1820374273 0.0022724618 0.0079573279 0.0540270885;
frequency C40pi14 = 0.3594965940 0.0072407229 0.0033421456 0.0031484357 0.0251417178 0.0049014279 0.0064962700 0.1194682267 0.0022970448 0.0458766662 0.0468053893 0.0050168849 0.0215568816 0.0092020461 0.0443915884 0.0465270945 0.0477755293 0.0024540215 0.0046450361 0.1942162766;
frequency C40pi15 = 0.2015583874 0.0430161610 0.0425386444 0.0954149893 0.0032365302 0.0772010857 0.1534908791 0.0667291678 0.0155218808 0.0067740832 0.0165114429 0.0547322644 0.0060162992 0.0025643300 0.0091970560 0.1185981804 0.0625472744 0.0009565508 0.0031150007 0.0202797924;
frequency C40pi16 = 0.1042731047 0.0147062345 0.0621645800 0.2424069523 0.0022450116 0.0356498946 0.1774821588 0.1697819523 0.0132648834 0.0018929517 0.0042542620 0.0220651981 0.0016441234 0.0012570256 0.0317041583 0.0778636230 0.0288515782 0.0006930898 0.0017741945 0.0060250231;
frequency C40pi17 = 0.0781183281 0.0111498472 0.0159270309 0.0041541669 0.0194448667 0.0240151620 0.0116633921 0.0111524105 0.0063589385 0.1354530457 0.2457574952 0.0093729846 0.1087781166 0.0262793949 0.0055294038 0.0408518858 0.0860514305 0.0031547586 0.0085108496 0.1482764918;
frequency C40pi18 = 0.0856592432 0.0101233167 0.0441923073 0.0135061568 0.0136072878 0.0092590642 0.0078602552 0.0245400880 0.0055379075 0.0100591561 0.0103343559 0.0127318506 0.0080675803 0.0047153035 0.0175273997 0.3406479487 0.3573294650 0.0014243098 0.0035099810 0.0193670227;
frequency C40pi19 = 0.0674594695 0.1161734658 0.1163107783 0.0662588409 0.0021634231 0.0939360452 0.0865501280 0.0368556575 0.0381149118 0.0033238825 0.0093839985 0.1899736999 0.0039487389 0.0018212730 0.0151207830 0.0842204423 0.0565953680 0.0007187305 0.0046189437 0.0064514195;
frequency C40pi20 = 0.0572262322 0.0494723554 0.1083882793 0.1793932771 0.0015301521 0.0903668522 0.1992261265 0.0316472274 0.0291392067 0.0045804559 0.0100739563 0.1015624916 0.0040204606 0.0013701849 0.0063674130 0.0621142922 0.0496102162 0.0006669285 0.0046497641 0.0085941279;
frequency C40pi21 = 0.0036020163 0.0102712927 0.0013455508 0.0020871647 0.0045484804 0.0032718114 0.0017857730 0.0056391633 0.0064968790 0.0029292916 0.0232635081 0.0010419846 0.0044592278 0.0855714596 0.0024991984 0.0030671803 0.0025900250 0.7617821954 0.0678809532 0.0058668443;
frequency C40pi22 = 0.2032018418 0.0083895722 0.0143743754 0.0135011707 0.0098131618 0.0044514580 0.0083818173 0.6184886075 0.0027747899 0.0011828492 0.0039826789 0.0044598895 0.0020631785 0.0019619615 0.0085870399 0.0739919851 0.0108922273 0.0018606145 0.0015638674 0.0060769136;
frequency C40pi23 = 0.0050898779 0.0028740788 0.0057092962 0.0016126151 0.0061776450 0.0024693148 0.0012040415 0.0016334183 0.0393460780 0.0059088776 0.0249343597 0.0013713662 0.0049795162 0.3563126947 0.0014136424 0.0059527667 0.0036536770 0.0357987380 0.4853645852 0.0081934106;
frequency C40pi24 = 0.0403335679 0.0540186397 0.0216052457 0.0098218598 0.0081549541 0.0383639077 0.0375406578 0.0047934404 0.0176735565 0.1893424159 0.1051859862 0.0607377395 0.0305599836 0.0119140782 0.0077550551 0.0257110173 0.1009913165 0.0028780020 0.0115276935 0.2210908828;
frequency C40pi25 = 0.0790086293 0.1065441152 0.0309384274 0.0546012394 0.0024947877 0.1843375981 0.1997882784 0.0192655847 0.0270700474 0.0075667489 0.0254542392 0.1553108816 0.0098024439 0.0023773444 0.0056640684 0.0332370813 0.0359574739 0.0011682801 0.0048820809 0.0145306498;
frequency C40pi26 = 0.0722240672 0.0489728405 0.0678929607 0.1194883992 0.0064755348 0.0708969573 0.1345886574 0.0287815397 0.0699011334 0.0173588702 0.0519870084 0.0490341790 0.0154411043 0.0348233029 0.0145597486 0.0589579876 0.0425972780 0.0087913770 0.0554386705 0.0317883834;
frequency C40pi27 = 0.1085842431 0.0206450023 0.0441956285 0.1529666596 0.0012502570 0.0405398136 0.1664851192 0.0336098469 0.0134902179 0.0038821795 0.0089861440 0.0576227094 0.0024339036 0.0014553522 0.1990095021 0.0846749753 0.0454715217 0.0005902831 0.0027650162 0.0113416246;
frequency C40pi28 = 0.0309526387 0.3195887318 0.0301336637 0.0082352132 0.0065593963 0.0832608108 0.0291974083 0.0154206187 0.0310385092 0.0098251607 0.0237900204 0.3062634996 0.0097071728 0.0036891639 0.0095029109 0.0295285439 0.0303052301 0.0028125285 0.0068850639 0.0133037148;
frequency C40pi29 = 0.0098953741 0.0019604525 0.0007307935 0.0003748228 0.0028276741 0.0017337004 0.0009182100 0.0006997068 0.0010419482 0.3115040359 0.3750387796 0.0013960508 0.0474451070 0.0298607430 0.0025296256 0.0014628019 0.0075738968 0.0016799771 0.0040259930 0.1973003069;
frequency C40pi30 = 0.1163213921 0.0273321006 0.0250163656 0.0731917718 0.0034792282 0.0586677248 0.1380880502 0.0193193469 0.0160240740 0.0712243431 0.0771473538 0.0355120487 0.0242841072 0.0094117688 0.0508926833 0.0475560280 0.0726552233 0.0026892716 0.0076166020 0.1235705162;
frequency C40pi31 = 0.1285218235 0.0373073487 0.1179844215 0.0402749992 0.0172928883 0.0439706110 0.0250692272 0.1127033137 0.0606981059 0.0109350265 0.0258415767 0.0288749652 0.0167592956 0.0199118302 0.0180674983 0.1741489481 0.0648967655 0.0063574951 0.0321771650 0.0182066946;
frequency C40pi32 = 0.0372286941 0.0094528028 0.0053377315 0.0023703173 0.0144940088 0.0079097138 0.0048585146 0.0046433943 0.0186795102 0.1820459527 0.1780099317 0.0058198481 0.0371334296 0.1463772419 0.0048538601 0.0103570678 0.0284161577 0.0211293603 0.0958905187 0.1849919442;
frequency C40pi33 = 0.0535643726 0.1159797757 0.0239172676 0.0113537364 0.0096256227 0.0928585070 0.0391699080 0.0120279334 0.0384887950 0.0522748270 0.1892392595 0.0996037748 0.0712219098 0.0264213736 0.0083720574 0.0299114019 0.0389484845 0.0104232046 0.0265030050 0.0500947835;
frequency C40pi34 = 0.1332424803 0.0033147683 0.0022704992 0.0012739239 0.0246514263 0.0030843469 0.0040461524 0.0089139209 0.0015864680 0.1971284995 0.1251288442 0.0023713225 0.0286947200 0.0156995251 0.0118845743 0.0171461828 0.0563298009 0.0017341820 0.0048778410 0.3566205216;
frequency C40pi35 = 0.1498658185 0.0326607222 0.0176452820 0.0280354786 0.0035437399 0.0348151308 0.0435380704 0.0311112643 0.0140625707 0.0101953314 0.0251433928 0.0393124980 0.0051548319 0.0047533945 0.3923800449 0.0874496981 0.0473306717 0.0015215239 0.0043208299 0.0271597054;
frequency C40pi36 = 0.4214366359 0.0061425967 0.0121590498 0.0073305074 0.0187609694 0.0072748556 0.0086837775 0.0902333103 0.0030262044 0.0039362777 0.0047193320 0.0051508681 0.0038306586 0.0027156136 0.0208940236 0.2901188793 0.0651922314 0.0008108235 0.0023622848 0.0252211004;
frequency C40pi37 = 0.1770713890 0.1332782050 0.0311656783 0.0226500225 0.0078348946 0.0752471493 0.0509767242 0.0897389513 0.0220667143 0.0059519850 0.0205369728 0.1257689326 0.0092982479 0.0040514178 0.0264087912 0.1169591448 0.0565566955 0.0029947127 0.0049346701 0.0165087010;
frequency C40pi38 = 0.0293984032 0.0370901720 0.1483622633 0.1099709900 0.0031729093 0.0388688450 0.0464270335 0.4222420155 0.0272494642 0.0007997326 0.0037634298 0.0622314461 0.0016657052 0.0015039626 0.0056481827 0.0472252404 0.0086568982 0.0009176022 0.0027693124 0.0020363920;
frequency C40pi39 = 0.0265779317 0.0791104753 0.1318603134 0.0280314140 0.0101369144 0.0989710810 0.0269057233 0.0173376629 0.2815133703 0.0064646977 0.0268210053 0.0474749135 0.0072375268 0.0276960902 0.0083014995 0.0426276702 0.0259042511 0.0078528946 0.0891598394 0.0100147256;
frequency C40pi40 = 0.0096096503 0.0027136180 0.0013104432 0.0006331856 0.0077301682 0.0033899420 0.0010471898 0.0020227436 0.0039001415 0.0733098005 0.4451691588 0.0014931484 0.0732575295 0.2630171690 0.0042768091 0.0036117358 0.0057928403 0.0181275729 0.0370698053 0.0425173480;

model C40 = POISSON+G+FMIX{C40pi1:1:0.0223853788,C40pi2:1:0.0338891820,C40pi3:1:0.0577169375,C40pi4:1:0.0252416233,C40pi5:1:0.0108607921,C40pi6:1:0.0462373793,C40pi7:1:0.0102293175,C40pi8:1:0.0147523625,C40pi9:1:0.0143161352,C40pi10:1:0.0182302541,C40pi11:1:0.0204025079,C40pi12:1:0.0425505156,C40pi13:1:0.0248627269,C40pi14:1:0.0105892988,C40pi15:1:0.0188238725,C40pi16:1:0.0086663445,C40pi17:1:0.0148496147,C40pi18:1:0.0343037402,C40pi19:1:0.0225335203,C40pi20:1:0.0174068578,C40pi21:1:0.0112207827,C40pi22:1:0.0443532245,C40pi23:1:0.0392573370,C40pi24:1:0.0196756555,C40pi25:1:0.0287690328,C40pi26:1:0.0114441177,C40pi27:1:0.0112338740,C40pi28:1:0.0582694099,C40pi29:1:0.0444272279,C40pi30:1:0.0112010942,C40pi31:1:0.0145176111,C40pi32:1:0.0114629026,C40pi33:1:0.0239628061,C40pi34:1:0.0266266492,C40pi35:1:0.0481201159,C40pi36:1:0.0371147423,C40pi37:1:0.0160476688,C40pi38:1:0.0237249267,C40pi39:1:0.0235226203,C40pi40:1:0.0261998398};

CAT-C50 profile mixture model of Le, Gascuel & Lartillot (2008)
frequency C50pi1 = 0.1357566757 0.0328511938 0.0937692919 0.0757182069 0.0041887049 0.0448010470 0.0572805366 0.1210866186 0.0167465028 0.0049719235 0.0113823284 0.0458096069 0.0064563157 0.0029292810 0.0228705187 0.2060115780 0.1011347978 0.0012443033 0.0056104605 0.0093801079;
frequency C50pi2 = 0.0530862751 0.1905936010 0.0595772279 0.0320970468 0.0026608079 0.1152605895 0.0840617877 0.0196495178 0.0274729775 0.0064919200 0.0158709120 0.2635539775 0.0078171228 0.0017231166 0.0121639300 0.0449347664 0.0472425608 0.0008407188 0.0037608716 0.0111402722;
frequency C50pi3 = 0.0083279799 0.0007172026 0.0006359642 0.0003134388 0.0020547407 0.0007351595 0.0005373710 0.0005576905 0.0004858721 0.4370910601 0.1208722220 0.0006394909 0.0195499664 0.0090175268 0.0007265254 0.0007876194 0.0057076665 0.0006453449 0.0016797264 0.3889174318;
frequency C50pi4 = 0.2072868350 0.0166858699 0.0129177658 0.0020625574 0.0849982226 0.0151757635 0.0065903656 0.0472047575 0.0130289256 0.0345690755 0.1042722764 0.0075861385 0.0498042308 0.0572909747 0.0064928361 0.1183618036 0.0780339514 0.0128352368 0.0323576924 0.0924447209;
frequency C50pi5 = 0.0364181183 0.0076427099 0.0052725527 0.0020389950 0.0171009943 0.0064088232 0.0042399368 0.0053824238 0.0198596156 0.1361523026 0.1651892915 0.0045481616 0.0387479055 0.2025922657 0.0055053348 0.0121111950 0.0254621828 0.0327580458 0.1368025306 0.1357666147;
frequency C50pi6 = 0.0535489196 0.0099543365 0.0269073208 0.3076150732 0.0007101021 0.0574988641 0.4066173371 0.0204537673 0.0096286483 0.0025879708 0.0049721459 0.0280989086 0.0025143457 0.0010618006 0.0124317994 0.0247246015 0.0191107367 0.0006385967 0.0024132214 0.0085115039;
frequency C50pi7 = 0.0074733729 0.0025226602 0.0033967505 0.0005574007 0.0081158286 0.0037658904 0.0013610444 0.0022017759 0.0115142679 0.0195730439 0.1268878488 0.0018497296 0.0269141680 0.3821985941 0.0019970421 0.0057127939 0.0039692337 0.0553575998 0.3184099394 0.0162210153;
frequency C50pi8 = 0.2615592974 0.0027098854 0.0124908261 0.0020153852 0.2740228527 0.0017043893 0.0007667803 0.0463498030 0.0019474361 0.0082858275 0.0147048711 0.0010787235 0.0063051368 0.0062080862 0.0039442437 0.1940042648 0.0963699489 0.0016185483 0.0048431386 0.0590705550;
frequency C50pi9 = 0.1190557043 0.0956320251 0.0215995297 0.0378323341 0.0041536088 0.1151348174 0.1337084452 0.0179375220 0.0216767047 0.0336228770 0.0557402194 0.1132452331 0.0178407325 0.0063405927 0.0147606946 0.0478666925 0.0712091035 0.0022867238 0.0075728630 0.0627835766;
frequency C50pi10 = 0.0505010344 0.0281381134 0.0341872191 0.0178157543 0.0183140005 0.0271729546 0.0212018661 0.0176052654 0.1190104107 0.0161645217 0.0561232531 0.0203908848 0.0146521042 0.1553484132 0.0135251600 0.0478959652 0.0292963208 0.0376058633 0.2477283800 0.0273225153;
frequency C50pi11 = 0.1239446910 0.0355525870 0.0409769096 0.1479953346 0.0011563976 0.0908869312 0.2700270273 0.0283589709 0.0126760201 0.0064825033 0.0122101302 0.0787433823 0.0042467440 0.0016540857 0.0205717500 0.0552940245 0.0474239965 0.0008596621 0.0027823209 0.0181565313;
frequency C50pi12 = 0.0160542063 0.0027359185 0.0014708079 0.0007004900 0.0034820152 0.0061470051 0.0016359686 0.0022137927 0.0013207229 0.1640035117 0.4616043506 0.0021342205 0.2174099502 0.0143751693 0.0013694259 0.0037614383 0.0172651408 0.0011454338 0.0019438536 0.0792265779;
frequency C50pi13 = 0.1548192401 0.0131324559 0.0280584102 0.0095301620 0.0166267416 0.0175228950 0.0170969133 0.0179616718 0.0078385586 0.0865181208 0.0523369910 0.0132802182 0.0326348210 0.0083511229 0.0145594414 0.1096327081 0.2218108602 0.0015829972 0.0062173360 0.1704883347;
frequency C50pi14 = 0.2950313592 0.0027580697 0.0021616268 0.0015364190 0.0375439186 0.0028808733 0.0042976283 0.0261726702 0.0008294969 0.0834938143 0.0553606311 0.0022642314 0.0181259911 0.0074433078 0.0126794048 0.0382913338 0.0783205173 0.0010015148 0.0034016419 0.3264055498;
frequency C50pi15 = 0.1683177099 0.0820396152 0.0526048706 0.0822517150 0.0023029997 0.0969341246 0.1488943001 0.0535291188 0.0179803231 0.0032503636 0.0114941086 0.1156402642 0.0039439899 0.0015002945 0.0066854154 0.0924511658 0.0480769504 0.0006152103 0.0025022919 0.0089851683;
frequency C50pi16 = 0.0334088176 0.0134485791 0.1590918150 0.3657542471 0.0025127086 0.0327665151 0.1820739351 0.0740807194 0.0202010901 0.0016650025 0.0036700956 0.0295517886 0.0017087810 0.0011422805 0.0073155123 0.0426788071 0.0211162106 0.0005931485 0.0034724580 0.0037474882;
frequency C50pi17 = 0.0777586977 0.0174438357 0.0053423343 0.0043431532 0.0062523949 0.0220851281 0.0161769285 0.0053903202 0.0080675581 0.1052945216 0.1617365895 0.0148319919 0.0288253912 0.0168985297 0.2565426868 0.0202089662 0.0542929694 0.0060146095 0.0078109966 0.1646823969;
frequency C50pi18 = 0.0727013979 0.0048977192 0.0026095383 0.0011420120 0.0198747408 0.0066949336 0.0030401434 0.0079074845 0.0026492900 0.1685788878 0.3185489163 0.0026024909 0.0735597038 0.0490419983 0.0051699104 0.0128630830 0.0305356924 0.0050857840 0.0095279173 0.2029683559;
frequency C50pi19 = 0.0658153836 0.0833432992 0.0224582275 0.0107735824 0.0092974677 0.0745951987 0.0299754097 0.0146336557 0.0148026634 0.0671888719 0.2198675990 0.0868172087 0.1084156835 0.0155812696 0.0071132147 0.0381451947 0.0562948237 0.0056421684 0.0102813038 0.0589577740;
frequency C50pi20 = 0.0525278351 0.0364897390 0.0903013988 0.1854660991 0.0037795400 0.0776857292 0.1789287290 0.0232011648 0.0687702011 0.0135825419 0.0337350646 0.0458143770 0.0108457797 0.0191020037 0.0088729983 0.0495289201 0.0389358438 0.0046292762 0.0354195947 0.0223831639;
frequency C50pi21 = 0.0026515970 0.0080885204 0.0010572021 0.0016052142 0.0036540307 0.0022979498 0.0014681767 0.0046230912 0.0043887616 0.0020669456 0.0172444871 0.0006593575 0.0034691503 0.0658351447 0.0019185467 0.0022498420 0.0021278866 0.8183345006 0.0515918357 0.0046677595;
frequency C50pi22 = 0.0548133174 0.0692044159 0.0211265710 0.0207779125 0.0072646572 0.0567865657 0.0738456579 0.0051797705 0.0168408457 0.1386104888 0.0713795154 0.0896393340 0.0201205491 0.0082150393 0.0104049016 0.0282344422 0.0995597110 0.0019722093 0.0074054035 0.1986186919;
frequency C50pi23 = 0.0047955268 0.0028033787 0.0050506238 0.0014080516 0.0061671241 0.0019350126 0.0009861551 0.0014396818 0.0389623239 0.0048950388 0.0151748150 0.0012306644 0.0032520404 0.3601993060 0.0011266316 0.0054509935 0.0034763921 0.0362899931 0.4980200998 0.0073361467;
frequency C50pi24 = 0.0365462996 0.0280070630 0.0183606115 0.0070525803 0.0093251684 0.0300239431 0.0221812842 0.0047778642 0.0178840316 0.2025947306 0.1973012130 0.0250209750 0.0557862640 0.0258067541 0.0042772210 0.0209374223 0.0731398943 0.0049738166 0.0200601168 0.1959427463;
frequency C50pi25 = 0.0684197684 0.0111619750 0.0544764241 0.0224313301 0.0106958312 0.0091799953 0.0097436799 0.0255871619 0.0055558006 0.0059416697 0.0076746853 0.0144198991 0.0056892166 0.0037356845 0.0172554137 0.3527301149 0.3586913194 0.0012501907 0.0028636710 0.0124961682;
frequency C50pi26 = 0.0495330775 0.1060064564 0.1511923969 0.0483471288 0.0080946362 0.0886108407 0.0449556763 0.0331436148 0.1447288287 0.0061850770 0.0190407203 0.0948075276 0.0063418871 0.0126162987 0.0100869563 0.0799801169 0.0445418973 0.0044765096 0.0363930724 0.0109172804;
frequency C50pi27 = 0.0702411901 0.0642050323 0.0779553908 0.0510328304 0.0042438849 0.0723300485 0.0883747710 0.0177347101 0.0233800891 0.0198779320 0.0183537117 0.1051267065 0.0107865869 0.0037987118 0.0112811107 0.1345081583 0.1805543234 0.0014252764 0.0055089381 0.0392805971;
frequency C50pi28 = 0.1207399152 0.1741788075 0.0385528120 0.0162689581 0.0118494185 0.0760068404 0.0337935391 0.0653431008 0.0342783806 0.0085426053 0.0256788075 0.1434443984 0.0112347894 0.0061270793 0.0294493558 0.1091415488 0.0634181251 0.0046156419 0.0085374279 0.0187984481;
frequency C50pi29 = 0.0064521696 0.0021817337 0.0005939658 0.0003904032 0.0021538307 0.0019099968 0.0008007758 0.0005208471 0.0011374294 0.2850758996 0.4278536740 0.0013920239 0.0561988528 0.0449501501 0.0026289702 0.0011053664 0.0055157148 0.0022753671 0.0059612583 0.1509015707;
frequency C50pi30 = 0.0969092741 0.0359723370 0.0633194168 0.0411020773 0.0145578946 0.0466661704 0.0469223767 0.0374614202 0.0537149580 0.0394603009 0.0856256544 0.0283577862 0.0346435320 0.0507298072 0.0167177549 0.0990945318 0.0806503833 0.0128373826 0.0598972198 0.0553597218;
frequency C50pi31 = 0.0840212010 0.0214242172 0.2240668646 0.0354684798 0.0265031681 0.0235675678 0.0076026464 0.1173325117 0.0516019781 0.0048917455 0.0067211727 0.0173653354 0.0079342101 0.0087501486 0.0093276105 0.2637097946 0.0630157977 0.0022314593 0.0170994247 0.0073646661;
frequency C50pi32 = 0.0055061507 0.0012508737 0.0004824961 0.0004530173 0.0054435931 0.0011315076 0.0004150379 0.0012285001 0.0019884532 0.0617431901 0.4342418135 0.0008161868 0.0554628445 0.3289659386 0.0025814794 0.0021197505 0.0029510440 0.0172981374 0.0412097497 0.0347102358;
frequency C50pi33 = 0.0442014612 0.1295816316 0.0258622052 0.0148900471 0.0076165815 0.1301765579 0.0636708052 0.0105339122 0.0662542863 0.0423977240 0.1434197528 0.1040381429 0.0403363621 0.0260540342 0.0089335090 0.0242573966 0.0317938092 0.0077831996 0.0309973779 0.0472012033;
frequency C50pi34 = 0.0571984155 0.0034929878 0.0031324721 0.0012472712 0.0113230439 0.0025279922 0.0040737817 0.0030647398 0.0020494153 0.3131200932 0.0901750144 0.0034699557 0.0242565205 0.0112345295 0.0048197020 0.0095675953 0.0529842025 0.0010645104 0.0041851135 0.3970126433;
frequency C50pi35 = 0.1141963934 0.0102229903 0.0178644126 0.0172307307 0.0056978908 0.0039055039 0.0085974326 0.7425714921 0.0026414175 0.0005602022 0.0019872568 0.0055400059 0.0004739977 0.0010663175 0.0054302447 0.0508318204 0.0055408544 0.0018890811 0.0012409205 0.0025110348;
frequency C50pi36 = 0.3531758625 0.0043402857 0.0031812423 0.0030024877 0.0165711581 0.0029126214 0.0042077690 0.4520896100 0.0021366362 0.0063692579 0.0120143269 0.0022586970 0.0080260130 0.0043865828 0.0111462027 0.0658344033 0.0182952730 0.0010872878 0.0023330172 0.0266312657;
frequency C50pi37 = 0.0310798708 0.0234519814 0.1273669012 0.1197925100 0.0031216960 0.0295858842 0.0470763446 0.4883046368 0.0193412101 0.0008855622 0.0032808220 0.0408430573 0.0014984226 0.0016298596 0.0063229464 0.0423452622 0.0082797260 0.0007718998 0.0024996877 0.0025217188;
frequency C50pi38 = 0.0370340667 0.0689410214 0.1704407181 0.1041817082 0.0018108784 0.0715495095 0.0659866718 0.2159298358 0.0443591808 0.0008668888 0.0064679416 0.1275300877 0.0027248464 0.0014178323 0.0060253154 0.0534574556 0.0147073432 0.0007999410 0.0037708147 0.0019979426;
frequency C50pi39 = 0.0160398536 0.0526622999 0.1051167149 0.0187352256 0.0085330116 0.0922616498 0.0154450839 0.0076235155 0.3848449137 0.0057129406 0.0277195224 0.0219347380 0.0071078308 0.0376358992 0.0072201969 0.0209969653 0.0142198783 0.0096946226 0.1384243143 0.0080708232;
frequency C50pi40 = 0.0165549167 0.0085856833 0.0049441851 0.0016567380 0.0086529073 0.0184087838 0.0033759867 0.0033844413 0.0084695063 0.0483923758 0.4963073963 0.0056997331 0.1949377866 0.0999527140 0.0060271256 0.0084289585 0.0122619536 0.0114013282 0.0192314834 0.0233259964;
frequency C50pi41 = 0.0227379959 0.0137060298 0.3162561805 0.2932103363 0.0037073869 0.0169119273 0.0380984220 0.0550224760 0.0319886436 0.0039219190 0.0041582288 0.0312539900 0.0019467591 0.0022276545 0.0059660826 0.0998736999 0.0462336456 0.0007310446 0.0069012376 0.0051463400;
frequency C50pi42 = 0.2406936002 0.0197081082 0.0462578641 0.0206379264 0.0186726798 0.0189843646 0.0129785315 0.1749109142 0.0118714342 0.0049349532 0.0126237761 0.0127876711 0.0095642661 0.0083606873 0.0326283314 0.2101300187 0.1130042042 0.0041951500 0.0069210515 0.0201344675;
frequency C50pi43 = 0.0214325714 0.3730744306 0.0220674626 0.0037495290 0.0069038342 0.0670391950 0.0159298773 0.0126211348 0.0284477629 0.0102051798 0.0242954287 0.3272456489 0.0093147452 0.0036403029 0.0070138928 0.0216860624 0.0232259733 0.0030422478 0.0065368590 0.0125278613;
frequency C50pi44 = 0.1567707052 0.0258059606 0.0161658338 0.0223946414 0.0074382689 0.0274455582 0.0410010574 0.0360501033 0.0159972680 0.0640941463 0.0944756654 0.0192586366 0.0312789234 0.0227728534 0.1653169011 0.0640177954 0.0549103568 0.0050980224 0.0138248643 0.1158824381;
frequency C50pi45 = 0.4345912387 0.0061142999 0.0097660767 0.0060102195 0.0197377879 0.0069062805 0.0082800652 0.0829075516 0.0029125126 0.0047747098 0.0054182241 0.0049974525 0.0039676868 0.0029052002 0.0193588692 0.2795854727 0.0677816788 0.0008196092 0.0025196339 0.0306454302;
frequency C50pi46 = 0.0296734965 0.1443250343 0.0128668160 0.0059561454 0.0129805897 0.0492311054 0.0262726056 0.0069437743 0.0676183913 0.0452364160 0.1374511139 0.0907089722 0.0308070846 0.0816441785 0.0060701025 0.0197130339 0.0299715868 0.0461468661 0.1119414237 0.0444412635;
frequency C50pi47 = 0.1089911217 0.0159187676 0.0643054232 0.2086425054 0.0016540963 0.0375565797 0.1791004993 0.0610564917 0.0144660242 0.0038322948 0.0067778708 0.0372270242 0.0022817918 0.0012634818 0.0851792013 0.1065821239 0.0524401536 0.0005901255 0.0027836060 0.0093508169;
frequency C50pi48 = 0.1429463629 0.0304191716 0.0191145368 0.0351867799 0.0031493079 0.0341248336 0.0508492526 0.0305914291 0.0134276644 0.0070227247 0.0197257013 0.0421442438 0.0038904796 0.0040697467 0.4052202085 0.0874406009 0.0445304918 0.0012842531 0.0039485525 0.0209136585;
frequency C50pi49 = 0.0580116857 0.0903213669 0.0369245281 0.0613603988 0.0022829951 0.2073851382 0.2225853236 0.0159476910 0.0311816018 0.0068543753 0.0217092509 0.1504781849 0.0084841006 0.0020581132 0.0046206107 0.0276754451 0.0321477211 0.0011651089 0.0051889637 0.0136173964;
frequency C50pi50 = 0.2153540940 0.0359173007 0.0219927944 0.0735128474 0.0037017294 0.0566408566 0.1350375818 0.0662986417 0.0157121780 0.0138456188 0.0266922211 0.0474338339 0.0088042600 0.0035035311 0.0739583083 0.0921989198 0.0575687235 0.0019306896 0.0044520833 0.0454437865;

model C50 = POISSON+G+FMIX{C50pi1:1:0.0164297003,C50pi2:1:0.0273175755,C50pi3:1:0.0460247610,C50pi4:1:0.0084864734,C50pi5:1:0.0125389252,C50pi6:1:0.0343549036,C50pi7:1:0.0130241102,C50pi8:1:0.0094755681,C50pi9:1:0.0190040551,C50pi10:1:0.0151902354,C50pi11:1:0.0320534760,C50pi12:1:0.0210059850,C50pi13:1:0.0237408547,C50pi14:1:0.0239841203,C50pi15:1:0.0213748021,C50pi16:1:0.0210717705,C50pi17:1:0.0050241805,C50pi18:1:0.0166262276,C50pi19:1:0.0143945956,C50pi20:1:0.0104391130,C50pi21:1:0.0107628277,C50pi22:1:0.0148818171,C50pi23:1:0.0321480239,C50pi24:1:0.0145477978,C50pi25:1:0.0332355807,C50pi26:1:0.0143190281,C50pi27:1:0.0234478734,C50pi28:1:0.0183044983,C50pi29:1:0.0403269452,C50pi30:1:0.0135629530,C50pi31:1:0.0091880799,C50pi32:1:0.0158270022,C50pi33:1:0.0121019379,C50pi34:1:0.0353560982,C50pi35:1:0.0404495617,C50pi36:1:0.0104569232,C50pi37:1:0.0146187792,C50pi38:1:0.0093984095,C50pi39:1:0.0146773809,C50pi40:1:0.0201635562,C50pi41:1:0.0255640273,C50pi42:1:0.0039486842,C50pi43:1:0.0393652608,C50pi44:1:0.0056415419,C50pi45:1:0.0382833580,C50pi46:1:0.0039735086,C50pi47:1:0.0140269355,C50pi48:1:0.0476703673,C50pi49:1:0.0204062788,C50pi50:1:0.0117835304};

CAT-C60 profile mixture model of Le, Gascuel & Lartillot (2008)
frequency C60pi1 = 0.1534363248 0.0444389067 0.0796726990 0.0546757288 0.0047306596 0.0514333025 0.0529324359 0.1103775749 0.0174480218 0.0050343887 0.0130294160 0.0603928711 0.0075550589 0.0035554315 0.0249523704 0.2029625968 0.0957668473 0.0014444483 0.0059800307 0.0101808864;
frequency C60pi2 = 0.0281984692 0.3031055487 0.0312954609 0.0091549350 0.0019503463 0.0939884393 0.0388530140 0.0084028325 0.0155384715 0.0107872879 0.0217786594 0.3476042929 0.0109904917 0.0015919288 0.0071539896 0.0197479052 0.0328352333 0.0009209994 0.0025714024 0.0135302919;
frequency C60pi3 = 0.0083680740 0.0007319768 0.0006123446 0.0002228366 0.0020433870 0.0009498685 0.0004731544 0.0004825748 0.0005189995 0.3768453098 0.2608334606 0.0006296168 0.0315700586 0.0123984358 0.0009595916 0.0009746383 0.0049990761 0.0008657759 0.0017132332 0.2938075872;
frequency C60pi4 = 0.2227229348 0.0064846074 0.0061206496 0.0007997588 0.1640285908 0.0051051888 0.0027280806 0.0202702520 0.0037183875 0.0455406072 0.0883350071 0.0022832871 0.0348094559 0.0228667054 0.0035471579 0.0850040072 0.1012848285 0.0048424833 0.0096500033 0.1698580069;
frequency C60pi5 = 0.0412139519 0.0067627055 0.0051067690 0.0017434391 0.0204715649 0.0057538477 0.0037263409 0.0069107492 0.0180293946 0.1154281623 0.1693562458 0.0042900270 0.0414066566 0.2239001858 0.0058416410 0.0149106129 0.0239548406 0.0332237129 0.1379349474 0.1200342049;
frequency C60pi6 = 0.0480550249 0.0308438053 0.0940628721 0.2084606133 0.0037801787 0.0747676701 0.1855184661 0.0191402239 0.0872162350 0.0094685435 0.0277340828 0.0375741243 0.0088308358 0.0196000958 0.0081267777 0.0439680761 0.0324588883 0.0034665720 0.0387499964 0.0181769181;
frequency C60pi7 = 0.0062848745 0.0026246919 0.0030342510 0.0005324147 0.0073027627 0.0034409089 0.0009741492 0.0019578159 0.0102225186 0.0180592309 0.1179064681 0.0016205916 0.0234721825 0.3974552519 0.0020165583 0.0056903327 0.0037091821 0.0598639097 0.3185565304 0.0152753744;
frequency C60pi8 = 0.1815005560 0.0026845411 0.0148484537 0.0025145485 0.4205633920 0.0014097001 0.0007088144 0.0461854175 0.0014374605 0.0041745536 0.0098310464 0.0006474254 0.0041611385 0.0068976432 0.0038767247 0.1864537050 0.0687189855 0.0027083549 0.0061033012 0.0345742379;
frequency C60pi9 = 0.0600740822 0.0367642654 0.0134869242 0.0170572285 0.0070719770 0.0142469806 0.0127486975 0.0343564471 0.0305859029 0.0204571345 0.0994551128 0.0212367087 0.0318165939 0.1140907926 0.0297628218 0.0505792699 0.0339368402 0.2312808862 0.1192491702 0.0217421638;
frequency C60pi10 = 0.0708394513 0.0474098489 0.0416822304 0.0324482918 0.0131641265 0.0494874703 0.0508264389 0.0183309196 0.0567272697 0.0650369079 0.1282255556 0.0343618389 0.0390362930 0.0594359563 0.0135608209 0.0551343199 0.0642260358 0.0137118382 0.0673934289 0.0789609573;
frequency C60pi11 = 0.0617689371 0.0076332888 0.0303081645 0.3430234188 0.0007199837 0.0307856241 0.3792509407 0.0284658686 0.0079592120 0.0016999627 0.0039945339 0.0216076877 0.0019734329 0.0009814186 0.0174791407 0.0337831940 0.0203426591 0.0006130268 0.0017102752 0.0058992300;
frequency C60pi12 = 0.0421559537 0.1042068314 0.0286980872 0.0164385240 0.0044450330 0.1393690851 0.0531949072 0.0134711207 0.0177764997 0.0267727728 0.1967237776 0.1323735242 0.1182827521 0.0086728324 0.0051837880 0.0255852718 0.0333292020 0.0045852327 0.0070281498 0.0217066546;
frequency C60pi13 = 0.2814809927 0.0100367066 0.0172867775 0.0064385734 0.0258337508 0.0133101925 0.0115046410 0.0270054934 0.0054629657 0.0188216093 0.0190993462 0.0098712843 0.0158719589 0.0050481705 0.0129510033 0.1886808600 0.2427104979 0.0012274627 0.0036052922 0.0837524211;
frequency C60pi14 = 0.2769188320 0.0017226995 0.0021315271 0.0011672545 0.0318292645 0.0018216251 0.0024752467 0.0199646887 0.0005170863 0.0983109006 0.0489264326 0.0016232163 0.0173414948 0.0070843906 0.0070179705 0.0336348952 0.0814141404 0.0007118144 0.0032942319 0.3620922883;
frequency C60pi15 = 0.1577797792 0.1112140270 0.0570403237 0.0648290471 0.0053318076 0.1065373681 0.0913586945 0.0906209718 0.0533809635 0.0029171632 0.0156225571 0.0782148712 0.0045758969 0.0025047816 0.0067077844 0.0929310045 0.0393122597 0.0028575821 0.0077590269 0.0085040899;
frequency C60pi16 = 0.0593735135 0.0354740772 0.1151175314 0.2189482708 0.0015332173 0.0688752402 0.1819422913 0.0813707101 0.0220478285 0.0020993577 0.0056191259 0.0750172075 0.0021871739 0.0010838321 0.0109737422 0.0726449461 0.0380238271 0.0007346460 0.0026664883 0.0042669729;
frequency C60pi17 = 0.0978066326 0.0265576438 0.0101843505 0.0120781428 0.0064138404 0.0307876446 0.0291282947 0.0128912798 0.0128036716 0.0723904209 0.1279438950 0.0245630658 0.0303267312 0.0198963719 0.2723524069 0.0350549441 0.0484557340 0.0046842467 0.0104773833 0.1152032995;
frequency C60pi18 = 0.0124023388 0.0030680354 0.0009239105 0.0006037316 0.0041885695 0.0032957441 0.0012524000 0.0011306791 0.0013542104 0.2344167852 0.4550557697 0.0016718177 0.0667307666 0.0610615367 0.0037076169 0.0019420934 0.0067612939 0.0038937184 0.0074911765 0.1290478057;
frequency C60pi19 = 0.0794230623 0.1294739355 0.0662792725 0.0587236242 0.0019919499 0.1143880588 0.1246900644 0.0325432311 0.0238605372 0.0036277150 0.0097987961 0.2147597316 0.0041846209 0.0012869951 0.0142410239 0.0615807386 0.0477333594 0.0006525371 0.0029420233 0.0078187231;
frequency C60pi20 = 0.0248148778 0.0083552910 0.1888915388 0.4278832998 0.0027839717 0.0210777725 0.1432386297 0.0643968435 0.0185736870 0.0022506941 0.0034558626 0.0179274104 0.0015714503 0.0014680353 0.0073768035 0.0377003132 0.0187767966 0.0005891859 0.0042602708 0.0046072655;
frequency C60pi21 = 0.0017003427 0.0060674330 0.0004222900 0.0010711490 0.0029059420 0.0016424179 0.0011731741 0.0035579609 0.0027630465 0.0012291190 0.0127420810 0.0004273804 0.0025671348 0.0513377024 0.0013536738 0.0011871674 0.0014033068 0.8640436936 0.0390912582 0.0033137266;
frequency C60pi22 = 0.0468360682 0.0639796924 0.0205603686 0.0185615516 0.0059954138 0.0557030821 0.0705436036 0.0045435329 0.0152062773 0.1550613356 0.0824253382 0.0866248354 0.0245854443 0.0080177192 0.0081485616 0.0237025617 0.0962054496 0.0018368673 0.0067131723 0.2047491243;
frequency C60pi23 = 0.0258764792 0.0201097124 0.0298384107 0.0107037437 0.0142503909 0.0158529432 0.0105649532 0.0073064999 0.1411078834 0.0114777629 0.0407992414 0.0119179202 0.0098798997 0.1876429961 0.0051228805 0.0275699644 0.0170764901 0.0405124999 0.3536390834 0.0187502449;
frequency C60pi24 = 0.0296285022 0.0046400334 0.0034944393 0.0008851024 0.0090046468 0.0055481111 0.0033046518 0.0027969482 0.0050701500 0.2583397750 0.2668085481 0.0046690936 0.0770825277 0.0408798247 0.0026918193 0.0068538089 0.0322265673 0.0035506055 0.0153353414 0.2271895033;
frequency C60pi25 = 0.0555725806 0.0098447861 0.0409064430 0.0140389597 0.0097418602 0.0068727710 0.0069443190 0.0157956555 0.0041631258 0.0069826497 0.0075271247 0.0139224817 0.0058762687 0.0034496730 0.0119733364 0.3482466393 0.4213655981 0.0010061491 0.0026576772 0.0131119012;
frequency C60pi26 = 0.0682671212 0.0615207091 0.0530661192 0.0360278709 0.0141433148 0.0612274332 0.0497415394 0.0268696520 0.1127674983 0.0132646615 0.0544493838 0.0482609047 0.0170033964 0.0803375967 0.0191949850 0.0671839752 0.0443995774 0.0199957919 0.1255070748 0.0267713947;
frequency C60pi27 = 0.0792618808 0.0638377192 0.0635289371 0.0436646174 0.0049503302 0.0666365188 0.0829639117 0.0183428565 0.0233169239 0.0249427251 0.0221483402 0.0932577596 0.0120893380 0.0049131149 0.0126360122 0.1334848656 0.1916745928 0.0018040086 0.0062353115 0.0503102360;
frequency C60pi28 = 0.0731759112 0.2105335985 0.0324200854 0.0110007149 0.0123458504 0.0858951989 0.0349942684 0.0224509173 0.0386903280 0.0246226304 0.0508307349 0.1783344831 0.0185740720 0.0093148787 0.0148722772 0.0603181436 0.0649574934 0.0051046395 0.0130597421 0.0385040321;
frequency C60pi29 = 0.0878402710 0.0110331750 0.0060801213 0.0032803903 0.0171147088 0.0109831614 0.0101465790 0.0087090941 0.0054902234 0.1987761871 0.1756460821 0.0082096925 0.0417232903 0.0191954435 0.0111283542 0.0209862621 0.0697718709 0.0031744014 0.0081905473 0.2825201446;
frequency C60pi30 = 0.0990215820 0.0349351987 0.0211149501 0.0118797946 0.0108995677 0.0557710676 0.0278999992 0.0240250097 0.0123445071 0.0776564721 0.2354511299 0.0322817789 0.1207665429 0.0214442058 0.0075655541 0.0524170141 0.0649785115 0.0047075806 0.0077328724 0.0771066610;
frequency C60pi31 = 0.0601641168 0.0161995226 0.2783522747 0.0337188808 0.0315066987 0.0210645987 0.0059839451 0.0543080710 0.0531523512 0.0070650825 0.0070698142 0.0139598368 0.0088298653 0.0069525877 0.0075834331 0.2829802556 0.0860317092 0.0014966551 0.0134849454 0.0100953553;
frequency C60pi32 = 0.0049781737 0.0018412331 0.0007012207 0.0005315368 0.0052978737 0.0024089907 0.0007630546 0.0015051317 0.0041575221 0.0443828633 0.4417417476 0.0011615060 0.0602807417 0.3351117140 0.0027847686 0.0025795769 0.0030288544 0.0171302592 0.0458455751 0.0237676560;
frequency C60pi33 = 0.0251996593 0.1114468110 0.0142031925 0.0041012288 0.0097099500 0.0620070749 0.0262571641 0.0038067269 0.0431938935 0.0974043253 0.2447197423 0.0824312856 0.0539323021 0.0429091639 0.0052658505 0.0096093107 0.0251183002 0.0146571900 0.0456965140 0.0783303143;
frequency C60pi34 = 0.0230361648 0.0014748749 0.0013534390 0.0006264439 0.0048580122 0.0009870046 0.0015762583 0.0011565336 0.0008899238 0.3952895890 0.0576537208 0.0014663528 0.0140986541 0.0072127040 0.0020177885 0.0028770237 0.0205580852 0.0005477695 0.0019539080 0.4603657493;
frequency C60pi35 = 0.1408776963 0.0297808449 0.0171297613 0.0285076933 0.0032213718 0.0320632225 0.0423838922 0.0299558472 0.0131321477 0.0066914481 0.0195120028 0.0383781635 0.0036276863 0.0041231064 0.4383466229 0.0851400095 0.0422765692 0.0013236871 0.0037087638 0.0198194632;
frequency C60pi36 = 0.4442491220 0.0050216551 0.0102305117 0.0057193038 0.0235405374 0.0055997640 0.0064889886 0.0822687710 0.0025505743 0.0033615104 0.0040990063 0.0038097073 0.0028683069 0.0024413211 0.0162890960 0.2999969708 0.0559664935 0.0007735426 0.0020639824 0.0226608347;
frequency C60pi37 = 0.0898717958 0.0070958305 0.0130067619 0.0129166888 0.0044131479 0.0023806547 0.0058957027 0.8087563021 0.0016517855 0.0004339282 0.0015564455 0.0033939025 0.0004253422 0.0008073572 0.0034128140 0.0362876891 0.0032887534 0.0015223902 0.0008537454 0.0020289624;
frequency C60pi38 = 0.0550840246 0.0472254260 0.1877829604 0.1273796123 0.0035824944 0.0527969268 0.0655884730 0.0637607521 0.0404883483 0.0075574152 0.0136304510 0.0867682792 0.0081684229 0.0040375032 0.0110681809 0.1263380956 0.0752544318 0.0013563681 0.0118590434 0.0102727908;
frequency C60pi39 = 0.0117681394 0.0442558806 0.0844144627 0.0144712108 0.0070388254 0.1038342049 0.0110901161 0.0049626578 0.4337194047 0.0061337038 0.0298794939 0.0137928558 0.0076237551 0.0338266335 0.0081346096 0.0140571089 0.0108276801 0.0080683065 0.1437251732 0.0083757773;
frequency C60pi40 = 0.0159285638 0.0048098656 0.0032692643 0.0010966937 0.0080519916 0.0134552459 0.0021324215 0.0025086365 0.0049192147 0.0501543893 0.5307634291 0.0035599431 0.2160085187 0.0743650717 0.0045247350 0.0066922196 0.0119092283 0.0070928134 0.0106565111 0.0281012433;
frequency C60pi41 = 0.0195973253 0.0105142992 0.3289103336 0.3099848991 0.0034539049 0.0116196758 0.0250777800 0.0627528956 0.0295961112 0.0032650434 0.0028246884 0.0240963907 0.0008425062 0.0019706550 0.0049062781 0.1064984500 0.0438053705 0.0006333959 0.0056197958 0.0040302013;
frequency C60pi42 = 0.0833804360 0.0125871438 0.0969824220 0.0686820704 0.0081981143 0.0121520930 0.0227415415 0.0982291876 0.0073954898 0.0017471177 0.0039653113 0.0129342146 0.0019557975 0.0024132583 0.0355924232 0.3115606483 0.2113368612 0.0016329034 0.0017991083 0.0047138579;
frequency C60pi43 = 0.0181409133 0.4129662563 0.0233205154 0.0033333547 0.0085143598 0.0526694251 0.0096531879 0.0224552642 0.0375238929 0.0035090482 0.0149146621 0.3208065790 0.0046098856 0.0035426859 0.0087197469 0.0262309419 0.0131791136 0.0034766995 0.0079588201 0.0044746474;
frequency C60pi44 = 0.2494227404 0.0185481724 0.0164119567 0.0169234299 0.0122862654 0.0228501981 0.0370491083 0.0347467705 0.0087069587 0.0595718359 0.0451065029 0.0177064733 0.0204556127 0.0077360919 0.0686403544 0.0889295672 0.0986017356 0.0028603862 0.0061938477 0.1672519917;
frequency C60pi45 = 0.1419737638 0.0373945961 0.0576296888 0.0537452477 0.0068856658 0.0286239972 0.0407540287 0.3988107872 0.0152895617 0.0016627616 0.0092348297 0.0314273807 0.0055425500 0.0040286132 0.0180328866 0.1123731997 0.0242478202 0.0025909098 0.0049054208 0.0048462908;
frequency C60pi46 = 0.0178903305 0.1958843646 0.0155853897 0.0031054277 0.0290304227 0.1051819261 0.0040503389 0.0100480293 0.1252696215 0.0016708003 0.0722356645 0.0233340169 0.0116142354 0.0238913260 0.0009938415 0.0181675536 0.0186260222 0.2260554691 0.0859787232 0.0113864962;
frequency C60pi47 = 0.1454758367 0.0420979067 0.0400419720 0.1294249748 0.0014186329 0.0906469055 0.2471353458 0.0319650773 0.0130426183 0.0058525371 0.0123593139 0.0818154090 0.0044178939 0.0017552077 0.0151135525 0.0656688174 0.0511289472 0.0007731441 0.0029258438 0.0169400635;
frequency C60pi48 = 0.0169799462 0.0242346701 0.1318047919 0.1043655101 0.0022087215 0.0269349684 0.0376379591 0.5404470183 0.0181137053 0.0007459679 0.0021146994 0.0508617611 0.0009473769 0.0006780593 0.0038754401 0.0297030159 0.0045836180 0.0006031889 0.0015704090 0.0015891728;
frequency C60pi49 = 0.0402646249 0.1152022601 0.0323829165 0.0293968352 0.0039388655 0.2497008043 0.1603524245 0.0129260411 0.0617967839 0.0098491259 0.0354918823 0.1448804422 0.0124818865 0.0041153375 0.0043374229 0.0243246958 0.0305645368 0.0026676598 0.0097227847 0.0156026694;
frequency C60pi50 = 0.2256914610 0.0523417493 0.0244308734 0.0637125217 0.0043390149 0.0578159236 0.1154830640 0.0867335173 0.0131066949 0.0085086217 0.0193314218 0.0660468804 0.0064877206 0.0027440054 0.0611149102 0.1070877179 0.0507677144 0.0013695913 0.0028982948 0.0299883012;
frequency C60pi51 = 0.0033164209 0.0015310773 0.0030830171 0.0008266472 0.0051890730 0.0011024889 0.0005134130 0.0010432830 0.0278451262 0.0041895268 0.0111212494 0.0007149922 0.0023621780 0.3801761447 0.0008365077 0.0035876698 0.0023608948 0.0333346985 0.5107889643 0.0060766272;
frequency C60pi52 = 0.1995014012 0.0236078675 0.0392254543 0.0094955104 0.0584590451 0.0254265363 0.0125535371 0.0939787338 0.0341857201 0.0140209879 0.0449387571 0.0118723304 0.0246990633 0.0634433944 0.0145385320 0.1663920640 0.0533159207 0.0129802666 0.0606346163 0.0367302614;
frequency C60pi53 = 0.0319448994 0.1011667268 0.2084709220 0.0378074649 0.0066040348 0.0766372935 0.0279488190 0.0365541130 0.2088643258 0.0047542347 0.0156545731 0.0868664783 0.0043253317 0.0108915768 0.0060899575 0.0577656939 0.0302051160 0.0026001883 0.0387897304 0.0060585202;
frequency C60pi54 = 0.0776799515 0.0142518583 0.0403216692 0.0080651725 0.0140092962 0.0179995517 0.0112622427 0.0136868237 0.0133729897 0.1239635380 0.0724670993 0.0129144967 0.0420745442 0.0173584908 0.0117084432 0.0922723571 0.2316899445 0.0028153633 0.0141726542 0.1679135132;
frequency C60pi55 = 0.1183662657 0.0805192606 0.0259524932 0.0495595439 0.0035624835 0.1204924917 0.1537589210 0.0194993426 0.0229373171 0.0302661211 0.0571250629 0.0982304112 0.0171727472 0.0068665705 0.0175153030 0.0486588400 0.0635796210 0.0023008307 0.0083027431 0.0553336300;
frequency C60pi56 = 0.0528559899 0.0193569043 0.0264743774 0.2092761515 0.0008625883 0.1212409715 0.4024189781 0.0155838458 0.0124148798 0.0054864832 0.0090256472 0.0497017031 0.0042357114 0.0012650715 0.0063185636 0.0197262901 0.0235463735 0.0008381610 0.0033948741 0.0159764347;
frequency C60pi57 = 0.0344366215 0.0426221820 0.1636716191 0.1139007491 0.0020985982 0.0605413987 0.0541780220 0.3361639671 0.0461776737 0.0003463416 0.0048355678 0.0667552967 0.0019704509 0.0031557619 0.0040369775 0.0481173332 0.0089148085 0.0006510101 0.0054145649 0.0020110555;
frequency C60pi58 = 0.1153088951 0.0151278638 0.0458476603 0.1755516676 0.0014962362 0.0366731222 0.1749410045 0.0394181311 0.0132401530 0.0056912974 0.0101409559 0.0433118387 0.0030332064 0.0015700232 0.1665802563 0.0871536033 0.0468260603 0.0007515702 0.0031432715 0.0141931831;
frequency C60pi59 = 0.3865149348 0.0037579334 0.0030420497 0.0022366810 0.0218928357 0.0021464743 0.0031387843 0.3694353983 0.0014672902 0.0085376076 0.0127257242 0.0018840458 0.0080581695 0.0039281367 0.0158688291 0.0808877279 0.0305195935 0.0009922880 0.0019020345 0.0410634615;
frequency C60pi60 = 0.0146570745 0.0028841333 0.0012998335 0.0005210575 0.0024317913 0.0049362750 0.0014874369 0.0020953252 0.0010181940 0.1913901476 0.4432797758 0.0022898369 0.2217427062 0.0091637503 0.0007685153 0.0027251487 0.0170997497 0.0008779380 0.0014756028 0.0778557075;

model C60 = POISSON+G+FMIX{C60pi1:1:0.0169698865,C60pi2:1:0.0211683374,C60pi3:1:0.0276589079,C60pi4:1:0.0065675964,C60pi5:1:0.0141221416,C60pi6:1:0.0068774834,C60pi7:1:0.0146909701,C60pi8:1:0.0067225777,C60pi9:1:0.0018396660,C60pi10:1:0.0102547197,C60pi11:1:0.0230896163,C60pi12:1:0.0057941033,C60pi13:1:0.0125394534,C60pi14:1:0.0204526478,C60pi15:1:0.0070629602,C60pi16:1:0.0117982741,C60pi17:1:0.0068334668,C60pi18:1:0.0433775839,C60pi19:1:0.0318278731,C60pi20:1:0.0222546108,C60pi21:1:0.0102264969,C60pi22:1:0.0150545891,C60pi23:1:0.0134159878,C60pi24:1:0.0148552065,C60pi25:1:0.0239111516,C60pi26:1:0.0128776278,C60pi27:1:0.0222318842,C60pi28:1:0.0247444742,C60pi29:1:0.0214274810,C60pi30:1:0.0115001882,C60pi31:1:0.0076017389,C60pi32:1:0.0130258568,C60pi33:1:0.0093701965,C60pi34:1:0.0467194264,C60pi35:1:0.0441940314,C60pi36:1:0.0322263154,C60pi37:1:0.0402999891,C60pi38:1:0.0150234227,C60pi39:1:0.0104589903,C60pi40:1:0.0214742395,C60pi41:1:0.0154957836,C60pi42:1:0.0101789953,C60pi43:1:0.0227980379,C60pi44:1:0.0123204539,C60pi45:1:0.0066777583,C60pi46:1:0.0004150083,C60pi47:1:0.0344385130,C60pi48:1:0.0113663379,C60pi49:1:0.0127143049,C60pi50:1:0.0124323741,C60pi51:1:0.0262124415,C60pi52:1:0.0064994957,C60pi53:1:0.0103203293,C60pi54:1:0.0142463512,C60pi55:1:0.0215600067,C60pi56:1:0.0199150700,C60pi57:1:0.0038964200,C60pi58:1:0.0113448855,C60pi59:1:0.0128595846,C60pi60:1:0.0117656776};

end;
"""

# execute program
if __name__ == "__main__":
    start_time = time.time()
    main(arguments)
    end_time = time.time()
    logging.info(f"Runtime: {end_time - start_time} seconds")
    # shutil.copy('fundi_wrapper.log', arguments.outdir + '/fundi_wrapper.log')
