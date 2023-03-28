# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 22:11:30 2023

@author: Genglin Guo
"""

import argparse
import sys
import pathlib
import multiprocessing
import subprocess
import time

__version__ = '2.0, 2023-01-11'

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'BacMLST', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-o', '--outdir', required = False, type = str, default = 'Bac_MLST_result.txt',
                              help = 'Output directory')
    parser_group_1.add_argument('-r', '--refs', required = False, type = str, default = 'MLST_database',
                              help = 'fasta file contain all 16S rRNA sequence from NCBI')

    # Parameters
    parser_group_2.add_argument('--species', required = False, type = str, 
                                help = 'Determine the species of input')
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('--outfmt', required = False, type = str, default = 'mix', 
                                help = 'Determine the output format, two types, "mix" and "same" were provided')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'BacSpecies v' + __version__, 
                        help = 'Show version number and exit')
    return parser

def check_dependencies():
    # Checks dependencies are available
    dependencies = ['makeblastdb', 'blastn', 'tblastn']
    for i in dependencies:
        try:
            subprocess.check_call(['which', i], stdout = subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('Error: could not find %s tool' % i, file=sys.stderr)
            sys.exit(1)

def process_16s_reference():
    # Generate a dict, genbank_number : bacteria_name
    species_database = {}
    with open('reference_database', 'rt') as file:
    	for line in file:
    	    if line.startswith('>'):
    	        line = line[1:]
    	        info = line.split(' ')
    	        id = info[0]
    	        name = info[1] + ' ' + info[2]
    	        species_database[id] = name
    	    else:
    	        continue
    return species_database

def solve_path(refs, species, inputfile):
    # get the absolute path
    abrivate = ['Achromobacter', 'Aeromonas', 'Arcobacter', 'Borrelia', 'Brucella', 'Citrobacter', 
                'Cronobacter', 'Edwardsiella', 'Escherichia', 'Neisseria', 'Rhodococcus', 'Shewanella', 
                'Sinorhizobium', 'Streptomyces', 'Taylorella', 'Tenacibaculum', 'Ureaplasma', 'Wolbachia',
                'Geotrichum', 'Leptospira', 'Mycobacteria']
    if species.split('\t')[0] in abrivate:
        species = species.split('\t')[0] + ' spp'
    refpath = pathlib.Path(refs) / species
    if not pathlib.Path(refpath).is_dir():
        print('Error: Sequence file {} --- no {} database provided'.format(inputfile, species))
        sys.exit(1)
    profilepath = pathlib.Path(refpath).resolve() / 'MLST_profiles'
    return refpath, profilepath

def process_mlst_reference(profilepath):
    # Generate a dict, ST : alleles
    database = {}
    with open(profilepath, 'rt') as file:
    	for line in file.readlines()[1:]:
            profile = line.strip('\n')
            profile = profile.split('\t')
            ST = profile[0]
            alleles = profile[1:]
            database[ST] = alleles
    return database

def get_labels(profilepath):
    labels = []
    with open(profilepath, 'rt') as file:
        line = file.readlines()[0]
        line = line.strip('\n')
        line = line.split('\t')
        for i in line[1:]:
            labels.append(i)
    return labels

def get_best_species_result(inputfile, threads, species_database):
    # Perform blast, find the best match in reference_database
    repa = pathlib.Path('reference_database').resolve()
    blast_hits = run_blast(inputfile, repa, threads)
    best_match = []
    best_cov = 0.0
    best_pident = 0.0
    for hit in blast_hits:
        if hit.length < 1000:
            continue
        elif hit.pident >= best_pident and hit.query_cov >= best_cov:
            best_pident = hit.pident
            best_cov = hit.query_cov
            if species_database[hit.sqseqid] in best_match:
                continue
            else:
                best_match.append(species_database[hit.sqseqid])
    return best_match

def get_best_result(inputfile, refpath, threads, database, labels):
    # Perform blast, find the best match in reference_database
    best_match = []
    best_ST = ''
    for i in labels:
        filename = i + '.fas'
        repath = refpath / filename
        repa = pathlib.Path(repath).resolve()
        blast_hits = run_blast(inputfile, repa, threads)
        best_allele = ''
        best_cov = 0.0
        best_pident = 0.0
        for hit in blast_hits:
            if hit.length < 100:
                continue
            if hit.pident >= best_pident and hit.query_cov >= best_cov:
                best_pident = hit.pident
                best_cov = hit.query_cov
                best_allele = hit.qseqid
        if not best_cov == 100.0 and best_pident == 100:
            best_allele += '?'
        best_match.append(best_allele)
    for i in database.keys():
        if database[i] == best_match:
            best_ST = i
    return best_match, best_ST

def run_blast(inputfile, repa, threads):
    inpa = pathlib.Path(inputfile).resolve()
    blast_hits = []
    command = ['blastn', '-query', repa, '-subject', inpa, '-num_threads', str(threads), '-evalue', '0.00001', '-perc_identity', '99.5',
		'-outfmt', '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    if out == 0:
        print('The sequence file {} has no blast hit, please check the fasta file'.format(inputfile))
    for line in line_iterator(out):
        blast_hits.append(BlastResult(line))
    return blast_hits
        
def line_iterator(line_breaks):
    # Handle the BLAST output and remove the line breaks 
    line = -1
    while True:
        nextline = line_breaks.find('\n', line + 1)
        if nextline < 0:
            break
        yield line_breaks[line + 1:nextline]
        line = nextline

class BlastResult(object):
    # Handle the BLAST output
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0].split('_')[-1]
        self.sqseqid = parts[0]
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])

def generate_output(outdir, labels, outfmt):
    # Generate a blank output table file
    if pathlib.Path(outdir).is_file():
        return
    headers = ['Sequence', 'Species', 'ST']
    if outfmt == 'mix':
        sur = 'allele'
        for i in range(len(labels)):
            con_num = str(i + 1)
            headers.append(sur + con_num)
    else:
        for i in labels:
            headers.append(i)
    with open(outdir, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def output(outdir, best_match, best_ST, inputfile, species, outfmt, labels):
    # Generate output
    if best_ST:
        simple_output = str(inputfile) + ' : ' + species + ' ' + 'ST' + best_ST
        line = [str(inputfile), species, best_ST]
    else:
        simple_output = str(inputfile) + ' : ' + species + ' ' + 'Unavailable'
        line = [str(inputfile), species, 'NA']
    if outfmt == 'mix':
        for i in range(len(labels)):
            item = labels[i] + '(' + best_match[i] + ')'
            line.append(item)
    else:
        line += best_match
    print(simple_output)
    with open(outdir, 'at') as table:
        table.write('\t'.join(line))
        table.write('\n')   

def main():
    print('If you have any quesion or suggestion in using BacMLST, please contact Genglin Guo, e-mail: 2019207025@njau.edu.cn')
    starttime = time.perf_counter()
    # Initialize
    args = get_argument().parse_args()
    check_dependencies()
    # Get species of input file
    if not args.species:
        species_database = process_16s_reference()
    for inputfile in args.input:
        best_species = ''
        if not args.species:
            best_species_match = get_best_species_result(inputfile, args.threads, species_database)
            if best_species_match:
                best_species = best_species_match[0]
            else:
                print('no species match for file {}'.format(inputfile))
                sys.exit(1)
                with open(args.outdir, 'at') as table:
                    table.write(inputfile + '\t' + 'no species match' + '\n')
                continue
        else:
            best_species = args.species
        #Fix path and run blast
        refpath, profilepath = solve_path(args.refs, best_species, inputfile)
        labels = get_labels(profilepath)
        database = process_mlst_reference(profilepath)
        best_match, best_ST = get_best_result(inputfile, refpath, args.threads, database, labels) 
    #Generate output
        generate_output(args.outdir, labels, args.outfmt)
        output(args.outdir, best_match, best_ST, inputfile, best_species, args.outfmt, labels)
    endtime = time.perf_counter() - starttime
    print('Total time consumed : {:.1f}h{:.1f}m{:.1f}s'.format(endtime // 3600, endtime % 3600 // 60, endtime % 60))
   
main()
