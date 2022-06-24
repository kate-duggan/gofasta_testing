#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from Bio import SeqIO

    
# process reference genome from fasta file
def process_reference(genome):
    records = SeqIO.parse(genome, 'fasta')
    for record in records:
        seq = str(record.seq)
    return(seq)

seq = process_reference('path/to/reference')

# sort ambiguities so internal ambiguities are acccessible as int, rather than range stored as string 'amb_start-amb_end'
def process_ambiguities(ambiguities):
    ambs_all = []
    # loop through ambiguities list
    for amb in ambiguities:
        # if ambiguity is a range, store start and end
        if '-' in str(amb):
            loc = int(amb.find('-'))
            lower_lim = int(amb[:loc])
            upper_lim = int(amb[loc+1:])
            # add bases between start and end
            for i in range(lower_lim, upper_lim+1):
                ambs_all.append(f'{i}')
    # otherwise, ambiguity is already in base format, no processing necessary
        else:
            if amb != '':
                ambs_all.append(int(amb))

    return(ambs_all)

# process list of SNPs to give lists of ancestral bases, SNP locations and derived bases
def process_snps(snps_unprocessed):
    anc_list = []
    mut_list = []
    bases_list = []
    for site in snps_unprocessed:
        ancestral = site[0]
        anc_list.append(ancestral)
        mut = site[len(site)-1]
        mut_list.append(mut)
        bases = site[1:len(site)-1]
        bases_list.append(bases)
    # return separate lists
    return(anc_list, bases_list, mut_list)  

# randomly generate specified percentage of genome-wide mutations
def generate_mutations(ambiguities, snps, perc):
    # ensure ambiguities are in correct format to check if mutation overlaps
    gofasta_ambs = process_ambiguities(ambiguities)
    # get SNP coordinates
    gofasta_snps = process_snps(snps)
        
    pos_list = []
    pos_counter = 0
    while pos_counter <= int((29903/100)*perc):
        pos = random.randint(0, len(seq)-1)
        
        # avoid duplicates, sites that are ambiguous in either query or reference
        if pos not in pos_list and seq[pos] != 'N' and seq[pos] != '-' and str(pos) not in gofasta_ambs:
            pos_list.append(pos)
            pos_counter += 1

     # loop through selected positions
    for base in pos_list:
        # if site is already mutated in query, store current mutation and ancestral base
        if str(base) in gofasta_snps[1]:
            index = gofasta_snps[1].index(str(base))
            current = gofasta_snps[2][index]
            ancestral = gofasta_snps[0][index]
            snps.pop(index)
        # if base is not currently a SNP, ancestral is taken from the reference
        else:
            ancestral = seq[base].upper()
            current = ''

        new = random.choice(['A', 'T', 'G', 'C'])
        
        while new == current or new  == ancestral:
            new = random.choice(['A', 'T', 'G', 'C'])
        
        # create SNP of new base at specified site
        snp = f'{ancestral}{base}{new}'
        snps.append(snp)
    return(snps)

# take input gofasta list, mutate up to 10% in steps of 0.2%
def mutations(file):
    # define percentage mutation
    for p in range(2, 102, 2):
        percent = p/10
        f = open(file)
        new_file = open(f'path_to_file/mutation_{percent}.csv','a')
    
    
        for line in f:
            # write header to output
            if line.startswith('query'):
                new_file.write(f'{line}')
                continue
            else:
                line = line.strip('\n').split(',')
                # convert SNP, ambs fields into lists
                snps = line[1].split('|')
                ambs = line[2].split('|')
                # process ambiguities
                ambs_processed = process_ambiguities(ambs)
                # add specified percentage of random mutations
                snps = generate_mutations(ambs_processed, snps, percent)
            
            # rejoin, write to output
            line[3] = str(len(snps))
            line[1] = '|'.join(snps)
            line = ','.join(line)
            new_file.write(f'{line}\n')
    
        new_file.close()
            
    return

mutations('path_to_file')
