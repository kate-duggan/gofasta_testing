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


def generate_ambiguities(ambiguities, snp, perc):
    gofasta_ambs = process_ambiguities(ambiguities)
    gofasta_snps = process_snps(snp)
    
    pos_list = []
    pos_counter = 0
    while pos_counter <= int((29903/100)*percent):
        pos = random.randint(0, len(seq)-1)
        if pos not in pos_list and seq[pos] != 'N' and seq[pos] != '-' and str(pos) not in gofasta_ambs:
            pos_list.append(pos)
            pos_counter += 1

    for base in pos_list:
        
        # if ambiguity site is listed as a SNP, remove from SNP list
        if str(base) in gofasta_snps[1]:
            index = gofasta_snps[1].index(str(base))
            snp.pop(index)
            
            # add to ambiguities
            ambiguities.append(base)
        
        return(ambiguities, snp)


def ambs_general(file):
    
    for p in range(20, 1020, 20):
        percent = p/100
        f = open(file)
        new_file = open(f'/Users/s2261824/Downloads/Data/ambiguity_{percent}.csv',"a")
        
        
        for line in f:
            if line.startswith('query'):
                new_file.write(f'{line}')
                continue
            else:
                line = line.strip('\n').split(',')
                snps = line[1].split('|')
                print(snps)
                ambs = line[2].split('|')
                
                # process SNPs, ambiguities
                gofasta_ambs = process_ambiguities(ambs)
                gofasta_snps = process_snps(snps)
                
                pos_list = []
                pos_counter = 0
                indices = []
                # select sites at random throughout the genome (excluding existing ambiguities in either query or reference)
                while pos_counter <= int((29903/100)*percent):
                    pos = random.randint(0, len(seq)-1)
                    if pos not in pos_list and seq[pos] != 'N' and seq[pos] != '-' and pos not in gofasta_ambs:
                        pos_list.append(pos)
                        pos_counter += 1
            
                for base in pos_list:
                    # if ambiguity site is listed as a SNP, remove from SNP list
                    if str(base) in gofasta_snps[1]:
                        index = gofasta_snps[1].index(str(base))
                        indices.append(index)
                # add ambiguity to ambs list
                    ambs.append(str(base))
    
    # remove ambiguities from SNP lisr
                for i in range(0, len(indices)):
                    snps.pop(i)
                        
                # rejoin, write to output
                line[4] =str(int(line[4])+ pos_counter)
                line[3] = str(len(snps))
                line[1] = '|'.join(snps)
                line[2] = '|'.join(ambs)
                line = ','.join(line)
                new_file.write(f'{line}\n')
                            
        new_file.close()

return

ambs_general('path_to_file')
