#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 11:33:36 2022

@author: s2261824
"""
import random
from Bio import SeqIO

# function takes subset of SNPs, mutates to new (non-ancestral) base, returns final SNP list
def novel_mutation(targeted_snps, all_snps):
    # pull location, current ancestral and derived bases from SNPs to be mutated
    for snp in targeted_snps:
        ancestral = snp[0]
        base = int(snp[1:len(snp)-1])
        current_snp = snp[-1:]
        new = random.choice(['A', 'T', 'G', 'C'])
        # if new choice is the same as the current derived or ancestral base, choose again
        while new == current_snp or new  == ancestral:
            new = random.choice(['A', 'T', 'G', 'C'])

        # define new SNP at same location
        new_snp = f'{ancestral}{base}{new}'
        # remove old snp from final list
        all_snps.remove(snp)
        # add new snp to final list
        all_snps.append(new_snp)
    return(all_snps)

# function takes gofasta list file, mutation type (novel or reversion)
def targeted_denovo(input_file, mutation_type):
    
    # loop through max number of SNPs (most mutated sequence has 91)
    for i in range(1,91):
        file = open(f'path_to_output/{input_file}')
        output = open(f'path_to_output/denovo_targeted_{i}.csv', 'a')
        for line in file:
            # write header line to file
            if line.startswith('q'):
                output.write(line)
            else:
                line = line.strip('\n').split(',')
                snps = line[1].split('|')
                ambs = line[2].split('|')

                # only keep mutating if no. of SNPs to be mutated is less than number of query SNPs
                if i <= len(snps):
                    sites = random.sample(snps, i)
                    
                    # if call specifies novel mutation, call novel_mutation function
                    if mutation_type == 'novel':
                        
                        snps = novel_mutation(sites, snps)
                    
                    elif mutation_type == 'reversion':
      
                    # if type is reversion, SNP must be removed from list
                        for snp in sites:
                            snps.remove(snp)
                            
                            # find location of SNP, add to ambiguities list
                            base = snp[1:len(snp)-1]
                            ambs.append(base)

                # if no. of SNPs is >= number of query SNPs, remove SNPs from file
                else:
                    snps = ''
                # rejoin, write to output
                line[3] = str(len(snps))
                line[1] = '|'.join(snps)
                line[2] = '|'.join(ambs)
                line[4] = str(len(ambs))
                line = ','.join(line)
                output.write(f'{line}\n')
        output.close()
    return

                            
targeted_denovo(input_file, 'novel')
targeted_denovo(input_file, 'reversion')

        


    
    



        
                           
                            
