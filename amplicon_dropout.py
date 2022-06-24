#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# extracts coordinates of amplicon from ARTIC v3 primer scheme
def get_coordinates(file):
    f = open(file)       
    start = []
    end = []
 
    for line in f:
        line = line.strip('\n')
        line = line.split('\t')
        if 'alt' not in line[3]:
            start.append(line[1])
            end.append(line[2])
    return(start,end)

coords_all = get_coordinates('/Users/s2261824/Downloads/amplicons.csv')


all_starts = []
all_ends = []
# add start and end of each amplicon from primer scheme
for i in range(0,196,2):
    all_starts.append(coords_all[0][i])
    all_ends.append(coords_all[1][i+1])

 
# takes base listfile as input, produces a list file with each 3-amplicon window missing across the 98 amplicons in the ARTIC v3 primer scheme
def amplicon_dropout_gofasta(listfile):
    # move in sliding window of three amplicons across the genome
    for n in range(0, 97):
        amplicons = [n, n+1, n+2]
        # create list of start, end for each amplicon in 3-amplicon window
        starts = []
        ends = []
        # loop through amplicon, pull start, end coordinates from primer list
        for amps in amplicons:
            coord_start = int(all_starts[amps])
            coord_end = int(all_ends[amps])

            # add primer coordinates to list
            starts.append(coord_start)
            ends.append(coord_end)
        
        # find start, end coordinates of 3-amplicon window
        proper_start = int(starts[0])
        proper_end = int(ends[-1])
        
        # create new ambiguity range covering dropout
        new_amb = f'{proper_start}-{proper_end}'
        
        # open list file, create output file
        output = open(f'path_to_output/dropout_{n+1}_to_{n+3}.csv', 'a')
        file = open(f'path_to_file/listfile')
            
        # add dropout for each query in the list file
        for line in file:
            ambs_final = []
            snps_final = []
            used_ambs = []
            if line.startswith('query'):
                output.write(line)
                continue
            else:
                line = line.strip('\n').split(',')
                ambs = line[2].split('|')
                snps = line[1].split('|')
                # process ambiguities
                ambs_processed_before = process_ambiguities(ambs)
                for amb in ambs_processed_before:
                    # only include ambiguities that do not overlap with dropout window in final list
                    if int(amb) not in range(proper_start, proper_end+1):
                        ambs_final.append(str(amb))
                # calculate length of window
                length = str(len(ambs_final)+(proper_end-proper_start)+1)
                ambs_final.append(new_amb)
                    
                # process SNPs, only include ones not in window in final list
                snps_processed = process_snps(snps)
                for snp in snps_processed[1]:
                    if int(snp) not in range(proper_start, proper_end+1):
                        polymorphism = snps[snps_processed[1].index(snp)]
                        snps_final.append(polymorphism)
                
                # removes empty strings in list so they don't affect final totals
                ambs_final = list(filter(None, ambs_final))
                snps_final = list(filter(None, snps_final))
                
                # rejoin line, write to output
                line[2] = '|'.join(ambs_final)
                line[1] = '|'.join(snps_final)
                line[3] = str(len(snps_final))
                line[4] = length
                line = ','.join(line)
                output.write(f'{line}\n')
        output.close()
    return

amplicon_dropout_gofasta(input_listfile)
