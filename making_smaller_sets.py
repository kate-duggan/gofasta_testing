#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# takes lineages.csv from PANGO designation GitHub as input, processes to give lineage, taxa lists
def get_lineages_taxa(pango_file):
    lineage_file = open(pango_file)
    lineage_list = []
    taxa_list = []
    for line in lineage_file:
        line = line.rstrip('\n').split(',')
        taxon = line[0]
        lineage = line[1]
        lineage_list.append(lineage)
        taxa_list.append(taxon)
    return(lineage_list, taxa_list)

lineages = get_lineages_taxa('path_to_lineages/lineages.csv')

# find unique lineages
used_lins = []
for lineage in lineages[0]:
    if lineage not in used_lins:
        used_lins.append(lineage)
        

# pull least ambiguous sequence for each lineage from file
def pull_best_from_processed(file, lineage):
    f = open(file)
    indices = []
    inner_counter = 0
    ambs_list = []
    for line in f:

        if line.startswith('query'):
        
            continue
        else:
            # increase counter to match index of query
            inner_counter +=1
        
            line = line.strip('\n').split(',')
           
           # add index to list for given input lineage
            if line[0] == lineage:
                indices.append(inner_counter)
                # store number of ambiguities
                ambs_list.append(int(line[4]))
                
    # for non-empty lists, pull least ambiguous representative
    if len(ambs_list) != 0:
        best = min(ambs_list)
        chosen = ambs_list.index(best)
        
        chosen_taxa = indices[chosen]
    
    # otherwise, lineage is not represented in the input dataset
    else:
        chosen_taxa = 'NA'

    return(chosen_taxa)

# create gofasta list of least-ambiguous n representatives of each lineage (output_size = n)
def create_dataset(input_dataset, output_size):
    
    # loop through n times to pull least ambiguous queries after each round
    for n in range(0,int(output_size)):
        chosen_lins = []
        counter = 0
        # loop through list of unique lineage names
        for lineage in used_lins:
            counter += 1
            least_ambiguous = pull_best_from_processed(input_dataset, lineage)
            chosen_lins.append(least_ambiguous)
        
        
        f = open(f'path_to_input_dataset/{input_dataset}')
        output = open(f'path_to_output/one_of_each_round{n}.csv', 'a')
        output2 = open(f'path_to_output/dataset_after_round{n}.csv', 'a')
        
        counter = 0
        for line in f:
            # write header to both output files
            if line.startswith('query'):
                output.write(f'{line}')
                output2.write(f'{line}')
                continue
            # count index of each line
            else:
                counter += 1
    
                # if index of query is in list of chosen indices, write to output database
                if counter in chosen_lins:
                    
                    output.write(line)
                
                # if query is not chosen as least-ambiguous at this round, write to second database (search database for output_size = 1, otherwise input database for future rounds)
                else:
    
                    output2.write(f'{line}')

                        
            f.close()
            output.close()
            output2.close()
