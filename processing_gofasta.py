#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# function takes gofasta bin (in list form) for a given query as input (level = up, down, side, same)
def find_branches(input_list, level, query):
    branch_list = []
    perc_list = []
    # if no matches in this bin, return message
    if input_list == ['']:
        print(f'No matches for query {query} at level {level}.')

    else:
        for element in input_list:
            # ensure no repetition- per lineage not per hit
            if element not in branch_list:
                # determine percentage of hits represented by each lineage
                branch_perc = (input_list.count(element)/len(input_list))*100
                perc_list.append(branch_perc)
                branch_list.append(element)

    # find most abundant hit, return hit lineage, percentage
        most_common = max(perc_list)
        index = perc_list.index(most_common)

    return(most_common, branch_list[index])

# takes as input gofasta output CSV file and txt file of push distances, returns dictionary of distance pushed up for each query (default = 1)
def push_distances(input_push_file, input_gofasta_file)

    file = open(input_push_file)
    directions = []
    queries = []
    push_dist = []
    for line in file:
        line= line.strip('\n')
        directions.append(line[5])
        start = 1+line.find('(')
        end = line.find(')')
        
        
        push_start = 2+ line.find(':')
        # pull only ranges pushed up
        if line[5] == 'u':
            queries.append(line[start:end])
            push_dist.append(line[push_start:])
    file.close()
    
    gofasta_file = open(input_gofasta_file)
    up_dict = {}
    for line in gofasta_file:
        if line.startswith('q'):
            continue
        else:
            line = line.strip('\n').split(',')
            start = line[0].find('(')+1
            end = line[0].find(')')
            query_name = line[0][start:end]
            if query_name in queries:
                index = queries.index(query_name)
                up_dict[query_name] = push_dist[index]
            # if name is not in ????
            else:
                up_dict[query_name] = 1
        
    gofasta_file.close()
    return(up_dict)
            

# function takes as input gofasta output plus text file of push distances
def process_gofasta(gofasta,push)
    # create dictionary of distances pushed for each query to fill up bin
    up_dict = push_distances(push,gofasta)
    f = open(gofasta)
    # create output file, write header
    output = open('path_to_output/output.csv', 'a')
    output.write('query,lineage,same_call,percentage,closest_same,up_call,percentage,closest_up,push_up,parent,assessment,is_same,n_lins\n')
    
    # alreadt have ancestry dict????????????????
    for line in f:
        if line.startswith('query'):
            continue
        else:
            line = line.strip('\n').split(',')
            query_name = line[0]
            same = line[1]
            up = line[2]
            down = line[3]
            side = line[4] 

            if same != '':
                calls = same.split(';')
                branches = find_branches(calls, 'same', query_name)
                call_same = branches[1]
                perc_same = branches[0]
                closest_same = calls[0]
                is_same = 'TRUE'

            elif same == '':
                call_same = ''
                perc_same = ''
                closest_same = ''
                is_same = 'FALSE'
                
            if up != '':
                up = up.split(';')
                up_closest = up[0]
                start = query_name.find('(')+1
                end = query_name.find(')')
                query = query_name[start:end]
                up_counter += 1
                push = int(up_dict[query])
                branches = find_branches(up, 'up', query_name)
                call_up = branches[1]
                perc_up = branches[0]
                
            # adds additional column for testing, sorting call by whether it is the correct, incorrect, parental or child lineage
                if call_up == query:
                    call_type = 'correct'
                elif call_up == ancestry[query]:
                    call_type = 'parent'
                elif ancestry[call_up] == query:
                    call_type = 'child'
                else:
                    call_type = 'incorrect'

        # overwrites up call if same call present
            if call_same == query:
                call_type = 'correct'
            output.write(f'{query_name},{query},{call_same},{perc_same},{closest_same},{call_up},{perc_up},{up_closest},{push},{ancestry[query]},{call_type},{is_same},{lins_dict[query]}\n')
    output.close()
    f.close()
    return
