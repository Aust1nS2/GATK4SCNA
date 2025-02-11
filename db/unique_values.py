#!/usr/bin/env python3

import sys
import os


filepath = sys.argv[1]
file = open(filepath, "r")
header = file.readline()
id_dict = dict()
newID_list = list()
id_list = list()
c = 0
for line in file:
    line_list = line.strip('\n').split('\t')
    #print(c)
    #print(line_list)
    ID = line_list[3]
    #print(ID)
    c += 1
    #ID = line_list[0]+line_list[1] #this should be unique
    if ID in id_list:
        print("current_line:", line_list)
        print("overlap: ", ID, id_dict[ID])
        print()
    else: 
        id_list.append(ID)
        id_dict[ID] = line_list
file.close()
