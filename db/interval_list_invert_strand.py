import sys
import os

input_file_path = sys.argv[1]

output_file_path = sys.argv[2]

input_interval_list = open(input_file_path, 'r')
output_interval_list = open(output_file_path, 'w')

for line in input_interval_list:
    if "@" in line:
        output_interval_list.write(line)
    if line[0:3] == "chr":
        line_list = line.strip().split("\t")
        if line_list[3] == "+":
            output_interval_list.write(line_list[0]+"\t"+line_list[1]+"\t"+line_list[2]+"\t"+"-"+"\t"+line_list[4]+"\n")
        elif line_list[3] == "-":
            output_interval_list.write(line_list[0]+"\t"+line_list[1]+"\t"+line_list[2]+"\t"+"+"+"\t"+line_list[4]+"\n")
        else:
            output_interval_list.write(line_list[0]+"\t"+line_list[1]+"\t"+line_list[2]+"\t"+line_list[3]+"\t"+line_list[4]+"\n")

input_interval_list.close()
output_interval_list.close()
