#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
#
# Usage ram_usage.py <input_file> <output_file>
#
# Computes a table with RAM-Usage from a file containing output of `time -v`.
import argparse
import os
import pandas

parser = argparse.ArgumentParser(description='Compute RAM-Usage', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', type=str, help='File containing all outputs of `time -v`.')
parser.add_argument('output', type=str, help='File to write output to.')
arguments = parser.parse_args()

file_names = []
ram_usages = []

with open(arguments.input, 'r') as input_file:
    parsing_ram_usage = False
    for line_number, line in enumerate(input_file):
        if line_number % 23 == 0:
            index_of_unit = line.rfind('unit')
            if index_of_unit != - 1:
                parsing_ram_usage = True
                file_names.append(line[index_of_unit:][:-2])
            else:
                parsing_ram_usage = False
        if parsing_ram_usage and ((line_number - 9) % 23) == 0:
            ram_usages.append(int(line.split(' ')[-1]) // 1024)

with open(arguments.output, 'w') as output_file:
    df = pandas.DataFrame({'File' : file_names, 'RAM in MiB' : ram_usages})
    df = df.sort_values(by=['RAM in MiB', 'File'], ascending=False)
    df.to_csv(output_file, index=False)
