#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# EXERCISE BREAK: Give all values of Skewi (GAGCCACCGCGATA) for i ranging from 0 to 14.

# Sample Input:
#      CATGGGCATCGGCCATACGCC

# Sample Output:
#      0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

# input_lines = IO.readlines(ARGV[0])
puts BioInfoAlgos.new.skew("CATTCCAGTACTTCATGATGGCGTGAAGA").join(" ")
