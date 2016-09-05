#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Solve the Minimum Skew Problem.

# Sample Input:
#      TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT

# Sample Output:
#      11 24

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.minimum_skew(input_lines[0].chomp).join(" ")
puts BioInfoAlgos.new.minimum_skew("CATTCCAGTACTTCATGATGGCGTGAAGA").join(" ")
