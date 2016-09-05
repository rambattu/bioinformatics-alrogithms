#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Pattern Matching Problem: Find all occurrences of a pattern in a string.
#      Input: Two strings, Pattern and Genome.
#      Output: All starting positions where Pattern appears as a substring of Genome.

# Sample Input:
# AGT
# Sample Output:
# 11

# puts BioInfoAlgos.new.pattern_to_number("AGT")
puts BioInfoAlgos.new.pattern_to_number("ATGCAA")

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.pattern_to_number(input_lines[0].chomp)
