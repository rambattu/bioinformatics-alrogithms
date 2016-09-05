#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Pattern Matching Problem: Find all occurrences of a pattern in a string.
#      Input: Two strings, Pattern and Genome.
#      Output: All starting positions where Pattern appears as a substring of Genome.

# Sample Input:
#      ATAT
#      GATATATGCATATACTT

# Sample Output:
#      1 3 9

input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.pattern_matching("ATAT", "GATATATGCATATACTT").join(" ")
puts BioInfoAlgos.new.pattern_matching(input_lines[0].chomp, input_lines[1].chomp).join(" ")
