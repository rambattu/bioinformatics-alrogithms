#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Solve the Approximate Pattern Matching Problem.

# Sample Input:
#      ATTCTGGA
#      CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
#      3

# Sample Output:
#      6 7 26 27

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.approx_pattern_matching(input_lines[0].chomp, input_lines[1].chomp, input_lines[2].to_i ).join(" ")
# puts BioInfoAlgos.new.approx_pattern_matching("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3).join(" ")
puts BioInfoAlgos.new.approx_pattern_matching("AAAAA", "AACAAGCTGATAAACATTTAAAGAG", 2).length
