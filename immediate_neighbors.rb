#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Input: ACG
# Output: CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG

# # input_lines = IO.readlines(ARGV[0])
# # puts BioInfoAlgos.new.approx_pattern_count(input_lines[0].chomp, input_lines[1].chomp, input_lines[2].to_i )
# puts BioInfoAlgos.new.approx_pattern_count("CGTGACAGTGTATGGGCATCTTT", "TGT", 1)

puts BioInfoAlgos.new.immediate_neighbors("ACG").join(" ")
