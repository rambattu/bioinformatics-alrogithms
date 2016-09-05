#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Sample Input:
# NQEL
# 0 99 113 114 128 227 257 299 355 356 370 371 484
# Sample Output:
# 11

input_lines = IO.readlines(ARGV[0])
peptide = input_lines[0].chomp
spectrum_a = input_lines[1].chomp.split(" ")
spectrum_i_a = spectrum_a.map {|i| i.to_i}
puts BioInfoAlgos.new.cyclopeptide_score(peptide, spectrum_i_a)
# puts BioInfoAlgos.new.cyclopeptide_score("NQEL",[0,99,113,114,128,227,257,299,355,356,370,371,484])
