#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Solve the Peptide Encoding Problem.

# Sample Input:
#      ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
#      MA

# Sample Output:
#      ATGGCC
#      GGCCAT
#      ATGGCC

input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.encode_peptide(input_lines[0].chomp, input_lines[1].chomp).join(" ")
puts BioInfoAlgos.new.encode_peptide(input_lines[0].chomp, input_lines[1].chomp).length
# puts BioInfoAlgos.new.encode_peptide("ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA", "MA").join(" ")
