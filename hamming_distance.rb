#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Sample Input:
# GGGCCGTTGGT
# GGACCGTTGAC
# Sample Output:
# 3

input_lines = IO.readlines(ARGV[0])
puts BioInfoAlgos.new.hamming_distance(input_lines[0].chomp, input_lines[1].chomp)
# puts BioInfoAlgos.new.hamming_distance("CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT","CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG")
