#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Implement ComputingFrequencies to generate a frequency array.
#     Input: A DNA string Text followed by an integer k.
#     Output: FrequencyArray(Text).

# Sample Input:
# ACGCGGCTCTGAAA
# 2
# Sample Output:
# 2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

# puts BioInfoAlgos.new.frequency_array("ACGCGGCTCTGAAA", 2).join(" ")

input_lines = IO.readlines(ARGV[0])
puts BioInfoAlgos.new.frequency_array(input_lines[0].chomp.to_s,  input_lines[1].chomp.to_i).join(" ")
