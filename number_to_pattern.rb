#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Implement NumberToPattern.
#      Input: Integers index and k.
#      Output: The string NumberToPattern(index, k). 

# Sample Input:
# 45
# 4
# Sample Output:
# AGTC

# puts BioInfoAlgos.new.number_to_pattern(45, 4)
puts BioInfoAlgos.new.number_to_pattern(5437, 8)

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.number_to_pattern(input_lines[0].chomp.to_i,  input_lines[1].chomp.to_i)
