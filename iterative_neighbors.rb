#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Input: ACG 
#        1
# Output: CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG


input_lines = IO.readlines(ARGV[0])
puts BioInfoAlgos.new.iterative_neighbors(input_lines[0].chomp, input_lines[1].chomp.to_i)

# puts BioInfoAlgos.new.iterative_neighbors("ACG", 1).join(" ")
