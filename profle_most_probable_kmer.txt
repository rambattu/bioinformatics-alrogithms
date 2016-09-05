#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

# Sample Input:
#      ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
#      5
#      0.2 0.2 0.3 0.2 0.3
#      0.4 0.3 0.1 0.5 0.1
#      0.3 0.3 0.5 0.2 0.4
#      0.1 0.2 0.1 0.1 0.2

# Sample Output:
#      CCGAG

input_lines = IO.readlines(ARGV[0])
text = input_lines[0].chomp
k = input_lines[1].chomp.to_i
matrix = []
(2..(input_lines.input_lines-1)).each do |i|
  row_str = input_lines[i].chomp.to_s
  row = []
  row_str.split(" ").each {|val| row << val.to_f}
  matrix << row
end

puts BioInfoAlgos.new.profile_most_probable_kmer(text, k, Matrix[matrix])
