#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Clump Finding Problem: Find patterns forming clumps in a string.
#      Input: A string Genome, and integers k, L, and t.
#      Output: All distinct k-mers forming (L, t)-clumps in Genome.

# Sample Input:
#      CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
#      5 50 4

# Sample Output:
#      CGACA GAAGA   

# puts BioInfoAlgos.new.find_clump("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4).join(" ")

# input_lines = IO.readlines(ARGV[0])
# genome = input_lines[0].chomp.to_s
# next_line = input_lines[1].chomp.split(" ")
# k = next_line[0].to_i
# l = next_line[1].to_i
# t = next_line[2].to_i
# puts BioInfoAlgos.new.find_clump(genome, k, l, t).join(" ")

# input_lines = IO.readlines(ARGV[0])
# genome = input_lines[0].chomp.to_s
# next_line = input_lines[1].chomp.split(" ")
genome = "GCACAAGGCCGACAATAGGACGTAGCCTTGAAGACGACGTAGCGTGGTCGCATAAGTACAGTAGATAGTACCTCCCCCGCGCATCCTATTATTAAGTTAATT"
k = 4
l = 30
t = 3
puts BioInfoAlgos.new.find_clump(genome, k, l, t).join(" ")
# puts BioInfoAlgos.new.find_clump(genome, k, l, t).uniq.count
