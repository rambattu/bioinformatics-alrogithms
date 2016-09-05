#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
#      Input: A collection of strings Dna, and integers k and d.
#      Output: All (k, d)-motifs in Dna.
# 3 1
# ATTTGGC
# TGCCTTA
# CGGTATC
# GAAAATT


input_lines = IO.readlines(ARGV[0])
k = input_lines[0].chomp.split(" ")[0].to_i
d = input_lines[0].chomp.split(" ")[1].to_i
dna = []
(1..(input_lines.length-1)).each do |i|
    dna << input_lines[i].chomp
end
puts BioInfoAlgos.new.motif_enumeration(dna, k, d).join(" ")

# puts BioInfoAlgos.new.motif_enumeration(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"], 3, 1).join(" ")
