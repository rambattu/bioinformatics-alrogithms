#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

# min_distance_pattern_text(GATTCTCA, GCAAAGACGCTGACCAA) = 3

# puts BioInfoAlgos.new.min_distance_pattern_text("GATTCTCA", "GCAAAGACGCTGACCAA") 


# ttaccttAAc  1
# gAtAtctgtc  1
# Acggcgttcg  2
# ccctAAAgag  0
# cgtcAgAggt  1

# For example, for the strings Dna shown below, d(AAA, Dna) = 1 + 1 + 2 + 0 + 1 = 5.

# puts BioInfoAlgos.new.min_distance_pattern_dna("AAA", ["ttaccttAAc", "gAtAtctgtc", "Acggcgttcg", "ccctAAAgag", "cgtcAgAggt"]) 


# Sample Input:
#      3
#      AAATTGACGCAT
#      GACGACCACGTT
#      CGTCAGCGCCTG
#      GCTGAGCACCGG
#      AGTACGGGACAG

# Sample Output:
#      GAC

# puts BioInfoAlgos.new.median_string(["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC", "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"], 7) 
puts BioInfoAlgos.new.median_strings(["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC", "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"], 7) 

# input_lines = IO.readlines(ARGV[0])

# k = input_lines[0].chomp.to_i

# dna = []
# (1..(input_lines.length-1)).each {|i| dna << input_lines[i].chomp}

# puts BioInfoAlgos.new.median_string(dna, k)
