#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

# Sample Input:
#      3 5
#      GGCGTTCAGGCA
#      AAGAATCAGTCA
#      CAAGGAGTTCGC
#      CACGTCAATCAC
#      CAATAATATTCG

# Sample Output:
#      CAG
#      CAG
#      CAA
#      CAA
#      CAA

# puts BioInfoAlgos.new.greedy_motif_search(["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5) 

input_lines = IO.readlines(ARGV[0])

k = input_lines[0].chomp.split(" ")[0].to_i
t = input_lines[0].chomp.split(" ")[1].to_i

dna = []
(1..(input_lines.length-1)).each {|i| dna << input_lines[i].chomp}

puts BioInfoAlgos.new.greedy_motif_search(dna, k, t)
