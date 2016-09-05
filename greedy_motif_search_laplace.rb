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
     # TTC
     # ATC
     # TTC
     # ATC
     # TTC

# puts BioInfoAlgos.new.greedy_motif_search_laplace(["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5) 

input_lines = IO.readlines(ARGV[0])

k = input_lines[0].chomp.split(" ")[0].to_i
t = input_lines[0].chomp.split(" ")[1].to_i

dna = []
(1..(input_lines.length-1)).each {|i| dna << input_lines[i].chomp}

puts BioInfoAlgos.new.greedy_motif_search_laplace(dna, k, t)
