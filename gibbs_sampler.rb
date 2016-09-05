#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

# Sample Input:
#      8 5 100
#      CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
#      GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
#      TAGTACCGAGACCGAAAGAAGTATACAGGCGT
#      TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
#      AATCCACCAGCTCCACGTGCAATGTTGGCCTA

# Sample Output:
#      TCTCGGGG
#      CCAAGGTG
#      TACAGGCG
#      TTCAGGTG
#      TCCACGTG

# puts BioInfoAlgos.new.gibbs_sampler(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 8, 5, 1000) 

input_lines = IO.readlines(ARGV[0])

k = input_lines[0].chomp.split(" ")[0].to_i
t = input_lines[0].chomp.split(" ")[1].to_i
n = input_lines[0].chomp.split(" ")[2].to_i

dna = []
(1..(input_lines.length-1)).each {|i| dna << input_lines[i].chomp}

puts BioInfoAlgos.new.gibbs_sampler(dna, k, t, n)
