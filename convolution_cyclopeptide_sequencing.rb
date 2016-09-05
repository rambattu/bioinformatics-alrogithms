#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Implement ConvolutionCyclopeptideSequencing.
#      Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
#      Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements
#      (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size
#      of Leaderboard is restricted to the top N (and ties).

# Sample Input:
#      20
#      60
#      57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493

# Sample Output:
#      99-71-137-57-72-57

input_lines = IO.readlines(ARGV[0])
m = input_lines[0].chomp.to_i
n = input_lines[1].chomp.to_i
spectrum_a = input_lines[2].chomp.split(" ")
spectrum_i_a = spectrum_a.map {|i| i.to_i}
puts BioInfoAlgos.new.convolution_cyclopeptide_sequencing(m,n,spectrum_i_a)
# puts BioInfoAlgos.new.convolution_cyclopeptide_sequencing(20,60, [0,71,113,129,147,200,218,260,313,331,347,389,460], 10)
