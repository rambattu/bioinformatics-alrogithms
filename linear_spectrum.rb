#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Implement LinearSpectrum.
#      Input: An amino acid string Peptide.
#      Output: The linear spectrum of Peptide.

# Extra Dataset

# Sample Input:
# NQEL
# Sample Output:
# 0 113 114 128 129 242 242 257 370 371 484

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.linear_spectrum(input_lines[0].chomp).join(" ")
# puts BioInfoAlgos.new.linear_spectrum("NQEL").join(" ")
puts BioInfoAlgos.new.linear_spectrum("N").join(" ")
