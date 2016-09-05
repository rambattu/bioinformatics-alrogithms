#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
#      Input: An amino acid string Peptide.
#      Output: Cyclospectrum(Peptide).

# CODE CHALLENGE: Solve the Generating Theoretical Spectrum Problem.

# Sample Input:
#      LEQN

# Sample Output:
#      0 113 114 128 129 227 242 242 257 355 356 370 371 484

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.cyclic_spectrum(input_lines[0].chomp).join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("LQEN").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("TMLA").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("IAMT").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("TALM").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("TLAM").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("ALTM").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("TMIA").join(" ")
puts BioInfoAlgos.new.linear_spectrum("VAQ").join(" ")
puts BioInfoAlgos.new.linear_spectrum("ETC").join(" ")
puts BioInfoAlgos.new.linear_spectrum("TCQ").join(" ")
puts BioInfoAlgos.new.linear_spectrum("CTV").join(" ")
puts BioInfoAlgos.new.linear_spectrum("CET").join(" ")
puts BioInfoAlgos.new.linear_spectrum("AVQ").join(" ")
# puts BioInfoAlgos.new.cyclic_spectrum("NQEL").join(" ")
