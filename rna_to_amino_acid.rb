#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Protein Translation Problem: Translate an RNA string into an amino acid string.
#      Input: An RNA string Pattern and the array GeneticCode.
#      Output: The translation of Pattern into an amino acid string Peptide.

# CODE CHALLENGE: Solve the Protein Translation Problem.

# Sample Input:
#      AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

# Sample Output:
#      MAMAPRTEINSTRING

# input_lines = IO.readlines(ARGV[0])
# puts BioInfoAlgos.new.rna_to_amino_acid(input_lines[0].chomp)
puts BioInfoAlgos.new.rna_to_amino_acid("CCAAGAACAGAUAUCAAU")
puts BioInfoAlgos.new.rna_to_amino_acid("CCUCGUACAGAAAUCAAC")
puts BioInfoAlgos.new.rna_to_amino_acid("CCUCGUACUGAUAUUAAU")
puts BioInfoAlgos.new.rna_to_amino_acid("CCCAGGACUGAGAUCAAU")
