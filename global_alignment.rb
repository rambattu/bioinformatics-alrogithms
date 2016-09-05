#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Solve the Global Alignment Problem.
     Input: Two protein strings written in the single-letter amino acid alphabet.
     Output: The maximum alignment score of these strings followed by an alignment achieving this
     maximum score. Use the BLOSUM62 scoring matrix and indel penalty σ = 5.

Sample Input:
     PLEASANTLY
     MEANLY

Sample Output:
     8
     PLEASANTLY
     -MEA--N-LY
=end

input_lines = IO.readlines(ARGV[0])
v = input_lines[0].chomp
w = input_lines[1].chomp

# puts BioInfoAlgos.new.global_alignment(v, w, 5)
BioInfoAlgos.new.global_alignment(v, w, 5)
