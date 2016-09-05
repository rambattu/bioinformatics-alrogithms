#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Solve the Local Alignment Problem.
     Input: Two protein strings written in the single-letter amino acid alphabet.
     Output: The maximum score of a local alignment of the strings, followed by a local alignment of these
     strings achieving the maximum score. Use the PAM250 scoring matrix and indel penalty σ = 5.

Download PAM250 scoring matrix

Sample Input:
     MEANLY
     PENALTY

Sample Output:
     15
     EANL-Y
     ENALTY
=end

input_lines = IO.readlines(ARGV[0])
v = input_lines[0].chomp
w = input_lines[1].chomp

BioInfoAlgos.new.local_alignment(v, w, 5)
