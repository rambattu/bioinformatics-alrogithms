#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
CODE CHALLENGE: Solve the String Composition Problem.
     Input: An integer k and a string Text.
     Output: Compositionk(Text), where the k-mers are written in lexicographic order.

Sample Input:
     5
     CAATCCAAC

Sample Output:
     AATCC
     ATCCA
     CAATC
     CCAAC
     TCCAA
=end

input_lines = IO.readlines(ARGV[0])
k = input_lines[0].chomp.to_i
text = input_lines[1].chomp
puts BioInfoAlgos.new.string_composition(k, text)

# puts BioInfoAlgos.new.string_composition(5, "CAATCCAAC")
