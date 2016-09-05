#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
Edit Distance Problem: Find the edit distance between two strings.
     Input: Two strings.
     Output: The edit distance between these strings.

CODE CHALLENGE: Solve the Edit Distance Problem.

Sample Input:
     PLEASANTLY
     MEANLY

Sample Output:
     5
=end

input_lines = IO.readlines(ARGV[0])
v = input_lines[0].chomp
w = input_lines[1].chomp

BioInfoAlgos.new.edit_distance(v, w)
