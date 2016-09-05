#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Use OUTPUTLCS (reproduced below) to solve the Longest Common Subsequence Problem.
     Input: Two strings s and t.
     Output: A longest common subsequence of s and t. (Note: more than one solution may exist,
     in which case you may output any one.)

Sample Input:
     AACCTTGG
     ACACTGTGA

Sample Output:
     AACTGG
=end

input_lines = IO.readlines(ARGV[0])
v = input_lines[0].chomp
w = input_lines[1].chomp

backtrack = BioInfoAlgos.new.lcs_backtrack(v, w)
# backtrack.each do |a|
#     puts a.join(" ")
#     puts ""
# end
puts BioInfoAlgos.new.output_lcs(backtrack, v, (v.length) , (w.length) )
