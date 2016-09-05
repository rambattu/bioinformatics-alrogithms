#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
Number of Breakpoints Problem: Find the number of breakpoints in a permutation.
     Input: A permutation.
     Output: The number of breakpoints in this permutation.

CODE CHALLENGE: Solve the Number of Breakpoints Problem.

Sample Input:
     (+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14)

Sample Output:
     8
=end

input_lines = IO.readlines(ARGV[0])[0]
input_lines.gsub!(/\(/,"")
input_lines.gsub!(/\)/,"")
p = input_lines.split(" ").map {|i| i.to_i}
p.unshift(0)
p << (p.length)
# puts p
num = BioInfoAlgos.new.breakpoints(p)
puts num