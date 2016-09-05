#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
Shared k-mers Problem: Given two strings, find all their shared k-mers.
     Input: An integer k and two strings.
     Output: All k-mers shared by these strings, in the form of ordered pairs (x, y).

CODE CHALLENGE: Solve the Shared k-mers Problem.

Sample Input:
     3
     AAACTCATC
     TTTCAAATC

Sample Output:
     (0, 4)
     (0, 0)
     (4, 2)
     (6, 6)
=end

input_lines = IO.readlines(ARGV[0])
k = input_lines[0].chomp.to_i
x = input_lines[1].chomp
y = input_lines[2].chomp

pairs = BioInfoAlgos.new.shared_kmers(k, x, y)

pairs.each {|pair| puts "(" + pair.join(", ")  + ")" }
