#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
     4

Sample Output:
     0000110010111101
=end

# input_lines = IO.readlines(ARGV[0])
# k = input_lines.shift.chomp.to_i
k = 8
puts BioInfoAlgos.new.universal_string(k)

