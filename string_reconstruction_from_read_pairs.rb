#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
4 2
GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT
Sample Output:
GTGGTCGTGAGATGTTGA
=end

input_lines = IO.readlines(ARGV[0])
k_d = input_lines.shift.chomp.split(" ")
k = k_d[0].to_i
d = k_d[1].to_i
input_lines.each {|l| l.chomp!}

puts BioInfoAlgos.new.string_reconstruction_from_read_pairs(k, d, input_lines)

