#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
     ATG
     ATG
     TGT
     TGG
     CAT
     GGA
     GAT
     AGA

Sample Output:
     AGA
     ATG
     ATG
     CAT
     GAT
     TGGA
     TGT
=end

input_lines = IO.readlines(ARGV[0])
input_lines.each {|l| l.chomp!}

contigs = BioInfoAlgos.new.contig_generation(input_lines)
puts contigs