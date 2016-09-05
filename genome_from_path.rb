#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT
Sample Output:
ACCGAAGCT
=end

input_lines = IO.readlines(ARGV[0])
input_lines.each {|str| str.chomp!}
puts BioInfoAlgos.new.genome_from_path(input_lines)

# puts BioInfoAlgos.new.genome_from_path(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"])
