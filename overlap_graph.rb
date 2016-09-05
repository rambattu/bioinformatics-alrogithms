#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
ATGCG
GCATG
CATGC
AGGCA
GGCAT
Sample Output:
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG
=end

input_lines = IO.readlines(ARGV[0])
input_lines.each {|str| str.chomp!}
graph = BioInfoAlgos.new.overlap_graph(input_lines)

graph.keys.each do |k|
    list = graph[k].join(",")
    puts k + " -> " + list
end

# puts BioInfoAlgos.new.overlap_graph(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"])
