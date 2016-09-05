#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG
Sample Output:
AGG -> GGG
CAG -> AGG,AGG
GAG -> AGG
GGA -> GAG
GGG -> GGA,GGG
=end

input_lines = IO.readlines(ARGV[0])
input_lines.each {|l| l.chomp!}

graph = BioInfoAlgos.new.debruijn_graph_from_paired_kmers(input_lines, 4, 2)

graph.keys.sort.each do |k|
    # puts k
    # puts graph[k].length
    list = graph[k].sort.join(",")

    puts k + " -> " + list
end

