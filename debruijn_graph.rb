#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
4
AAGATTCTCTAAGA
Sample Output:
AAG -> AGA,AGA
AGA -> GAT
ATT -> TTC
CTA -> TAA
CTC -> TCT
GAT -> ATT
TAA -> AAG
TCT -> CTA,CTC
TTC -> TCT
=end

input_lines = IO.readlines(ARGV[0])
k = input_lines[0].chomp.to_i
text = input_lines[1].chomp
graph = BioInfoAlgos.new.debruijn_graph(k, text)

graph.keys.sort.each do |k|
    # puts k
    # puts graph[k].length
    list = graph[k].sort.join(",")

    puts k + " -> " + list
end

# puts BioInfoAlgos.new.overlap_graph(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"])
