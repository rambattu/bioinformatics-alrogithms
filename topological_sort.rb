#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Implement TOPOLOGICALORDERING.
     Input: The adjacency list of a graph (with nodes represented by integers).
     Output: A topological ordering of this graph.

Sample Input: 
0 -> 1
1 -> 2
3 -> 1
4 -> 2

Sample Output: 
0, 3, 4, 1, 2
=end

input_lines = IO.readlines(ARGV[0])
graph = {}
input_lines.each do |line|
    line.chomp!
    k_v_a = line.split("->")
    k_v_a[0].gsub!(/\s+/,"")
    k_v_a[1].gsub!(/\s+/,"")
    graph[k_v_a[0]] = k_v_a[1].split(",")
end

puts graph
topo_order = BioInfoAlgos.new.topological_ordering(graph)
puts topo_order.join(", ")
