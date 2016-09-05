#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Solve the Longest Path in a DAG Problem.
     Input: An integer representing the source node of a graph, followed by an integer representing the
     sink node of the graph, followed by a list of edges in the graph. The edge notation 0->1:7 indicates
     that an edge connects node 0 to node 1 with weight 7. 
     Output: The length of a longest path in the graph, followed by a longest path.

Sample Input:
     0
     4
     0->1:7
     0->2:4
     2->3:2
     1->4:1
     3->4:3

Sample Output:
     9
     0->2->3->4

=end

input_lines = IO.readlines(ARGV[0])
source = input_lines.shift.chomp
sink = input_lines.shift.chomp
graph = {}
graph_with_weights = []
# puts input_lines
input_lines.each do |line|
    graph_with_weights << line
    line.chomp!
    k_v_a = line.split("->")
    k_v_a[0].gsub!(/\s+/,"")
    k_v_a[1].gsub!(/\s+/,"")
    k_v_a[1].gsub!(/:(\d+)/,"")
    graph[k_v_a[0]] = [] unless graph[k_v_a[0]]
    graph[k_v_a[0]] << k_v_a[1]
end
# puts graph_with_weights
# puts graph
(length, path) = BioInfoAlgos.new.dag_longpath(source, sink, graph_with_weights, graph)
puts length
puts path.join("->")
# puts topo_order.join(", ")
