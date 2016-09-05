#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
     0 -> 2
     1 -> 3
     2 -> 1
     3 -> 0,4
     6 -> 3,7
     7 -> 8
     8 -> 9
     9 -> 6

Sample Output:
     6->7->8->9->6->3->0->2->1->3->4
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

graph_path = BioInfoAlgos.new.eulerian_path(graph)
puts graph_path.join("->")
