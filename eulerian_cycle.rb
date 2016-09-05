#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
     0 -> 3
     1 -> 0
     2 -> 1,6
     3 -> 2
     4 -> 2
     5 -> 4
     6 -> 5,8
     7 -> 9
     8 -> 7
     9 -> 6

Sample Output:
     6->8->7->9->6->5->4->2->1->0->3->2->6
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
# input_lines.each {|l| l.chomp!}

# graph = BioInfoAlgos.new.debruijn_graph_from_kmers(input_lines)

# graph.keys.sort.each do |k|
#     # puts k
#     # puts graph[k].length
#     list = graph[k].sort.join(",")

#     puts k + " -> " + list
# end
# graph = {
#     "0" => ["3"],
#     "1" => ["0"],
#     "2" => ["1","6"],
#     "3" => ["2"],
#     "4" => ["2"],
#     "5" => ["4"],
#     "6" => ["5","8"],
#     "7" => ["9"],
#     "8" => ["7"],
#     "9" => ["6"],

# }

graph_path = BioInfoAlgos.new.eulerian_cycle(graph)
puts graph_path
