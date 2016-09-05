#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
Sample Input:
PairedComposition3,1(TAATGCCATGGGATGTT)

Sample Output:
    (AAT|CCA) (ATG|CAT) (ATG|GAT) (CAT|GGA) (CCA|GGG) (GCC|TGG) (GGA|GTT) (GGG|TGT) (TAA|GCC) (TGC|ATG) (TGG|ATG)
=end

# input_lines = IO.readlines(ARGV[0])
# graph = {}
# input_lines.each do |line|
#     line.chomp!
#     k_v_a = line.split("->")
#     k_v_a[0].gsub!(/\s+/,"")
#     k_v_a[1].gsub!(/\s+/,"")
#     graph[k_v_a[0]] = k_v_a[1].split(",")
# end

# graph_path = BioInfoAlgos.new.eulerian_path(graph)
# puts graph_path.join("->")

pairs =  BioInfoAlgos.new.paired_composition("TAATGCCATGGGATGTT", 3, 2)

pairs.each do |pair|
    print "(" + pair + ") "
end
puts ""