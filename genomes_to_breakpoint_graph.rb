#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Implement ChromosomeToCycle.
     Input: A chromosome Chromosome containing n synteny blocks.
     Output: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle
     to Chromosome.

Extra Dataset

Sample Input:
(+1 -2 -3 +4)
Sample Output:
(1 2 4 3 6 5 7 8)
=end

# input_lines = IO.readlines(ARGV[0])[0]
# input_lines = "(+1 -2 -3 +4)"
# input_lines = "(+1 +2 +3 +4 +5 +6)"
# input_lines.gsub!(/\(/,"")
# input_lines.gsub!(/\)/,"")
# p = input_lines.split(" ").map {|i| i.to_i}
# # puts p
# seq = BioInfoAlgos.new.chromosome_to_cyle(p)
# puts seq.join(" ")

# input_lines = "(+1 -3 -6 -5)"
# input_lines.gsub!(/\(/,"")
# input_lines.gsub!(/\)/,"")
# p = input_lines.split(" ").map {|i| i.to_i}
# # puts p
# seq = BioInfoAlgos.new.chromosome_to_cyle(p)
# puts seq.join(" ")


# input_lines = "(+2 -4)"
# input_lines.gsub!(/\(/,"")
# input_lines.gsub!(/\)/,"")
# p = input_lines.split(" ").map {|i| i.to_i}
# # puts p
# seq = BioInfoAlgos.new.chromosome_to_cyle(p)
# puts seq.join(" ")

=begin

CODE CHALLENGE: Implement CycleToChromosome.
     Input: A sequence Nodes of integers between 1 and 2n.
     Output: The chromosome Chromosome containing n synteny blocks resulting from applying
     CycleToChromosome to Nodes.

Sample Input:
(1 2 4 3 6 5 7 8)
Sample Output:
(+1 -2 -3 +4)

=end

# input_lines = "(1 2 4 3 6 5 8 7 9 10 12 11)"
# input_lines.gsub!(/\(/,"")
# input_lines.gsub!(/\)/,"")
# p = input_lines.split(" ").map {|i| i.to_i}
# # puts p
# seq = BioInfoAlgos.new.cycle_to_chromosome(p)
# puts seq.join(" ")


=begin
CODE CHALLENGE: Implement ColoredEdges.
     Input: A genome P.
     Output: The collection of colored edges in the genome graph of P in the form (x, y).

Sample Input:
(+1 -2 -3)(+4 +5 -6)
Sample Output:
(2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)
=end


# # input_lines = "(+1 -2 -3)(+4 +5 -6)"
# input_lines = "(+1 +2 +3 +4 +5 +6)"
# chromosomes = input_lines.split(")(")
# chrs = []
# chromosomes.each do |chrome|
#     chrome.gsub!(/\(/,"")
#     chrome.gsub!(/\)/,"")
#     p = chrome.split(" ").map {|i| i.to_i}
#     chrs << p.dup
# end
# # puts p
# seq = BioInfoAlgos.new.colored_edges(chrs)
# puts seq.join(" ")

# input_lines = "(+1 -3 -6 -5)(+2 -4)"
# chromosomes = input_lines.split(")(")
# chrs = []
# chromosomes.each do |chrome|
#     chrome.gsub!(/\(/,"")
#     chrome.gsub!(/\)/,"")
#     p = chrome.split(" ").map {|i| i.to_i}
#     chrs << p.dup
# end
# # puts p
# seq = BioInfoAlgos.new.colored_edges(chrs)
# puts seq.join(" ")


=begin
CODE CHALLENGE: Implement 2-BreakOnGenomeGraph.
     Input: The colored edges of a genome graph GenomeGraph, followed by indices i, i', j, and j'.
     Output: The colored edges of the genome graph resulting from applying the 2-break operation
     2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′).

Extra Dataset

Sample Input:
(2, 4), (3, 8), (7, 5), (6, 1)
1, 6, 3, 8
Sample Output:
(2, 4), (3, 1), (7, 5), (6, 8)    
=end

# genome_graph = [[2,4], [3,8], [7,5], [6,1]]
# out_graph = BioInfoAlgos.new.two_break_on_genome_graph(genome_graph, 1, 6, 3, 8)
# out_graph.each do |node|
#     print node.join(", ") + " ; "
# end
# puts ""



=begin
CODE CHALLENGE: Implement 2-BreakOnGenome.
     Input: A genome P, followed by indices i, i', j, and j'.
     Output: The genome P' resulting from applying the 2-break operation
     2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′).

Sample Input:
(+1 -2 -4 +3)
1, 6, 3, 8
Sample Output:
(+1 -2)(-3 +4)
=end

# input_lines = IO.readlines(ARGV[0])[0]
# input_lines.gsub!(/\(/,"")
# input_lines.gsub!(/\)/,"")
# p = input_lines.split(" ").map {|i| i.to_i}
# puts p
# p = [[+1, -2, -4, +3]]
# genome = BioInfoAlgos.new.two_break_on_genome([p], 116, 118, 73, 72)
# genome.each {|p| print "(" + p.join(" ") +")" }
# puts ""


=begin
CODE CHALLENGE: Implement GraphToGenome.
     Input: The colored edges ColoredEdges of a genome graph.
     Output: The genome P corresponding to this genome graph.

Extra Dataset

Sample Input:
(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)
Sample Output:
(+1 -2 -3)(-4 +5 -6)
=end
colored_edges = [[2,4], [3,6], [5,1], [7,9], [10,12], [11,8]]
genome = BioInfoAlgos.new.graph_to_genome(colored_edges)
genome.each do |chrome|
     print "(" + chrome.join(", ") + ")"
end
puts ""
