#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Solve the 2-Break Distance Problem.
     Input: Genomes P and Q.
     Output: The 2-break distance d(P, Q).

Sample Input:
     (+1 +2 +3 +4 +5 +6)
     (+1 -3 -6 -5)(+2 -4)

Sample Output:
     3
=end


input_lines = IO.readlines(ARGV[0])
# input_lines = "(+1 +2 +3 +4 +5 +6)"
chromosomes = input_lines[0].split(")(")
p_nodes = []
chromosomes.each do |chrome|
    chrome.gsub!(/\(/,"")
    chrome.gsub!(/\)/,"")
    p = chrome.split(" ").map {|i| i.to_i}
    p_nodes << p.dup
end

# input_lines = "(+1 -3 -6 -5)(+2 -4)"
chromosomes = input_lines[1].split(")(")
q_nodes = []
chromosomes.each do |chrome|
    chrome.gsub!(/\(/,"")
    chrome.gsub!(/\)/,"")
    p = chrome.split(" ").map {|i| i.to_i}
    q_nodes << p.dup
end

distance = BioInfoAlgos.new.two_break_distance(p_nodes,q_nodes)
puts distance
