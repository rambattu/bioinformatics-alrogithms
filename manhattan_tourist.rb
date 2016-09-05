#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
    
CODE CHALLENGE: Find the length of a longest path in the Manhattan Tourist Problem.
     Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right.
     The two matrices are separated by the - symbol.
     Output: The length of a longest path from source (0, 0) to sink (n, m) in the n × m rectangular grid
     whose edges are defined by the matrices Down and Right.

Sample Input:
     4 4
     1 0 2 4 3
     4 6 5 2 1
     4 4 5 2 1
     5 6 8 5 3
     -
     3 2 4 0
     3 2 4 2
     0 7 3 3
     3 3 0 2
     1 3 2 2

Sample Output:
     34
=end

ip_line_index = 0
input_lines = IO.readlines(ARGV[0])
(n,m) = input_lines[ip_line_index].chomp.split(" ").map {|i| i.to_i}
ip_line_index += 1
# puts m.to_s + " " + n.to_s

down = []
(1..n).each do |j|
    down << input_lines[ip_line_index].chomp.split(" ").map {|i| i.to_i}
    ip_line_index += 1
end

ip_line_index += 1 # ignore the hypen line

right = []
(1..(n+1)).each do |j|
    right << input_lines[ip_line_index].chomp.split(" ").map {|i| i.to_i}
    ip_line_index += 1
end

# puts down.join(" ")
# puts right.join(" ")
puts BioInfoAlgos.new.manhattan_tourist(n, m, down, right)
