#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
CODE CHALLENGE: Implement GREEDYSORTING.
     Input: A permutation P.
     Output: The sequence of permutations corresponding to applying GREEDYSORTING to P, ending with
     the identity permutation.

Sample Input:
     (-3 +4 +1 +5 -2)

Sample Output:
     (-1 -4 +3 +5 -2)
     (+1 -4 +3 +5 -2)
     (+1 +2 -5 -3 +4)
     (+1 +2 +3 +5 +4)
     (+1 +2 +3 -4 -5)
     (+1 +2 +3 +4 -5)
     (+1 +2 +3 +4 +5)
=end

input_lines = IO.readlines(ARGV[0])[0]
input_lines.gsub!(/\(/,"")
input_lines.gsub!(/\)/,"")
p = input_lines.split(" ").map {|i| i.to_i}
# puts p
sorted_lines = BioInfoAlgos.new.greedy_sorting(p)
# puts sorted_lines.join("   ")
sorted_lines.each do |line|
    str = ""
    str += "("
    line.each do |i|
        # puts i
        if i > 0
            str += "+" + i.to_s
        else
            str += i.to_s
        end
        str += " "
    end 
    str.chop!
    str += ")"
    puts str           
end