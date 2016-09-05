#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

#  Input: A string Text and an integer k.
#  Output: All most frequent k-mers in Text.

# Sample Input:
#      ACGTTGCATGTCGCATGATGCATGAGAGCT
#      4

# Sample Output:
#      CATG GCAT

# puts BioInfoAlgos.new.faster_frequent_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4).join(" ")

input_lines = IO.readlines(ARGV[0])
puts BioInfoAlgos.new.faster_frequent_words(input_lines[0].chomp.to_s,  input_lines[1].chomp.to_i).join(" ")
