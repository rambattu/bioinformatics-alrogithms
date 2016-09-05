#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Sample Input:
#      0 113 128 186 241 299 314 427

# Sample Output:
#      186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

input_lines = IO.readlines(ARGV[0])
spectrum_a = input_lines[0].chomp.split(" ")
spectrum_i_a = spectrum_a.map {|i| i.to_i}
puts BioInfoAlgos.new.cyclopeptide_sequencing(spectrum_i_a).join(" ")
# puts BioInfoAlgos.new.cyclopeptide_sequencing([0,113,128,186,241,299,314,427]).join(" ")
