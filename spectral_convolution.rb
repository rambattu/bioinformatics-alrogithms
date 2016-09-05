#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# Sample Input:
#      0 137 186 323

# Sample Output:
#      137 137 186 186 323 49

# input_lines = IO.readlines(ARGV[0])
# spectrum_a = input_lines[0].chomp.split(" ")
# spectrum_i_a = spectrum_a.map {|i| i.to_i}
# puts BioInfoAlgos.new.spectral_convolution(spectrum_i_a).join(" ")
puts BioInfoAlgos.new.spectral_convolution([0,57,118,179,236,240,301]).sort.join(" ")
