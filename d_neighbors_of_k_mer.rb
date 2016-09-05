#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# The d-neighborhood of the k-mer Pattern is the collection of all k-mers 
# that are at most Hamming distance d from Pattern.
#  How many 4-mers are in the 3-neighborhood of Pattern = TAGC? 
# Note that the d-neighborhood of Pattern includes Pattern.

puts BioInfoAlgos.new.d_neighbors_of_k_mer("TGCA",4,3)
