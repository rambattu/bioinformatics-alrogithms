#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

# Motifs        
#    T   C   G   G   G   G   g   T   T   T   t   t           
#    c   C   G   G   t   G   A   c   T   T   a   C
#    a   C   G   G   G   G   A   T   T   T   t   C
#    T   t   G   G   G   G   A   c   T   T   t   t
#    a   a   G   G   G   G   A   c   T   T   C   C
#    T   t   G   G   G   G   A   c   T   T   C   C
#    T   C   G   G   G   G   A   T   T   c   a   t
#    T   C   G   G   G   G   A   T   T   c   C   t
#    T   a   G   G   G   G   A   a   c   T   a   C
#    T   C   G   G   G   t   A   T   a   a   C   C


# Score
#    3 + 4 + 0 + 0 + 1 + 1 + 1 + 5 + 2 + 3 + 6 + 4 = 30     


# Count
# A:   2   2   0   0   0   0   9   1   1   1   3   0          
# C:   1   6   0   0   0   0   0   4   1   2   4   6  
# G:   0   0  10  10   9   9   1   0   0   0   0   0  
# T:   7   2   0   0   1   1   0   5   8   7   3   4 

# Profile
# A:  .2  .2   0   0   0   0  .9  .1  .1  .1  .3   0            
# C:  .1  .6   0   0   0   0   0  .4  .1  .2  .4  .6  
# G:   0   0   1   1  .9  .9  .1   0   0   0   0   0  
# T:  .7  .2   0   0  .1  .1   0  .5  .8  .7  .3  .4  

# input_lines = IO.readlines(ARGV[0])

# Entropy is a measure of the uncertainty of a probability distribution (p1, …, pN), and is defined as follows:
# H(p1,…,pN)= − for i = 1 to N { pi·log2pi }
# For example, the entropy of the probability distribution (0.2, 0.6, 0.0, 0.2) corresponding to the 2nd column of the NF-κB profile matrix is

# −(0.2log(2)0.2 + 0.6log(2)0.6 + 0.0log(2)0.0 + 0.2log(2)0.2) ~ 1.371

# motifs = Matrix[
#            ["T",   "C",   "G",   "G",   "G",   "G",   "g",   "T",   "T",   "T",   "t",   "t"],           
#            ["c",   "C",   "G",   "G",   "t",   "G",   "A",   "c",   "T",   "T",   "a",   "C"],
#            ["a",   "C",   "G",   "G",   "G",   "G",   "A",   "T",   "T",   "T",   "t",   "C"],
#            ["T",   "t",   "G",   "G",   "G",   "G",   "A",   "c",   "T",   "T",   "t",   "t"],
#            ["a",   "a",   "G",   "G",   "G",   "G",   "A",   "c",   "T",   "T",   "C",   "C"],
#            ["T",   "t",   "G",   "G",   "G",   "G",   "A",   "c",   "T",   "T",   "C",   "C"],
#            ["T",   "C",   "G",   "G",   "G",   "G",   "A",   "T",   "T",   "c",   "a",   "t"],
#            ["T",   "C",   "G",   "G",   "G",   "G",   "A",   "T",   "T",   "c",   "C",   "t"],
#            ["T",   "a",   "G",   "G",   "G",   "G",   "A",   "a",   "c",   "T",   "a",   "C"],
#            ["T",   "C",   "G",   "G",   "G",   "t",   "A",   "T",   "a",   "a",   "C",   "C"],
#     ]
# puts BioInfoAlgos.new.motifs_entropy(motifs)



puts BioInfoAlgos.new.probability_entropy([0.5, 0, 0, 0.5])
puts BioInfoAlgos.new.probability_entropy([0.25, 0.25, 0.25, 0.25])
puts BioInfoAlgos.new.probability_entropy([0, 0, 0, 1])
puts BioInfoAlgos.new.probability_entropy([0.25, 0, 0.5, 0.25])
