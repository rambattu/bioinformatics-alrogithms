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
# # puts BioInfoAlgos.new.profile_motifs(motifs)
# profile =  BioInfoAlgos.new.profile_motifs(motifs)


# #  Pr(TCGGGGATTTCC | Profile) = 0.7 · 0.6 · 1.0 · 1.0 · 0.9 · 0.9 · 0.9 · 0.5 · 0.8 · 0.7 · 0.4 · 0.6 = 0.0205753  

# puts BioInfoAlgos.new.dna_profile_probability("TCGTGGATTTCC", profile)


profile = Matrix[
  [0.4,  0.3,  0.0,  0.1,  0.0,  0.9],
  [0.2,  0.3,  0.0,  0.4,  0.0,  0.1],
  [0.1,  0.3,  1.0,  0.1,  0.5,  0.0],
  [0.3,  0.1,  0.0,  0.4,  0.5,  0.0]
]
puts BioInfoAlgos.new.dna_profile_probability("GAGCTA", profile)
