#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

# Function calcuate the score of the given set of motifs
# This is the score of the number of unpopular lower case letters in the motif matrix
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

# input_lines = IO.readlines(ARGV[0])

motifs = Matrix[
           ["T",   "C",   "G",   "G",   "G",   "G",   "g",   "T",   "T",   "T",   "t",   "t"],           
           ["c",   "C",   "G",   "G",   "t",   "G",   "A",   "c",   "T",   "T",   "a",   "C"],
           ["a",   "C",   "G",   "G",   "G",   "G",   "A",   "T",   "T",   "T",   "t",   "C"],
           ["T",   "t",   "G",   "G",   "G",   "G",   "A",   "c",   "T",   "T",   "t",   "t"],
           ["a",   "a",   "G",   "G",   "G",   "G",   "A",   "c",   "T",   "T",   "C",   "C"],
           ["T",   "t",   "G",   "G",   "G",   "G",   "A",   "c",   "T",   "T",   "C",   "C"],
           ["T",   "C",   "G",   "G",   "G",   "G",   "A",   "T",   "T",   "c",   "a",   "t"],
           ["T",   "C",   "G",   "G",   "G",   "G",   "A",   "T",   "T",   "c",   "C",   "t"],
           ["T",   "a",   "G",   "G",   "G",   "G",   "A",   "a",   "c",   "T",   "a",   "C"],
           ["T",   "C",   "G",   "G",   "G",   "t",   "A",   "T",   "a",   "a",   "C",   "C"],
    ]
puts BioInfoAlgos.new.score_motifs(motifs)
