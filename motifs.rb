#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

#                       A:  4/5  0    0  1/5               ttaccttaac
#             Profile   C:   0  3/5  1/5  0          Dna   gatgtctgtc
#                       G:  1/5 1/5  4/5  0                acggcgttag
#                       T:   0  1/5   0  4/5               ccctaacgag
#                                                          cgtcagaggt

# Taking the Profile-most probable 4-mer from each row of Dna produces the following 4-mers (shown in red):

#                                              ttACCTtaac
#                                              gATGTctgtc
#                        Motifs(Profile,Dna)   acgGCGTtag
#                                              ccctaACGAg
#                                              cgtcagAGGT

profile = Matrix[
  [0.8, 0, 0, 0.2], 
  [0, 0.6, 0.2, 0], 
  [0.2, 0.2, 0.8, 0], 
  [0, 0.2, 0, 0.8], 
]
puts BioInfoAlgos.new.motifs(profile, ["ttaccttaac", "gatgtctgtc", "acggcgttag", "ccctaacgag", "cgtcagaggt"], 4)
