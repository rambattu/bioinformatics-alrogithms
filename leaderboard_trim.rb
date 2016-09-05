#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Implement Trim (reproduced below).
#      Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
#      Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.

# Sample Input:
# LAST ALST TLLT TQAS
# 0 71 87 101 113 158 184 188 259 271 372
# 2
# Sample Output:
# LAST ALST

input_lines = IO.readlines(ARGV[0])
lb = input_lines[0].chomp.split(" ")
spectrum_a = input_lines[1].chomp.split(" ")
spectrum_i_a = spectrum_a.map {|i| i.to_i}
n = input_lines[2].chomp.to_i
# puts lb.join(" ")
# puts spectrum_i_a.join(" ")
# puts n.to_s
puts BioInfoAlgos.new.leaderboard_trim(lb, spectrum_i_a, n).join(" ")
# puts BioInfoAlgos.new.leaderboard_trim(["LAST","ALST","TLLT","TQAS"],[0,71,87,101,113,158,184,188,259,271,372],2)
