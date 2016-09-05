#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

# CODE CHALLENGE: Implement LEADERBOARDCYCLOPEPTIDESEQUENCING.
#      Input: An integer N and a collection of integers Spectrum.
#      Output: LeaderPeptide after running LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N).

# Sample Input:
#      10
#      0 71 113 129 147 200 218 260 313 331 347 389 460

# Sample Output:
#      113-147-71-129

input_lines = IO.readlines(ARGV[0])
n = input_lines[0].chomp.to_i
spectrum_a = input_lines[1].chomp.split(" ")
spectrum_i_a = spectrum_a.map {|i| i.to_i}
puts BioInfoAlgos.new.leaderboard_cyclopeptide_sequencing(spectrum_i_a, n)
# puts BioInfoAlgos.new.leaderboard_cyclopeptide_sequencing([0,71,113,129,147,200,218,260,313,331,347,389,460], 10)
