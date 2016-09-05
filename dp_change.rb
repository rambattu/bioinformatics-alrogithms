#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
    
CODE CHALLENGE: Solve the Change Problem. The DPCHANGE pseudocode is reproduced below for your convenience.
     Input: An integer money and an array Coins = (coin1, ..., coind).
     Output: The minimum number of coins with denominations Coins that changes money.

Sample Input:
     40
     50,25,20,10,5,1

Sample Output:
     2

=end

input_lines = IO.readlines(ARGV[0])
money = input_lines[0].chomp.to_i
coins = input_lines[1].chomp.split(",").map {|i| i.to_i}
puts BioInfoAlgos.new.dp_change(money,coins)
