#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'
require 'matrix'

=begin
Sample Input:
     4
     CTTA
     ACCA
     TACC
     GGCT
     GCTT
     TTAC

Sample Output:
     GGCTTACCA
=end

input_lines = IO.readlines(ARGV[0])
k= input_lines.shift.chomp
input_lines.each {|l| l.chomp!}

puts BioInfoAlgos.new.string_reconstruction(k, input_lines)

