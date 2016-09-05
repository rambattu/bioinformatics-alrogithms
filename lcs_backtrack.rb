#!/usr/bin/env ruby 
require_relative 'BioInfoAlgos.rb'

=begin
Sample Input:
     AACCTTGG
     ACACTGTGA
=end

backtrack = BioInfoAlgos.new.lcs_backtrack("AACCTTGG", "ACACTGTGA")
puts backtrack.join(" ")
