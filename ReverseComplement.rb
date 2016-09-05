#!/usr/bin/env ruby -w

# Reverse Complement Problem: Find the reverse complement of a DNA string.
#      Input: A DNA string Pattern.
#      Output: Pattern, the reverse complement of Pattern.

# Sample Input:
#      AAAACCCGGT
# Sample Output:
#      ACCGGGTTTT

def reverse_complement(txt)
    comp = []
    txt.chars.each do |ch|
        comp << "A" if ch == 'T'
        comp << "T" if ch == 'A'
        comp << "G" if ch == 'C'
        comp << "C" if ch == 'G'
    end
    comp.join.reverse
end

# input_lines = IO.readlines(ARGV[0])
# puts reverse_complement(input_lines[0].chomp)
puts reverse_complement("TTGTGTC")