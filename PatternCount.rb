#!/usr/bin/env ruby -w
# CODE CHALLENGE: Implement PatternCount 
#      Input: Strings Text and Pattern (File containing text and pattern)
#      Output: Count(Text, Pattern).


#         PatternCount(Text, Pattern)
#         count ← 0
#         for i ← 0 to |Text| − |Pattern|
#             if Text(i, |Pattern|) = Pattern
#                 count ← count + 1
#         return count

# Input text file format (PatternCount1.txt, PatternCount.txt)
# 1 Text line
# 2 Pattern

def pattern_count(txt, pattern)
    count = 0
    for i in 0..(txt.length - pattern.length)
        count += 1 if (txt.slice(i,pattern.length) == pattern)
    end
    return count
end

# Read the lines
input_lines = IO.readlines(ARGV[0])
# Pass in the Text line and the Pattern by removing the new lines from the input file
puts pattern_count(input_lines[0].chomp, input_lines[1].chomp)
