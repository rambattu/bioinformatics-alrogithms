#!/usr/bin/env ruby -w

#  Input: A string Text and an integer k.
#  Output: All most frequent k-mers in Text.

# Sample Input:
#      ACGTTGCATGTCGCATGATGCATGAGAGCT
#      4

# Sample Output:
#      CATG GCAT

# FrequentWords(Text, k)
#     FrequentPatterns ← an empty set
#     for i ← 0 to |Text| − k
#         Pattern ← the k-mer Text(i, k)
#         Count(i) ← PatternCount(Text, Pattern)
#     maxCount ← maximum value in array Count
#     for i ← 0 to |Text| − k
#         if Count(i) = maxCount
#             add Text(i, k) to FrequentPatterns
#     remove duplicates from FrequentPatterns  
#     return FrequentPatterns

# Copied from PatternCount.rb, should find a better way to share code
def pattern_count(txt, pattern)
    count = 0
    for i in 0..(txt.length - pattern.length)
        count += 1 if (txt.slice(i,pattern.length) == pattern)
    end
    return count
end

def frequent_words(txt, k)
    max_count = 0
    count_a = []
    for i in 0..(txt.length - k)
        count_a[i] = pattern_count(txt, txt.slice(i, k))
        max_count = count_a[i] if max_count < count_a[i]
    end
    freq_patterns = []
    for i in 0..(txt.length - k)
        freq_patterns << txt.slice(i,k) if count_a[i] == max_count
    end
    return freq_patterns.uniq.reverse
end

input_lines = IO.readlines(ARGV[0])
puts frequent_words(input_lines[0].chomp, input_lines[1].chomp.to_i).join(" ")
