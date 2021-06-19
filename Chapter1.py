def Max(arr):
    max = 0
    for i in arr:
        if i > max:
            max = i
    return max


def Min(arr):
    min = arr[0]
    for i in arr:
        if i < min:
            min = i
    return min


def counting_words(text, pattern):
    text = text.upper()
    pattern = pattern.upper()
    counter = 0
    k = len(pattern)
    for i in range(0, len(text) - (k - 1)):
        if pattern == text[i:i + k]:
            counter += 1
    return counter


def frequent_words(text, k):
    freq_patterns = []
    n = len(text)
    count = []
    for i in range(0, n - k):
        pattern = text[i:i + k]
        count.append(counting_words(text, pattern))
    max_element = Max(count)  # finding max element
    for i in range(0, n - k):
        if count[i] == max_element:
            pattern = text[i:i + k]
            freq_patterns.append(pattern)
    freq_patterns = list(dict.fromkeys(freq_patterns))  # remove duplicates
    return freq_patterns


def better_frequent_words(text, k):
    text = text.upper()
    freq_patterns = []
    freq_map = frequency_map(text, k)
    max_count = Max(freq_map.values())
    for pattern in freq_map:
        if freq_map[pattern] == max_count:
            freq_patterns.append(pattern)
    return freq_patterns


def frequency_map(text, k):
    freq_map = {}
    n = len(text)
    for i in range(0, n - k):
        pattern = text[i:i + k]
        if pattern not in freq_map.keys():
            freq_map[pattern] = 1
        else:
            freq_map[pattern] = freq_map[pattern] + 1
    return freq_map


def complement(text):
    dna = ""
    for c in text.upper():
        if c == 'A':
            dna += 'T'
        elif c == 'T':
            dna += 'A'
        elif c == 'G':
            dna += 'C'
        if c == 'C':
            dna += 'G'
    return dna


def reverse(text):
    return text[::-1]


def reverse_complement(text):
    return reverse(complement(text))


def pattern_match(text, pattern):
    pattern = pattern.lower()
    text = text.lower()
    indexes = []
    k = len(pattern)
    for i in range(0, len(text) - (k - 1)):
        if pattern == text[i:i + k]:
            indexes.append(i)
    return indexes


def find_clumps(text, k, L, t):
    patterns = []
    n = len(text)
    for i in range(0, n - L):
        window = text[i:i + L]
        freq_map = frequency_map(window, k)
        for pattern in freq_map:
            if freq_map[pattern] >= t:
                patterns.append(pattern)
    patterns = list(dict.fromkeys(patterns))  # remove duplicates
    return patterns


def skew_array(dna):
    n = len(dna)
    arr = [0]
    for i in range(1, n):
        arr.append(arr[i - 1] + skew(dna[i - 1]))
    return arr


def skew(symbol):
    if symbol == 'G':
        return 1
    elif symbol == 'C':
        return -1
    else:
        return 0


def min_skew(dna):
    indexes = []
    n = len(dna)
    arr = skew_array(dna)
    m = Min(arr)
    for i in range(0, n):
        if arr[i] == m:
            indexes.append(i)
    return indexes


vibrio_cholerae = """
atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagata
ggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaa
tgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaa
gcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagcca
aagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgat
ccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctc
tattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc 
"""

test_input = """
CGGACTCGACAGATGTGAAGAACGACAATGTGAAGAC
TCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
"""

print("Counting Words: ")
print(counting_words("CGATATATCCATAG", "ATA"))
print("--------------------------------------")
print("Frequent Words: ")
print(frequent_words("ACGTTTCACGTTTTACGG", 3))
print("--------------------------------------")
print("Better Frequent Words: ")
print(better_frequent_words("ACGTTTCACGTTTTACGG", 3))
print("--------------------------------------")
print("Vibrio Cholerae 9-mer: ")
print(better_frequent_words(vibrio_cholerae, 9))
print("--------------------------------------")
print("Reverse Complement Problem: ")
print(reverse_complement("AGTCGCATAGT"))
print("--------------------------------------")
print("Pattern Matching Problem: ")
print(pattern_match(vibrio_cholerae, "ATGATCAAG"))
print("--------------------------------------")
print("Clump Finding Problem: ")
print(find_clumps(test_input, 5, 50, 4))
print("--------------------------------------")
print("Minimum Skew Problem: ")
print(min_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))

