import random


def create_random_kmer(k):
    nucleotids = ['A', 'T', 'G', 'C']
    kmer = ""
    for _ in range(0, k):
        kmer += random.choice(nucleotids)

    return kmer


def d_func(pattern, dna):
    total_score = 0
    k = len(pattern)
    for d in dna:
        score_list = []
        for i in range(0, len(d) - k):
            score = 0
            index = 0
            for j in d[i:i + k]:
                if j != pattern[index]:
                    score += 1
                index += 1
            score_list.append(score)
        total_score += min(score_list)
    return total_score


def base4_convert(i, k):
    result = []
    index = 0
    while i > 0:
        result.append(i % 4)
        i = i // 4
        index += 1
    if len(result) < k:
        for _ in range(index, k):
            result.append(0)
    return result[::-1]


def median_string_problem(dna, k):
    t = len(dna)
    L = len(dna[0])
    kmers = []
    distance = 999999
    median = ""
    for i in range(0, L ** t):
        if i < 4 ** k:
            # print(i)
            base_list = base4_convert(i, k)
            # print(base_list)
            kmer = ""
            for n in base_list:
                if n == 0:
                    kmer += "A"
                elif n == 1:
                    kmer += "C"
                elif n == 2:
                    kmer += "G"
                elif n == 3:
                    kmer += "T"
                else:
                    print("There is an error check line 60!")
            kmers.append(kmer)
    for pattern in kmers:
        if d_func(pattern, dna) <= distance:
            median = pattern
            distance = d_func(pattern, dna)
    return median


def most_frequent(list):
    dict = {}
    counter = 0
    tmp_item = ''
    for item in reversed(list):
        dict[item] = dict.get(item, 0) + 1
        if dict[item] >= counter:
            counter, tmp_item = dict[item], item
    return tmp_item


def calculate_score_helper(frequent, all):
    score = 0
    for c in all:
        if frequent != c:
            score += 1
    return score


def calculate_score(dna):
    l = []
    score = 0
    for i in range(0, len(dna[0])):
        for d in dna:
            l.append(d[i])
        score += calculate_score_helper(most_frequent(l), l)
        l.clear()
    return score


def find_profile_score_helper(dna_list, p_score):
    l = len(dna_list)
    a = 0
    t = 0
    c = 0
    g = 0
    for col in dna_list:
        if col == "A":
            a += 1
        elif col == "T":
            t += 1
        elif col == "C":
            c += 1
        elif col == "G":
            g += 1
        else:
            print("There is an error occurred please check line 99")
    p_score['A'].append(a / l)
    p_score['T'].append(t / l)
    p_score['C'].append(c / l)
    p_score['G'].append(g / l)
    # print(f"A: {a / l} -- T: {t / l} -- C: {c / l} -- G: {g / l}")


def make_profile(motif):
    p_score = {
        'A': [],
        'T': [],
        'C': [],
        'G': []
    }
    columns = {}
    l = []
    for i in range(0, len(motif[0])):
        for m in motif:
            l.append(m[i])
        find_profile_score_helper(l, p_score)
        l.clear()
    return p_score


def finding_motif(dna, pattern):
    result = []
    p_len = len(pattern)
    d_len = len(dna)
    for i in range(0, d_len - p_len + 1):
        if dna[i:i + p_len] == pattern:
            result.append(i + 1)
    return result


def base_convert(i, t, L):
    result = []
    index = 0
    while i > 0:
        result.append(i % t)
        i = i // t
        index += 1
    if len(result) < L:
        for _ in range(index, L):
            result.append(0)
    return result[::-1]


def motif_finding_problem(dna, k):
    best_score = 999999
    t = len(dna)
    L = len(dna[0]) - k
    motifs = []
    best_motif = None
    for i in range(1, L ** t):
        motifs.clear()
        index = 0
        for j in base_convert(i, L, t):
            #print(j)
            motifs.append(dna[index][j:j + k])
            index += 1
        #print(motifs)
        if calculate_score(motifs) < best_score:
            best_score = calculate_score(motifs)
            best_motif = motifs.copy()


    return best_motif


def profile_most_probable_k(text, k, matrix_profile):
    dna = [text[i:i + k] for i in range(0, len(text), k)]
    best_score = 0
    profile = ""
    for line in dna:
        score = 1
        for i in range(0, len(line)):
            if line[i] == "A":
                score *= matrix_profile['A'][i]
            elif line[i] == "C":
                score *= matrix_profile['C'][i]
            elif line[i] == "G":
                score *= matrix_profile['G'][i]
            elif line[i] == "T":
                score *= matrix_profile['T'][i]
            else:
                break
        if score > best_score:
            best_score = score
            profile = line
    return profile


def profile(motifs):
    profile = []
    columns = []
    for i in range(len(motifs[0])):
        for j in range(len(motifs)):
            columns.append(motifs[j][i])
        col = ''.join(columns)
        profile.append([float(col.count(n)) / float(len(col)) for n in 'ACGT'])
    return profile


def profile_most_probable_kmer(dna, k, profile):
    nucs = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    max_prob = [-1, None]
    for i in range(len(dna) - k + 1):
        current_prob = 1
        j = 0
        for n in (dna[i:i + k]):
            current_prob *= profile[j][nucs[n]]
            j += 1
        if current_prob > max_prob[0]:
            max_prob = [current_prob, dna[i:i + k]]

    return max_prob[1]


def greedy_motif_search(dna_list, k, t):
    best_score = 999999
    best_motif = None

    for i in range(len(dna_list[0]) - k + 1):
        motifs = [dna_list[0][i:i + k]]
        c_profile = profile(motifs)

        for j in range(1, len(dna_list)):
            motifs.append(profile_most_probable_kmer(dna_list[j], k, c_profile))
            c_profile = profile(motifs)

        score = calculate_score(motifs)
        if score <= best_score:
            best_score = score
            best_motif = motifs

    return best_motif


def motifs_func(dna, profile):
    k = (len(profile['A']))
    best_motifs = []
    motifs = []
    for line in dna:
        best_score = 0
        kmers = []
        for i in range(0, len(line) - k + 1):
            kmers.append(line[i:i + k])
        for kmer in kmers:
            # print(kmer)
            score = 1
            for j in range(0, len(kmer)):
                if kmer[j] == "A":
                    score *= profile['A'][j]
                elif kmer[j] == "C":
                    score *= profile['C'][j]
                elif kmer[j] == "G":
                    score *= profile['G'][j]
                elif kmer[j] == "T":
                    score *= profile['T'][j]
                else:
                    break
            if score > best_score:
                best_score = score
                motifs = kmer
        best_motifs.append(motifs)
    return best_motifs


def pseudocounts_profile_score_helper(dna_list, p_score):
    l = len(dna_list) + 1
    a = 1
    t = 1
    c = 1
    g = 1
    for col in dna_list:
        if col == "A":
            a += 1
        elif col == "T":
            t += 1
        elif col == "C":
            c += 1
        elif col == "G":
            g += 1
        else:
            print("There is an error occurred please check line 297")
    p_score['A'].append(a / l)
    p_score['T'].append(t / l)
    p_score['C'].append(c / l)
    p_score['G'].append(g / l)
    # print(f"A: {a / l} -- T: {t / l} -- C: {c / l} -- G: {g / l}")


def make_profile_pseudocounts(motif):
    p_score = {
        'A': [],
        'T': [],
        'C': [],
        'G': []
    }
    columns = {}
    l = []
    for i in range(0, len(motif[0])):
        for m in motif:
            l.append(m[i])
        pseudocounts_profile_score_helper(l, p_score)
        l.clear()
    return p_score


def randomized_motif_search(dna, k, t):
    # for seq in dna:
    #     rand_number = random.randint(0, len(seq) - k)
    #     motifs.append(seq[rand_number:rand_number + k])
    # best_motif = motifs
    best_motif = None
    best_score = 999999
    rand_ints = [random.randint(0, len(dna[0]) - k) for a in range(len(dna))]
    motifs = [dna[i][r:r + k] for i, r in enumerate(rand_ints)]
    best_score = calculate_score(motifs)
    best_motif = motifs
    while True:
        profile = make_profile_pseudocounts(motifs)
        motifs = motifs_func(dna, profile)
        current_score = calculate_score(motifs)
        if current_score <= best_score:
            best_score = current_score
            best_motif = motifs

        return best_score, best_motif


def random_iterator(dna, k, t):
    best_motif = None
    best_score = 999999
    for _ in range(1000):
        c_score, c_motif = randomized_motif_search(dna, k, t)
        if c_score <= best_score:
            best_motif = c_motif
            best_score = c_score
    return best_motif


def gibbs_sampler(dna, k, t, N):
    best_motif = None
    best_score = 99999
    rand_nums = []
    for a in range(t):
        rand_nums.append(random.randint(0, len(dna[0]) - k))
    motifs = [dna[i][r:r + k] for i, r in enumerate(rand_nums)]
    best_score = calculate_score(motifs)
    best_motif = motifs
    for j in range(N):
        r = random.randint(0, t - 1)
        c_profile = make_profile_pseudocounts([motif for index, motif in enumerate(motifs) if index != r])
        motifs = [profile_most_probable_k(dna[index], k, c_profile) if index == r else motif for index, motif in
                  enumerate(motifs)]
        c_score = calculate_score(motifs)
        if c_score <= best_score:
            best_score = c_score
            best_motif = motifs

    return best_score, best_motif


def gibbs_iterator(dna, k, t, N):
    best_score = 9999999
    best_motif = None
    for _ in range(500):
        c_score, c_motif = gibbs_sampler(dna, k, t, N)
        if c_score <= best_score:
            best_score = c_score
            best_motif = c_motif

    return best_motif


# motif_finding_problem(input, 10)
# print(median_string_problem(input, 10))


# print(motif_finding_problem(dna_sample, 8))
# print(calculate_score(string,8))


dna_motif = ["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC",
             "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"]

# make_profile(dna_motif)


matrix_profile = [[0.2, 0.2, 0.3, 0.2, 0.3], [0.4, 0.3, 0.5, 0.2, 0.4],
                  [0.3, 0.3, 0.5, 0.2, 0.4], [0.1, 0.2, 0.1, 0.1, 0.2]]
test_input = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"

# print(profile_most_probable_kmer(test_input, 5, matrix_profile))
# print(finding_motif(test_input,"CCGAG"))

test_2 = "AGCAGCTTTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATCTGAACTGGTTACCTGCCGTGAGTAAAT"

# matrix_profile2 = {'A': [0.7, 0.2, 0.1, 0.5, 0.4, 0.3, 0.2, 0.1],
# 'C': [0.2, 0.2, 0.5, 0.4, 0.2, 0.3, 0.1, 0.6],
# 'G': [0.1, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.4, 0.2],
# 'T': [0.0, 0.3, 0.2, 0.0, 0.2, 0.3, 0.3, 0.1]}

# print(profile_most_probable_kmer(test_2, 8, matrix_profile2))


# print(motifs_func(["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"], {'A': [0.8, 0, 0, 0.2],
#                                                                                            'C': [0, 0.6, 0.2, 0],
#                                                                                            'G': [0.2, 0.2, 0.8, 0],
#                                                                                            'T': [0, 0.2, 0, 0.8]}))
print("Motif Finding Problem: ")
print(motif_finding_problem(["TTACCTTAAC", "GATATCTGTC", "ACGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"], 3))

print("Median String Problem: ")
print(median_string_problem(["AAATTGACGCAT",
                             "GACGACCACGTT",
                             "CGTCAGCGCCTG",
                             "GCTGAGCACCGG",
                             "AGTACGGGACAG"], 3))


print("Greedy Motif Search: ")
print(greedy_motif_search(["GGCGTTCAGGCA",
                           "AAGAATCAGTCA",
                           "CAAGGAGTTCGC",
                           "CACGTCAATCAC",
                           "CAATAATATTCG"], 3, 5))
print("Randomized Motif Search: ")
print(random_iterator(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
                         "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                         "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
                         "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                         "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 8, 5))

print("Gibbs Sampling Motif Search: ")

print(gibbs_iterator(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 8, 5, 100))


# print(d_func("AAA", ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"]))
