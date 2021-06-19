def ksorting_reversal(P):
    P.reverse()
    P = map(lambda x: -x, P)
    return P


def greedy_sorting(P):
    result = []
    approved_reversal_dist = 0

    for i in range(len(P)):
        if P[i] != i + 1:
            if i + 1 in P:
                id_i = P.index(abs(i + 1))
            else:
                id_i = P.index(-abs(i + 1))
            P[i: id_i + 1] = ksorting_reversal(P[i:id_i + 1])
            approved_reversal_dist += 1
            result.append(P[:])
        if P[i] == -(i + 1):
            P[i] = i + 1
            approved_reversal_dist += 1
            result.append(P[:])

    return result, approved_reversal_dist


def number_of_breakpoints(P):
    result = 0
    p_list = list()
    p_list.append(0)
    for p in P:
        p_list.append(p)
    p_list.append(len(P) + 1)
    # print(p_list)
    for i in range(0, len(p_list) - 1):
        if p_list[i + 1] - p_list[i] == 1:
            result += 1
    return len(p_list) - 1 - result


def reverse_complement(text):
    text.upper()
    rev_nucs = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    text = text[::-1]
    result = ""
    for t in text:
        result += rev_nucs[t]
    return result


def kmers_dict(k, text):
    kmers = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        kmers[kmer] = kmers.setdefault(kmer, []) + [i]
        kmers[reverse_complement(kmer)] = kmers[kmer]
    # print(kmers)
    return kmers


def shared_kmers(k, str1, str2):
    result = list()
    kmer2 = kmers_dict(k, str2)
    for i in range(len(str1) - k + 1):
        kmer1 = str1[i: i + k]
        if kmer1 in kmer2:
            for j in kmer2[kmer1]:
                result.append((i, j))
    result.sort()
    return result


print("Greedy Sorting: ")
result, ard = greedy_sorting([-3, +4, +1, +5, -2])
print(result)
print("Approved Reversal Distance:", ard)
print("-----------------------------------------")
print("Number Of Breakpoints Problem: ")
print(number_of_breakpoints([+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]))
print("Shared k-mers Problem: ")
print(shared_kmers(3, "AAACTCATC", "TTTCAAATC"))
