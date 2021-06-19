# CONSTANTS
AMINO_ACIDS = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T", "AGA": "R",
    "AGC": "S",
    "AGG": "R", "AGU": "S", "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I", "CAA": "Q", "CAC": "H", "CAG": "Q",
    "CAU": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P", "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R", "CUA": "L",
    "CUC": "L",
    "CUG": "L", "CUU": "L", "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D", "GCA": "A", "GCC": "A", "GCG": "A",
    "GCU": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G", "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V", "UAA": "*",
    "UAC": "Y",
    "UAG": "*", "UAU": "Y", "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S", "UGA": "*", "UGC": "C", "UGG": "W",
    "UGU": "C",
    "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F"}
AMINO_ACIDS_MASS = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103, "I": 113, "L": 113,
                    "N": 114, "D": 115, "K": 128, "Q": 128, "E": 129, "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163,
                    "W": 186}


def pattern_translation_problem(rna):
    rna = rna.upper()
    aa = ""
    for i in range(0, len(rna) - 3, 3):
        if AMINO_ACIDS[rna[i:i + 3]] == "*":  # STOP CODONS
            break
        else:
            aa += AMINO_ACIDS[rna[i:i + 3]]
    return aa


def generating_theoretical_spectrum(peptide):
    new_peptide = peptide + peptide  # for getting rid of indexing problem in peptide string
    kmers = set()
    score_list = []
    for i in range(0, len(peptide)):
        for j in range(0, len(peptide)):
            comb_str = new_peptide[i:(i + j % len(peptide))]
            if len(comb_str) <= len(peptide):
                kmers.add(comb_str)
    kmers.add(peptide)

    for kmer in kmers:
        score = 0
        for k in kmer:
            score += AMINO_ACIDS_MASS[k]
        score_list.append(score)
    score_list.sort()
    # print(kmers)
    return score_list


def cyclopeptide_scoring(peptide, spectrum):
    theoretical_spectrum = generating_theoretical_spectrum(peptide)
    score = 0
    for ts in theoretical_spectrum:
        for s in spectrum:
            if ts == s:
                score += 1
    return score


def cyclospectrum_mass_peptide(peptide):
    spectrum = [0, sum(peptide)]
    tmp = peptide + peptide  # for getting rid of indexing problem in peptide string
    for k in range(1, len(peptide)):
        for i in range(len(peptide)):
            subpeptide = tmp[i:i + k]
            spectrum.append(sum(subpeptide))
    spectrum.sort()
    return spectrum


def score_func(peptide, spectrum):
    peptide_spectrum = cyclospectrum_mass_peptide(peptide)
    result = 0
    mass_set = set(peptide_spectrum + spectrum)
    for mass in mass_set:
        result += min(peptide_spectrum.count(mass), spectrum.count(mass))
    return result


def trim(leader_board, spectrum, N):
    if len(leader_board) <= N:
        return leader_board

    scores = {}
    for i, peptide in enumerate(leader_board):
        scores[i] = score_func(peptide, spectrum)

    sorted_scores = sorted(scores.values())
    sorted_scores.reverse()
    threshold = sorted_scores[N - 1]

    result = []
    for index, score in scores.items():
        if score >= threshold:
            result.append(leader_board[index])
    return result


def expand(peptides):
    result = []
    aa_masses = list(AMINO_ACIDS_MASS.values())
    for p in peptides:
        for m in aa_masses:
            result.append(p + [m])
    return result


def leaderboard_cylopeptide_sequencing(N, spectrum):
    leader_board = [[]]
    leader_peptide = []

    while len(leader_board) > 0:
        leader_board = expand(leader_board)
        for peptide in leader_board:
            if sum(peptide) == spectrum[-1]:
                if score_func(peptide, spectrum) > score_func(leader_peptide, spectrum):
                    leader_peptide = peptide
            elif sum(peptide) > spectrum[-1]:
                leader_board = [pep for pep in leader_board if pep != peptide]
        leader_board = trim(leader_board, spectrum, N)
    return leader_peptide


def spectral_convolution(spectrum):
    convolution = []
    for i in spectrum:
        for j in spectrum:
            if i - j > 0:
                convolution.append(i - j)
    convolution.sort()
    result = [i for i in convolution[1:]]
    result.append(convolution[0])
    return result


print("Pattern Translation Problem: ")
print(pattern_translation_problem("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
print("Generating Theoretical Spectrum Problem: ")
print(generating_theoretical_spectrum("NQEL"))
print("CycloPeptide Scoring Problem: ")
print(cyclopeptide_scoring("NQEL", [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]))
print("Leaderboard CycloPeptide Sequencing: ")
print(leaderboard_cylopeptide_sequencing(10, [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]))
print("Spectral Convolution Problem: ")
print(spectral_convolution([0, 137, 186, 323]))
