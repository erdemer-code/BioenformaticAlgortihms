def change_problem(money, coins):
    min_num_coins = [0]
    for m in range(1, money + 1):
        tmp = 9999999
        for i in range(0, len(coins)):
            if m >= coins[i]:
                if min_num_coins[m - coins[i]] + 1 < tmp:
                    tmp = min_num_coins[m - coins[i]] + 1
        min_num_coins.append(tmp)
        # print(min_num_coins)
    return min_num_coins[money]


def longest_common_subsequence_problem(string1, string2):
    m = len(string1)
    n = len(string2)
    lcs = [[0 for i in range(n + 1)] for i in range(m + 1)]
    # print(lcs)

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                lcs[i][j] = 0
            elif string1[i - 1] == string2[j - 1]:
                lcs[i][j] = lcs[i - 1][j - 1] + 1
            else:
                lcs[i][j] = max(lcs[i - 1][j], lcs[i][j - 1])
    idx = lcs[m][n]
    result = ["" for _ in range(idx + 1)]

    i = m
    j = n
    while i > 0 and j > 0:
        if string1[i - 1] == string2[j - 1]:
            result[idx - 1] = string1[i - 1]
            i -= 1
            j -= 1
            idx -= 1

        elif lcs[i - 1][j] > lcs[i][j - 1]:
            i -= 1
        else:
            j -= 1
    return ''.join(result)


print("The Change Problem: ")
print(change_problem(40, [1, 5, 10, 20, 25, 50]))
print("Longest Common Subsequence Problem: ")
print(longest_common_subsequence_problem("AACCTTGG", "ACACTGTGA"))
