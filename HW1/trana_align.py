from collections import deque
import sys

purines = set({"a", "g"})
pyrimidines = set({"c", "t"})
GAP_PENALTY = -10

def score(tide1, tide2):
    if tide1 == tide2:
        return 5
    elif (tide1 in purines and tide2 in purines) or (tide1 in pyrimidines and tide2 in pyrimidines):
        return -2
    else:
        return -3
    
def generate_dp(seqs):
    dp = [[0] * (len(seqs[1]) + 1) for _ in range(len(seqs[0]) + 1)]
    for i in range(1, len(dp)): 
        dp[i][0] = dp[i-1][0] + GAP_PENALTY
    for j in range(1, len(dp[0])): 
        dp[0][j] = dp[0][j-1] + GAP_PENALTY

    for i in range(1, len(dp)):
        for j in range(1, len(dp[0])):
            dp[i][j] = max(dp[i-1][j] + GAP_PENALTY, dp[i][j-1] + GAP_PENALTY, dp[i-1][j-1] + score(seqs[0][i-1],seqs[1][j-1]))
    return dp

def traceback(seqs, dp):
    i = len(dp) - 1
    j = len(dp[0]) - 1
    ans = [deque(), deque()]
    while i > 0 or j > 0:
        # print("i:",i,"j:",j)
        if i == 0:
            ans[1].appendleft(seqs[1][j-1])
            ans[0].appendleft("-")
            j -= 1
        elif j == 0:
            ans[0].appendleft(seqs[0][i-1])
            ans[1].appendleft("-")
            i -= 1 
        else:
            if dp[i][j] == dp[i-1][j] + GAP_PENALTY:
                ans[0].appendleft(seqs[0][i-1])
                ans[1].appendleft("-")
                i -= 1 
            elif dp[i][j] == dp[i][j-1] + GAP_PENALTY:
                ans[1].appendleft(seqs[1][j-1])
                ans[0].appendleft("-")
                j -= 1
            else:
                ans[0].appendleft(seqs[0][i-1])
                ans[1].appendleft(seqs[1][j-1])
                i -=1
                j -=1

    ans[0] = ''.join(ans[0])
    ans[1] = ''.join(ans[1])
    return ans

def main(file):
    print("Reading ", file)
    file1 = open(file, 'r')
    Lines = file1.readlines()
    seqs = ["",""]
    idx = -1
    for line in Lines:
        line = line.strip()
        line = line.lower()
        if line == "":
            continue
        elif line[0] == ">":
            idx += 1
            if idx > 1:
                print("Fasta file should only contain two sequences to align!")
                return
        else:
            seqs[idx] += line
    
    dp = generate_dp(seqs)
    
    # for i in range(1, len(dp)):
    #     print(dp[i])

    ans = traceback(seqs, dp)

    print("Best alignment\n{}\n{}\nBest Score: {}".format(ans[0], ans[1], dp[-1][-1]))



if __name__ == "__main__":
    if (len(sys.argv) < 2):
        print("please give filename in argument.")
    else:
        main(sys.argv[1]);