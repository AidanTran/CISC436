import math
import random
import sys


class Hmm:
    def __init__(self):
        trans_0 = random.randint(40, 60) / 100
        trans_1 = random.randint(40, 60) / 100
        percentages_0 = [random.randint(40, 60) for _ in range(4)]
        percentages_0 = [p/sum(percentages_0) for p in percentages_0]
        percentages_1 = [random.randint(40, 60) for _ in range(4)]
        percentages_1 = [p/sum(percentages_1) for p in percentages_1]

        # [0][1] is transition from neg to positive, [1][1] is transition from positive to positive.
        self.transition = [[trans_0, 1 - trans_0], [trans_1, 1 - trans_1]]
        # self.transition = [[.75413, .24586], [.2984, .70517]]

        # Initialize the emissions randomly
        # self.emissions = [{"A": .3519, "T": .31454, "G": .2298, "C": .103}, {
        #     "A": .12078, "T": .14475, "G": .35588, "C": .37858}]
        self.emissions = [{"A": percentages_0[0], "T": percentages_0[1], "G": percentages_0[2], "C": percentages_0[3]}, {
            "A": percentages_1[0], "T": percentages_1[1], "G": percentages_1[2], "C": percentages_1[3]}]


def generate_dp(hmm: Hmm, seq, log_transition, log_emissions):
    viterbi = [[0, 0] for _ in range(len(seq))]
    traceback = [[0, 0] for _ in range(len(seq))]

    viterbi[0][0] = math.log(.5 * hmm.emissions[0][seq[0]])
    viterbi[0][1] = math.log(.5 * hmm.emissions[1][seq[0]])

    for i in range(1, len(viterbi)):
        # For cell 0
        neg_prev_neg = viterbi[i-1][0] + log_transition[0][0]
        neg_prev_pos = viterbi[i-1][1] + log_transition[1][0]
        viterbi[i][0] = log_emissions[0][seq[i]] + max(neg_prev_neg,
                                                       neg_prev_pos)
        traceback[i][0] = 0 if neg_prev_neg > neg_prev_pos else 1

        # For cell 1
        pos_prev_neg = viterbi[i-1][0] + log_transition[0][1]
        pos_prev_pos = viterbi[i-1][1] + log_transition[1][1]
        viterbi[i][1] = log_emissions[1][seq[i]] + max(pos_prev_neg,
                                                       pos_prev_pos)
        traceback[i][1] = 0 if pos_prev_neg > pos_prev_pos else 1

    # print(viterbi)
    # print(traceback)
    return (viterbi, traceback)


def generate_path(viterbi, traceback):
    n = len(traceback)
    path = [0] * n
    # print(viterbi[0], viterbi[1])
    path[n-1] = 0 if viterbi[0] > viterbi[1] else 1
    for i in range(n-1, 0, -1):
        path[i-1] = traceback[i][path[i]]
    return path


def count(seq, path, transition_freq, counting_freq):
    for i in range(0, len(path) - 1):
        # print(path[i], seq[i])
        counting_freq[path[i]][seq[i]] += 1
        transition_freq[path[i]][path[i+1]] += 1

    counting_freq[path[len(path)-1]][seq[len(path)-1]] += 1


def train(hmm: Hmm, seqs):
    while (True):
        log_transition = [[math.log(hmm.transition[0][0]), math.log(hmm.transition[0][1])], [
            math.log(hmm.transition[1][0]), math.log(hmm.transition[1][1])]]
        log_emissions = [{"A": math.log(hmm.emissions[0]["A"]), "T": math.log(hmm.emissions[0]["T"]), "G": math.log(hmm.emissions[0]["G"]), "C": math.log(hmm.emissions[0]["C"])}, {
            "A": math.log(hmm.emissions[1]["A"]), "T": math.log(hmm.emissions[1]["T"]), "G": math.log(hmm.emissions[1]["G"]), "C": math.log(hmm.emissions[1]["C"])}]

        transition_freq = [[0, 0], [0, 0]]
        counting_freq = [{"A": 0, "T": 0, "G": 0, "C": 0},
                         {"A": 0, "T": 0, "G": 0, "C": 0}]
        for seq in seqs:
            (viterbi, traceback) = generate_dp(
                hmm, seq, log_transition, log_emissions)
            path = generate_path(viterbi[-1], traceback)
            # print(path)
            count(seq, path, transition_freq, counting_freq)

        # print(transition_freq)
        # print(counting_freq)

        # To prevent divide by 0 errors
        for i in range(2):
            for j in range(2):
                if transition_freq[i][j] == 0:
                    transition_freq[i][j] += 1
            for j in ["A", "T", "G", "C"]:
                if counting_freq[i][j] == 0:
                    counting_freq[i][j] += 1

        sum_0 = sum(counting_freq[0].values())
        sum_1 = sum(counting_freq[1].values())
        new_emissions = [{"A": counting_freq[0]["A"]/sum_0, "T": counting_freq[0]["T"]/sum_0, "G": counting_freq[0]["G"]/sum_0, "C": counting_freq[0]["C"]/sum_0},
                         {"A": counting_freq[1]["A"]/sum_1, "T": counting_freq[1]["T"]/sum_1, "G": counting_freq[1]["G"]/sum_1, "C": counting_freq[1]["C"]/sum_1}]
        new_transition = [[p/sum(transition_freq[0]) for p in transition_freq[0]],
                          [p/sum(transition_freq[1]) for p in transition_freq[1]]]

        trans_diff = abs(hmm.transition[0][0] - new_transition[0][0]) + abs(hmm.transition[0][1] - new_transition[0][1]) + abs(
            hmm.transition[1][0] - new_transition[1][0]) + abs(hmm.transition[1][1] - new_transition[1][1])
        emiss_diff = abs(hmm.emissions[0]["A"] - new_emissions[0]["A"]) + abs(hmm.emissions[0]["T"] - new_emissions[0]["T"]) + abs(
            hmm.emissions[0]["G"] - new_emissions[0]["G"]) + abs(hmm.emissions[0]["C"] - new_emissions[0]["C"]) + abs(hmm.emissions[1]["A"] - new_emissions[1]["A"]) + abs(hmm.emissions[1]["T"] - new_emissions[1]["T"]) + abs(
            hmm.emissions[1]["G"] - new_emissions[1]["G"]) + abs(hmm.emissions[1]["C"] - new_emissions[1]["C"])

        diff = trans_diff + emiss_diff
        # print("Updating model, difference: ", diff)
        hmm.transition = new_transition
        hmm.emissions = new_emissions

        if (diff < .01):
            break


def main(training_data, testing_data):

    file1 = open(training_data, 'r')
    Lines = file1.readlines()
    seqs = []
    idx = -1
    for line in Lines:
        line = line.strip()
        line = line.upper()
        if line == "":
            continue
        elif line[0] == ">":
            seqs.append("")
            idx += 1
        else:
            seqs[idx] += line

    hmm = Hmm()
    # train(hmm, ["CAGCAATAGTGAACTACAGATCTGATAATTAATATCCGGGCGGACCTACGTAGACCAAGCAACGTTGGGCGACTGCATGGCCGCAGGAAGTATCAGCCGAC"])
    train(hmm, seqs)
    print("Finish Testing, hmm with values:", hmm.transition, hmm.emissions)

    file2 = open(testing_data, 'r')
    Lines = file2.readlines()
    seqs = []
    soln = []
    idx = -1
    for line in Lines:
        line = line.strip()
        line = line.upper()
        if line == "":
            continue
        elif line[0] == ">":
            seqs.append("")
            soln.append("")
            idx += 1
        elif line[0] == "-" or line[0] == "+":
            soln[idx] += line
        else:
            seqs[idx] += line

    test(hmm, seqs, soln)


def test(hmm, seqs, soln):
    true_pos = 0
    true_neg = 0
    false_pos = 0
    false_neg = 0

    log_transition = [[math.log(hmm.transition[0][0]), math.log(hmm.transition[0][1])], [
        math.log(hmm.transition[1][0]), math.log(hmm.transition[1][1])]]
    log_emissions = [{"A": math.log(hmm.emissions[0]["A"]), "T": math.log(hmm.emissions[0]["T"]), "G": math.log(hmm.emissions[0]["G"]), "C": math.log(hmm.emissions[0]["C"])}, {
        "A": math.log(hmm.emissions[1]["A"]), "T": math.log(hmm.emissions[1]["T"]), "G": math.log(hmm.emissions[1]["G"]), "C": math.log(hmm.emissions[1]["C"])}]
    for i in range(len(seqs)):
        (viterbi, traceback) = generate_dp(
            hmm, seqs[i], log_transition, log_emissions)
        path = generate_path(viterbi[-1], traceback)
        true_path = soln[i]
        for j in range(len(path)):
            if (path[j] == 0 and true_path[j] == "-"):
                true_neg += 1
            elif (path[j] == 1 and true_path[j] == "+"):
                true_pos += 1
            elif (path[j] == 0 and true_path[j] == "+"):
                false_neg += 1
            elif (path[j] == 1 and true_path[j] == "-"):
                false_pos += 1

    recall = true_pos/(true_pos + false_neg)
    precision = true_pos/(true_pos + false_pos)
    print("Recall: ", true_pos/(true_pos + false_neg))
    print("Precision: ", true_pos/(true_pos + false_pos))
    print("F1: ", 2 * (recall * precision) / (recall + precision))
    print("Accuracy: ", (true_pos + true_neg) /
          (true_neg + true_pos + false_neg + false_pos))


if __name__ == "__main__":
    if (len(sys.argv) < 3):
        print("Please give <train_data> followed by <test_data> in the arguments.")
    else:
        main(sys.argv[1], sys.argv[2])
