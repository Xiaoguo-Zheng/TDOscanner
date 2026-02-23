import numpy as np


def extract_features(sequence):
    features = []

    # Feature 1: Sum of G and C counts at positions 2-6
    features.append(sequence[1:6].count('G') + sequence[1:6].count('C'))

    # Feature 2: Sum of G and C counts at positions 9-12
    features.append(sequence[8:12].count('G') + sequence[8:12].count('C'))

    # Feature 3: Sum of G and C counts at positions 13-17
    features.append(sequence[12:17].count('G') + sequence[12:17].count('C'))

    # Feature 4: Sum of G and C counts at positions 9-17
    features.append(sequence[8:17].count('G') + sequence[8:17].count('C'))

    # Feature 5: Annealing temperature of the combined subsequence (positions 2-6 and 9-17),
    # calculated by T = 2×(A+T) + 4×(G+C)
    sub_sequence_2_6_9_17 = sequence[1:6] + sequence[8:17]
    T = 2 * (sub_sequence_2_6_9_17.count('A') + sub_sequence_2_6_9_17.count('T')) + 4 * (
                sub_sequence_2_6_9_17.count('G') + sub_sequence_2_6_9_17.count('C'))
    features.append(T)

    # Feature 6: Count of subsequence "TT" in positions 2-6 and 9-17
    features.append(sub_sequence_2_6_9_17.count('TT'))

    # Feature 7: Count of subsequence "AA" in positions 2-6 and 9-17
    features.append(sub_sequence_2_6_9_17.count('AA'))

    # Feature 8: Count of subsequence "ACTT" in positions 2-6
    features.append(sequence[1:6].count('ACTT'))

    # Feature 9: Count of subsequence "CTT" in positions 9-12
    features.append(sequence[8:12].count('CTT'))

    # Feature 10: Count of subsequence "TC" in positions 27-31
    features.append(sequence[26:31].count('TC'))

    # Feature 11: Whether position 35 is C (1 if yes, else 0)
    features.append(1 if sequence[34] == 'C' else 0)

    # Feature 12: Whether position 36 is C (1 if yes, else 0)
    features.append(1 if sequence[35] == 'C' else 0)

    # Feature 13: Count of subsequence "TTC" in positions 22-26
    features.append(sequence[21:26].count('TTC'))

    # Feature 14: Count of subsequence "TC" in positions 22-26
    features.append(sequence[21:26].count('TC'))

    # Feature 15: Count of subsequence "TTT" in positions 13-17
    features.append(sequence[12:17].count('TTT'))

    # Feature 16: Count of subsequence "TT" in positions 13-17
    features.append(sequence[12:17].count('TT'))

    # Feature 17: Whether position 6 is G or C (1 if yes, else 0)
    features.append(1 if sequence[5] in ['G', 'C'] else 0)

    # Feature 18: Whether position 9 is G or C (1 if yes, else 0)
    features.append(1 if sequence[8] in ['G', 'C'] else 0)

    # Feature 19: Whether position 17 is G or C (1 if yes, else 0)
    features.append(1 if sequence[16] in ['G', 'C'] else 0)

    # Feature 20: Whether position 22 is G or C (1 if yes, else 0)
    features.append(1 if sequence[21] in ['G', 'C'] else 0)

    # Feature 21: Whether position 31 is G or C (1 if yes, else 0)
    features.append(1 if sequence[30] in ['G', 'C'] else 0)

    # Feature 22: Whether position 35 is G or C (1 if yes, else 0)
    features.append(1 if sequence[34] in ['G', 'C'] else 0)

    return np.array(features)


# Testing the function with a sample sequence
# sequence = "AGCTGACTGACTGACTGACTGACTGACTGACTGACTGACT"
# features = extract_features(sequence)
# print(features)
