import random
import pandas
import numpy as np


def index_to_nucleotide(i):
    return "ACGT"[i]


def nucleotide_to_index(n):
    return "ACGT".index(n)


def get_motif_counts(motifs, pseudocounts):
    length = len(motifs[0])
    counts = []
    for i in range(length):
        if pseudocounts:
            col_count = np.ones(4)
        else:
            col_count = np.zeros(4)
        for dna_sequence in motifs:
            col_count[nucleotide_to_index(dna_sequence[i])] += 1
        counts.append(col_count)
    return counts


def score_motifs(motifs):
    score = 0
    for count in get_motif_counts(motifs, False):
        score += np.sum(count) - count[np.argmax(count)]
    return score


def motif_to_profile(motifs):
    counts = get_motif_counts(motifs, True)
    profile = []
    for count in counts:
        total = float(np.sum(count))
        probability = count / total
        profile.append(probability)
    return profile


def most_probable_kmer(profile, kmer):
    prob = 1.0
    for i in range(len(kmer)):
        prob *= profile[i][nucleotide_to_index(kmer[i])]
    return prob


def random_motif(dna, k):
    return [get_random_kmer(sequence, k) for sequence in dna]


def get_random_kmer(sequence, k):
    rand = random.randint(0, (len(sequence) - k))
    return sequence[rand:rand+k]


def gibbs_random(prob):
    total = float(sum(prob))
    rand = random.random()
    p_sum = 0.0
    for i in range(len(prob)):
        p_sum += prob[i]
        if p_sum/total >= rand:
            return i
    return -1


def gibbs_sampler(dna, k, t, n):
    motifs = random_motif(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for j in range(n):
        i = random.randint(0, t-1)
        profile = motif_to_profile(motifs[:i] + motifs[i+1:])
        probs = [most_probable_kmer(profile, dna[i][x:x+k]) for x in range(len(dna[i]) - (k+1))]
        loc = gibbs_random(probs)
        motifs[i] = dna[i][loc:loc+k]
        score = score_motifs(motifs)
        if score < best_score:
            best_motifs = motifs
            best_score = score
    return best_motifs, best_score


def repeat_gibbs_sampler(dna, k, t, n):
    best_motifs = random_motif(dna, k)
    best_score = score_motifs(best_motifs)
    for i in range(0, 2000):
        motifs, score = gibbs_sampler(dna, k, t, n)
        if score < best_score:
            best_motifs = motifs
            best_score = score
    return best_motifs, best_score


def consensus(motifs):
    consensus_string = ""
    counts = get_motif_counts(motifs, False)
    for count in counts:
        consensus_string += index_to_nucleotide(np.argmax(count))
    return consensus_string


def get_data(filename):
    with open(filename) as file:
        dna = []
        for line in file:
            dna.append(line.strip())
    return dna


dna = get_data('randomData.txt')

ans, score = repeat_gibbs_sampler(dna, 7, len(dna), 200)
motif = consensus(ans)
print(motif)
print(score)

