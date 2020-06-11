from itertools import product
from collections import Counter

import numpy as np

NUCLEOTIDES = "ACGT"
nucl2ind = {nucleotide: i for i, nucleotide in enumerate(NUCLEOTIDES)}

def generate_all_kmers(k):
  for pattern_tuple in product(NUCLEOTIDES, repeat=k):
    yield "".join(pattern_tuple)

def generate_kmers(dna_string, k):
    for i in range(len(dna_string)-k+1):
        yield dna_string[i:i+k]

def hamming_distance(kmer1, kmer2):
    distance = 0
    for c1, c2 in zip(kmer1,kmer2):
        if c1!=c2:
            distance += 1
    return distance

def approximate_matching(pattern, genome, distance):
    indices = []
    for i in range(len(genome)-len(pattern)+1):
        kmer = genome[i:i+len(pattern)]
        if hamming_distance(pattern,kmer) <= distance:
            indices.append(i)
    return indices

def kmer_neighbors(kmer, distance):
    if distance == 0:
        return {kmer}
    
    if len(kmer) == 1:
        return {"A", "T", "G", "C"}
    
    neighbors = dict()
    suffix = kmer[1:]
    suffix_neighbors = kmer_neighbors(suffix, distance)
    for neighbor in suffix_neighbors:
        if hamming_distance(suffix, neighbor) < distance:
            for nucleotide in NUCLEOTIDES:
                neighbors[nucleotide + neighbor] = 0
        else:
             neighbors[kmer[0] + neighbor] = 0
    
    return neighbors

def create_profile(motifs):
  profile = np.zeros((4,len(motifs[0])))
  for coli in range(len(motifs[0])):
    for rowi in range(len(motifs)):
      i = nucl2ind[motifs[rowi][coli]]
      profile[i, coli] += 1/len(motifs)
  return profile

def create_laplace_profile(motifs):
  profile = np.ones((4,len(motifs[0])))
  for coli in range(len(motifs[0])):
    for rowi in range(len(motifs)):
      i = nucl2ind[motifs[rowi][coli]]
      profile[i, coli] += 1/len(motifs)
  return profile

def profile_most_probable_kmer(dna_string, k, profile):
  max_prob = 0
  most_prob_kmer = dna_string[0:k]
  for kmer in generate_kmers(dna_string, k):
    prob = np.prod([profile[nucl2ind[nucl],i] for i, nucl in enumerate(kmer)])
    if prob > max_prob:
      max_prob = prob
      most_prob_kmer = kmer
  return most_prob_kmer

def score_motifs(motifs):
  score = 0
  for coli in range(len(motifs[0])):
    c = Counter()
    for rowi in range(len(motifs)):    
      c[motifs[rowi][coli]] += 1
    most_frequent_nucl = max(c.items(), key=lambda x: x[1])[0]
    c.pop(most_frequent_nucl)
    score += sum(c.values())
  return score
