from itertools import product

NUCLEOTIDES = "ACGT"

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
