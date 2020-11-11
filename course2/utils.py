from collections import defaultdict

complement_map = {
    "A": "T",
    "G": "C",
    "T": "A",
    "C": "G"
}

def reverse_complement(dna_string):
    reverse = list()
    for nucleotide in dna_string[::-1]:
        reverse.append(complement_map[nucleotide])
    return "".join(reverse)

def composition(k, dna_string):
    return set(composition_iterator(k, dna_string))

def composition_iterator(k, dna_string):
    for i in range(0,len(dna_string)-k+1):
        yield dna_string[i:i+k]

def path_to_genome(kmers):
    k = len(kmers[0])
    genome = kmers[0]
    for i, kmer in enumerate(kmers[1:]):
        prev_kmer = kmers[i]
        assert kmer[0:k-1] == prev_kmer[1:k]
        genome += kmer[-1]
    return genome

def debruijn(genome, k):
    graph = defaultdict(list)
    for kmer in composition_iterator(k, genome):
        graph[kmer[0:k-1]].append(kmer[1:])
    return graph

def debruijn_from_patterns(patterns):
    k = len(patterns[0])
    graph = defaultdict(list)
    for pattern in patterns:
        prefix = pattern[0:k-1]
        suffix = pattern[1:]
        graph[prefix].append(suffix)
    return graph

def eulerian_cycle(a_list, location):
    stack = list()
    circuit = list()
    
    while a_list[location] or stack:
        neighbors = a_list[location]
        if neighbors:
            stack.append(location)
            location = neighbors.pop()
        else:
            circuit.append(location)
            location = stack.pop()
    
    circuit.append(location)      
    
    return circuit[::-1] 

def find_start_node(a_list):
    outgoing = dict()
    incoming = defaultdict(int)
    
    for h, l in a_list.items():
        outgoing[h] = len(l)
        for t in l:
            incoming[t] += 1
            
    for h, cnt in outgoing.items():
        if incoming[h] < cnt:
            return h

        
def eulerian_path(a_list):
    start_node = find_start_node(a_list)
    path = eulerian_cycle(a_list, start_node)
    return path

def assemble_genome(patterns):
    a_list = debruijn_from_patterns(patterns)
    path = eulerian_path(a_list)
    return path_to_genome(path)

def universal_string(patterns, start_node):
    a_list = debruijn_from_patterns(patterns)
    path = eulerian_cycle(a_list, start_node)
    path = path[:-(len(start_node))]
    return path_to_genome(path)

def debruijn_from_readpairs(patterns, k, d):
    graph = defaultdict(list)
    for pattern in patterns:
        kmer1, kmer2 = pattern.strip().split('|')
        prefix = f"{kmer1[0:k-1]}|{kmer2[0:k-1]}"
        suffix = f"{kmer1[1:k]}|{kmer2[1:k]}"
        graph[prefix].append(suffix)
    return graph

def pairpath_to_genome(pairs, k, d):
    k = k -1
    genome = pairs[0].split("|")[0]
    for i, pair in enumerate(pairs[1:]):
        prefix = pair.split("|")[0]
        prev_pair = pairs[i].split("|")[0]
        assert prefix[0:k-1] == prev_pair[1:k]
        genome += prefix[-1]
    prev_pair = pairs[-1].split("|")[0]
    for i in range(len(pairs)-k-d-1,len(pairs)):
        prefix = pairs[i].split("|")[1]
        print(prefix, prev_pair)
        assert prefix[0:k-1] == prev_pair[1:k]
        genome += prefix[-1]
        prev_pair = prefix
    return genome

def load_mass():
    aa2mass = dict()
    with open("../data/integer_mass_table.txt") as fin:
        for line in fin:
            aa, mass = line.split()
            aa2mass[aa.strip()] = int(mass.strip())
    return aa2mass

def cyclic_spectrum(peptide, aa2mass):
    prefix_mass = [0]
    for i, aa in enumerate(peptide):
        prefix_mass.append(prefix_mass[i]+aa2mass[aa])
    spectrum = [0]    
    for i in range(len(peptide)):
        for j in range(i+1,len(peptide)+1):
            m = prefix_mass[j]-prefix_mass[i]
            spectrum.append(m)
            if i>0 and j<len(peptide):
                spectrum.append(prefix_mass[-1]-m)
    return sorted(spectrum)

def linear_spectrum(peptide, aa2mass):
    prefix_mass = [0]
    for i, aa in enumerate(peptide):
        prefix_mass.append(prefix_mass[i]+aa2mass[aa])
    spectrum = [0]    
    for i in range(len(peptide)):
        for j in range(i+1,len(peptide)+1):
            m = prefix_mass[j]-prefix_mass[i]
            spectrum.append(m)
    return sorted(spectrum)
