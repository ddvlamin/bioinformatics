from collections import defaultdict
import sys
import numpy as np

def viterbi(sequence, state_transitions, emission_matrix, state_probabilities, backtrack):
    #print(state_probabilities)
    #sorted_states = sorted(state_probabilities.items(), key=lambda tpl: tpl[1], reverse=True)
    #yield sorted_states[0][0]
    
    if len(sequence) > 0:
        symbol = sequence[0]
        new_state_probabilities = dict()
        backtrack.append({})
        for to_state in state_transitions:
            max_prob = -sys.float_info.max
            max_state = None
            for from_state, prob in state_probabilities.items():
                p = prob + np.log(state_transitions[from_state][to_state] *  emission_matrix[to_state][symbol])
                if p > max_prob:
                    max_prob = p
                    max_state = from_state
            new_state_probabilities[to_state] = max_prob
            backtrack[-1][to_state] = max_state
        return viterbi(sequence[1:], state_transitions, emission_matrix, new_state_probabilities, backtrack)
    else:
        sorted_states = sorted(state_probabilities.items(), key=lambda tpl: tpl[1], reverse=True)
        return sorted_states[0][0]
    
def backtrack_path(backtrack, start_state):
    states = [start_state]
    next_state = backtrack[-1][start_state]
    for s in backtrack[-2::-1]:
        states.append(next_state)
        next_state = s[next_state]
    states.append(next_state)
    return "".join(states[::-1])

def get_seed_alignment(alignment, threshold):
    counts = [defaultdict(int) for i in range(len(alignment[0]))]
    columns2remove = []
    for row in alignment:
        for i, v in enumerate(row):
            counts[i][v] += 1
    
    for col, col_counts in enumerate(counts):
        freq = col_counts["-"]/len(alignment)
        if freq >= threshold:
            columns2remove.append(col)
    return set(columns2remove)

def init_transition_emission(n, alphabet):
    states = ["S", "I0"]
    for i in range(n):
        states.append(f"M{i+1}")
        states.append(f"D{i+1}")
        states.append(f"I{i+1}")
    states.append("E")
    
    transition = dict()
    emission = dict()
    for state in states:
        transition[state] = {s: 0 for s in states}
        emission[state] = {a: 0 for a in sorted(alphabet)}
        
    return transition, emission

def compute_transition(transition, emission, alignment, columns2remove, alphabet, sigma):
    compute_emission_transition_counts(transition, emission, alignment, columns2remove, alphabet)
    normalize_transition(transition)
    normalize_emission(emission)
    add_pseudocount_transition(transition, sigma=sigma)
    normalize_transition(transition)
    add_pseudocount_emission(emission, sigma=sigma)
    normalize_emission(emission)

def compute_emission_transition_counts(transition, emission, alignment, columns2remove, alphabet):
    col = 0
    prev_states = ["S"]*len(alignment)
    for i in range(len(alignment[0])):
        for row, seq in enumerate(alignment):
            if seq[i] in alphabet:
                if i not in columns2remove:
                    state_label = f"M{col+1}"
                else:
                    state_label = f"I{col}"
            else:
                if i not in columns2remove:
                    state_label = f"D{col+1}"
                else:
                    continue
            #print(f"i: {i}, col: {col}, letter: {seq[col]}, state_label: {state_label}, prev_state: {prev_states[row]}")
            transition[prev_states[row]][state_label] += 1
            prev_states[row] = state_label
            if seq[i] != "-":
                emission[state_label][seq[i]] += 1
        if i not in columns2remove:
            col += 1            
    for prev_state in prev_states:
        transition[prev_state]["E"] += 1
    
def normalize_transition(transition):
    for prev_state, next_states in transition.items():
        total = sum((v for v in next_states.values()))
        if total != 0:
            for k in next_states.keys():
                next_states[k] /= total    

def normalize_emission(emission):
    for state, state_emissions in emission.items():
        total = sum((v for v in state_emissions.values()))
        if total != 0:
            for k in state_emissions.keys():
                state_emissions[k] /= total                

def add_pseudocount_transition(transition, sigma=0.01):
    i2label = {i: k for i, k in enumerate(transition.keys())}
    
    for i in range(-1,len(transition)-1,3):
        col_start = i + 2
        col_stop = col_start + (1 if (col_start+2)==len(transition) else 2) 
        for row in range(3):
            if i+row > -1:
                for j in range(col_start,col_stop+1):
                    transition[i2label[i+row]][i2label[j]] += sigma
  
def add_pseudocount_emission(emission_matrix, sigma=0.01):
    for state, emissions in emission_matrix.items():
        if state[0] in {"I", "M"}:
            for symbol in emissions.keys():
                emissions[symbol] += sigma

def print_output(transition, emissions, alphabet):
    header = " " + " ".join(transition.keys())
    print(header)
    for k, states in transition.items():
        row = k + " " + " ".join((str(s) for s in states.values()))
        print(row)
    print("--------")
    header = " " + " ".join(sorted(alphabet))
    print(header)
    for state, state_emissions in emissions.items():
        row = state + " " + " ".join((str(s) for s in state_emissions.values()))
        print(row)   
        
class HiddenNode():
    def __init__(self, observation_it, state_label):
        self.observation_it = observation_it
        self.state_label = state_label
        self.parents = []
        self.children = []
        self.children_edges = []
        self.parent_edges = []
        self.prob = 0.0
        self.backtrack_node = None
        
    def add_edge(self, node, edge_weight):
        self.children.append(node)
        self.children_edges.append(edge_weight)
        node.parents.append(self)
        node.parent_edges.append(edge_weight)
        
        node_prob = self.prob*edge_weight
        if node_prob > node.prob:
            node.prob = node_prob
            node.backtrack_node = self
        
    def __str__(self):
        return f"{self.state_label}-{self.observation_it}"
       
class ViterbiGraph():
    def __init__(self, observations, transition, emission):
        self.root = HiddenNode(-1, "S")
        self.root.prob = 1.0
        self.transition = transition
        self.emission = emission
        self.observations = observations
        self.nodes = dict()
        self.queue = list()
        self.end_node = None

        self.queue.append(self.root)
        
        self._create()
        
    def hidden_path(self):
        path = [str(node).split("-")[0] for node in self._hidden_path()]
        return path[::-1]

    def _hidden_path(self):
        node = self.end_node
        while str(node.backtrack_node) != "S--1":
            yield node.backtrack_node
            node = node.backtrack_node
        
    def _create(self):
        while len(self.queue) > 0:
            node = self.queue.pop(0)
            print(f"node in set: {node}")
            self._add_child(node)
    
    def compute_node_probabilities(self, node):
        ski = 0
        max_ski_node = ""
        for parent, weight in zip(node.parents, node.parents_edges):
            w = parent.prob*weight
            if w > ski:
                ski = w
                max_ski_node = parent
        node.backtrack_node = max_ski_node
        #for child in node.children
        
    def _lookup_node(self, node):
        """
        adds the node if not present yet
        """
        if str(node) not in self.nodes:
            print(f"adding {node} to nodes")
            self.queue.append(node)
            self.nodes[str(node)] = node
        return self.nodes[str(node)]        
    
    def _add_child(self, node):
        observation_it = node.observation_it
        label_suffix = None
        h1 = None
        h2 = None
        h3 = None
        if node.state_label[0] == "S":
            h1label = "I0"
            h1 = HiddenNode(0, h1label)
            h2label = "M1"
            h2 = HiddenNode(0, h2label)
            h3label = "D1"
            h3 = HiddenNode(-1, h3label)
        elif node.state_label[0] != "E":
            label_suffix = int(node.state_label[1:])
            h1label = f"I{label_suffix}"
            h1 = HiddenNode(observation_it+1, h1label)
            h2label = f"M{label_suffix+1}"
            h2 = HiddenNode(observation_it+1, h2label)
            h3label = f"D{label_suffix+1}"
            h3 = HiddenNode(observation_it, h3label) 
            
        if observation_it + 1 == len(self.observations):
            if label_suffix == (len(self.transition)/3)-1:
                he = HiddenNode(observation_it+1, "E")
                if str(he) not in self.nodes:
                    self.nodes[str(he)] = he
                print(f"{str(node)} adds end child {self.nodes[str(he)]}")
                node.add_edge(self.nodes[str(he)], 1)
                self.end_node = self.nodes[str(he)]
                return
            h3 = self._lookup_node(h3)
            weight = self.transition[node.state_label][h3label]
            node.add_edge(h3, weight)
            print(f"{str(node)} - {str(h3)}: {weight}")
        elif label_suffix == (len(self.transition)/3)-1:
            h1 = self._lookup_node(h1)
            observation = self.observations[observation_it+1]
            weight = self.transition[node.state_label][h1label]*self.emission[h1label][observation]
            node.add_edge(h1, weight)  
            print(f"{str(node)} - {str(h1)}: {weight}")
        else:    
            h3 = self._lookup_node(h3)
            h1 = self._lookup_node(h1)
            h2 = self._lookup_node(h2)

            observation = self.observations[observation_it+1]
            weight = self.transition[node.state_label][h1label]*self.emission[h1label][observation]
            node.add_edge(h1, weight)   
            print(f"{str(node)} - {str(h1)}: {weight}")
            observation = self.observations[observation_it+1]
            weight = self.transition[node.state_label][h2label]*self.emission[h2label][observation]
            node.add_edge(h2, weight)
            print(f"{str(node)} - {str(h2)}: {weight}")
            weight = self.transition[node.state_label][h3label]
            node.add_edge(h3, weight)
            print(f"{str(node)} - {str(h3)}: {weight}")

    
    
        