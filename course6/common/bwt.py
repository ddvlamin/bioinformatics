#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict


symbol = lambda c, i: f"{c}{i}"

def bwt(text):
    rotations = []
    for i in range(len(text)):
        rotation = text[len(text)-i:] + text[0:len(text)-i]
        rotations.append(rotation)
    rotations = sorted(rotations)
    
    last = [k[-1] for k in rotations]
    
    counts = defaultdict(int)
    new_last = list()
    for s in last:
        counts[s] += 1
        new_last.append(f"{s}{counts[s]}")
    return new_last


def get_first_occurence(last):
    first_column = sorted(last)
    first = dict()
    for i, c in enumerate(first_column):
        if c[0] not in first:
            first[c[0]] = i
    return first
            
def get_partial_suffix_array(last, first_occurence, C=5):
    suffix_array = dict()
    j = 0
    c = last[j]
    for i in range(0, len(last)):
        suffix_array[c] = len(last)-i-2
        j = first_occurence[c[0]] + int(c[1:]) - 1
        #j = last2first[c]
        c = last[j]
    suffix_array[symbol("$",1)] = len(last)-1
    
    suffix_array = sorted(suffix_array.items())
    print(f"suffix array {suffix_array}")
    
    partial = dict()
    for c, i in suffix_array:
        if i % C == 0:
            partial[c] = i
    
    
    return partial

def get_counts(column, k):
    counts = dict()
    for c in sorted({c[0] for c in column}):
        counts[c] = [0]
        
    for c, l in counts.items():
        for i, s in enumerate(column):
            if s[0] == c:
                l.append(l[i]+1)
            else:
                l.append(l[i])

    selected_indices = [i for i in range(0,len(column)+1) if i%k==0]
    
    for c, l in counts.items():
        counts[c] = [l[i] for i in selected_indices]
            
    return counts

def contains_symbol(last, top, bottom, symbol):
    for s in last[top:bottom+1]:
        if s[0] == symbol:
            return True
    return False

def get_count_symbol(symbol, pos, last, partial_count, C):
    start_pos = int(pos / C)
    symbol_counts = partial_count[symbol]
    count = symbol_counts[start_pos]
    for s in last[pos - (pos % C): pos]:
        if s[0] == symbol:
            count += 1
    return count

def better_bw_matching(first_occurence, last, pattern, partial_count, C):
    top = 0
    bottom = len(last)-1
    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if contains_symbol(last, top, bottom, symbol):
                top = first_occurence[symbol] + get_count_symbol(symbol, top, last, partial_count, C)
                bottom = first_occurence[symbol] + get_count_symbol(symbol, bottom+1, last, partial_count, C) -1
                #print(f"symbol {symbol}: top {top} - bottom {bottom}")
            else:
                #print(f"no symbol {symbol} between top {top} and bottom {bottom}")
                return 0
        else:
            return bottom - top + 1
        
def bwt_matching_positions(suffix_array, first_occurence, last, pattern, partial_count, C):
    top = 0
    bottom = len(last)-1
    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if contains_symbol(last, top, bottom, symbol):
                top = first_occurence[symbol] + get_count_symbol(symbol, top, last, partial_count, C)
                bottom = first_occurence[symbol] + get_count_symbol(symbol, bottom+1, last, partial_count, C) -1
                #print(f"symbol {symbol}: top {top} - bottom {bottom}")
            else:
                #print(f"no symbol {symbol} between top {top} and bottom {bottom}")
                return []
        else:
            #print(f"found pattern with last first symbol {symbol} from {top} to {bottom}")
            return lookup_suffixes(symbol, suffix_array, last, first_occurence, top, bottom)    
        
def lookup_suffixes(symbol, suffix_array, last, first_occurence, top, bottom):
    start_positions = []
    for i in range(top, bottom+1):
        symbol_suffix = i-first_occurence[symbol]+1
        symbol_with_suffix = f"{symbol}{symbol_suffix}"
        start_position = lookup_suffix(symbol_with_suffix, suffix_array, last, first_occurence, i)
        start_positions.append(start_position)
    return start_positions
    
def lookup_suffix(symbol, suffix_array, last, first_occurence, i):
    start_position = 0
    offset = 0
    while symbol not in suffix_array:
        symbol = last[i]
        i = first_occurence[symbol[0]] + int(symbol[1:]) - 1
        offset += 1
    start_position = suffix_array[symbol]+offset
    return start_position

def index_symbols(last):
    symbol_count = defaultdict(int)
    new_last = []
    for s in last:
        symbol_count[s] += 1
        new_last.append(f"{s}{symbol_count[s]}")
    return new_last


def bwt_matching(text, patterns, C=5):
    text += "$"
    bwt_text = bwt(text)
    partial_counts = get_counts(bwt_text, C)
    first_occurence = get_first_occurence(bwt_text)
    partial_suffix_array = get_partial_suffix_array(bwt_text, first_occurence, C)
    found = defaultdict(list)
    p5 = int(0.05*len(patterns))+1
    for i, pattern in enumerate(patterns):
        c = bwt_matching_positions(partial_suffix_array, first_occurence, bwt_text, pattern, partial_counts, C)
        found[pattern].extend(c)
        if i % p5 == 0:
            print(f"pattern {i} of {len(patterns)} processed")
    return found
        
def print_matches(matches):
    for k, v in matches.items():
        pos = " ".join([str(i) for i in sorted(v)])
        print(f"{k}: {pos}")


