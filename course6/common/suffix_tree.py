BLUE = 0
RED = 1
PURPLE = 2

class Node(object):
    def __init__(self, children, position, label, color=-1):
        self.children = children
        self.position = position
        self.length = len(label)
        self.label = label
        self.color = color
        
    def __repr__(self):
        children = "-".join(self.children.keys())
        return f"{self.label}|{self.color}|{self.position}|{self.length}:{children}"
              
    def get_color(self):
        if self.color != -1:
            return self.color
        else:
            colors = [child.get_color() for child in self.children.values()]
            if len(set(colors)) > 1:
                self.color = PURPLE
            elif len(colors)>0:
                color_sum = sum(colors)
                if color_sum == 0:
                    self.color = BLUE
                elif color_sum == len(colors):
                    self.color = RED
                elif color_sum == 2*len(colors):
                    self.color = PURPLE
                else:
                    raise ValueError(f"all colors should be same value {colors}")            
            elif self.node.label[-1] == "$": #leaf
                pass
            else:
                raise ValueError(f"{self} should be leaf")
                
            return self.color

class SuffixTree(object):
    #position of leaf node is the start position of the suffix 
    
    #Node = namedtuple('Node', ["children","position","label"])
    stop_symbol = "$"
    comparison_symbol = "#"
    
    def __init__(self, text):
        self.text = text + "$"
        self.separating_position = -1
        
        self._root = Node({}, -1, "")
        
        for i, s in enumerate(self.text):
            self._add_suffix(self.text[i:], i)
            if s == self.comparison_symbol:
                self.separating_position = i
            
        self.compact()
        if self.separating_position != -1:
            #print("coloring nodes")
            self.color_nodes()
        
    def _color_leafs(self, node):
        for child in node.children.values():
            if child.label[-1] == self.stop_symbol:
                #print(f"coloring leaf {child}")
                if child.position <= self.separating_position:
                    child.color = BLUE
                else:
                    child.color = RED
            else:
                self._color_leafs(child)
                
    def color_nodes(self):
        self._color_leafs(self._root)
        return self._root.get_color()
        
    def _add_suffix(self, suffix, start_position):
        self._current_node = self._root
        for i, symbol in enumerate(suffix):
            if symbol in self._current_node.children:
                self._current_node = self._current_node.children[symbol]
            else:
                if symbol == self.stop_symbol:
                    self._add_leaf(start_position)
                else:
                    newNode = self._add_symbol(symbol, start_position+i)
                    self._current_node = newNode
        
    def _add_symbol(self, symbol, position):
        newNode = Node({}, position, symbol)
        self._current_node.children[symbol] = newNode
        return newNode
    
    def _add_leaf(self, start_position):
        newNode = Node({}, start_position, self.stop_symbol)
        self._current_node.children[self.stop_symbol] = newNode
         
    def thread(self, pattern):
        self._current_node = self._root
        return self._thread(pattern, 0)
            
    def _thread(self, pattern, position):
        if position == len(pattern):
            return position
            
        symbol = pattern[position]
        self._current_node = self._current_node.children.get(symbol)
        if self._current_node is None:
            return position
        
        step = min(len(self._current_node.label), len(pattern)-position)
        if pattern[position:position+step] == self._current_node.label[:step]:
            position += step
            if position == len(pattern):
                return position
        else:
            p = 0
            for (s1, s2) in zip(pattern[position:position+step], self._current_node.label[:step]):
                if s1 == s2:
                    p += 1
                else:
                    break
            return position + p 
        
        return self._thread(pattern, position)
    
    def compact(self):
        self._compact(self._root)
    
    def _compact(self, node):
        if len(node.children) == 1:
            child = next(iter(node.children.values()))
            node.label += child.label
            node.length = len(node.label)
            node.children = child.children
            if child.label == self.stop_symbol:
                node.position = child.position
            self._compact(node=node)
        else:
            for label, child in node.children.items():
                self._compact(node=child)
            
    def labels(self):
        node_labels = []
        self._labels(self._root, node_labels)
        return node_labels
    
    def _labels(self, node, node_labels=[]):
        for label, child in node.children.items():
            node_labels.append(child.label)
            self._labels(child, node_labels)
            
    def get_leafs(self, node):
        if node.label!="" and node.label[-1] == self.stop_symbol:
            return 1
        else:
            nleafs = 0
            for child in node.children.values():
                nleafs += self.get_leafs(child)
            return nleafs
            
    def find_longest_repeat(self):
        substring_count = Counter()
        self._find_longest_repeat(self._root, substring_count)
        return sorted(filter(lambda tpl: tpl[1]>=2, substring_count.items()), 
               key=lambda tpl: len(tpl[0]), reverse=True)[0][0]
    
    def _find_longest_repeat(self, node, substring_count, substring=""):
        for child in node.children.values():
            nleafs = self.get_leafs(child)
            substring_count[substring + child.label] = nleafs
            self._find_longest_repeat(child, substring_count, substring + child.label)
    
    def _find_common_substrings(self, node, substring="", common_strings=set()):
        if node.color != PURPLE:
            return
        
        stop = True
        for child in node.children.values():
            if child.color == PURPLE:
                self._find_common_substrings(child, substring + child.label, common_strings)
                stop = False
                
        if stop:
            common_strings.add(substring)
            
    def find_common_substrings(self):
        common_strings = set()
        self._find_common_substrings(self._root, "", common_strings)
        return common_strings
        
    def find_longest_substring(self):
        #todo remove special symbols from string if they appear in there
        substrings = self.find_common_substrings()
        return sorted([(s, len(s)) for s in substrings], key=lambda tpl: tpl[1])[-1][0]
    
    def find_shortest_nonshared_substring(self, text):
        nonshared = set()
        for i, _ in enumerate(text):
            s = text[i:]
            pos = self.thread(s)
            if pos != len(s):
                nonshared.add(s[:pos+1])
        return sorted([(s, len(s)) for s in nonshared], key=lambda tpl: tpl[1])[0][0]
    
    
    def _suffix_array(self, node, suffix_array=list()):
        for _, child in sorted(node.children.items(), key = lambda tpl: tpl[0]):
            if child.label[-1] == self.stop_symbol:
                suffix_array.append(child.position)
            else:
                self._suffix_array(child, suffix_array)
                
    def suffix_array(self):
        suffix_array = list()
        self._suffix_array(self._root, suffix_array)
        return suffix_array
