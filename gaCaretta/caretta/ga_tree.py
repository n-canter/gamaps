from __future__ import annotations
import copy
import numpy as np
import time
from multiprocess import Pool
import math

from typing import Type, Optional, List, Tuple

class Node:
    val: int

    left_edges: int
    right_edges: int
    
    left: Optional[Node]
    right: Optional[Node]

    traverse_index: int = 0
    
    def __init__(self, val: int):
        self.left = None
        self.right = None
        self.left_edges = 0
        self.right_edges = 0
        self.traverse_index = 0
        self.val = val

    def split_given_edge(self, src: Node, dst: Node, val: int) -> None:
        node = Node(-1)
        leaf = Node(val)

        if src.left == dst:
            src.left = node
            src.left_edges += 2
        else:
            src.right = node
            src.right_edges += 2

        node.left = dst
        node.right = leaf
        node.left_edges = dst.left_edges + dst.right_edges + 1
        node.right_edges = 1
        
    def split_edge(self, edge: int, val: int) -> None:
        if edge < self.left_edges:
            if edge == self.left_edges - 1:
                self.split_given_edge(self, self.left, val)
            else:
                self.left_edges += 2
                self.left.split_edge(edge, val)
        else:
            if edge == self.left_edges + self.right_edges - 1:
                self.split_given_edge(self, self.right, val)
            else:
                self.right_edges += 2
                self.right.split_edge(edge - self.left_edges, val)

    def traverse(self, result: List[int]) -> None:
        if self.left is not None:
            self.left.traverse(result)

        if self.val != -1:
            result.append(self.val)
            
        if self.right is not None:
            self.right.traverse(result)

    def permutate(self, src: List[int], dst: List[int]) -> None:
        Node.traverse_index = 0
        self.permutation_apply(src, dst)

    def permutation_apply(self, src: List[int], dst: List[int]) -> None:
        if self.left is not None:
            self.left.permutation_apply(src, dst)

        if self.val != -1:
            if src[Node.traverse_index] != dst[Node.traverse_index]:
                self.val = dst[Node.traverse_index]
            Node.traverse_index += 1
            
        if self.right is not None:
            self.right.permutation_apply(src, dst)

    def align(self, size: int) -> List[List[int]]:
        Node.traverse_index = size
        result: List[List[int]] = []
        self.align_do(result)
        return result

    def align_do(self, result: List[List[int]]) -> int:
        if self.left is None and self.right is None:
            return self.val
        
        if self.left is not None and self.right is not None:
            cur = Node.traverse_index
            
            Node.traverse_index += 1
            left = self.left.align_do(result)
            right = self.right.align_do(result)
            
            result.append([left, cur])
            result.append([right, cur])
            return cur

        return -2

    def display(self):
        lines, *_ = self._display_aux()
        for line in lines:
            print(line)
        print()

    def _display_aux(self):
        """Returns list of strings, width, height, and horizontal coordinate of the root."""
        # No child.
        if self.right is None and self.left is None:
            line = '%s' % self.val
            width = len(line)
            height = 1
            middle = width // 2
            return [line], width, height, middle

        # Only left child.
        if self.right is None:
            lines, n, p, x = self.left._display_aux()
            s = '%s' % self.val
            u = len(s)
            first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s
            second_line = x * ' ' + '/' + (n - x - 1 + u) * ' '
            shifted_lines = [line + u * ' ' for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, n + u // 2

        # Only right child.
        if self.left is None:
            lines, n, p, x = self.right._display_aux()
            s = '%s' % self.val
            u = len(s)
            first_line = s + x * '_' + (n - x) * ' '
            second_line = (u + x) * ' ' + '\\' + (n - x - 1) * ' '
            shifted_lines = [u * ' ' + line for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, u // 2

        # Two children.
        left, n, p, x = self.left._display_aux()
        right, m, q, y = self.right._display_aux()
        s = '%s' % self.val
        u = len(s)
        first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s + y * '_' + (m - y) * ' '
        second_line = x * ' ' + '/' + (n - x - 1 + u + y) * ' ' + '\\' + (m - y - 1) * ' '
        if p < q:
            left += [n * ' '] * (q - p)
        elif q < p:
            right += [m * ' '] * (p - q)
        zipped_lines = zip(left, right)
        lines = [first_line, second_line] + [a + u * ' ' + b for a, b in zipped_lines]
        return lines, n + m + u, max(p, q) + 2, n + u // 2
                

class Tree:
    size: int
    score: Optional[float]
    root: Optional[Node]
    
    def __init__(self, size, perm):
        self.score = None
        self.size = size
        self.root = None

        if perm:
            self.from_perm(perm)
            return
        for i in range(0, size):
            self.insert(i)
   

    def insert(self, val: int) -> None:
        if self.root is None:
            self.root = Node(val)
            return

        edges = self.root.left_edges + self.root.right_edges
        roll = np.random.randint(0, edges + 1)
        if roll == edges:
            new_root = Node(-1)

            new_root.left = self.root
            new_root.right = Node(val)
            new_root.left_edges = edges + 1
            new_root.right_edges = 1
            self.root = new_root
            return

        self.root.split_edge(roll, val)
        
    def from_perm(self, perm: List[int]):
        if self.root is None:
            self.root = Node(0)
            
        for i, roll in enumerate(perm[1:]):
            edges = self.root.left_edges + self.root.right_edges
            if roll == edges:
                new_root = Node(-1)

                new_root.left = self.root
                new_root.right = Node(i+1)
                new_root.left_edges = edges + 1
                new_root.right_edges = 1
                self.root = new_root
                continue
            self.root.split_edge(roll, i+1)
        
    def mutate(self) -> None:
        src = self.traverse()
        dst = list(src)

        a, b = random_pair(self.size)
        dst[a], dst[b] = src[b], src[a]
        
        self.score = None
        self.root.permutate(src, dst)

    def align_order(self) -> List[List[int]]:
        order = self.root.align(self.size)
        # order.sort(key=lambda t: t[1])
        return order

    @staticmethod
    def crossover(a: Tree, b: Tree) -> Tuple[Tree, Tree]:
        child_a: Tree = copy.deepcopy(a)
        child_b: Tree = copy.deepcopy(b)

        child_a.score = None
        child_b.score = None
        
        src_a = a.traverse()
        src_b = b.traverse()
        
        dst_a, dst_b = pmx(src_a, src_b)
        child_a.root.permutate(src_a, dst_a)
        child_b.root.permutate(src_b, dst_b)
        
        return (child_a, child_b)

    def traverse(self) -> List[int]:
        ret: List[int] = []
        self.root.traverse(ret)
        return ret

    def display(self):
        print(f"({self.score})")
        self.root.display()

class Ga:
    n: int
    size: int

    max_score: float
    avg_score: float
    best: Optional[Tree]
            
    lim_time: int
    lim_iter: int
    lim_same: int

    elite: int
    convergence_threshold: float
    
    population: List[Tree]
    
    def __init__(self, n: int, score, score_ctx,
                 lim_time: int = 0, lim_iter: int = 0, lim_same: int = 0,
                 pop_size: int = 0,
                 elite: int = 10,
                 convergence_threshold: float = 0.05):

        if pop_size != 0:
            self.n = pop_size
        else:
            self.n = self.population_size(n)
            
        self.size = n
        self.score = score
        self.score_ctx = score_ctx
        
        self.lim_time = lim_time
        self.lim_iter = lim_iter
        self.lim_same = lim_same
        
        self.elite = elite
        self.convergence_threshold = convergence_threshold

        self.max_score = 0
        self.avg_score = 0
        self.best = None

        self.population = [Tree(n, None) for i in range(self.n)]

    def population_size(self, size: int) -> int:
        if size > 300:
            return 2 * size
        
        x1 = 1.0
        y1 = 10.0
        x2 = 300.0
        y2 = 2

        y = ((y2 - y1) * size + x2 * y1 - x1 * y2) / (x2 - x1);
        return math.ceil(y * size)

    def selection(self):
        self.population.sort(key=lambda t: t.score, reverse=True)
        
        offset = self.size // self.elite
        if offset == 0:
            offset = 1
            
        ret = self.population[:offset]
           
        candidates = self.population[offset:]
        while len(ret) < self.n:
            if len(candidates) == 1:
                ret.append(candidates[0])
                break
            
            perm = np.random.permutation(len(candidates))
            loosers: List[Tree] = []
            for i in range(len(candidates) // 2):
                a = candidates[2 * i]
                b = candidates[2 * i + 1]
                if (a.score / (a.score + b.score)) < np.random.uniform():
                    ret.append(a)
                    loosers.append(b)
                else:
                    ret.append(b)
                    loosers.append(a)
                    

                if len(ret) == self.n:
                    break

            if len(candidates) % 2 == 1:
                loosers.append(candidates[-1])
            candidates = loosers
                    
        self.population = ret
        self.stats()

    def fitness(self):
        tasks = []
        for i in range(len(self.population)):
            v = self.population[i]
            if v.score is None:
                tasks.append((np.array(v.align_order()), *self.score_ctx))
            else:
                tasks.append((None, *self.score_ctx))

        with Pool(processes=8) as pool:
            scores = pool.starmap(self.score, tasks)
            
            for i in range(len(self.population)):
                if self.population[i].score is None:
                    self.population[i].score = scores[i]
        
        self.stats()

    def converged(self) -> bool:
        if self.max_score == 0:
            return True
        return ((self.max_score - self.avg_score)/ self.max_score) < self.convergence_threshold

    def iteration(self, iteration: int):
        if iteration > 0 and self.converged():
            for i in range(self.n):
                self.population.append(Tree(self.size, None))

        if iteration == 0 or self.converged():
            self.fitness()

        self.selection()
        self.crossover()
        self.fitness()
        self.mutation()
        self.fitness()

        
    def stats(self):
        max_score: float = 0
        sum_score: float = 0

        cur_best: float = 0
        cur_best_idx = -1
        if self.best is not None:
            cur_best = self.best.score
            
        for i in range(len(self.population)):
            v = self.population[i].score
            max_score = max(max_score, v)
            sum_score += v
            if cur_best < v:
                cur_best_idx = i
                cur_best = v

        if cur_best_idx >= 0:
            self.best = copy.deepcopy(self.population[i])
        
        self.max_score = max_score
        self.avg_score = sum_score / len(self.population)
        
        
    def mutation(self):
        for v in self.population:
            if np.random.uniform() < self.mutation_prob(v.score):
                v.mutate()

    def crossover(self):
        perm = np.random.permutation(self.n)
        for i in range(self.n // 2):
            a = self.population[perm[2 * i]]
            b = self.population[perm[2 * i + 1]]
            if np.random.uniform() >= self.crossover_prob(a.score, b.score):
                continue

            child_a, child_b = Tree.crossover(a, b)
            self.population.append(child_a)
            self.population.append(child_b)
            
    def mutation_prob(self, score: float) -> float:
        if self.max_score == 0:
            return 1

        if self.max_score == self.avg_score:
            return 0.01
        
        if score < self.avg_score:
            return 0.5
        
        p = 0.5 * (self.max_score - score) / (self.max_score - self.avg_score)
        if p < 0.01:
            p = 0.01
        return p

    def crossover_prob(self, a_score: float, b_score: float) -> float:
        if self.max_score == 0:
            return 1

        if self.max_score == self.avg_score:
            return 0.5
        
        m = (a_score + b_score) / 2;
        if m < self.avg_score:
            return 0.05

        return (m - self.avg_score) / (self.max_score - self.avg_score)

    def terminate(self, start, n_iter: int, unchanged: int) -> bool:
        if self.lim_time != 0 and time.time() - start >= self.lim_time:
            return True
        
        if self.lim_iter != 0 and n_iter >= self.lim_iter:
            return True

        if self.lim_same != 0 and unchanged >= self.lim_same:
            return True

        return False


    def run(self):
        n_iter: int = 0
        unchanged: int = 0
        start = time.time()
        best = self.best
        while not self.terminate(start, n_iter, unchanged):
            self.iteration(n_iter)
            n_iter += 1
            if best == self.best:
                unchanged += 1
            else:
                unchanged = 0
            best = self.best
            if best is not None:
                best.display()
            print(f"iteration {n_iter}: elapsed {time.time() - start}, avg:{self.avg_score} max:{self.max_score}")
        return best.align_order()
                
def random_pair(n: int) -> Tuple[int, int]:
    a = np.random.randint(n)
    b = np.random.randint(n)
    while a == b:
        b = np.random.randint(n)
    if a > b:
        a, b = b, a
    return (a, b)

def pmx_do(a: List[int], b: List[int], m: int, n: int) -> List[int]:
    ret = [0] * len(a)
    
    s = {}
    mapping = {}

    for i in range(m, n):
        mapping[b[i]] = i
        s[b[i]] = True
        ret[i] = b[i]

    for l, r in [(0, m), (n, len(a))]:
        for i in range(l, r):
            if a[i] not in s:
                ret[i] = a[i]
            else:
                v = a[mapping[a[i]]]
                while v in s:
                    v = a[mapping[v]]
                ret[i] = v
    return ret

def pmx(a: List[int], b: List[int]) -> Tuple[List[int], List[int]]:
    m, n = random_pair(len(a))
    return (pmx_do(a, b, m, n), pmx_do(b, a, m, n))




