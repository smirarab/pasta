from collections import deque


class DSet:
    def __init__(self,elements):
        self.elements = deque(elements)
        self.count = len(elements)

    def extend(self,S):
        self.elements.extend(S.elements)
        self.count += S.count

class DisjointSets:
    def __init__(self,n):
        self.groupings = [DSet([i]) for i in range(n)]
        self.n_groups = n

    def find(self,k):
        return self.groupings[k]

    def join(self,p,q):
        # merge the smaller set to the larger
        S_q = self.groupings[q]
        S_p = self.groupings[p]

        if S_q.count > S_p.count:
            S = S_q
            S_q = S_p
            S_p = S

        S_p.extend(S_q)
        for r in S_q.elements:
            self.groupings[r] = S_p 
        self.n_groups -= 1

    def same_set(self,p,q):
        return self.find(p) is self.find(q)

    def is_single(self):
        return self.n_groups == 1

