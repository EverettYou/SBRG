import random
from operator import attrgetter
from itertools import combinations
""" Mat: tensor product of Pauli matrices
Mat.Xs     :: set : collection of sites of X gates
Mat.Zs     :: set : collection of sites of Z gates
Mat.ipower :: int : number of overlap between Xs and Zs (num of Y gates)
"""
class Mat:
    def __init__(self, *arg):
        l_arg = len(arg)
        if l_arg == 2: 
            self.Xs = set(arg[0])
            self.Zs = set(arg[1])
        elif l_arg == 1:
            inds = arg[0]
            self.Xs = set()
            self.Zs = set()
            if isinstance(inds, dict): # dict of inds rules
                for (i, mu) in inds.items():
                    if mu == 1:
                        self.Xs.add(i)
                    elif mu == 3:
                        self.Zs.add(i)
                    elif mu == 2:
                        self.Xs.add(i)
                        self.Zs.add(i)
            elif isinstance(inds, (tuple, list)): # list of inds
                for (i, mu) in enumerate(inds):
                    if mu == 0:
                        continue
                    elif mu == 1:
                        self.Xs.add(i)
                    elif mu == 3:
                        self.Zs.add(i)
                    elif mu == 2:
                        self.Xs.add(i)
                        self.Zs.add(i)
        elif l_arg == 0:
            self.Xs = set()
            self.Zs = set()
        self.ipower = len(self.Xs & self.Zs)
        self._key = None
    def __repr__(self):
        return "<Xs:%s Zs:%s>" % (sorted(list(self.Xs)), sorted(list(self.Zs)))
    def __hash__(self):
        if self._key is None:
            self._key = hash((tuple(self.Xs), tuple(self.Zs)))
        return self._key
    def __eq__(self, other):
        return self.Xs == other.Xs and self.Zs == other.Zs
    def __neq__(self, other):
        return self.Xs != other.Xs or self.Zs != other.Zs
# commutativity check
def is_commute(mat1, mat2):
    return (len(mat1.Xs & mat2.Zs) - len(mat1.Zs & mat2.Xs))%2 == 0
# merging Pauli indices (coefficient not determined here)
def pdot(mat1, mat2):
    return Mat(mat1.Xs ^ mat2.Xs, mat1.Zs ^ mat2.Zs)
""" Term: a Mat with coefficient and position
Term.mat :: Mat : matrix of Pauli operator
Term.val :: numeric : coefficient
Term.pos :: int : my position in Hamiltonian.terms
"""
class Term:
    def __init__(self, *arg):
        l_arg = len(arg)
        if l_arg == 2:
            self.mat, self.val = arg
        elif l_arg == 1:
            self.mat = arg[0]
            self.val = 1.        
        elif l_arg == 0:
            self.mat = Mat()
            self.val = 1.
        self.pos = 0
    def __repr__(self):
        return "%s %s" % (self.val, self.mat)
# dot product of two terms
def dot(term1, term2):
    mat1 = term1.mat
    mat2 = term2.mat
    mat = pdot(mat1, mat2)
    n = mat1.ipower + mat2.ipower - mat.ipower
    n = n + 2*len(mat1.Zs & mat2.Xs)
    s = (-1)**(n/2)
    return Term(mat, s*term1.val*term2.val)
# dot product of two terms (times additional i)
def idot(term1, term2):
    mat1 = term1.mat
    mat2 = term2.mat
    mat = pdot(mat1, mat2)
    n = mat1.ipower + mat2.ipower - mat.ipower
    n = n + 2*len(mat1.Zs & mat2.Xs) + 1
    s = (-1)**(n/2)
    return Term(mat, s*term1.val*term2.val)
""" Ham: a collection of Terms
Ham.terms :: list : terms stored in binary heap structure
Ham.mats  :: dict : mapping tat to term
Ham.imap  :: dict : mapping site to covering terms
"""
class Ham:
    def __init__(self, *arg):
        self.terms = []
        self.mats = {}
        self.imap = {}
        if len(arg) == 1:
            self.extend(arg[0])
    def __repr__(self):
        return "%s" % self.terms
    def terms_push(self, term):
        pos = len(self.terms)
        term.pos = pos
        self.terms.append(term)
        self.terms_shiftUV(pos)
    def terms_adjust(self, term):
        pos = term.pos
        self.terms_shiftUV(pos)
        self.terms_shiftIR(pos)
    def terms_shiftUV(self, pos):
        terms = self.terms
        this_term = terms[pos]
        # Follow the path to the root, moving parents down until fits.
        while pos > 0:
            parent_pos = (pos - 1) >> 1
            parent_term = terms[parent_pos]
            if abs(this_term.val) > abs(parent_term.val):
                parent_term.pos = pos
                terms[pos] = parent_term
                pos = parent_pos
                continue
            break
        if pos != this_term.pos: # if pos is new
            this_term.pos = pos
            terms[pos] = this_term
    def terms_shiftIR(self, pos):
        terms = self.terms
        end_pos = len(terms)
        this_term = terms[pos]
        child_pos = 2*pos + 1 # left child position
        while child_pos < end_pos:
            # Set child_pos to index of larger child.
            rchild_pos = child_pos + 1 # right child position
            if rchild_pos < end_pos and abs(terms[child_pos].val) < abs(terms[rchild_pos].val):
                child_pos = rchild_pos
            # Move the larger child up.
            child_term = terms[child_pos]
            if abs(this_term.val) < abs(child_term.val):
                child_term.pos = pos
                terms[pos] = child_term
                pos = child_pos
                child_pos = 2*pos + 1 # left child position
                continue
            break
        if pos != this_term.pos: # if pos is new
            this_term.pos = pos
            terms[pos] = this_term
    def imap_add(self, term):
        mat = term.mat
        for i in mat.Xs | mat.Zs:
            try:
                self.imap[i].add(term)
            except:
                self.imap[i] = {term}
    def imap_del(self, term):
        mat = term.mat
        for i in mat.Xs | mat.Zs:
            self.imap[i].remove(term)
    def push(self, term):
        if term.mat in self.mats: # if mat already exist
            old_term = self.mats[term.mat]
            old_term.val += term.val
            self.terms_adjust(old_term)
        else: # if mat is new
            self.terms_push(term)
            self.mats[term.mat] = term
            self.imap_add(term)
    def extend(self, terms):
        for term in terms:
            self.push(term)
    def remove(self, term):
        terms = self.terms
        end_pos = len(terms)
        pos = term.pos
        if pos == end_pos - 1:
            del terms[pos]
        elif 0 <= pos <= end_pos - 1:
            last_term = terms.pop()
            last_term.pos = pos
            terms[pos] = last_term
            self.terms_adjust(last_term)
        self.imap_del(term)
        del self.mats[term.mat]
    """
    def pop(self):
        top_term = self.terms[0]
        self.remove(top_term)
        return top_term"""
    def C4(self, gen, sgn = +1):
        mats = self.mats
        imap = self.imap
        gen_mat = gen.mat
        # collect terms to be transformed
        relevant_terms = set() # start with empty set
        for i in gen_mat.Xs | gen_mat.Zs: # supporting sites of gen
            if i in imap: # if i registered in imap
                relevant_terms.update(imap[i])
        relevant_terms = [term for term in relevant_terms if not is_commute(term.mat, gen_mat)]
        for term in relevant_terms:
            # remove mat
            self.imap_del(term)
            del self.mats[term.mat]
            # C4 by idot with gen
            new_term = idot(term, gen)
            # update mat & val only
            term.mat = new_term.mat
            term.val = sgn * new_term.val
        # add new mats, NOT COMBINE TO ABOVE LOOP
        for term in relevant_terms:
            mats[term.mat] = term
            self.imap_add(term)
    def forward(self, Rs):
        for R in Rs:
            self.C4(R)
    def backward(self, Rs):
        for R in reversed(Rs):
            self.C4(R,-1)
""" SBRG:
"""
class SBRG:
    def __init__(self, model):
        self.tol = 1.e-8
        self.max_rate = 2.
        self.bits = model.bits
        self.H = Ham(model.terms)
        self.Hbdy = model.terms
        self.Hblk = []
        self.Heff = []
        self.EHM = []
        self.phybits = set(range(self.bits))
        self.trash = []
    def findRs(self, mat):
        if len(mat.Xs) > 0: # if X or Y exists, use it to pivot the rotation
            pbit = min(mat.Xs) # take first off-diag qubit
            return ([idot(Term(Mat([],[pbit])), Term(mat))], pbit)
        else: # if only Z
            if len(mat.Zs) > 1:
                for pbit in sorted(list(mat.Zs)): # find first Z in phybits
                    if (pbit in self.phybits):
                        tmp = Term(Mat([pbit],[])) # set intermediate term
                        return ([idot(tmp, Term(mat)), idot(Term(Mat([],[pbit])), tmp)], pbit)
            elif len(mat.Zs) == 1:
                pbit = min(mat.Zs)
        return ([], pbit)
    def perturbation(self, H0, offdiag):
        h0 = H0.val # set h0
        min_prod = abs(h0)**2*self.tol # set minimal product
        # SiSj for commuting terms whose product val > min_prod
        SiSj = [dot(term1, term2) for (term1, term2) in combinations(offdiag, 2)
                if is_commute(term1.mat,term2.mat) and abs(term1.val*term2.val) > min_prod]
        SiSj.sort(key=attrgetter('val')) # sort by val
        # term number truncation
        max_len = round(self.max_rate*len(offdiag))
        if len(SiSj) > max_len:
            self.trash.extend([term.val/h0 for term in SiSj[:-max_len]])
            SiSj = SiSj[-max_len:]
        # multiply by H0 inverse
        H0inv = Term(H0.mat,1/h0)
        pert = [dot(H0inv,term) for term in SiSj]
        # add backward correction
        pert.append(Term(H0.mat, sum((term.val)**2 for term in offdiag)/(2*h0)))
        return pert
    def nextstep(self):
        if len(self.phybits) == 0: # return if no physical bits
            return self
        # get leading energy scale
        H0 = self.H.terms[0]
        if abs(H0.val) == 0.: # if leading scale vanishes
            self.phybits = set() # quench physical space
            return self
        # find Clifford rotations
        Rs, pbit = self.findRs(H0.mat)
        self.EHM.extend(Rs) # add to EHM
        self.H.forward(Rs) # apply to H
        # pick out offdiag terms
        offdiag = [term for term in self.H.imap[pbit] if pbit in term.mat.Xs]
        pert = self.perturbation(H0, offdiag) # 2nd order perturbation
        for term in offdiag:
            self.H.remove(term) # remove off-diagonal terms
        self.H.extend(pert) # add perturbation to H
        self.phybits.remove(pbit) # reduce physical bits
        # remove identity terms in physical space
        for term in list(self.H.imap[pbit]): # NOT REMOVE list(...)
            if len((term.mat.Xs | term.mat.Zs) & self.phybits) == 0:
                self.Heff.append(term)
                self.H.remove(term)
    def flow(self, step = float('inf')):
        step = min(step, len(self.phybits)) # adjust RG steps
        # carry out RG flow
        stp_count = 0
        while (len(self.phybits) > 0 and stp_count < step):
            self.nextstep()
            stp_count += 1
""" Model: defines Hilbert space and Hamiltonian
Model.bits  :: int : num of bits
Model.terms :: list : terms in the Hamiltonian
"""
class Model:
    def __init__(self):
        self.bits = 0
        self.terms = []
def TFIsing(L, **para):
    # L - number of sites (assuming PBC)
    # model - a dict of model parameters
    try: # set parameter alpha
        alpha = para['alpha']
        alpha_J = alpha
        alpha_K = alpha
        alpha_h = alpha
    except:
        alpha_J = para.get('alpha_J',1)
        alpha_K = para.get('alpha_K',1)
        alpha_h = para.get('alpha_h',1)
    model = Model()
    model.bits = L
    # translate over the lattice by deque rotation
    H_append = model.terms.append
    rnd_beta = random.betavariate
    for i in range(L):
        H_append(Term(Mat({i: 1, (i+1)%L: 1}), para['J']*rnd_beta(alpha_J, 1)))
        H_append(Term(Mat({i: 3, (i+1)%L: 3}), para['K']*rnd_beta(alpha_K, 1)))
        H_append(Term(Mat({i: 3}), para['h']*rnd_beta(alpha_h, 1)))
    model.terms = [term for term in model.terms if abs(term.val) > 0]
    return model