import random
from operator import attrgetter
from itertools import combinations
from copy import deepcopy
''' Mat: tensor product of Pauli matrices
Mat.Xs :: frozenset : collection of sites of X gates
Mat.Zs :: frozenset : collection of sites of Z gates
'''
class Mat:
    def __init__(self, Xs, Zs):
        self.Xs = Xs
        self.Zs = Zs
        self._ipower = None
        self._key = None
    def __repr__(self):
        return '<Xs:%s Zs:%s>' % (sorted(list(self.Xs)), sorted(list(self.Zs)))
    def __hash__(self):
        if self._key is None:
            self._key = hash((self.Xs, self.Zs))
        return self._key
    def __eq__(self, other):
        return self.Xs == other.Xs and self.Zs == other.Zs
    def __neq__(self, other):
        return self.Xs != other.Xs or self.Zs != other.Zs
    def ipower(self): # number of overlap between Xs and Zs (num of Y gates)
        if self._ipower is None: # if ipower has not been calculated
            self._ipower = len(self.Xs & self.Zs)
            # once calculated the result is stored in self._ipower
        return self._ipower
# use mkMat to construct Mat
def mkMat(*arg):
    l_arg = len(arg)
    if l_arg == 2:
        return Mat(frozenset(arg[0]),frozenset(arg[1]))
    elif l_arg == 1:
        inds = arg[0]
        Xs = set()
        Zs = set()
        if isinstance(inds, dict): # dict of inds rules
        # example: mkMat({i:mu, ...})
            for (i, mu) in inds.items():
                if mu == 1:
                    Xs.add(i)
                elif mu == 3:
                    Zs.add(i)
                elif mu == 2:
                    Xs.add(i)
                    Zs.add(i)
        elif isinstance(inds, (tuple, list)): # list of inds
        # example: mkMat([mu0, mu1, mu2, ...])
            for (i, mu) in enumerate(inds):
                if mu == 0:
                    continue
                elif mu == 1:
                    Xs.add(i)
                elif mu == 3:
                    Zs.add(i)
                elif mu == 2:
                    Xs.add(i)
                    Zs.add(i)
        return Mat(frozenset(Xs), frozenset(Zs))
    elif l_arg == 0: # empty Mat by mkMat()
        return Mat(frozenset(), frozenset())
    else:
        raise TypeError('mkMat expected at most 2 arguments, got %s.' % l_arg)
# commutativity check
def is_commute(mat1, mat2):
    return (len(mat1.Xs & mat2.Zs) - len(mat1.Zs & mat2.Xs))%2 == 0
# merging Pauli indices (coefficient not determined here)
def pdot(mat1, mat2):
    return Mat(mat1.Xs ^ mat2.Xs, mat1.Zs ^ mat2.Zs)
''' Term: a Mat with coefficient and position
Term.mat :: Mat : matrix of Pauli operator
Term.val :: numeric : coefficient
Term.pos :: int : my position in Ham.terms
'''
class Term:
    def __init__(self, *arg):
        l_arg = len(arg)
        if l_arg == 2:
            self.mat, self.val = arg
        elif l_arg == 1:
            self.mat = arg[0]
            self.val = 1.        
        elif l_arg == 0:
            self.mat = mkMat()
            self.val = 1.
        self.pos = 0
    def __repr__(self):
        return '%s %s' % (self.val, self.mat)
# dot product of two terms
def dot(term1, term2):
    mat1 = term1.mat
    mat2 = term2.mat
    mat = pdot(mat1, mat2)
    n = mat1.ipower() + mat2.ipower() - mat.ipower()
    n = n + 2*len(mat1.Zs & mat2.Xs)
    s = (-1)**(n/2)
    term = Term(mat, s*term1.val*term2.val)
    return term
# dot product of two terms (times additional i)
def idot(term1, term2):
    mat1 = term1.mat
    mat2 = term2.mat
    mat = pdot(mat1, mat2)
    n = mat1.ipower() + mat2.ipower() - mat.ipower()
    n = n + 2*len(mat1.Zs & mat2.Xs) + 1
    s = (-1)**(n/2)
    return Term(mat, s*term1.val*term2.val)
''' Ham: a collection of Terms
Ham.terms :: list : terms stored in binary heap structure
Ham.mats  :: dict : mapping mat to term
Ham.imap  :: dict : mapping site to covering terms
'''
class Ham:
    def __init__(self, *arg):
        self.terms = []
        self.mats = {}
        self.imap = {}
        if len(arg) == 1:
            self.extend(arg[0])
    def __repr__(self):
        return '%s' % self.terms
    def __len__(self):
        return len(self.terms)
    def __bool__(self):
        return bool(self.terms)
    def __iter__(self):
        return iter(self.terms)
    # add a term to the heap tree (self.terms)
    def terms_push(self, term):
        pos = len(self.terms) # set pos to the end of self.terms
        term.pos = pos
        self.terms.append(term) # append from IR end
        self.terms_shiftUV(pos) # shifted to UV
    # adjust the position of a term in the heap tree
    def terms_adjust(self, term):
        pos = term.pos
        self.terms_shiftUV(pos)
        self.terms_shiftIR(pos)
    # shifting a term indexed by pos in the heap tree towards UV (upward)
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
    # shifting a term indexed by pos in the heap tree towards IR (downward)
    def terms_shiftIR(self, pos):
        terms = self.terms
        end_pos = len(terms) - 1
        this_term = terms[pos]
        child_pos = 2*pos + 1 # left child position
        while child_pos <= end_pos:
            # Set child_pos to index of larger child.
            rchild_pos = child_pos + 1 # right child position
            if rchild_pos <= end_pos and abs(terms[child_pos].val) < abs(terms[rchild_pos].val):
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
    # push a term into the Hamiltonian
    def push(self, term):
        if term.mat in self.mats: # if mat already exist
            old_term = self.mats[term.mat]
            old_term.val += term.val
            self.terms_adjust(old_term)
        else: # if mat is new
            self.terms_push(term)
            self.mats[term.mat] = term
            self.imap_add(term)
    # extend Hamiltonian by adding terms (given by iterator)
    def extend(self, terms):
        for term in terms:
            self.push(term)
    # remove a term from the Hamiltonian
    def remove(self, term):
        terms = self.terms
        end_pos = len(terms) - 1
        pos = term.pos
        del self.mats[term.mat]
        self.imap_del(term)
        if pos == end_pos:
            del terms[pos]
        elif 0 <= pos < end_pos:
            last_term = terms.pop()
            last_term.pos = pos
            terms[pos] = last_term
            self.terms_adjust(last_term)
    # perform C4 rotation generated by sgn*gen to Hamiltonian
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
            del mats[term.mat]
            self.imap_del(term)
            # C4 by idot with gen
            new_term = idot(term, gen)
            # update mat & val only
            term.mat = new_term.mat
            term.val = sgn * new_term.val
        # add new mats, NOT COMBINE TO ABOVE LOOP
        for term in relevant_terms:
            mats[term.mat] = term
            self.imap_add(term)
    # perform a series of C4 rotations Rs forward
    def forward(self, Rs):
        for R in Rs:
            self.C4(R)
    # perform a series of C4 rotations Rs backward
    def backward(self, Rs):
        for R in reversed(Rs):
            self.C4(R,-1)
''' Ent: calculate entanglement entropy of stablizers
Ent.mat2is :: dict : mapping from mat to the supporting sites
Ent.i2mats :: dict : mapping from site to the covering mat
Ent.subsys :: set  : entanglement subsystem (a set of sites)
Ent.shared :: set  : a set of mats shared between region and its complement
'''
import numpy as np
from fortran_ext import z2rank
class Ent:
    def __init__(self, taus):
        self.mat2is = {}
        self.i2mats = {}
        for term in taus:
            mat = term.mat
            sites = mat.Xs | mat.Zs
            self.mat2is[term.mat] = sites
            for i in sites:
                try:
                    self.i2mats[i].add(mat)
                except:
                    self.i2mats[i] = {mat}
        self.clear()
    def is_shared(self, mat):
        sites = self.mat2is[mat]
        return 0 < len(sites & self.subsys) < len(sites)
    def update_shared(self, sites):
        mats = set() # prepare to collect relevant mats
        for i in sites: # scan over relevant sites
            mats.update(self.i2mats[i]) # union into mats
        for mat in mats:
            if self.is_shared(mat): # if shared
                self.shared.add(mat) # add to shared
            else: # if not shared, discard if present in shared
                self.shared.discard(mat)
    # include sites to entanglement region
    def include(self, sites):
        self.subsys.update(sites)
        self.update_shared(sites)
    # exclude sites from entanglement region
    def exclude(self, sites):
        self.subsys.difference_update(sites)
        self.update_shared(sites)
    # clear
    def clear(self):
        self.subsys = set()
        self.shared = set()
    # return entropy of the entanglement region
    def entropy(self):
        mats = [Mat(mat.Xs & self.subsys, mat.Zs & self.subsys) for mat in self.shared]
        # mats is a list of Pauli monomials as generators
        n = len(mats) # get num of projected stablizers
        adj = np.zeros((n, n), dtype=int) # prepare empty adj mat
        # construct adj mat
        for k1 in range(n):
            for k2 in range(k1 + 1, n):
                if not is_commute(mats[k1], mats[k2]):
                    adj[k1, k2] = adj[k2, k1] = 1
        return z2rank(adj)/2
# half-system-size bipartite entropy (averaged over translation)
def bipartite_entropy(system):
    ent = Ent(system.taus)
    l_cut = 0
    L = int(system.size/2)
    S = 0
    ent.include(range(l_cut, l_cut + L))
    for l_cut in range(0, system.size):
        S += ent.entropy()
        ent.exclude({l_cut})
        ent.include({(l_cut + L) % system.size})
    return S/system.size
''' SBRG: doing RG, holding RG data and performing data analysis
SBRG.tol      :: float : terms with energy < leading energy * tol will be truncated
SBRG.max_rate :: float : each RG step allows at most (max_rate * num of off-diagonal terms) amount of new terms
SBRG.size     :: int : num of bits in the Hilbert space
SBRG.phybits  :: set : a collection of physical bits
SBRG.H        :: Ham : where the Hamiltonian is held and processed
SBRG.Hbdy     :: list : keep the original terms passed in with the model
SBRG.Hblk     :: list : holographic bulk Hamiltonian transformed by RCC
SBRG.Heff     :: list : terms in the effective Hamiltonian
SBRG.RCC      :: list : C4 transformations from beginning to end
SBRG.taus     :: Ham : stabilizers
SBRG.trash    :: list : hold the energy scales that has been truncated
'''
class SBRG:
    tol = 1.e-8
    max_rate = 2.
    def __init__(self, model):
        self.size = model.size
        self.phybits = set(range(self.size))
        self.H = Ham(deepcopy(model.terms))
        self.Hbdy = model.terms
        self.Hblk = None
        self.Heff = []
        self.RCC = []
        self.taus = None
        self.trash = []
    def findRs(self, mat):
        if len(mat.Xs) > 0: # if X or Y exists, use it to pivot the rotation
            pbit = min(mat.Xs) # take first off-diag qubit
            return ([idot(Term(mkMat(set(),{pbit})), Term(mat))], pbit)
        else: # if only Z
            if len(mat.Zs) > 1:
                for pbit in sorted(list(mat.Zs)): # find first Z in phybits
                    if (pbit in self.phybits):
                        tmp = Term(mkMat({pbit},set())) # set intermediate term
                        return ([idot(tmp, Term(mat)), idot(Term(mkMat(set(),{pbit})), tmp)], pbit)
            elif len(mat.Zs) == 1:
                pbit = min(mat.Zs)
        return ([], pbit)
    def perturbation(self, H0, offdiag):
        h0 = H0.val # set h0
        min_prod = abs(h0)**2*SBRG.tol # set minimal product
        # SiSj for commuting terms whose product val > min_prod
        SiSj = [dot(term1, term2) for (term1, term2) in combinations(offdiag, 2)
                if is_commute(term1.mat,term2.mat) and abs(term1.val*term2.val) > min_prod]
        SiSj.sort(key=attrgetter('val')) # sort by val
        # term number truncation
        max_len = round(SBRG.max_rate*len(offdiag))
        if len(SiSj) > max_len:
            self.trash.extend([term.val/h0 for term in SiSj[:-max_len]])
            SiSj = SiSj[-max_len:]
        # multiply by H0 inverse
        H0inv = Term(H0.mat,1/h0)
        pert = [dot(H0inv,term) for term in SiSj]
        # add backward correction
        var = sum((term.val)**2 for term in offdiag) # also used in error estimate
        pert.append(Term(H0.mat, var/(2*h0)))
        return pert
    def nextstep(self):
        if not (self.phybits and self.H): # return if no physical bits or no H
            self.phybits = set() # clear physical bits
            return self
        # get leading energy scale
        H0 = self.H.terms[0]
        h0 = H0.val
        if not abs(h0): # if leading scale vanishes
            self.phybits = set() # quench physical space
            return self
        # find Clifford rotations
        Rs, pbit = self.findRs(H0.mat)
        self.RCC.extend(Rs) # add to RCC
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
            if not ((term.mat.Xs | term.mat.Zs) & self.phybits):
                self.Heff.append(term)
                self.H.remove(term)
        return (Term(H0.mat,h0), Rs, offdiag)
    def flow(self, step = float('inf')):
        step = min(step, len(self.phybits)) # adjust RG steps
        # carry out RG flow
        stp_count = 0
        while self.phybits and stp_count < step:
            self.nextstep()
            stp_count += 1
    def make(self):
        # reconstruct stabilizers
        stabilizers = []
        blkbits = set(range(self.size))
        for term in self.Heff:
            if len(term.mat.Zs) == 1:
                stabilizers.append(deepcopy(term))
                blkbits -= term.mat.Zs
        stabilizers.extend(Term(mkMat(set(),{i}),0) for i in blkbits)
        self.taus = Ham(stabilizers)
        self.taus.backward(self.RCC)
        # reconstruct holographic bulk Hamiltonian
        self.Hblk = Ham(deepcopy(self.Hbdy))
        self.Hblk.forward(self.RCC)
    def run(self):
        self.flow()
        self.make()
        return self
    # calculate Anderson correlator between pairs in terms
    def correlate(self, terms):
        ops = Ham(terms)
        ops.forward(self.RCC)
        cor = {}
        L = self.size
        for (i,j) in combinations(range(len(ops)),2):
            if len(ops.terms[i].mat.Xs ^ ops.terms[j].mat.Xs) == 0:
                d = int(abs((j - i + L/2)%L - L/2))
                cor[d] = cor.get(d,0) + 1
        return cor
''' Model: defines Hilbert space and Hamiltonian
Model.size  :: int : num of bits
Model.terms :: list : terms in the Hamiltonian
'''
class Model:
    def __init__(self):
        self.size = 0
        self.terms = []
# quantum Ising model
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
    model.size = L
    # translate over the lattice by deque rotation
    H_append = model.terms.append
    rnd_beta = random.betavariate
    for i in range(L):
        H_append(Term(mkMat({i: 1, (i+1)%L: 1}), para['J']*rnd_beta(alpha_J, 1)))
        H_append(Term(mkMat({i: 3, (i+1)%L: 3}), para['K']*rnd_beta(alpha_K, 1)))
        H_append(Term(mkMat({i: 3}), para['h']*rnd_beta(alpha_h, 1)))
    model.terms = [term for term in model.terms if abs(term.val) > 0]
    return model
# XYZ model
def XYZ(L, **para):
    # L - number of sites (assuming PBC)
    # model - a dict of model parameters
    try: # set parameter alpha
        alpha = para['alpha']
        alpha_X = alpha
        alpha_Y = alpha
        alpha_Z = alpha
    except:
        alpha_X = para.get('alpha_x',1)
        alpha_Y = para.get('alpha_y',1)
        alpha_Z = para.get('alpha_z',1)
    model = Model()
    model.size = L
    # translate over the lattice by deque rotation
    H_append = model.terms.append
    rnd_beta = random.betavariate
    for i in range(L):
        H_append(Term(mkMat({i: 1, (i+1)%L: 1}), para['Jx']*rnd_beta(alpha_X, 1)))
        H_append(Term(mkMat({i: 2, (i+1)%L: 2}), para['Jy']*rnd_beta(alpha_Y, 1)))
        H_append(Term(mkMat({i: 3, (i+1)%L: 3}), para['Jz']*rnd_beta(alpha_Z, 1)))
    model.terms = [term for term in model.terms if abs(term.val) > 0]
    return model

# Toolbox 
# I/O 
# JSON pickle: export to communicate with Mathematica 
import jsonpickle
def export(filename, obj):
    with open(filename + '.json', 'w') as outfile:
        outfile.write(jsonpickle.encode(obj))
def export_Ham(filename, ham):
    export(filename, [[term.val,[list(term.mat.Xs),list(term.mat.Zs)]] for term in ham])
import pickle
# pickle: binary dump and load for python.
def dump(filename, obj):
    with open(filename + '.dat', 'bw') as outfile:
        pickle.dump(obj, outfile)
def load(filename):
    with open(filename + '.dat', 'br') as infile:
        return pickle.load(infile)