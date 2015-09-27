from numpy import array
import random
from operator import attrgetter
from itertools import combinations

### Fortran extensions ###
# from fortran_ext import z2rank, ipu

### Pauli Operator ###
class Pauli:
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
        self.key = hash((tuple(self.Xs), tuple(self.Zs)))
    def __repr__(self):
        return "<Xs:%s Zs:%s>" % (self.Xs, self.Zs)
# commutativity check
def is_commute(mat1, mat2):
    return (len(mat1.Xs & mat2.Zs) - len(mat1.Zs & mat2.Xs))%2 == 0
# merging Pauli indices (coefficient not determined here)
def pdot(mat1, mat2):
    return Pauli(mat1.Xs ^ mat2.Xs, mat1.Zs ^ mat2.Zs)

### Pauli Monomial ###
class Term:
    def __init__(self, *arg):
        l_arg = len(arg)
        if l_arg == 1:
            self.mat = arg[0]
            self.val = 1.
        elif l_arg == 2:
            self.mat, self.val = arg
        elif l_arg == 0:
            self.mat = Pauli()
            self.val = 1.
        self.UV = None
        self.IR = None
    def __repr__(self):
        return "%s %s" % (self.val, self.mat)
    def __hash__(self):
        return self.mat.key
    def __cmp__(self, other):
        mk0 = self.mat.key
        mk1 = other.mat.key
        return (mk0 > mk1) - (mk0 < mk1)
    def __eq__(self, other):
        if isinstance(other, Term):
            return self.mat.key == other.mat.key
        else:
            return False
    def __ne__(self, other):
        if isinstance(other, Term):
            return self.mat.key != other.mat.key
        else:
            return True
    def __lt__(self, other): return self.mat.key < other.mat.key
    def __gt__(self, other): return self.mat.key > other.mat.key
    def __le__(self, other): return self.mat.key <= other.mat.key
    def __ge__(self, other): return self.mat.key >= other.mat.key
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

### Pauli Polynomial ###
class Poly:
    def __init__(self, *arg):
        self.isempty = True
        self.UVscale = None
        self.IRscale = None
        self.terms = []
        self.imap = dict()
        if len(arg) == 1:
            self.extend(arg[0])
    def __repr__(self):
        return "%s" % self.terms
    def _registerUV(self, C):
        Cval = abs(C.val)
        if self.UVscale == None and self.IRscale == None: # scale not established yet
            self.UVscale = C
            self.IRscale = C
            C.UV = None
            C.IR = None
        elif Cval >= abs(self.UVscale.val): # over UVscale
            C.UV = None
            C.IR = self.UVscale
            C.IR.UV = C
            self.UVscale = C
        elif Cval <= abs(self.IRscale.val): # under IRscale
            C.IR = None
            C.UV = self.IRscale
            C.UV.IR = C
            self.IRscale = C
        else: # within scales, search for position
            CUV = self.IRscale
            while Cval > abs(CUV.val):
                CUV = CUV.UV
            C.UV = CUV
            C.IR = CUV.IR
            C.UV.IR = C
            C.IR.UV = C
    def _updateUV(self, C):
        Cval = abs(C.val)
        CUV = C.UV
        if CUV != None and Cval > abs(C.UV.val): # if over UV
            if C.IR == None:
                C.UV.IR = None
                self.IRscale = CUV
            else:
                C.UV.IR, C.IR.UV = C.IR, C.UV
            while CUV != None and Cval > abs(CUV.val):
                CUV = CUV.UV
            if CUV == None:
                C.UV = None
                C.IR = self.UVscale
                C.IR.UV = C
                self.UVscale = C
            else:
                C.UV = CUV
                C.IR = CUV.IR
                C.IR.UV = C
                C.UV.IR = C
        CIR = C.IR
        if CIR != None and Cval < abs(C.IR.val): # if under IR
            if C.UV == None:
                C.IR.UV = None
                self.UVscale = CIR
            else:
                C.IR.UV, C.UV.IR = C.UV, C.IR
            while CIR != None and Cval < abs(CIR.val):
                CIR = CIR.IR
            if CIR == None:
                C.IR = None
                C.UV = self.IRscale
                C.UV.IR = C
                self.IRscale = C
            else:
                C.IR = CIR
                C.UV = CIR.UV
                C.UV.IR = C
                C.IR.UV = C
    def _imap_add(self, term):
        mat = term.mat
        sites = mat.Xs | mat.Zs
        for i in sites:
            try:
                self.imap[i].add(term)
            except:
                self.imap[i] = {term}
    def _imap_remove(self, term):
        mat = term.mat
        for i in mat.Xs | mat.Zs:
            self.imap[i].remove(term)
    def extend(self, other):
        # orther - a list of terms
        As = self.terms
        Bs = sorted(other)
        nA = len(As)
        nB = len(Bs)
        iA = 0
        iB = 0
        Cs = []
        while iA < nA and iB < nB:
            A = As[iA]
            B = Bs[iB]
            if A < B:
                Cs.append(A)
                C = A
                iA += 1
            elif A == B:
                A.val += B.val
                Cs.append(A)
                C = A
                self._updateUV(C)
                iA += 1
                iB += 1
            else: # A > B
                if len(Cs) == 0 or C != B: # if B is new
                    Cs.append(B)
                    C = B
                    self._imap_add(C)
                    self._registerUV(C)
                else: # if B existed as C
                    C.val += B.val
                    self._updateUV(C)
                iB += 1
        # As or Bs has been exhausted
        if iB >= nB: # B exhausted
            Cs.extend(As[iA:nA]) # dump A
        if iA >= nA: # A exhausted
            for B in Bs[iB:nB]:
                if len(Cs) == 0 or C != B: # if B is new
                    Cs.append(B)
                    C = B
                    self._imap_add(B)
                    self._registerUV(B)
                else: # if B existed as C
                    C.val += B.val
                    self._updateUV(C)
        self.terms = Cs
    def _terms_remove(self, term):
        mk0 = term.mat.key
        n = len(self.terms)
        a = 0
        b = n-1
        if self.terms[a].mat.key == mk0:
            del self.terms[a]
            return
        if self.terms[b].mat.key == mk0:
            del self.terms[b]
            return
        c = (a+b)//2
        while a < b:
            mkc = self.terms[c].mat.key
            if mkc > mk0:
                b = c
                c = (a+b)//2
            elif mkc < mk0:
                a = c
                c = (a+b)//2
            else:
                del self.terms[c]
                return
        del self.terms[a]
    def remove(self, term):
        # fix UV-IR chain
        if term.IR == None: # term is IRscale
            self.IRscale = term.UV
        else: # term has IR
            term.IR.UV = term.UV
        if term.UV == None: # term is UVscale
            self.UVscale = term.IR
        else: # term has UV
            term.UV.IR = term.IR
        self._imap_remove(term) # remove term from imap
        self._terms_remove(term) # remove from terms
    def ascending(self):
        terms = []
        app = terms.append
        term = self.IRscale
        while term != None:
            app(term)
            term = term.UV
        return terms
    def descending(self):
        terms = []
        app = terms.append
        term = self.UVscale
        while term != None:
            app(term)
            term = term.IR
        return terms
    def _C4(self, gen, sgn = +1):
        # collect terms to be transformed
        terms = set() # start with empty set
        for i in gen.mat.Xs | gen.mat.Zs: # supporting sites of gen
            if i in self.imap: # if i registered in imap
                terms.update(self.imap[i]) # add the related terms
        # filter out non-commuting terms
        terms = [term for term in terms if not is_commute(term.mat, gen.mat)]
        for term in terms: # now the term will be transformed
            self._imap_remove(term) # first remove from imap
            # C4 by idot with gen
            new_term = idot(term, gen)
            # update mat & val only, without affecting UV-IR
            term.mat = new_term.mat
            term.val = sgn * new_term.val
        # add new terms to imap
        # NOTE: be done after all terms are transformed, otherwise
        # the untransformed terms may block new terms to be added
        # if they are the same
        for term in terms:
            self._imap_add(term)
    def forward(self, Rs):
        for R in Rs:
            self._C4(R)
        if len(Rs) > 0:
            self.terms.sort()
    def backward(self, Rs):
        for R in reversed(Rs):
            self._C4(R,-1)
        if len(Rs) > 0:
            self.terms.sort()

### SBRG ###
class SBRG:
    def __init__(self, model):
        self.tol = 1.e-8
        self.max_rate = 2.
        self.L = model.L
        self.H = Poly(model.terms)
        self.Hbdy = model.terms
        self.Hblk = []
        self.Heff = []
        self.EHM = []
        self.phybits = set(range(self.L))
        self.bit = 0
        self.trash = []
    def _findR(self, mat):
        if len(mat.Xs) > 0: # if X or Y exists, use it to pivot the rotation
            self.bit = min(mat.Xs) # take first off-diag qubit
            return [idot(Term(Pauli([],[self.bit])), Term(mat))]
        else: # if only Z
            if len(mat.Zs) > 1:
                for self.bit in sorted(list(mat.Zs)): # find first Z in phybits
                    if (self.bit in self.phybits):
                        tmp = Term(Pauli([self.bit],[])) # set intermediate term
                        return [idot(tmp, Term(mat)), idot(Term(Pauli([],[self.bit])), tmp)]
            elif len(mat.Zs) == 1:
                self.bit = min(mat.Zs)
        return []
    def _perturbation(self, H0, offdiag):
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
    def _nextstep(self):
        if len(self.phybits) == 0: # return if no physical bits
            return self
        # get leading energy scale
        H0 = self.H.UVscale
        if abs(H0.val) == 0.: # if leading scale vanishes
            self.phybits = set() # quench physical space
            return self
        # find Clifford rotations
        Rs = self._findR(H0.mat)
        self.EHM.extend(Rs) # add to EHM
        self.H.forward(Rs) # apply to H
        # pick out offdiag terms
        offdiag = [term for term in self.H.imap[self.bit] if self.bit in term.mat.Xs]
        pert = self._perturbation(H0, offdiag) # 2nd order perturbation
        for term in offdiag:
            self.H.remove(term) # remove off-diagonal terms
        self.H.extend(pert) # add perturbation to H
        self.phybits.remove(self.bit) # reduce physical bits
        # remove identity terms in physical space
        for term in list(self.H.imap[self.bit]):
            if len(term.mat.Xs & self.phybits) + len(term.mat.Zs & self.phybits) == 0:
                self.Heff.append(term)
                self.H.remove(term)
    def flow(self, step = float('inf')):
        step = min(step, len(self.phybits)) # adjust RG steps
        # carry out RG flow
        stp_count = 0
        while (len(self.phybits) > 0 and stp_count < step):
            self._nextstep()
            stp_count += 1

### Model ###
class Model:
    def __init__(self):
        self.L = 0
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
    model.L = L
    # translate over the lattice by deque rotation
    H_append = model.terms.append
    rnd_beta = random.betavariate
    for i in range(L):
        H_append(Term(Pauli({i: 1, (i+1)%L: 1}), para['J']*rnd_beta(alpha_J, 1)))
        H_append(Term(Pauli({i: 3, (i+1)%L: 3}), para['K']*rnd_beta(alpha_K, 1)))
        H_append(Term(Pauli({i: 3}), para['h']*rnd_beta(alpha_h, 1)))
    model.terms = [term for term in model.terms if abs(term.val) > 0]
    return model