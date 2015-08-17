from numpy import array
import random
from operator import methodcaller

### Fortran extensions ###
from fortran_ext import z2rank, ipu
def overlap(lst1, lst2):
    ipu.setab(array(lst1), array(lst2))
    return ipu.overlap()
def merge(lst1, lst2):
    ipu.setab(array(lst1), array(lst2))
    ipu.merge()
    return list(ipu.c[:ipu.nc])
def combine(lst1, lst2):
    ipu.setab(array(lst1), array(lst2))
    ipu.combine()
    return list(ipu.c[:ipu.nc])
def intersect(lst1, lst2):
    ipu.setab(array(lst1), array(lst2))
    ipu.intersect()
    return list(ipu.c[:ipu.nc])

### Pauli Operator ###
class Pauli:
    def __init__(self, *arg):
        l_arg = len(arg)
        if l_arg == 2:
            (self.Xs, self.Zs) = arg
        elif l_arg == 1:
            inds = arg[0]
            self.Xs = []
            self.Zs = []
            if isinstance(inds, dict): # dict of inds rules
                for (i, mu) in inds.items():
                    if mu == 1:
                        self.Xs.append(i)
                    elif mu == 3:
                        self.Zs.append(i)
                    elif mu == 2:
                        self.Xs.append(i)
                        self.Zs.append(i)
                # dict not ordered, need to sort
                self.Xs.sort()
                self.Zs.sort()
            elif isinstance(inds, (tuple, list)): # list of inds
                for (i, mu) in enumerate(inds):
                    if mu == 0:
                        continue
                    elif mu == 1:
                        self.Xs.append(i)
                    elif mu == 3:
                        self.Zs.append(i)
                    elif mu == 2:
                        self.Xs.append(i)
                        self.Zs.append(i)
        self.ipower = overlap(self.Xs, self.Zs)
    def __repr__(self):
        return "<Xs:%s Zs:%s>" % (self.Xs, self.Zs)
# commutativity check
def is_commute(mat1, mat2):
    return (overlap(mat1.Xs, mat2.Zs) - overlap(mat1.Zs, mat2.Xs))%2 == 0
# merging Pauli indices (coefficient not determined here)
def pdot(mat1, mat2):
    return Pauli(merge(mat1.Xs, mat2.Xs), merge(mat1.Zs, mat2.Zs))

### Pauli Monomial ###
class Term:
    def __init__(self, *arg):
        l_arg = len(arg)
        if (l_arg == 1):
            self.mat = arg[0]
            self.val = 1.
        elif (l_arg == 2):
            self.mat, self.val = arg
        self.UV = None
        self.IR = None
    def __repr__(self):
        return "%s %s" % (self.val, self.mat)
    def __neg__(self):
        self.val = - self.val
        return self
    def __hash__(self):
        return hash(self.matkey())
    def __cmp__(self, other):
        mk0 = self.matkey()
        mk1 = other.matkey()
        if mk0 == mk1:
            return 0
        elif mk0 > mk1:
            return 1
        else:
            return -1
    def __eq__(self, other):
        if isinstance(other, Term):
            return self.__cmp__(other) == 0
        else:
            return False
    def __ne__(self, other):
        if isinstance(other, Term):
            return not self.__cmp__(other) == 0
        else:
            return True
    def __lt__(self, other): return self.__cmp__(other) < 0
    def __gt__(self, other): return self.__cmp__(other) > 0
    def __le__(self, other): return self.__cmp__(other) <= 0
    def __ge__(self, other): return self.__cmp__(other) >= 0
    def matkey(self):
        return (tuple(self.mat.Xs), tuple(self.mat.Zs))
# dot product of two terms
def dot(term1, term2):
    mat1 = term1.mat
    mat2 = term2.mat
    mat = pdot(mat1, mat2)
    n = mat1.ipower + mat2.ipower - mat.ipower
    n = n + 2*overlap(mat1.Zs, mat2.Xs)
    s = (-1)**(n/2)
    return Term(mat, s*term1.val*term2.val)
# dot product of two terms (times additional i)
def idot(term1, term2):
    mat1 = term1.mat
    mat2 = term2.mat
    mat = pdot(mat1, mat2)
    n = mat1.ipower + mat2.ipower - mat.ipower
    n = n + 2*overlap(mat1.Zs, mat2.Xs) + 1
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
            self.__add__(arg[0])
    def __repr__(self):
        return "%s" % self.terms
    def __add__(self, other):
        As = self.terms
        Bs = sorted(other)
        Cs = []
        Capp = Cs.append
        nA = len(As)
        nB = len(Bs)
        iA = 0
        iB = 0
        while iA < nA and iB < nB:
            A = As[iA]
            B = Bs[iB]
            mkA = A.matkey()
            mkB = B.matkey()
            print(mkA, mkB)
            if mkA < mkB:
                Capp(A)
                iA += 1
            elif mkA == mkB:
                A.val += B.val
                Capp(A)
                iA += 1
                iB += 1
                Aval = abs(A.val)
                # update UV/IR relation
                AUV = A.UV
                if AUV != None and Aval > abs(A.UV.val): # if over UV
                    if A.IR == None:
                        A.UV.IR = None
                        self.IRscale = AUV
                    else:
                        (A.UV.IR, A.IR.UV) = (A.IR, A.UV)
                    while AUV != None and Aval > abs(AUV.val):
                        AUV = AUV.UV
                    if AUV == None:
                        A.UV = None
                        A.IR = self.UVscale
                        A.IR.UV = A
                        self.UVscale = A
                    else:
                        A.UV = AUV
                        A.IR = AUV.IR
                        A.IR.UV = A
                        A.UV.IR = A
                AIR = A.IR
                if AIR != None and Aval < abs(A.IR.val): # if under IR
                    if A.UV == None:
                        A.IR.UV = None
                        self.UVscale = AIR
                    else:
                        (A.IR.UV, A.UV.IR) = (A.UV, A.IR)
                    while AIR != None and Aval < abs(AIR.val):
                        AIR = AIR.IR
                    if AIR == None:
                        A.IR = None
                        A.UV = self.IRscale
                        A.UV.IR = A
                        self.IRscale = A
                    else:
                        A.IR = AIR
                        A.UV = AIR.UV
                        A.UV.IR = A
                        A.IR.UV = A
            else: # mkA > mkB
                Capp(B)
                self._imap_add(B)
                iB += 1
                Bval = abs(B.val)
                if Bval >= abs(self.UVscale.val):
                    B.UV = None
                    B.IR = self.UVscale
                    B.IR.UV = B
                    self.UVscale = B
                elif Bval <= abs(self.IRscale.val):
                    B.IR = None
                    B.UV = self.IRscale
                    B.UV.IR = B
                    self.IRscale = B
                else:
                    BUV = self.IRscale
                    while Bval > abs(BUV.val):
                        BUV = BUV.UV
                    B.UV = BUV
                    B.IR = BUV.IR
                    B.UV.IR = B
                    B.IR.UV = B
        # As or Bs has been exhausted
        if iB >= nB: # B exhausted
            Cs.extend(As[iA:nA]) # dump A
        if iA >= nA: # A exhausted
            for B in Bs[iB:nB]:
                Capp(B)
                self._imap_add(B)
                Bval = abs(B.val)
                if self.isempty:
                    # initialize
                    self.UVscale = B
                    self.IRscale = B
                    B.UV = None
                    B.IR = None
                    self.isempty = False
                elif Bval >= abs(self.UVscale.val):
                    B.UV = None
                    B.IR = self.UVscale
                    B.IR.UV = B
                    self.UVscale = B
                elif Bval <= abs(self.IRscale.val):
                    B.IR = None
                    B.UV = self.IRscale
                    B.UV.IR = B
                    self.IRscale = B
                else:
                    BUV = self.IRscale
                    while Bval > abs(BUV.val):
                        BUV = BUV.UV
                    B.UV = BUV
                    B.IR = BUV.IR
                    B.UV.IR = B
                    B.IR.UV = B
        self.terms = Cs
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
    def _imap_add(self, term):
        mat = term.mat
        sites = combine(mat.Xs, mat.Zs)
        for i in sites:
            try:
                self.imap[i].add(term)
            except:
                self.imap[i] = {term}
    def _C4(self, gen):
        terms = set()
        for i in combine(gen.mat.Xs, gen.mat.Zs):
            if i in self.imap:
                terms.update(self.imap[i])
        terms = [term for term in terms if not is_commute(term.mat, gen.mat)]
        for term in terms:
            for i in combine(term.mat.Xs, term.mat.Zs):
                self.imap[i].remove(term)
            new_term = idot(term, gen)
            term.mat = new_term.mat
            term.val = new_term.val
        for term in terms:
            self._imap_add(term)
    def forward(self, Rs):
        for R in Rs:
            self._C4(R)
        if len(Rs) > 0:
            self.terms.sort()
    def backward(self, Rs):
        for R in reversed(Rs):
            self._C4(-R)
        if len(Rs) > 0:
            self.terms.sort()


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