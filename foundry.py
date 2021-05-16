#!/usr/bin/env python

import os, sys, logging, codecs, numpy as np
from opticks.sysrap.OpticksCSG import CSG_ as CSG 
log = logging.getLogger(__name__)

class Node(object):
    def __init__(self, item, fd):
        self.fd = fd 
        self.tc = item[3,2].view(np.int32)
        self.ty = CSG.desc(self.tc)

        self.tref = item[3,3].view(np.int32) & 0x7fffffff 
        self.cm = bool(item[3,3].view(np.int32) & 0x8000000) 
        self.bb = item[2:4].ravel()[:6]
        self.pa = item[0:2].ravel()[:6]
        self.tr = self.fd.tran[self.tref-1] if self.tref > 0 else None
        self.it = self.fd.itra[self.tref-1] if self.tref > 0 else None

    def __repr__(self):
        return "Node %10s tref:%5d cm:%1d bb %30s pa %30s " % (self.ty, self.tref, self.cm, str(self.bb), str(self.pa))   


class Prim(object):
    def __init__(self, item, fd):
        self.fd = fd 
        self.numNode = item[0,0].view(np.int32)
        self.nodeOffset = item[0,1].view(np.int32)
        self.node = list(map(lambda item:Node(item, fd), self.fd.node[self.nodeOffset:self.nodeOffset+self.numNode]))
        self.bb = item[2:4].ravel()[:6]

    def __repr__(self):
        return "Prim %3d %5d %30s " % (self.numNode, self.nodeOffset, str(self.bb))   

    def __str__(self):
        return "\n".join([repr(self)]+list(map(repr, self.node)))

    def __getitem__(self, nodeIdx):
        return self.node[nodeIdx]

class Solid(object):
    @classmethod
    def DecodeLabel(cls, item):
        return item.tobytes()[:8].split(b'\0')[0].decode("utf-8")     # TODO: test full 8 char label  

    def __init__(self, item, fd):
        self.item = item  
        self.fd = fd 
        self.label = self.DecodeLabel(item)
        self.numPrim = item[0,2].view(np.int32)
        self.primOffset = item[0,3].view(np.int32)
        self.prim = list(map(lambda item:Prim(item, fd), self.fd.prim[self.primOffset:self.primOffset+self.numPrim]))
        self.ce = item[1].view(np.float32)

    def __repr__(self):
        return "Solid %10s : %4d %5d  ce %35s  "  % (self.label, self.numPrim, self.primOffset, str(self.ce))

    def __str__(self):
        return "\n".join([repr(self)]+list(map(repr, self.prim)))

    def __getitem__(self, primIdx):
        return self.prim[primIdx]


class Foundry(object):
    BASE = "$TMP/CSG_GGeo/CSGFoundry"
    NAMES = "solid prim node tran itra inst".split()
    def __init__(self, base=None):
        if base is None: base = self.BASE  
        base = os.path.expandvars(base) 
        for n in self.NAMES+["plan"]:
            path = os.path.join(base, "%s.npy" % n)
            if not os.path.exists(path): continue
            setattr( self, n, np.load(path))
        pass   
        self.name = np.loadtxt(os.path.join(base, "name.txt"), dtype=np.object)
        self.label = list(map(Solid.DecodeLabel, self.solid))
        self.solids = list(map(lambda item:Solid(item,self), self.solid))

    def __repr__(self):
        return "\n".join(["Foundry"] + list(map(lambda n:"%10s : %s " % (n,repr(getattr(self, n).shape)), self.NAMES))) 

    def __str__(self):
        return "\n".join(["Foundry"]+ list(map(repr, self.solids)))

    def index(self, solid_label):
        return self.label.index(solid_label) 

    def __getitem__(self, arg):
        """
        Access solids nodes prims via string specification::

             fd["d1"]        # solid with label "d1"
             fd["d1/0"]      # 1st prim
             fd["d1/0/0"]    # 1st node of 1st prim

        """
        ret = None
        elem = arg.split("/")

        solid_label = elem[0] if len(elem) > 0 else None
        primIdx = int(elem[1]) if len(elem) > 1 else None
        nodeIdx = int(elem[2]) if len(elem) > 2 else None

        solidIdx = self.index(solid_label) if not solid_label is None else None

        so = None
        pr = None
        nd = None

        if not solidIdx is None:
            so = self.solids[solidIdx]
        pass  
        if not so is None and not primIdx is None:
            pr = so[primIdx]
        pass
        if not pr is None and not nodeIdx is None:
            nd = pr[nodeIdx]
        pass

        if not nd is None:
            return nd
        elif not pr is None:
            return pr
        elif not so is None:
            return so
        else:
            return None
        pass
        return None



if __name__ == '__main__':
    fd = Foundry("$TMP/CSG_GGeo/CSGFoundry")
    print(repr(fd))
    #print(str(fd))
    args = sys.argv[1:] if len(sys.argv) > 1 else "d1 d2 d3 d4".split()
    for arg in args: 
        obj = fd[arg] 
        print(arg)
        print(obj)
    pass


