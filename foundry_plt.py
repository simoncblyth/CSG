#!/usr/bin/env python

import os, sys, logging, numpy as np
log = logging.getLogger(__name__)
from foundry import Foundry 

try:
    import matplotlib.pyplot as plt 
    import matplotlib.lines as mlines
except ImportError:
    plt = None
pass


class BB(object):
    """

                +--------+       
                |        | 
      z1  +--------+     | 
          |     |  |     |
          |     +--|-----+   y1  
          |        |
      z0  +--------+     y0
         x0        x1
    """
    def __init__(self, bb=[-0.,-0.,-0.,0.,0.,0.]):
        self.x0 = bb[0]
        self.y0 = bb[1]
        self.z0 = bb[2]
        self.x1 = bb[3]
        self.y1 = bb[4]
        self.z1 = bb[5]

    def plot_xz(self, ax, c="r", l="-"):
        ax.add_line(mlines.Line2D([self.x0, self.x1], [self.z0, self.z0], c=c, linestyle=l))  # bottom horizonal
        ax.add_line(mlines.Line2D([self.x0, self.x0], [self.z0, self.z1], c=c, linestyle=l))  # left vertical
        ax.add_line(mlines.Line2D([self.x1, self.x1], [self.z0, self.z1], c=c, linestyle=l))  # right vertical
        ax.add_line(mlines.Line2D([self.x0, self.x1], [self.z1, self.z1], c=c, linestyle=l))  # top horizonal

    def adjust_xz(self, ax, scale=1.2):
        ax.set_aspect('equal')
        ax.set_xlim(scale*self.x0, scale*self.x1) 
        ax.set_ylim(scale*self.z0, scale*self.z1) 

    def include(self, other):
        if other.x0 < self.x0: self.x0 = other.x0
        if other.y0 < self.y0: self.y0 = other.y0
        if other.z0 < self.z0: self.z0 = other.z0

        if other.x1 > self.x1: self.x1 = other.x1
        if other.y1 > self.y1: self.y1 = other.y1
        if other.z1 > self.z1: self.z1 = other.z1


def plot_solid_bb(s, pnbb = True):
    numPrim = len(s.prim)

    if numPrim <= 5:
       layout = (2,3)
    else:
       layout = (2,4)
    pass
    plt.ion()

    fig, axs = plt.subplots(*layout)
    if not type(axs) is np.ndarray: axs = [axs]
    plt.suptitle("foundry_plt.py %s " % s.label)

    axs = axs.ravel() 
    color = "rgbcmyk"
    sbb = BB()
    for i in range(numPrim): 
        c = color[i % len(color)]
        p = s.prim[i]
        pbb = BB(p.bb)
        sbb.include(pbb)

        pbb.plot_xz(axs[0], c=c)     # collect all onto first 
        pbb.plot_xz(axs[1+i], c=c) 

        if pnbb:
            l = ":" # '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
            for n in p.node:
                nbb = BB(n.bb)
                nbb.plot_xz(axs[1+i], c=c, l=l)   
            pass 
        pass
    pass
    for ax in axs:
        sbb.adjust_xz(ax)
    pass
    fig.show()




if __name__ == '__main__':
    fd = Foundry()
    args = sys.argv[1:] if len(sys.argv) > 1 else "r1 d1".split() 
    for arg in args:
        s = fd[arg]
        plot_solid_bb(s)
    pass


