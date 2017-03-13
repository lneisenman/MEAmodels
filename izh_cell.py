# -*- coding: utf-8 -*-

from __future__ import (print_function, division, absolute_import,
                        unicode_literals)

from neuron import h, load_mechanisms


load_mechanisms('mod_files')    # load the model dll files


class IzhCell():
    """
    A python version of the icell template from Ted's Izhmodel with the synapse
    model from Brette et al 2007
    """
    def __init__(self):
        self.soma = h.Section(name="soma")
        h.pt3dclear()
        h.pt3dadd(0, 0, 0, 1)
        h.pt3dadd(15, 0, 0, 1)
        self.x = self.y = self.z = 0
        self.soma.L = 10
        self.soma.diam = 3.1831
        self.soma.Ra = 100
        self.soma.cm = 1
        self.soma.insert('pas')
        self.soma.g_pas = 5e-5
        self.soma.e_pas = -65
        self.izh = h.Izh(0.5, sec=self.soma)
        self.all = h.SectionList()
        self.all.append()
        self.make_synapses()

    def position(self, x, y, z):
        for i in range(h.n3d()):
            h.pt3dchange(i, x+self.x+h.x3d(i), y+self.y+h.y3d(i),
                         z+self.z+h.z3d(i), h.diam3d(i))
            self.x = x
            self.y = y
            self.z = z

    def is_art(self):
        return 0

    def make_synapses(self):
        self.synlist = list()
        syn = h.ExpSyn(0.5, sec=self.soma)
        syn.tau = 5
        self.synlist.append(syn)

        syn = h.ExpSyn(0.5, sec=self.soma)
        syn.tau = 10
        syn.e = -80
        self.synlist.append(syn)

    def connect2target(self, target):
        self.soma.push()
        nc = h.NetCon(self.soma(0.5)._ref_v, target)
        h.pop_section()
        return nc
    
    def class1(self):
        ''' set params for class 1 behavior '''
        self.izh.a = 0.02
        self.izh.b = -0.1
        self.izh.c = -55
        self.izh.d = 6
        self.izh.e = 0.04
        self.izh.f = 4.1
        self.izh.g = 108
        self.izh.delay = 10
        self.izh.dur = 1e9
        self.izh.amp = 0.34

    def class2(self):
        ''' set params for class 2 behavior '''
        self.izh.a = 0.2
        self.izh.b = 0.26
        self.izh.c = -65
        self.izh.d = 0
        self.izh.e = 0.04
        self.izh.f = 5
        self.izh.g = 140
        self.izh.delay = 20
        self.izh.dur = 1e9
        self.izh.amp = 0.3412

if __name__ == '__main__':
    print('hello')
    cell = IzhCell()
    print(cell.soma.diam, cell.soma.cm)
    h.psection()
