# -*- coding: utf-8 -*-

from __future__ import (print_function, division, absolute_import,
                        unicode_literals)

from neuron import h, load_mechanisms

import netutils as nu

from izh_cell import IzhCell
from ranstream import RandomStream


#load_mechanisms('mod_files')    # load the model dll files
#h.load_file("stdrun.hoc")       # load the standard run libraries
#h.load_file("hoc_files/cobacell.hoc")     # load the cell model
#h.load_file("hoc_files/ranstream.hoc")    # to give each cell its own sequence generator


class IzhNet(nu.BaseNet):
    """
    A random network based on the COBA network from Vogels and Abbott 2005
    from the HOC implementation of Brette et al 2007 usings Ted's Izh cells
    """
    def __init__(self, NCELL=4000, make_cell_fcn=IzhCell,
                 inhibitory=20, ex_con=2, inh_con=2, syn_delay=0, stim_delay=1,
                 ampa_gmax=0.006, gaba_gmax=0.067):

        super().__init__(NCELL, make_cell_fcn)
        self.INHIBITORY = inhibitory
        self.N_I = int(self.NCELL*inhibitory/100)       # Number of inibitory cells
        self.N_E = self.NCELL - self.N_I                # Number of excitatory cells
        self.EX_CON = ex_con                            # excitatory Connection probability
        self.INH_CON = inh_con                          # inhibitory Connection probability
        self.C_I = int(self.N_I*inh_con/100)            # Number inhib synapses per neuron
        self.C_E = int(self.N_E*ex_con/100)             # Number excit synapses per neuron
        self.AMPA_GMAX = ampa_gmax                      # (uS)
        self.GABA_GMAX = gaba_gmax                      # (uS)
        self.SYN_DELAY = syn_delay                      # (ms)
        self.STIM_DELAY = stim_delay                    # (ms)
        if self.C_E == self.N_E:
            self.C_E -= 1

        if self.C_I == self.N_I:
            self.C_I -= 1

        self.insert_ranstreams()
        self.connect_cells()
        self.make_stim()

    def insert_ranstreams(self):
        for i, cell in enumerate(self.cellList):
            cell.rs = RandomStream(i)

    def connect_cells(self):
        self.ncList = list()
        for i, cell in enumerate(self.cellList):
            # excitatory connections
            syn = cell.synlist[0]    # excitatory synapse
            cell.rs.discunif(0, self.N_E-1)  # return integer in range 0..N_E-1
            u = list()
            while (len(u) < self.C_E):
                r = int(cell.rs.repick())
                # no self-connection, & only one connection from any source
                if (r != i):
                    if r not in u:
                        # set up connection from source to target
                        src = self.cellList[r]
                        nc = src.connect2target(syn)
                        self.ncList.append(nc)
                        nc.delay = self.SYN_DELAY
                        nc.weight[0] = self.AMPA_GMAX
                        u.append(r)

            # inhibitory connections
            if self.C_I < 1:
                continue

            syn = cell.synlist[1]    # inhibitory synapse
            cell.rs.discunif(self.N_E, self.NCELL-1)  # return integer in range N_E..NCELL-1
            u = list()
            length = len(u)
            while (length < self.C_I):
                r = int(cell.rs.repick())
                # no self-connection, & only one connection from any source
                if (r != i):
                    if r not in u:
                        # set up connection from source to target
                        src = self.cellList[r]
                        nc = src.connect2target(syn)
                        self.ncList.append(nc)
                        nc.delay = self.SYN_DELAY
                        nc.weight[0] = self.GABA_GMAX
                        u.append(r)

                length = len(u)

    def make_stim(self):
        ''' get rid of stim on cell 0 '''
        pass
