# -*- coding: utf-8 -*-

from __future__ import (print_function, division, absolute_import,
                        unicode_literals)

import matplotlib.pyplot as plt
from neuron import h
import numpy as np

import netutils as nu
import netutils.raster as nur

from izh_cell import IzhCell
from ranstream import RandomStream


class IzhNet(nu.BaseNet):
    """
    A random network usings Ted's Izh cells and InGauss noise mod files
    to replicate Izhikevich 2003 Figure 3
    """
    def __init__(self, NCELL=1000, make_cell_fcn=IzhCell,
                 inhibitory=20, syn_delay=0, noise=0.01,
                 ampa_gmax=0.01, gaba_gmax=0.02):

        super().__init__(NCELL, make_cell_fcn)
        self.INHIBITORY = inhibitory
        self.N_I = int(self.NCELL*inhibitory/100)       # Number of inibitory cells
        self.N_E = self.NCELL - self.N_I                # Number of excitatory cells
        self.AMPA_GMAX = ampa_gmax                      # (uS)
        self.GABA_GMAX = gaba_gmax                      # (uS)
        self.SYN_DELAY = syn_delay                      # (ms)
        self.NOISE = noise

        self.insert_ranstreams()
        self.set_cell_properties()
        self.connect_cells()
        self.insert_noise()
        self.make_stim()

    def insert_ranstreams(self):
        for i, cell in enumerate(self.cellList):
            cell.rs = RandomStream(i)

    def set_cell_properties(self):
        for i in range(self.N_E):
            cell = self.cellList[i]
            cell.rs.uniform(0, 1)
            r = cell.rs.repick()
            (c, d) = np.asarray([-65., 8.]) + np.asarray([15., -6.])*r*r
            cell.set_params(a=0.02, b=0.2, c=c, d=d)

        for i in range(self.N_E, self.NCELL):
            cell = self.cellList[i]
            cell.rs.uniform(0, 1)
            r = cell.rs.repick()
            (a, b) = np.asarray([0.02, 0.25]) + np.asarray([0.08, -0.05])*r*r
            cell.set_params(a=a, b=b, c=-65, d=2)

    def insert_noise(self):
        for cell in self.cellList:
            cell.noise = h.InGauss(0.5, sec=cell.soma)
            cell.noise.mean = 0
            cell.noise.stdev = self.NOISE
            cell.noise.delay = 0
            cell.noise.dur = 1e9
            cell.rs.normal(0, 1)
            cell.noise.noiseFromRandom(cell.rs.r)

        for i in range(self.N_E, self.NCELL):
            self.cellList[i].noise.stdev = 0.6 * self.NOISE

    def connect_cells(self):
        self.ncList = list()
        for i in range(self.N_E):
            cell1 = self.cellList[i]
            cell1.rs.uniform(0, 1)
            for j in range(self.NCELL):
                if not (i == j):
                    nc = cell1.connect2target(self.cellList[j].synlist[0])
                    nc.delay = 0
                    nc.weight[0] = self.AMPA_GMAX*cell1.rs.repick()
                    self.ncList.append(nc)

        # inhibitory connections
        for i in range(self.N_E, self.NCELL):
            cell1 = self.cellList[i]
            cell1.rs.uniform(0, 1)
            for j in range(self.NCELL):
                if not (i == j):
                    nc = cell1.connect2target(self.cellList[j].synlist[1])
                    nc.delay = 0
                    nc.weight[0] = self.GABA_GMAX*cell1.rs.repick()
                    self.ncList.append(nc)

    def make_stim(self):
        ''' get rid of stim on cell 0 '''
        pass


if __name__ == '__main__':
    h.load_file("stdrun.hoc")

    net = IzhNet(ampa_gmax=0.0000125, gaba_gmax=0.0005, noise=0.02)
    recording = nu.SpikeRecorder(net.cell_list())
    sc = nu.SimulationController()
    tvec = h.Vector()
    tvec.record(h._ref_t)
    vvec = h.Vector()
    vvec.record(net.cell_list()[300].soma(0.5)._ref_v)
    sc.set_fixed_dt(0.1)
    sc.set_v_init(-70)
    sc.stdrun(10000)
    plt.plot(tvec, vvec)
    plt.figure()
    nur.draw_raster_plot(recording, display=True)
