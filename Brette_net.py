# -*- coding: utf-8 -*-

from __future__ import (print_function, division, absolute_import,
                        unicode_literals)

from neuron import h

import netutils as nu
import netutils.raster as nur

from COBA_net import COBANet


class BretteNet(COBANet):
    """
    A random network based on the COBA network from Vogels and Abbott 2005
    from the HOC implementation of Brette et al 2007
    """
    def __init__(self, NCELL=4000, make_cell_fcn=h.CobaCell,
                 inhibitory=20, ex_con=2, inh_con=2, syn_delay=0, stim_delay=1,
                 ampa_gmax=0.006, gaba_gmax=0.067, crls=1, rrls=2):

        super().__init__(NCELL, make_cell_fcn)

    def make_stim(self):
        # random external input for first 50 ms
        # The first N_STIM (excitatory) cells are stimulated.
        self.N_STIM = int(self.NCELL / 50)    # number of neurons stimulated
        self.STOPSTIM = 50                  # duration of stimulation (ms)
        self.NSYN_STIM = 20                 # nb of stim (exc) synapses per neuron
        self.STIM_INTERVAL = 70             # mean interval between stims (ms)
        self.stimList = list()
        self.ncStimList = list()
        self.stvec = h.Vector()
        self.svec = h.Vector()

        h.mcell_ran4_init(self.run_random_low_start)
        for i in range(self.N_STIM):
            # The ith cell and its corresponding RandomStream.
            cell = self.cellList[i]
            rs = self.ranList[i]            # The ith cell and its corresponding RandomStream.
            stim = h.NetStim()
            stim.interval = self.STIM_INTERVAL
            stim.number = 1000              # but will shut off after STOPSTIM
            stim.noise = 1
            stim.start = 0
            stim.noiseFromRandom(rs.r)
            rs.r.negexp(1)
            rs.start()
            nc = h.NetCon(stim, cell.synlist.object(0))
            nc.delay = self.STIM_DELAY      # =1 to be consistent with Brette
            nc.weight[0] = self.AMPA_GMAX
            nc.record(self.stvec, self.svec, i)
            self.ncStimList.append(nc)
            self.stimList.append(stim)

        stim = h.NetStim()                  # will turn off all the others
        stim.number = 1
        stim.start = self.STOPSTIM
        for st in self.stimList:
            nc = h.NetCon(stim, st)
            nc.delay = 1    # to be consistent with Brette
            nc.weight[0] = -1
            self.ncStimList.append(nc)

        self.stimList.append(stim)

    def print_stim(self):
        print("\nthere are ", len(self.stvec), "stim")
        print("\ntime\t nc")
        for (st, sid) in zip(self.stvec, self.svec):
            print("%.5f\t %d" % (st, int(sid)))


def simulate():
    """
    this test duplicates Brette et al 2007 COBA network
    Only 100 ms are simulated to save time
    The results below are taken from running ModelDB code
    """

    #make the net
    h('random_stream_offset_ = 500')    # to match the statement in Brette /commmon/init.hoc
                                        # set as 1000 in Brette /common/ranstream.hoc
    net=BretteNet()

    #instumentation
    recording = nu.SpikeRecorder(net.cell_list())
    sc = nu.SimulationController()

    #run the simulation
    sc.set_fixed_dt(0.1)
    #sc.set_nthread(2)
    sc.set_v_init(-60)
    sc.set_celsius(36)
    sc.stdrun(100)

    #output the results
    net.print_stim()
    recording.print_spikes()
    print("all done")
    nur.draw_raster_plot(recording, display=True)


if __name__ == "__main__":
    simulate()

""" simulate should generate the following output

...    ...
99.8	1061
99.8	1230
99.8	2390
99.8	2498
99.8	3082
99.8	3168
99.8	3665
99.9	148
99.9	334
99.9	1665
99.9	1870
99.9	1897
99.9	2183
99.9	2525
99.9	2793
99.9	2831
99.9	3171
all done
"""
