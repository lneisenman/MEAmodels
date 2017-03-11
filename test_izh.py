# -*- coding: utf-8 -*-

from __future__ import (print_function, division, absolute_import,
                        unicode_literals)

import matplotlib.pyplot as plt
from neuron import h
import pandas as pd

import netutils as nu
import netutils.raster as nur

from izh_cell import IzhCell
from Izh_net import IzhNet


def simulate():
    ''' create a neuron '''

    cell = IzhCell()
    cell.class2()

    tvec = h.Vector()
    tvec.record(h._ref_t)
    vvec = h.Vector()
    vvec.record(cell.soma(0.5)._ref_v)
    h.v_init = -60
    h.dt = 0.1
    h.tstop = 200
    h.run()
    plt.plot(tvec, vvec)
    plt.show()


def simulate_pair():

    cell1 = IzhCell()
    cell1.class2()
    cell2 = IzhCell()
    cell2.class1()
    nc1 = cell1.connect2target(cell2.synlist[0])
    nc1.weight[0] = 0.0002
    nc1.delay = 1

    tvec = h.Vector()
    tvec.record(h._ref_t)
    vvec = h.Vector()
    vvec.record(cell2.soma(0.5)._ref_v)
    h.v_init = -60
    h.dt = 0.1
    h.tstop = 1000
    h.run()
    plt.plot(tvec, vvec)
    plt.show()


def simulate_net(stim_percentage=5):

    net = IzhNet(500, inhibitory=20, ex_con=20, inh_con=50, ampa_gmax=0.008)

    # instumentation
    recording = nu.SpikeRecorder(net.cell_list())
    sc = nu.SimulationController()
#    for cell in net.cell_list():
#        cell.set_bias(0.1)
#    net.cell_list()[0].set_bias(0.1001)
#    net.cell_list()[-1].set_bias(0.100003)
#    net.cell_list()[5].set_bias(0.100001)
#    net.cell_list()[-5].set_bias(0.1000002)

    stims = list()
    ncstims = list()
    num_stim = int(net.NCELL*stim_percentage/100)
    for i in range(num_stim):
        stim = h.NetStim(0.5)
        stim.interval = 50000 / (num_stim - i)
        stim.noise = 0.5
        stims.append(stim)
        stim.number = 10e12
        ncstim = h.NetCon(stim, net.cell_list()[int(net.NCELL*i/num_stim)].synlist.object(0))
        ncstim.weight[0] = 0.0002
        ncstim.delay = 0
        ncstims.append(ncstim)

    tvec = h.Vector()
    tvec.record(h._ref_t)
    vvec = h.Vector()
    vvec.record(net.cell_list()[0].soma(0.5)._ref_v)

    # run the simulation
    sc.set_fixed_dt(0.1)
    sc.set_v_init(-60)
    sc.set_celsius(36)
    sc.stdrun(60000)

    # output the results
#    recording.print_spikes()
    print("all done")
#    plt.plot(tvec, vvec)
#    plt.figure()
#    recording.save_spikes('baseline.csv')
    nur.draw_raster_plot(recording)
    nur.draw_raster_plot(read_baseline(), color='r', display=True)


def read_baseline(file_name='baseline.csv'):
    df = pd.read_csv(file_name)
    df.rename(index=str, columns={'time': 'tvec', 'cell id': 'idvec'},
              inplace=True)
    return df


if __name__ == "__main__":
    ''' '''
    h.load_file("stdrun.hoc")
#    simulate()
#    simulate_pair()
    simulate_net()
#    read_baseline()
