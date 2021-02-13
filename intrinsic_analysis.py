#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 13:05:31 2021

@author: gemmagothard
"""


import pyabf
import matplotlib.pyplot as plt
import efel
import numpy as np
import pandas as pd
from scipy.signal import chirp, find_peaks, peak_widths
import statistics as st

from intrinsic_functions import get_membrane_potential, get_input_resistance, get_ISIs, get_spikes, get_spike_rates, get_spike_threshold, get_peak_widths

file = "file/path"

# stim start and end in ms
stim_start = 650
stim_end = 4650

# sampling rate is 10kHz
samples_per_s = 10000
samples_per_ms = samples_per_s / 1000

# stim duration in seconds
stim_duration_s = (stim_end-stim_start) / samples_per_s

# load the file
abf = pyabf.ABF(file)

# channel 0 is mV, channel 1 is pA
mV = 0
pA = 1

# set these variables to zero
data = []
voltage = []
command = []
time = []
current_step = []

# Make array of all sweep voltage traces
for count,sweep in enumerate(abf.sweepList):
    
    abf.setSweep(sweepNumber = sweep, channel = mV)
    
    voltage.append(abf.sweepY)
    command.append(abf.sweepC)
    time.append(abf.sweepX)
    
    this_command = abf.sweepC
    
    current_step.append(np.median(this_command[stim_start:stim_end]))


#%% Get variables we are interested in 

# find resting membrane potential
RMP = get_membrane_potential(voltage,stim_start)

# find input resistance (mOhm) 
Input_R = get_input_resistance(RMP,voltage,current_step)

# find spikes
spikes, spike_amplitude = get_spikes(voltage, stim_start, stim_end)

# find spike threshold
spike_threshold_mV, spike_threshold_ind = get_spike_threshold(spikes,voltage,time,stim_start,stim_end,samples_per_ms)

# find spike widths
spike_width50, spike_width_heights, left_ips, right_ips = get_peak_widths(voltage,spikes,samples_per_ms,spike_threshold_ind,stim_start,stim_end)

# find ISIs
ISIs = get_ISIs(spikes,samples_per_ms)

# find spike rates
spike_rates_Hz = get_spike_rates(spikes,stim_duration_s)


#%% put everything into a dataframe
data = ({'SweepNo':abf.sweepList,
             'current_step':current_step,
             'membrane_potential':RMP,
             'input_resistance':Input_R,
             'ISIs':ISIs,
             'spike_rate_Hz':spike_rates_Hz,
             'spike_threshold':spike_threshold_mV,
             'spike_amplitude':spike_amplitude,
             'spike_width50':spike_width50})

results = pd.DataFrame(data)


