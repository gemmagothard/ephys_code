#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 15:34:54 2021

This script calculates contains functions necessary to calculate intrinsic spiking properties of neurons recording in current clamp during a series of current injections

"voltage" is a list of numpy arrays, where each list is a new sweep/current step of a single file

@author: gemmagothard
"""

import pyabf
import matplotlib.pyplot as plt
import efel
import numpy as np
import pandas as pd
from scipy.signal import chirp, find_peaks, peak_widths, peak_prominences



def get_membrane_potential(voltage,stim_start):
    
    RMP = [np.median(x[0:stim_start]) for x in voltage]
    
    
    return RMP


def get_input_resistance(RMP,voltage,current_step):
    
    # find voltage during current step:
    step_voltage = [np.median(x[2000:4000]) for x in voltage]
        
    # find difference in voltage between resting membrane potential and current step
    delta_v = [x - RMP[c] for c,x in enumerate(step_voltage)]
    
    # find resistance by dividing voltage difference by current step value (R = V / C)
    resistance = [delta_v[c] / current_step[c] *1000 for c,x in enumerate(delta_v)]
    
    # take the average of the first 3 current injections 
    resistance = np.nanmean(resistance[0:3])
    

    return resistance


def get_spikes(voltage,stim_start,stim_end):
    
    all_peaks = []
    all_peak_heights = []
    
    for x in voltage:

        peaks,peak_info = find_peaks(x,height=10,distance=50) 
        
        all_peaks.append(peaks)
        all_peak_heights.append(peak_info['peak_heights'])

    
    return all_peaks,all_peak_heights
    


    
def get_ISIs(spikes,samples_per_ms):
    
    ISIs_ms = []
    for i in spikes:
        ISIs = np.diff(i)
        
        ISIs_ms.append(ISIs / samples_per_ms)
        
    
    return ISIs_ms



def get_spike_rates(spikes,stim_duration_s):
    
    spike_rate_Hz = []
    for s in spikes:
    
        spike_rate_Hz.append(len(s) / stim_duration_s)
    
    return spike_rate_Hz
    
    
    
def get_spike_threshold(spikes,voltage,time,stim_start,stim_end,samples_per_ms):

    
    # spike threshold is defined as when dVm/dt crosses 25 preceding a spike
    threshold_dVmdt = 20
    
    spike_threshold_mV = []
    all_idx = []
    
    for c,s in enumerate(spikes):
        
        v = voltage[c][stim_start+50:stim_end]
        
        # find first derivative of membrane potential
        dVmdt = np.diff(v) * samples_per_ms

        # find where derivative crosses 25 (e.g. the threshold we are setting)
        idx = np.argwhere(np.diff(np.sign(dVmdt - threshold_dVmdt))).flatten()

        # find the voltage values are these indices - note these will have *some* offset because they are the closest time point before dVmdt crosses 25
        idx_v = v[idx] 
        
        # take the values where voltage is less than zero - this will be the crossing that happens just before the spike, instead of above zero which is after
        spike_threshold_mV.append(idx_v[idx_v<0])
        
        all_idx.append(idx+stim_start+50)


    return spike_threshold_mV, all_idx
    
    
    
    
def get_peak_widths(voltage,spikes,samples_per_ms,spike_threshold_ind,stim_start,stim_end):
    
    allspikes_width50 = []
    allspikes_width_h = []
    all_leftips = []
    all_rightips = []
    
    for c,s in enumerate(spikes):
        
        if len(spikes[c]) == 0:
            allspikes_width50.append(np.nan)
            allspikes_width_h.append(np.nan)
            all_leftips.append(np.nan)
            all_rightips.append(np.nan)
            continue
        
        idx = spike_threshold_ind[c][0]

        s = spikes[c]-idx
        v = voltage[c][idx:stim_end]

        width50, width_heights, left_ips, right_ips = peak_widths(v,s,rel_height=0.5)
        
        allspikes_width50.append(width50/samples_per_ms) 
        allspikes_width_h.append(width_heights)
        
        all_leftips.append(left_ips+idx)
        all_rightips.append(right_ips+idx)
    
    
    return allspikes_width50, allspikes_width_h, all_leftips, all_rightips







