#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:51:43 2020

Loads file from 'raw data file'

Manual sweep selection

Creates pandas dataframe for each file where each row is a sweep, columns have the trace from each cell and quality control information

Applies quality control criteria and appends "true/false" columns to dataframe depending on whether the sweep passed QC

Saves pandas dataframe as a HDF file in a folder corresponding to that experiment (e.g. pair of cells)

Saves list of passed sweeps to excel file

@author: gemmagothard
"""


import matplotlib.pyplot as plt
import pyabf.plot
import os
import numpy as np
from openpyxl import load_workbook
import pandas as pd
import scipy.signal as spsig

#%%

animal_type = 'GTST LED'

wb = load_workbook(filename = '/Users/gemmagothard/Documents/DPhil/Ephys analysis 2020/Patch data input.xlsx')
ws = wb[animal_type]


#% USER INPUT

mV = 10 # voltage step mV used to check access resistance

# pyabf inputs:
cell1_pA = 0 
cell2_pA = 2
cell1_mV = 1

# samples per second
samples_p_s = 4.096

 # data points either side of maximum to take as ipsc window
ipsc_window = int(np.round(samples_p_s * 1,0))

full_path = os.path.join('/Volumes/PHARM_AKERMAN/GG/Patch_data',animal_type,'Raw data')
                

# list of folders where each folder is one date of experiment, containing multiple cell pairs
folder_list = list((folder for folder in os.listdir(full_path) if folder.startswith('2')))


for folder in folder_list:
    
    folderpath = os.path.join(full_path,folder)
    
    file_list = list((file for file in os.listdir(folderpath) if file.endswith('abf')))
    
    
    step_voltage   = []
    step_currentc1 = []
    step_currentc2 = []


    for file in file_list:

        abf = pyabf.ABF(os.path.join(folderpath,file))
    
        stim_type = []
    
        # define time windows relevant for events in stimuli
        pt0     = int(0    * abf.dataPointsPerMs)
        pt04    = int(40   * abf.dataPointsPerMs)
        pt1     = int(100 * abf.dataPointsPerMs) 
        pt150   = int(150 * abf.dataPointsPerMs)
        pt2     = int(200 * abf.dataPointsPerMs)
        pt25    = int(250 * abf.dataPointsPerMs)
        pt3     = int(300 * abf.dataPointsPerMs) 
        pt4     = int(400 * abf.dataPointsPerMs)
        pt5     = int(500 * abf.dataPointsPerMs) 
        pt6     = int(600 * abf.dataPointsPerMs)
        pt7     = int(750 * abf.dataPointsPerMs) 
        pt8     = int(850 * abf.dataPointsPerMs) 
        pt9     = int(900 * abf.dataPointsPerMs)
        
        
        RA_check = 'N' # unless "RA check" is in stimulus file name, assuming no RA check was done in each sweep
        RM_check = 'N'
        
        # for each stimulus type, define the time points when the RA check is done
        if 'test' in file:
            stim_type = 'LED'
            RA1 = pt7
            RA2 = pt8
            
        if 'start' in file:
            RA1 = pt0
            RA2 = pt150
            
        if 'plasticity' in file:
            stim_type = 'Plasticity'
            RA_check = 'Y'
            RA1 = pt0
            RA2 = pt150
            
            RM_check = 'Y'
            RM1 = pt04
            RM2 = pt1
            
        if 'voltage' in file:
            stim_type = 'Voltage steps'
            RA1 = pt2
            RA2 = pt3
            RA_check = 'Y'
            
        if 'spontaneous' in file:
            stim_type = 'Spontaneous'
            RA1 = pt7
            RA2 = pt8
            
        if 'RA' in file:
            RA_check = 'Y'
            
        if 'RM' in file:
            RM_check = 'Y'
            RM1 = pt04
            RM2 = pt1
        
        
        data   = []
        result = []
        
        keep_mask = [True for i in range(len(abf.sweepList))]
        
        # Make data frame of all sweeps
        for count,sweep in enumerate(abf.sweepList):
        
            #%% Cell 1:
            abf.setSweep(sweepNumber = sweep, channel = cell1_pA)
            
            holding_currentc1 = np.median(abf.sweepY) 
            cell1             = abf.sweepY-holding_currentc1 # subtract holding current 
            max_vec           = cell1[pt2:pt3] # part of trace that contains IPSC peak
            
            if stim_type == 'Voltage steps':
                max_vec                 = cell1[pt4:pt6]
                step_currentc1          = np.average(cell1[pt3:pt4])
                
                #get voltage info:
                abf.setSweep(sweepNumber = sweep, channel = cell1_mV)
                holding_voltage          = np.average(abf.sweepY[pt0:pt2])
                step_voltage             = np.average(abf.sweepY[pt3:pt4]) - holding_voltage
                mV                       = step_voltage
                
            
            cell1_maxind = np.argmax(abs((max_vec))) # find index of IPSC peak 
            
            # find the average around the IPSC peak
            cell1_max    = np.average(max_vec[cell1_maxind-ipsc_window:cell1_maxind+ipsc_window])   # IPSC peak value
            
            # calculate access resistance
            RA_c1 = abs(mV/np.max(np.abs(cell1[RA1:RA2]))*1000) 
            
            if RM_check == 'Y':
                RM_c1 = abs(mV / np.average(cell1[RM1:RM2])*1000)
            if RM_check == 'N':
                RM_c1 = np.nan
                RM_c2 = np.nan
    
            # do exactly the same for cell 2
            #%% Cell 2:
            abf.setSweep(sweepNumber = sweep, channel = cell2_pA)
            
            holding_currentc2   = np.median(abf.sweepY)
            cell2               = abf.sweepY-holding_currentc2
            max_vec2            = cell2[pt2:pt3]
            
            if stim_type == 'Voltage steps':
                max_vec2        = cell2[pt4:pt6]
                step_currentc2  = np.average(cell2[pt3:pt4])
                
            cell2_maxind = np.argmax(abs((max_vec2)))
            cell2_max    = np.average(max_vec2[cell2_maxind-ipsc_window:cell2_maxind+ipsc_window]) 
        
            RA_c2 = abs(mV/np.max(np.abs(cell2[RA1:RA2]))*1000)

            if RM_check == 'Y':
                RM_c2 = abs(mV / np.average(cell2[RM1:RM2])*1000)
                
                        
            ratio       = cell1_max / cell2_max  # ratio between two IPSC peaks
            RA_diff     = abs(RA_c1 - RA_c2)  # difference between access resistance
            RA_percent  = abs( ( RA_diff / np.max([RA_c1, RA_c2]) ) * 100 ) # percentage change between RAs
            
            
            #%% command sweep 
            abf.setSweep(sweepNumber = sweep, channel = cell1_mV)
            pipette_offset           = np.average(abf.sweepY[pt0:pt1])
            holding_voltage          = abf.sweepY - pipette_offset

            #Plot sweep:
            plt.plot(abf.sweepX,cell2,color='blue')
            plt.plot(abf.sweepX,cell1,color='orange')
            plt.fill_between(abf.sweepX[RA1:RA2],max(cell2)) # check the area that is being used to calculate RA
            plt.title(file)
            plt.show()
            
            #%% Manual sweep selection - sweep will be plotted, to keep sweep press enter. To delete sweep type 'N'
            
            print(str("sweep number:" + str(count)))
            print("keep this sweep?")
            keep_sweep = input()
            
            if keep_sweep == 'N':
                keep_mask[count] = bool(0)
            
            #%% put everything into a dataframe
            data.append({'SweepNo':sweep,
                          'Trace1':cell1,
                          'Trace2':cell2,
                          'Timevec':abf.sweepX,
                          'Holding voltage':holding_voltage,
                          'Stim type': stim_type,
                          'Ra info': RA_check,
                          'Ra1': RA_c1,
                          'Ra2': RA_c2,
                          'Rm check':RM_check,
                          'RM1': RM_c1,
                          'RM2': RM_c2,
                          'Ra_percent_change': RA_percent,
                          'IPSC1': cell1_max,
                          'IPSC2': cell2_max,
                          'Ratio': ratio,
                          'Holding1': holding_currentc1,
                          'Holding2': holding_currentc2,
                          'Voltage steps':step_voltage,
                          'Step current 1': step_currentc1,
                          'Step current 2': step_currentc2})
    
        df = pd.DataFrame(data)
        
        # delete sweeps:
        df = df[keep_mask]
    
        #%% Exclude sweeps with certain criteria
        
        # Criteria: 
        #           Holding current must be under +/-300pA
        #           Ratio of IPSCs must not change by more than 20% (if applicable / LED is used)
        #           Ra is below 50 mOhm
        
        Ratio_changepos = df.Ratio > np.median(df.Ratio)-(np.median(df.Ratio)*0.2)
        Ratio_changeneg = df.Ratio < np.median(df.Ratio)+(np.median(df.Ratio)*0.2)
        
        Ratio_change = Ratio_changepos & Ratio_changeneg
        
        Holding1_pos = df.Holding1<300  
        Holding1_neg = df.Holding1>-300
        
        Holding1_pass = Holding1_pos & Holding1_neg
        
        Holding2_pos = df.Holding2<300  
        Holding2_neg = df.Holding2>-300
        
        Holding2_pass = Holding2_pos & Holding2_neg
        
        criteria = {'Ra1_pass':df.Ra1<50,
                    'Ra2_pass':df.Ra2<50,
                    'Holding1_pass':Holding1_pass,
                    'Holding2_pass':Holding2_pass,
                    'Ratio_change':Ratio_change }
        
        
        df2 = pd.DataFrame(criteria)
        result = pd.concat([df, df2], axis=1)
        
     #%% Save dataframe
     
        # save dataframe as HDF in folder for the cell pair
        save_folder_name = 'P ' + file[:10]
        save_file_path = os.path.join('/Volumes/PHARM_AKERMAN/GG/Patch_data',animal_type,'Pair pandas',save_folder_name)
        save_file_name = os.path.join(save_file_path,file)
        
        # creates a new folder with the cell pair ID as the folder title
        if not os.path.exists(save_file_path):
                os.makedirs(save_file_path)
        
        result.to_hdf(save_file_name,key='hdf',mode='w')



        

