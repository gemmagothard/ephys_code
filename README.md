# ephys_code
Code used to analyse in vitro ephys data, from Akerman lab, Oxford. Code written by Gemma Gothard.

Used for a range of stimuli in paired whole-cell patch clamp recordings in voltage clamp including; opsin stimulation, spontaneous recordings, voltage steps

All stimuli are recorded with a 10mv voltage in every sweep to calculate access resistance offline

Preprocessing script:
Input: abf files of dual whole-cell recordings 
Preprocessing elements include; access resistance & membrane resistance calculation, opsin response calculation, manual sweep selection, 
Generation of quality control criteria 
Output: generates a pandas dataframe (saved as a hdf file) which contains all relevant information for each file and saves it inside a folder for each pair of cells

Timings of stimulations and voltage steps can be easily defined and manipulated within the preprocessing script






