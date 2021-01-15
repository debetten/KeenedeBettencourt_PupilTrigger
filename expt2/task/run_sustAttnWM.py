#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Runs the experiment.
Author: MdB 04/2017 s
'''

import time, pickle, csv, sys
import numpy as np
import os.path
from psychopy import core, gui, event, parallel,visual
from task_sustAttnWM import *
from datetime import datetime
from shutil import copyfile
import settings_sustAttnWM as settings
import helper_functions_sustAttnWM as hf
if sys.platform == 'win32':
    from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
    import pylink

#get basic specifications for this experiment from GUI
gui_specs = hf.get_specs()
print "Subject: ", gui_specs['subj_name']

#set up the experimental design
if not gui_specs['restarting'] or gui_specs['debug']:
    #create experiment design
    expt_design = settings.ExperimentalDesign(gui_specs)
else:
    #load existing design
    filename_expdat = gui_specs['save_dir'] + gui_specs['subj_name'] + '_expdat.p'
    filename_expdesign = gui_specs['save_dir']+ gui_specs['subj_name'] + '_expdesign.p'
    if os.path.isfile(filename_expdat):
        expt_design = pickle.load( open( filename_expdat, "rb" ) )
    elif os.path.isfile(filename_expdesign):
        expt_design = pickle.load( open( filename_expdesign, "rb" ) )
    else:
        expt_design = settings.ExperimentalDesign(gui_specs)
    
    #overwrite the specs to match what you specified in the gui for this time
    expt_design.restarting = gui_specs['restarting']
    expt_design.env = gui_specs['env']
    expt_design.exptDate = gui_specs['expt_date']
    expt_design.debug = gui_specs['debug']

#set up the experimental display
expt_display = settings.ExperimentalDisplay(expt_design,gui_specs)

#create a Task class variable with both the deisgn and display information
exp = Task(expt_design,expt_display)

#save files
filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_' + 'expdat.p'
pickle.dump( exp.dsgn, open( filename, "wb" ) )

#open csv files to save
f = hf.open_csv_data_file(gui_specs,gui_specs['subj_name'] + '_explog')
f_cd = hf.open_csv_data_file(gui_specs,gui_specs['subj_name'] + '_explog_cd')

#calculate the number of first trial
filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_trialnum.p'
if gui_specs['restarting'] and os.path.isfile(filename):
    trial0 = pickle.load( open( filename, "rb" ) )
    exp.dsgn.counter_trial = trial0

#start eyetracker
if exp.dsgn.eyetrack:
    tk = hf.eyetracker_settings(exp.disp,exp.dsgn)
    tk.doTrackerSetup()
else:
    tk = None
    
#if exp.dsgn.env=='eeg':
#    parallel.setPortAddress(address=53328)
#    parallel.setData(0)

if __name__ == '__main__':
    
    #initialize the clock
    clock = core.Clock()
    clock.reset()
    
    #set mouse to not be visible 
    if not gui_specs['debug']:
        exp.disp.win.setMouseVisible(False)

    #welcome screen
    exp.welcomeToExperiment(tk=tk)
    
    ######################################################################
    ## Encoding/WM Phase
    if gui_specs['sustattn']:
        
#        #working memory long-term memory instructions
#        if exp.dsgn.counter_trial[0]==0:
#            exp.instructionsSustAttn(clock,tk=tk)
#            exp.instructionsWholeReport(clock)
#            exp.instructionsBoth(clock,tk=tk)
        
        #block loop
        block0 = 0 #FIX THIS TO RESTART BLOCKS?
        for iblock in range(block0,exp.dsgn.nblocks):
            print iblock
            exp.startOfBlock(iblock,clock,f=f,tk=tk)
            if iblock == 0:
                hf.write_to_csv_data_file(f,0,0,exp.dsgn,headerLine=True)
                
            #set mouse to not be visible 
            if not gui_specs['debug']:
                exp.disp.win.setMouseVisible(False)
            
            #exp.dsgn.counter_trial[iblock]
            #trial loop
            for itrial in np.arange(exp.dsgn.counter_trial[iblock],exp.dsgn.ntrials_perblock):
                exp.encodingArray(iblock,itrial,clock,f=f,tk=tk)
                #print (str(itrial)+'\t'+str(exp.dsgn.freq_trials[iblock,itrial])+'\t'+str(exp.dsgn.acc[iblock,itrial])+'\t'+str(exp.dsgn.within_tol[iblock,itrial])+'\t'+str(exp.dsgn.rts_trailingavg[iblock,itrial])+'\t'+str(' fast ')+str(exp.dsgn.rts_runningfastthresh[iblock,itrial])+'\t'+str('slow ')+str(exp.dsgn.rts_runningslowthresh[iblock,itrial])+'\t'+str(np.nansum(exp.dsgn.probe_trials[iblock,:itrial])))+'\t'+str(np.nansum(exp.dsgn.rt_triggered_fast[iblock,:itrial]))+'\t'+str(np.nansum(exp.dsgn.rt_triggered_slow[iblock,:itrial]))
                if exp.dsgn.probe_trials[iblock,itrial]==1:
                    exp.retentionInterval(iblock,itrial,clock,tk=tk)
                    if exp.dsgn.within_tol[iblock,itrial]==1:
                        exp.memoryProbe(iblock,itrial,clock,f=f,tk=tk)
                    exp.itiWM(iblock,itrial,clock,tk=tk)
                    #print exp.dsgn.rt_triggered_fast+'\t'+exp.dsgn.rt_triggered_slow+'\t'+exp.dsgn.rt_triggered_med
                    
                #update counter
                exp.dsgn.counter_trial[iblock] = exp.dsgn.counter_trial[iblock]+1
               
                #save files
                filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_' + 'expdat.p'
                pickle.dump( exp.dsgn, open( filename, "wb" ) )
                hf.write_to_csv_data_file(f,iblock,itrial,exp.dsgn) 
            #exp.endOfBlock(iblock,clock,f=None,tk=None)
        
        
        #close and copy explog for this section of the experiment
        today = datetime.now()
        filename_old = gui_specs['save_dir'] + gui_specs['subj_name'] + '_explog.csv'
        filename_new = gui_specs['save_dir'] + gui_specs['subj_name'] + '_explog_' + today.strftime('%Y%m%d_%H%M%S') + '.csv'
        copyfile(filename_old,filename_new)
    
    
    ######################################################################
    ## Change Detection Phase
    if gui_specs['changedetect']:
        iblock = exp.dsgn.nblocks
        
        #open csv files to save this section of the experiment
        f = hf.open_csv_data_file(gui_specs,gui_specs['subj_name'] + '_explog_cd')
        
        if exp.dsgn.counter_trial[iblock]<exp.dsgn.cd_ntrials_perblock:
            exp.startOfBlock(exp.dsgn.nblocks,clock,f=None,tk=tk)
            
        if exp.dsgn.counter_trial[iblock]==0:
            exp.changeDetectionInstructionsIntro()
            exp.changeDetectionInstructionsPracticeSS1(tk=tk)
            exp.changeDetectionInstructionsPracticeSS2Diff(tk=tk)
            exp.changeDetectionInstructionsPracticeSS6Same(tk=tk)
            
        if exp.dsgn.counter_trial[iblock]<exp.dsgn.cd_ntrials_perblock:
            exp.startOfBlock(exp.dsgn.nblocks,clock,f=None,tk=tk)
            if iblock == 0:
                hf.write_to_csv_data_file_cd(f,0,0,exp.dsgn,headerLine=True)
                
        iblock = 0
        for itrial in range(exp.dsgn.counter_trial[-1],exp.dsgn.cd_ntrials_perblock):
            exp.changeDetectionEncodingArray(iblock,itrial,clock,tk=tk)       #encoding array
            exp.changeDetectionRetentionInterval(iblock,itrial,clock)   #retention interval
            exp.changeDetectionMemoryProbe(iblock,itrial,clock,f=f,tk=tk)     #memory probe
            
            #update counter
            exp.dsgn.counter_trial[-1] = exp.dsgn.counter_trial[-1]+1

            #save files
            filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_' + 'trialnum.p'
            pickle.dump( exp.dsgn.counter_trial, open( filename, "wb" ) )
            filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_' + 'expdat.p'
            pickle.dump( exp.dsgn, open( filename, "wb" ) )
            hf.write_to_csv_data_file_cd(f,iblock,itrial,exp.dsgn) 
            
            #iti
            exp.changeDetectionITI(iblock,itrial,clock,f=f,tk=tk)
            
        #save a copy of data with timestamp that cannot be overwritten
        today = datetime.now()
        filename_old = gui_specs['save_dir'] + gui_specs['subj_name'] + '_expdat.p'
        filename_new = gui_specs['save_dir'] + gui_specs['subj_name'] + '_expdat_' + today.strftime('%Y%m%d_%H%M%S') + '.p'
        copyfile(filename_old,filename_new)
        filename_old = gui_specs['save_dir'] + gui_specs['subj_name'] + '_explog_cd.csv'
        filename_new = gui_specs['save_dir'] + gui_specs['subj_name'] + '_explog_cd_' + today.strftime('%Y%m%d_%H%M%S') + '.csv'
        copyfile(filename_old,filename_new)
    
    #
    #end of experiment screen 
    exp.endOfExperiment(f=f,tk=tk)