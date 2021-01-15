from __future__ import division
import time, os, platform, math, pickle, sys, types
import numpy as np
from datetime import datetime
if 'psychopy' in sys.modules:
    from psychopy import visual, core, tools, monitors, gui, event
    import psychopy.event
if sys.platform == 'win32':
    from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
    import pylink
    
def get_specs(subj_name='mmddyyn_wmLoadMem02'):
    """Opens a GUI to set up the experiment
    """
    
    dictDlg = gui.DlgFromDict({'Participant number': '08dd18_eyeSustAttnWM01',     #subject name: month, day, year, number, and project name
                                'Environment':['booth','eeg','Dirk VU','imac','edslab'],    #which environment (aka monitor name in Tools > Monitor Center)
                                'Restarting experiment?':False,     #whether restarting an experiment that was aborted in the middle
                                'Debug':False,                    #whether instructions are in english
                                'Eye tracker':True,
                                'Picture game':True,
                                'Block game': True},
                                title='Welcome to the experiment',  #title for GUI
                                fixed=[' '],
                                order=['Participant number','Environment','Restarting experiment?','Debug','Eye tracker']) #choose order of elements
    if dictDlg.OK == True:
        gui_specs = {}
        gui_specs['subj_name']=str(dictDlg.data[0])
        gui_specs['env'] = dictDlg.data[1]
        gui_specs['restarting'] = dictDlg.data[2]
        gui_specs['expt_date'] = datetime.now().strftime("%m/%d/%Y %H:%M")
        gui_specs['debug'] = dictDlg.data[3]
        gui_specs['eyetrack'] = dictDlg.data[4]
        gui_specs['seed'] = int(time.time())
        gui_specs['sustattn'] = dictDlg.data[5]
        gui_specs['changedetect'] = dictDlg.data[6]
    else:
        core.quit()
        
    #check the specs
    assert isinstance(gui_specs['subj_name'],basestring), "subj_name is not an string: %r" % specs['subj_name']
    assert isinstance(gui_specs['env'],basestring), "env is not a string: %r" % specs['env']
    assert isinstance(gui_specs['expt_date'],basestring), "expt_date is not a string: %r" % specs['expt_date']
    assert isinstance(gui_specs['seed'],int), "seed is not a int: %r" % specs['seed']

    #which environment is being used to run the experiment
    if gui_specs['env']=='imac':
        gui_specs['project_dir']='/Users/megan/Documents/projects/eyeSustAttnWM01/'
    elif gui_specs['env']=='booth':
        gui_specs['project_dir']='C:/Users/AwhVogelLab/Desktop/Megan/eyeSustAttnWM01/'
    elif gui_specs['env']=='edslab':
        gui_specs['project_dir']='/Users/edslab/Documents/projects/eyeSustAttnWM01/'
    elif gui_specs['env']=='Dirk VU':
        gui_specs['project_dir']='C:/Users/Dirk VU/Documents/Github/eyeSustAttnWM01/'
    elif gui_specs['env']=='eeg':
        if os.path.isdir('C:/Users/awhvogellab/Desktop/Megan/eegSustAttnWM01/'):
            gui_specs['project_dir']='C:/Users/awhvogellab/Desktop/Megan/eyeSustAttnWM01/' #FIX THIS LATER
    else:
        pass
        #print 'ERROR: unknown environment selected in the GUI'
        
    #which directory to use to save the data
    gui_specs['save_dir'] = gui_specs['project_dir'] + 'subjects/' + gui_specs['subj_name'] + '/data/beh/'
    
    #if the directory where data will be saved does not exist yet
    if not os.path.isdir(gui_specs['save_dir']):
        #print "saving files to: ", gui_specs['save_dir']
        os.makedirs(gui_specs['save_dir']) #this command can make multiple subfolders
        
    return gui_specs

def angle_diff(t1, t2):
    t1r = np.radians(t1)
    t2r = np.radians(t2)
    
    x = np.cos(t1r-t2r)
    y = np.sin(t1r-t2r)
    dr = np.arctan2(y,x)
    return np.degrees(dr) 


def wait_for_response(win,keylist=None,f=None):
    event.clearEvents(eventType='keyboard')
    if keylist == None:
        keys = psychopy.event.waitKeys()
    else:
        keys = psychopy.event.waitKeys(keyList=keylist)
    if keys[0] == 'escape':
        if f is not None: f.close()
        win.close()
        core.quit()
    response = keys[0]
                    
    return response

def poll_for_response(win,keylist=None,f=None):
    if keylist is not None:
        keys = psychopy.event.getKeys(keyList=keylist)
    else:
        keys = psychopy.event.getKeys()
    if np.size(keys)>0:# is not None:
        if keys[0] == 'escape':
            if f is not None: f.close()
            win.close()
            core.quit()
        response = keys[0]
    else:
        response = None
                    
    return response
   
def drawpath(expt_display,n_circles,i_fill=None):
    
    # draw the circle that indicates the start of experiment
    expt_display.start_circle.draw()
    expt_display.start_text.draw()
    
    # draw the star that indicates the end of experiment
    expt_display.finish_star.draw()
    expt_display.finish_text.draw()
    
    # center of all the circles
    circle_centerX = np.linspace(-3.5,3.5,n_circles)
    
    # draw each circle
    for i in range(n_circles):
        x = circle_centerX[i]
        y = -6  #np.sin(2*np.pi*(i/nCircles))-6
        expt_display.block_outter.fillColor = (-1,-1+2*(i/n_circles),0)
        expt_display.block_outter.pos = (x,y)
        expt_display.block_outter.draw()
        expt_display.block_inner.pos = (x,y)
        if i_fill is not None:
            if i<=i_fill:
                expt_display.block_inner.fillColor= (-1,-1+2*(i/n_circles),0)
            else:
                expt_display.block_inner.fillColor= expt_display.scrcolor
        else:
            expt_display.block_inner.fillColor= expt_display.scrcolor
        expt_display.block_inner.draw()
    

def open_csv_data_file(gui_specs,data_filename,overwrite_ok=None):
    """Opens the csv file and writes the header.
    Parameters:
        - gui_specs: all the selections from the gui at the start of the experiment
        - filename: the name of the CSV file (without the .csv)
    """
    
    #if the filename has .csv delete it
    if data_filename[-4:] == '.csv':
        data_filename = data_filename[:-4]
    
    #add the path to the file name
    data_filename = gui_specs['save_dir'] + data_filename + '.csv'
    
    #open the csv file
    data_file = open(data_filename, 'a')
        
    #write all the header information
    for key, value in gui_specs.iteritems():
        data_file.write('"')
        data_file.write(key)
        data_file.write(',')
        data_file.write(str(value))
        data_file.write('"')
        data_file.write('\n')
    data_file.write('\n')
    
    return data_file

def write_to_csv_data_file(data_file,iblock,itrial,dsgn,headerLine=False):
    if headerLine:
        data_file.write('block,trial,time,setsize,freq_trials,'+
                        'colorind,correctresp,actualresp,acc,rts,'+
                        'probe_trials,rt_triggered_fast,rt_triggered_slow,rt_triggered_med,'+
                        'rt_trailingavg,rt_runningavg,rt_runningstd,rt_runningslowthresh,rt_runningfastthresh,rt_runningmedslowthresh,rt_runningmedfastthresh,'+
                        'wholereportresp,wholereportresporder,wholereportrespcolorind,'+
                        'wholereportrespacc,wholereportRT1,wholereportRT2'+
                        'timeBlockStart,timeEncArrayOnset,timeResponse,' + 
                        'timeWholeReportOnset,timeWholeReportOffset,\n')
    else:
        data_file.write(str(iblock) +','+ str(itrial) +','+ datetime.utcnow().strftime("%H:%M:%S:%f") +','+ 
                        str(dsgn.setsize) +','+ str(dsgn.freq_trials[iblock,itrial]) +','+ 
                        str(dsgn.trials_colorind[iblock,itrial]) +','+
                        str(dsgn.correct_response[iblock,itrial]) +','+ 
                        str(dsgn.actual_response[iblock,itrial]) +','+ str(dsgn.acc[iblock,itrial]) +','+ 
                        str(dsgn.rts[iblock,itrial]) +','+ str(dsgn.probe_trials[iblock,itrial]) +','+ 
                        str(dsgn.rt_triggered_fast[iblock,itrial]) +','+ str(dsgn.rt_triggered_slow[iblock,itrial]) +','+ 
                        str(dsgn.rt_triggered_med[iblock,itrial]) +','+
                        str(dsgn.rts_trailingavg[iblock,itrial]) +','+ str(dsgn.rts_runningavg[iblock,itrial]) +','+
                        str(dsgn.rts_runningstd[iblock,itrial]) +','+ str(dsgn.rts_runningslowthresh[iblock,itrial]) +','+ 
                        str(dsgn.rts_runningfastthresh[iblock,itrial]) +','+
                        str(dsgn.rts_runningmedslowthresh[iblock,itrial]) +','+ str(dsgn.rts_runningmedfastthresh[iblock,itrial]) +','+ 
                        str(dsgn.wholereport_resp[iblock,itrial]) +','+ str(dsgn.wholereport_resporder[iblock,itrial]) +','+ 
                        str(dsgn.wholereport_respcolorind[iblock,itrial]) +','+ 
                        str(np.nansum(dsgn.wholereport_respacc[iblock,itrial])) +','+ str(dsgn.wholereport_rts[iblock,itrial,0:3]) +','+ str(dsgn.wholereport_rts[iblock,itrial,3:6]) +','+ #the csv has a weird glitch if we dont break the wr_rts up.
                        str(dsgn.time_blockstart[iblock]) +','+ str(dsgn.time_encarray_onset[iblock,itrial]) +','+ str(dsgn.time_response[iblock,itrial]) +','+ 
                        str(dsgn.time_memprobe_onset[iblock,itrial]) +','+ str(dsgn.time_memprobe_offset[iblock,itrial]))
        data_file.flush()
        data_file.write('\n')


def write_to_csv_data_file_cd(data_file,iblock,itrial,dsgn,headerLine=False):
    if headerLine:
        data_file.write('blockNum,trialNum,time,setSize,sameProbe,encArrayColorInd,encArrayQuad,encArrayMinDist,' +
                        'memProbeInd,memProbeOrigColorInd,memProbeOrigColorRGB,memProbeProbeColorInd,memPRobeProbeColorRGB,memProbeQuad,memProbeX,memProbeY' +
                        'corrResp,actualResp,actualRespString,acc,RT,' +
                        'timeBlockStart,timeEncArrayOnset,timeEncArrayOffset,timeMemProbeOnset,timeMemProbeOffset\n')
    else:
        data_file.write(str(iblock) + ','+ str(itrial) + ',' + datetime.utcnow().strftime("%H:%M:%S:%f") + ',' + str(dsgn.cd_setsize) +',' + str(dsgn.cd_same_probe[itrial]) + ',' + str(dsgn.cd_encarray_colorind[itrial]) + ',' + str(dsgn.cd_encarray_quad[itrial]) + ',' +str(dsgn.cd_encarray_mindist[itrial]) + ',' +str(dsgn.cd_memprobe_ind[itrial])+','+str(dsgn.cd_memprobe_origcolorind[itrial]) + ',' +str(dsgn.cd_memprobe_probecolorind[itrial]) + ',' + str(dsgn.cd_memprobe_probecolorrgb[itrial]) + ',' + str(dsgn.cd_memprobe_quad[itrial]) + ',' + str(dsgn.cd_memprobe_x[itrial]) + ',' + str(dsgn.cd_memprobe_y[itrial]) + ',' +str(dsgn.cd_correct_response[itrial]) + ',' + str(dsgn.cd_actual_response[itrial]) + ',' + dsgn.cd_actual_response_string[itrial] + ',' +str(dsgn.cd_acc[itrial]) + ',' + str(dsgn.cd_rt[itrial]) + ',' + str(dsgn.cd_time_blockstart) + ',' + str(dsgn.cd_time_encarray_onset[itrial]) + ',' +str(dsgn.cd_time_encarray_offset[itrial]) + ',' + str(dsgn.cd_time_memprobe_onset[itrial]) + ',' + str(dsgn.cd_time_memprobe_offset[itrial]))
        data_file.flush()
        data_file.write('\n')

def fullfact(levels):
    """
        Create a general full-factorial design
        
        Parameters
        ----------
        levels : array-like
        An array of integers that indicate the number of levels of each input
        design factor.
        
        Returns
        -------
        mat : 2d-array
        The design matrix with coded levels 0 to k-1 for a k-level factor
        
        Example
        -------
        ::
        
        >>> fullfact([2, 4, 3])
        array([[ 0.,  0.,  0.],
        [ 1.,  0.,  0.],
        [ 0.,  1.,  0.],
        [ 1.,  1.,  0.],
        [ 0.,  2.,  0.],
        [ 1.,  2.,  0.],
        [ 0.,  3.,  0.],
        [ 1.,  3.,  0.],
        [ 0.,  0.,  1.],
        [ 1.,  0.,  1.],
        [ 0.,  1.,  1.],
        [ 1.,  1.,  1.],
        [ 0.,  2.,  1.],
        [ 1.,  2.,  1.],
        [ 0.,  3.,  1.],
        [ 1.,  3.,  1.],
        [ 0.,  0.,  2.],
        [ 1.,  0.,  2.],
        [ 0.,  1.,  2.],
        [ 1.,  1.,  2.],
        [ 0.,  2.,  2.],
        [ 1.,  2.,  2.],
        [ 0.,  3.,  2.],
        [ 1.,  3.,  2.]])
        
        """
    n = len(levels)  # number of factors
    nb_lines = np.prod(levels)  # number of trial conditions
    H = np.zeros((nb_lines, n))
    
    level_repeat = 1
    range_repeat = np.prod(levels)
    for i in range(n):
        range_repeat //= levels[i]
        lvl = []
        for j in range(levels[i]):
            lvl += [j]*level_repeat
        rng = lvl*range_repeat
        level_repeat *= levels[i]
        H[:, i] = rng
    
    return H

def eyetracker_settings(disp,dsgn):
    tk = pylink.EyeLink('100.1.1.1')
    genv = EyeLinkCoreGraphicsPsychoPy(tk, disp.win)
    pylink.openGraphicsEx(genv)
    #print "dsgn.eyetracking_filename"
    tk.openDataFile(dsgn.eyetracking_filename)
    pylink.flushGetkeyQueue()
    #pylink.setTargetSize(int(disp.scrw/70), int(disp.scrw/300))
    tk.setOfflineMode()
    tk.sendCommand("add_file_preamble_text = %s" %('eegSustAttnWM01'))
    tk.sendCommand('elcl_select_configuration = BTABLER')
    tk.sendCommand("screen_pixel_coords = 0 0 %d %d" %(disp.scrw, disp.scrh))
    tk.sendCommand('calibration_type = HV5')    
    tk.sendCommand('sample_rate 1000') # can set 250, 500, or 1000 depending on currently selected mode
    #getEYELINK().sendCommand('select_parser_configuration 0') # 0-> standard, 1-> sensitive
    tk.sendCommand('elcl_tt_power = 1')
    tk.sendCommand('use_ellipse_fitter = yes')
    #tk.sendCommand('add_file_preamble_text %s' %(dsgn.save_dir))
    tk.sendCommand("screen_pixel_coords =  0 0 %d %d" %(disp.scrw, disp.scrh))
    tk.sendCommand('file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT')#EDF file contents
    tk.sendCommand('file_sample_data = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT')
    tk.sendCommand('link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT') #Set link data
    tk.sendCommand('link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT')
    tk.sendCommand('calibration_area_proportion 0.5 0.5') #Adjust size of calibrated area
    tk.sendCommand('validation_area_proportion 0.5 0.5')
    tk.sendCommand('saccade_velocity_threshold = 30');
    tk.sendCommand('saccade_acceleration_threshold = 8000');
    tk.sendMessage("DISPLAY_COORDS 0 0 %d %d" %(disp.scrw, disp.scrh))
    return tk



def eyetracker_gazecontrol(disp,tk,fix_tol=10,fb=False):

    #figure out which eye is being tracked
    if tk is not None:
        eye = tk.eyeAvailable()
    #    if eye ==2: #if both eye's data present: use left eye only
    #        eye = 0

        #check for new sample update
        dt = tk.getNewestSample() 

        # Gets the gaze position of the latest sample
        if(dt != None):
            if eye == 0:
                gaze_pos= dt.getLeftEye().getGaze() 
                gaze_pos_2 = 0
            elif eye == 1:
                gaze_pos = dt.getRightEye().getGaze()
                gaze_pos_2 = 0
            elif eye == 2:
                gaze_pos = dt.getLeftEye().getGaze() 
                gaze_pos_2= dt.getRightEye().getGaze()
            #print dt, eye, gaze_pos
            
            if(gaze_pos != None):
                #calculate deviation from center point
                gaze_dist = np.sqrt((gaze_pos[0]-disp.scrw/2)**2+(gaze_pos[1]-disp.scrh/2)**2)
                if (eye == 2) and (gaze_pos_2 != None):
                    gaze_dist_2 = np.sqrt((gaze_pos_2[0]-disp.scrw/2)**2+(gaze_pos_2[1]-disp.scrh/2)**2)
                else:
                    gaze_dist_2 =100
            else:
                gaze_dist = [disp.scrw/2,disp.scrh/2]
                gaze_dist_2 = 100
            
            gaze_dist_deg = tools.monitorunittools.pix2deg(gaze_dist,disp.mon)
            gaze_dist_deg_2 = tools.monitorunittools.pix2deg(gaze_dist_2,disp.mon)
        else:
            gaze_dist_deg = 0
            gaze_dist_deg_2 = 0
        #print "gaze_pos, gaze_dist, gaze_dist_deg", gaze_pos, gaze_dist, gaze_dist_deg

        #compare to tolerance
        if (gaze_dist_deg < fix_tol ) or (gaze_dist_deg_2 < fix_tol):
            within_tol = 1
        else: 
            within_tol = 0
            #print "gaze_pos, gaze_dist, gaze_dist_deg", gaze_pos, gaze_dist, gaze_dist_deg
            
            #should give x, y coordinates in degrees
            #gaze_pos_deg = tools.monitorunittools.pix2deg((float(gaze_pos[0]-disp.scrw/2),float(gaze_pos[1]-disp.scrh/2)),disp.mon)
            
            if fb == True:
                #draw line
                if gaze_dist_deg<50:
                    line = visual.Line(disp.win,start=(0,0),end = (gaze_pos[0]-disp.scrw/2,-1*(gaze_pos[1]-disp.scrh/2)),units='pix')
                    line.draw()
                else:
                    text = visual.TextStim(disp.win, text='Blink',alignHoriz='center',wrapWidth=20,height=1,pos=(0,7),color=(-1,-1,-.5),fontFiles=[disp.project_dir + 'display/fonts/BowlbyOneSC-Regular.ttf'],font=['Bowlby One SC'])
                    text.draw()
                disp.fix.draw()
                disp.win.flip()
                core.wait(1)
    else:
        #print "ERROR IN GAZE CONTROL! ALERT! WHY IS TK NONE???"
        within_tol = 1
    
    return within_tol