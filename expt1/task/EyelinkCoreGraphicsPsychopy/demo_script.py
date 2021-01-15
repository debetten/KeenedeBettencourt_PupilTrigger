######## Psychopy demo with the custom PsychopyCoreGraphics
# If you need to use a screen units other than 'pix', which we do not recommend as the gaze coordinates 
# returned by pylink is in 'pix' units, please make sure to properly set the size of the cursor and the
# position of the gaze. However, the calibration/validation routine should work fine regardless of the 
# screen units.

# import libraries
import pylink, numpy, os, random
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
from psychopy import visual, core, event, monitors

#### Set a few task parameters
useGUI = True # whether use the Psychopy GUI module to collect subject information
dummyMode = False # If in Dummy Mode, press ESCAPE to skip calibration/validataion

#### STEP I: get subject info with GUI ########################################################
expInfo = {'SubjectNO':'00', 'SubjectInitials':'TEST'}

if useGUI:
    from psychopy import gui
    dlg = gui.DlgFromDict(dictionary=expInfo, title="GC Example", order=['SubjectNO', 'SubjectInitials'])
    if dlg.OK == False: core.quit()  # user pressed cancel
else:
    expInfo['SubjectNo'] = raw_input('Subject # (1-99): ')
    expInfo['SubjectInitials'] = raw_input('Subject Initials (e.g., WZ): ')

#### SETP II: established a link to the tracker ###############################################
if not dummyMode: tk = pylink.EyeLink('100.1.1.1')
else:             tk = pylink.EyeLink(None)

#### STEP III: Open an EDF data file EARLY ####################################################
dataFolder = os.getcwd() + '/edfData/'
if not os.path.exists(dataFolder): os.makedirs(dataFolder)
dataFileName = expInfo['SubjectNO'] + '_' + expInfo['SubjectInitials'] + '.EDF'

# Note that for Eyelink 1000/II, the file name cannot exceeds 8 characters
# we need to open eyelink data files early so as to record as much info as possible
tk.openDataFile(dataFileName)

# add personalized header (preamble text)
tk.sendCommand("add_file_preamble_text 'Psychopy GC demo'") 

#### STEP IV: Initialize custom graphics for camera setup & drift correction ##################
scnWidth, scnHeight = (1920, 1080)
# you MUST specify the physical properties of your monitor first, otherwise you won't be able to properly use
# different screen "units" in psychopy 
mon = monitors.Monitor('myMac15', width=53.0, distance=70.0)
mon.setSizePix((scnWidth, scnHeight))
win = visual.Window((scnWidth, scnHeight), fullscr=True, monitor=mon, color=[0,0,0], units='pix')

# this functional calls our custom calibration routin "EyeLinkCoreGraphicsPsychopy.py"
genv = EyeLinkCoreGraphicsPsychoPy(tk, win)
pylink.openGraphicsEx(genv)

#### STEP V: Set up the tracker ################################################################
# we need to put the tracker in offline mode before we change its configrations
tk.setOfflineMode()
# sampling rate, 250, 500, 1000, or 2000
tk.sendCommand('sample_rate 500')

# Online parser configuration: 0-> standard/coginitve, 1-> sensitive/psychophysiological
# [see Eyelink User Manual, Section 4.3: EyeLink Parser Configuration]
tk.sendCommand('select_parser_configuration 0')
# Set the tracker to record Event Data in "GAZE" (or "HREF") coordinates
tk.sendCommand("recording_parse_type = GAZE")

# inform the tracker the resolution of the subject display
# [see Eyelink Installation Guide, Section 8.4: Customizing Your PHYSICAL.INI Settings ]
tk.sendCommand("screen_pixel_coords = 0 0 %d %d" % (scnWidth-1, scnHeight-1))

# stamp display resolution in EDF data file for Data Viewer integration
# [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
tk.sendMessage("DISPLAY_COORDS = 0 0 %d %d" % (scnWidth-1, scnHeight-1))

# specify the calibration type, H3, HV3, HV5, HV13 (HV = horiztonal/vertical), 
tk.sendCommand("calibration_type = HV13") # tk.setCalibrationType('HV9') also works, see the Pylink manual
# specify the proportion of subject display to calibrate/validate
tk.sendCommand("calibration_area_proportion 0.85 0.85")
tk.sendCommand("validation_area_proportion  0.85 0.85")

# allow buttons on the gamepad to accept calibration/dirft check target, so you
# do not need to press keys on the keyboard to initiate/accept calibration
tk.sendCommand("button_function 1 'accept_target_fixation'")

# data stored in data file and passed over the link (online)
# [see Eyelink User Manual, Section 4.6: Settting File Contents]
eyelinkVer = tk.getTrackerVersion()
if eyelinkVer >=3: # Eyelink 1000/1000 plus
    tk.sendCommand("file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT")
    tk.sendCommand("link_event_filter = LEFT,RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,BUTTON,INPUT")
    tk.sendCommand("file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HTARGET")
    tk.sendCommand("link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HTARGET")
else: # Eyelink II
    tk.sendCommand("file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT")
    tk.sendCommand("link_event_filter = LEFT,RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,BUTTON,INPUT")
    tk.sendCommand("file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT")
    tk.sendCommand("link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT")

#### STEP VI: specify all possible experimental cells #################################################
# one may read in s spreadsheet that defines the experimentl cells; usually, a simple list-like the one below
# should also do the job; if we need tweenty trials, simple go with "new_list = trials[:]*10", then 
# random.shuffle(new_list) 
trials = [['CondA', 'sacrmeto.jpg'],
          ['CondB', 'sacrmeto.jpg']]

#### SETP VII: a helper to run a single trial #########################################################
def runTrial(pars):
    """ pars corresponds to a row in the trial list"""
    
    # retrieve paramters from the trial list
    cond, pic = pars 
    
    # load the image to display
    img = visual.ImageStim(win, image=pic, size=(scnWidth, scnHeight)) # stretch the image to fill full screen
    gazeCursor = visual.GratingStim(win, tex='none', mask='circle', size=100.0, color=[1.0,1.0,1.0])

    # take the tracker offline
    tk.setOfflineMode()
    pylink.pumpDelay(50)

    # send the standard "TRIALID" message to mark the start of a trial
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tk.sendMessage('TRIALID')
    
    # record_status_message : show some info on the host PC
    tk.sendCommand("record_status_message 'Cond %s'"% cond)
    
    # drift check
    try:
        err = tk.doDriftCorrect(win.size[0]/2, win.size[1]/2,1,1)
    except:
        tk.doTrackerSetup()
        
    # uncomment this line to read out calibration/drift-correction results
    #print tk.getCalibrationMessage() 

    # start recording
    tk.setOfflineMode()
    pylink.pumpDelay(50)
    error = tk.startRecording(1,1,1,1)
    pylink.pumpDelay(100) # wait for 100 ms to make sure data of interest is recorded
    
    # show the image 
    img.draw()  
    win.flip()
    # this message marks the time 0 of a trial
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tk.sendMessage('DISPLAY_SCREEN') 
    
    #determine which eye(s) are available
    eyeTracked = tk.eyeAvailable() 
    
    # show the image indefinitely until a key is pressed
    gazePos =  (scnWidth/2, scnHeight/2)
    terminate = False
    event.clearEvents() # clear cached (keyboard/mouse etc.) events, if there is any
    while not terminate:
        # check for keypress to terminate a trial
        if len(event.getKeys())>0: # KEYBOARD
            terminate = True
        if True in tk.getLastButtonPress(): # GamePad connected to the tracker HOST PC
            terminate = True
           
        # check for new samples
        dt = tk.getNewestSample()
        if (dt != None):
            if eyeTracked == 1 and dt.isRightSample():
                gPos = dt.getRightEye().getGaze()
            elif eyeTracked == 0 and dt.isLeftSample():
                gPos = dt.getLeftEye().getGaze()
            gazePos = (gPos[0]-scnWidth/2, scnHeight/2-gPos[1])

        gazeCursor.pos = gazePos
        img.draw()
        gazeCursor.draw()
        win.flip()
        
    # clear the subject display
    win.color=[0,0,0]
    win.flip()
    
    # clear the host display, this command is needed if you are backdropping images
    # to the host display (not demonstrated in this script)
    tk.sendCommand('clear_screen 0') 

    # send trial variables for Data Viewer integration
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tk.sendMessage('!V TRIAL_VAR cond %s' %cond)

    # send a message to mark the end of trial
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tk.sendMessage('TRIAL_RESULTS 0')
    pylink.pumpDelay(100)
    tk.stopRecording() # stop recording


#### STEP VIII: The real experiment starts here ##########################################

# show some instructions here.
msg = visual.TextStim(win, text = 'Calibration instructions\n\n\n\tENTER--Show Camera Image\n\n\tC--Calibration\n\n\tV--Validation\n\n\tO--Start Recording', 
                        color = 'black', units = 'pix')
msg.draw(); win.flip()
event.waitKeys()

# set up the camera and calibrate the tracker at the beginning of each block
tk.doTrackerSetup()

# run a block of trials
testList = trials[:]*1 # construct the trial list
random.shuffle(testList) # randomize the trial list
# Looping through the trial list
for t in testList: 
    runTrial(t)

# close the EDF data file
tk.setOfflineMode()
tk.closeDataFile()
pylink.pumpDelay(50)

# Get the EDF data and say goodbye
msg.text='Data transfering.....'
msg.draw(); win.flip()
tk.receiveDataFile(dataFileName, dataFolder + dataFileName)

#close the link to the tracker
tk.close()

# close the graphics
pylink.closeGraphics()
win.close()
core.quit()
