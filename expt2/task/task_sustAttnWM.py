import sys
if 'psychopy' in sys.modules:
    from psychopy import visual, core, tools, event, monitors, logging, parallel
from helper_functions_sustAttnWM import *
import numpy as np
import scipy
from pylink import *
from scipy.spatial.distance import pdist
if sys.platform == 'win32':
    from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
    import pylink

#all the data files 
class Task(object):

    def __init__(self,expt_design,expt_display):

        self.dsgn = expt_design
        self.disp = expt_display             
        
    def welcomeToExperiment(self,f=None,tk=None):

        self.disp.text_welcome.draw()           #welcome text at top of screen
        #self.disp.text_advance.draw()           #advance text at bottom of screen
        self.disp.text_record_advance.draw()         #advance text at bottom of screen
        self.disp.fix.draw()
        drawpath(self.disp,self.dsgn.nblocks+1)   #draw the path of the experiment
        self.disp.text_record_screen.draw()             #Check if saving data
        self.disp.win.flip()                    #flip the window
        
        #if not self.dsgn.debug:
        wait_for_response(self.disp.win,keylist=self.dsgn.record_key+self.dsgn.key_escape,f=f)
        self.disp.fix.draw()
        self.disp.win.flip()
        core.wait(.5)

    def instructionsSustAttn(self,clock,f=None,tk=None):

        self.disp.text_instructions = visual.TextStim(self.disp.win, 
            text="During this task, a new array of circles like this one will appear on the screen every second\n\n" + 
                "You should press the d key using the left pointer finger whenever new circles appear",
                alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
        self.disp.text_instructions.draw()     #instructions text at top of screen
        self.disp.fix.draw()
        self.disp.text_advance.pos = (0,-6)
        self.disp.text_advance.draw()           #advance text at bottom of screen

        npracticetrials = 1
        practicetrials_colorind = np.zeros((npracticetrials,self.dsgn.setsize))
        practicetrials_colorrgb = np.zeros((npracticetrials,self.dsgn.setsize,3))
        #color assigned to each position
        for itrial in range(npracticetrials):
            practicetrials_colorind[itrial] = np.random.choice(self.dsgn.ncolors,
                size=self.dsgn.setsize, replace=False)

        #set up the rgb vales for each trial
        for itrial in range(npracticetrials):
            
            for ii in range(self.dsgn.setsize):
                print practicetrials_colorind[itrial][ii]
                print self.dsgn.colors_rgb
                print self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]
                practicetrials_colorrgb[itrial][ii] = self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]

        #draw the circles
        for itrial in range(npracticetrials):
            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw()
            for i in range(self.dsgn.setsize):
                self.disp.circle.fillColor = practicetrials_colorrgb[itrial,i]
                self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                self.disp.circle.draw()           
        self.disp.win.flip()
        wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)

        self.disp.fix.draw()
        self.disp.win.flip()
        core.wait(.5)

        self.disp.text_instructions = visual.TextStim(self.disp.win, 
                text="Sometimes, squares will appear instead of circles. "+
                "Press the s key when squares appear using your left middle finger.",               
                alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
        self.disp.text_instructions.draw()     #instructions text at top of screen
        self.disp.fix.draw()
        self.disp.text_advance.pos = (0,-6)
        self.disp.text_advance.draw()           #advance text at bottom of screen
        npracticetrials = 1
        practicetrials_colorind = np.zeros((npracticetrials,self.dsgn.setsize))
        practicetrials_colorrgb = np.zeros((npracticetrials,self.dsgn.setsize,3))
        #color assigned to each position
        for itrial in range(npracticetrials):
            practicetrials_colorind[itrial] = np.random.choice(self.dsgn.ncolors,
                size=self.dsgn.setsize, replace=False)

        #set up the rgb vales for each trial
        for itrial in range(npracticetrials):
            
            for ii in range(self.dsgn.setsize):
                print practicetrials_colorind[itrial][ii]
                print self.dsgn.colors_rgb
                print self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]
                practicetrials_colorrgb[itrial][ii] = self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]

        #draw the squares
        for itrial in range(npracticetrials):
            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw()
            for i in range(self.dsgn.setsize):
                self.disp.square.fillColor = practicetrials_colorrgb[itrial,i]
                self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                self.disp.square.draw()         
        self.disp.win.flip()
        wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)   

        self.disp.fix.draw()
        self.disp.win.flip()
        core.wait(.5)
        
        self.disp.text_instructions = visual.TextStim(self.disp.win, 
            text="Press s for squares and d for circles. Let's practice this now.\n\nIt is important that you keep your eyes on the black dot in the middle the whole time",
            alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
            fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
        self.disp.text_instructions.draw()     #instructions text at top of screen
        self.disp.fix.draw()
        self.disp.text_advance.draw()           #advance text at bottom of screen
        self.disp.win.flip() 
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape+self.dsgn.key_skip,f=f)       

        self.disp.fix.draw()
        self.disp.win.flip()
        
        #start recording
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            
        core.wait(.5)

        npracticetrials = 10

        acc = np.zeros(npracticetrials)
        while True:

            if response==self.dsgn.key_skip[0]:
                break
            
            practicetrials_colorind = np.zeros((npracticetrials,self.dsgn.setsize))
            practicetrials_colorrgb = np.zeros((npracticetrials,self.dsgn.setsize,3))
            actual_response = np.zeros((npracticetrials))
            correct_response = [1,1,1,0,1,1,1,0,1,1]

            for itrial in range(npracticetrials):

                #color assigned to each position
                practicetrials_colorind[itrial] = np.random.choice(self.dsgn.ncolors,
                    size=self.dsgn.setsize, replace=False)

            #ensure that the same color doesn't appear in the same position back to back
            for ii in range(self.dsgn.setsize):
                diff_zero = np.diff(self.dsgn.trials_colorind[:,ii])
                for it in np.where(diff_zero==0)[0]:
                    if it > 0:
                        used_colors = np.append(np.append(practicetrials_colorind[it,:],practicetrials_colorind[it-1,ii]),practicetrials_colorind[it+1,ii])
                        avail_colors = np.setdiff1d(self.dsgn.colors_ind,used_colors)
                        practicetrials_colorind[it,ii] = np.random.choice(avail_colors,size=1)

            for itrial in range(npracticetrials):
                #set up the rgb vales for each trial
                for ii in range(self.dsgn.setsize):
                    print practicetrials_colorind[itrial][ii]
                    print self.dsgn.colors_rgb
                    print self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]
                    practicetrials_colorrgb[itrial][ii] = self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]

            for itrial in range(npracticetrials):
                self.disp.fix.fillColor = (-1,-1,-1)
                self.disp.fix.draw()
                if itrial==3 or itrial==7:
                    for i in range(self.dsgn.setsize):
                        self.disp.square.fillColor = practicetrials_colorrgb[itrial,i]
                        self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                        self.disp.square.draw() 
                else:
                    for i in range(self.dsgn.setsize):
                        self.disp.circle.fillColor = practicetrials_colorrgb[itrial,i]
                        self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                        self.disp.circle.draw()           
                self.disp.win.flip()

                #save timing information
                practicetrial_encarray_onset = clock.getTime()
                
                response = None
                event.clearEvents(eventType='keyboard')
                while clock.getTime()<(practicetrial_encarray_onset+self.dsgn.t_stim):
                    if response is None:
                        response = poll_for_response(self.disp.win,self.dsgn.key_freq+self.dsgn.key_infreq+self.dsgn.key_escape)
                        if response is not None:
                            
                            if response==self.dsgn.key_freq[0]:
                                actual_response[itrial]=1
                            elif response==self.dsgn.key_infreq[0]:
                                actual_response[itrial]=0

                            #see if the response was correct or not
                            acc[itrial]=(correct_response[itrial]==actual_response[itrial])*1

                            #redraw the display
                            self.disp.fix.fillColor = (1,1,1)
                            self.disp.fix.draw()
                            if itrial==3 or itrial==7:
                                for i in range(self.dsgn.setsize):
                                    self.disp.square.fillColor = practicetrials_colorrgb[itrial,i]
                                    self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                                    self.disp.square.draw() 
                            else:
                                for i in range(self.dsgn.setsize):
                                    self.disp.circle.fillColor = practicetrials_colorrgb[itrial,i]
                                    self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                                    self.disp.circle.draw()           
                            self.disp.win.flip()


            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw()
            self.disp.win.flip() 
            
            if tk is not None:
                pylink.pumpDelay(100)
                tk.stopRecording()
                pylink.pumpDelay(100)
                tk.setOfflineMode()
                
            core.wait(.5)

            if np.mean(acc)>=.9: 
                self.disp.text_instructions = visual.TextStim(self.disp.win, 
                    text="You got " + str(int(np.round(np.mean(acc),decimals=2)*100)) + " percent correct. "+
                    "Great job!",
                    alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                    fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
                self.disp.text_instructions.draw()     #instructions text at top of screen
                self.disp.fix.draw()
                self.disp.text_advance.draw()           #advance text at bottom of screen
                self.disp.win.flip()             
            else:
                self.disp.text_instructions = visual.TextStim(self.disp.win, 
                    text="You got " + str(int(np.round(np.mean(acc),decimals=2)*100)) + " percent correct. "+
                    "Let's try that again. Remember to press d for circles, s for squares",
                    alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                    fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
                self.disp.text_instructions.draw()     #instructions text at top of screen
                self.disp.fix.draw()
                self.disp.text_advance.draw()           #advance text at bottom of screen
                self.disp.win.flip()       
            
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape+self.dsgn.key_skip,f=f)     

            if np.mean(acc)>=.9:
                break
            else:
                acc = np.zeros(npracticetrials)
            self.disp.fix.draw()
            self.disp.win.flip()                 
            core.wait(.5)

    def instructionsWholeReport(self,clock,f=None):

        self.disp.text_instructions = visual.TextStim(self.disp.win, 
            text="Sometimes, the whole screen will go blank and grey. " +
            "Then, small objects that look like a face of a Rubik's cube will appear at each spot.\n\n" +
            "You will use your mouse to click the color that had been at each location.",
            alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
            fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
        self.disp.text_instructions.draw()     #instructions text at top of screen
        self.disp.fix.draw()
        self.disp.text_advance.draw()           #advance text at bottom of screen 
        for i in range(self.dsgn.setsize):
            for j in range(9):
                self.disp.wholereportsquare[i,j].draw()
        self.disp.win.flip() 
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)       

        self.disp.fix.draw()
        self.disp.win.flip() 
        core.wait(.5)

        while True:

            self.disp.text_instructions = visual.TextStim(self.disp.win, 
                text="Let's practice just this Rubik's cube task now. It will probably feel hard, but that's okay!\n\n"+
                "Please try to respond as quickly and as accurately as possible. Don't agonize too much over your decision, but try to go with gut feeling and just try your best!",
                alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
            self.disp.text_instructions.draw()      #instructions text at top of screen
            self.disp.fix.draw()
            self.disp.text_advance.draw()           #advance text at bottom of screen
            self.disp.win.flip() 
            wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)  
                  
            self.disp.fix.draw()
            self.disp.win.flip() 
            core.wait(.5)

            npracticetrials = 1
            practicetrials_colorind = np.zeros((npracticetrials,self.dsgn.setsize))
            practicetrials_colorrgb = np.zeros((npracticetrials,self.dsgn.setsize,3))
            #color assigned to each position
            for itrial in range(npracticetrials):
                practicetrials_colorind[itrial] = np.random.choice(self.dsgn.ncolors,
                    size=self.dsgn.setsize, replace=False)

            #set up the rgb vales for each trial
            for itrial in range(npracticetrials):
                
                for ii in range(self.dsgn.setsize):
                    print practicetrials_colorind[itrial][ii]
                    print self.dsgn.colors_rgb
                    print self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]
                    practicetrials_colorrgb[itrial][ii] = self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]

            #draw the squares
            for itrial in range(npracticetrials):
                self.disp.fix.fillColor = (-1,-1,-1)
                self.disp.fix.draw()
                for i in range(self.dsgn.setsize):
                    self.disp.circle.fillColor = practicetrials_colorrgb[itrial,i]
                    self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                    self.disp.circle.draw()         
            self.disp.win.flip()
            core.wait(self.dsgn.t_stim)
                
            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw() 
            self.disp.win.flip()
            core.wait(self.dsgn.t_poststim)               

            #draw response screen                   
            self.disp.fix.draw()
            for i in range(self.dsgn.setsize):
                for j in range(9):
                    self.disp.wholereportsquare[i,j].draw()
            self.disp.win.flip()

            #display mouse at fixation
            self.disp.mouse.setPos(newPos=(0,0))
            self.disp.mouse.setVisible(1)

            responded = np.zeros(self.dsgn.setsize)
            respondedind = np.zeros(self.dsgn.setsize)
            acc = np.zeros(self.dsgn.setsize)
            while np.any(responded==0):
                for i in range(6):
                    if responded[i]==0:
                        for j in range(9):
                            if self.disp.mouse.isPressedIn(self.disp.wholereportsquare[i,j],buttons=[0]):
                                #log the time that they responded to that item
                                responded[i] = 1
                                respondedind[i] = j
                                if respondedind[i]==practicetrials_colorind[0][i]:
                                    acc[i] = 1

                                #redraw the screen
                                self.disp.fix.draw()
                                for i in range(self.dsgn.setsize):
                                    for j in range(9):
                                        self.disp.wholereportsquare[i,j].draw()
                                        #add black lines around the responded item
                                        if responded[i]==1:
                                            self.disp.square_pos[i].draw()
                                self.disp.win.flip()

            #remove memory probe after last response
            self.disp.fix.draw()
            self.disp.win.flip()
            self.disp.mouse.setVisible(0)

            self.disp.text_instructions = visual.TextStim(self.disp.win, 
                text="You got " + str(int((np.sum(acc)))) + " correct.\n\n" +
                "This task is super hard, and on average people get 3 correct!\n\n" +
                "Press r to repeat this if you want more practice.",
                alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
            self.disp.text_instructions.draw()     #instructions text at top of screen
            self.disp.fix.draw()
            self.disp.text_advance.draw()           #advance text at bottom of screen
            self.disp.win.flip() 
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape+self.dsgn.key_repeat+self.dsgn.key_skip,f=f)  

            if response==self.dsgn.key_skip[0]:
                break

            if response=='space':
                break

            self.disp.fix.draw()
            self.disp.win.flip() 
            core.wait(.5)

    def instructionsBoth(self,clock,f=None,tk=None):

        while True:
            self.disp.text_instructions = visual.TextStim(self.disp.win, 
                text="Now let's try both tasks together. Remember, for the main task, press d for circles, s for squares.\n\n" + 
                "Then, when the Rubik's cubes appear, use the mouse to click what the colors of the squares had been",
                alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
            self.disp.text_instructions.draw()      #instructions text at top of screen
            self.disp.fix.draw()
            self.disp.text_advance.draw()           #advance text at bottom of screen
            self.disp.win.flip() 
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape+self.dsgn.key_skip,f=f)  
            
            if response == self.dsgn.key_skip[0]:
                break

            self.disp.fix.draw()
            self.disp.win.flip() 
            #start recording
            if tk is not None:
                tk.startRecording(1,1,1,1) #start 
                pylink.pumpDelay(100)
                
            core.wait(.5)        

            npracticetrials = 10
            acc = np.zeros(npracticetrials)       
            practicetrials_colorind = np.zeros((npracticetrials,self.dsgn.setsize))
            practicetrials_colorrgb = np.zeros((npracticetrials,self.dsgn.setsize,3))
            actual_response = np.zeros((npracticetrials))
            correct_response = [1,1,0,1,1,1,1,1,0,1]
            acc = np.zeros(npracticetrials)

            for itrial in range(npracticetrials):
                #color assigned to each position
                practicetrials_colorind[itrial] = np.random.choice(self.dsgn.ncolors,
                    size=self.dsgn.setsize, replace=False)

            #ensure that the same color doesn't appear in the same position back to back
            for ii in range(self.dsgn.setsize):
                diff_zero = np.diff(self.dsgn.trials_colorind[:,ii])
                for it in np.where(diff_zero==0)[0]:
                    if it > 0:
                        used_colors = np.append(np.append(practicetrials_colorind[it,:],practicetrials_colorind[it-1,ii]),practicetrials_colorind[it+1,ii])
                        avail_colors = np.setdiff1d(self.dsgn.colors_ind,used_colors)
                        practicetrials_colorind[it,ii] = np.random.choice(avail_colors,size=1)

            for itrial in range(npracticetrials):
                #set up the rgb vales for each trial
                for ii in range(self.dsgn.setsize):
                    practicetrials_colorrgb[itrial][ii] = self.dsgn.colors_rgb[int(practicetrials_colorind[itrial][ii])]

            for itrial in range(npracticetrials):
                self.disp.fix.fillColor = (-1,-1,-1)
                self.disp.fix.draw()
                if itrial==2 or itrial==8:
                    for i in range(self.dsgn.setsize):
                        self.disp.square.fillColor = practicetrials_colorrgb[itrial,i]
                        self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                        self.disp.square.draw() 
                else:
                    for i in range(self.dsgn.setsize):
                        self.disp.circle.fillColor = practicetrials_colorrgb[itrial,i]
                        self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                        self.disp.circle.draw()           
                self.disp.win.flip()

                #save timing information
                practicetrial_encarray_onset = clock.getTime()
                print practicetrial_encarray_onset
                
                response = None
                event.clearEvents(eventType='keyboard')
                while clock.getTime()<(practicetrial_encarray_onset+self.dsgn.t_stim):
                    if response is None:
                        response = poll_for_response(self.disp.win,self.dsgn.key_freq+self.dsgn.key_infreq+self.dsgn.key_escape)
                        if response is not None:
                            
                            if response==self.dsgn.key_freq[0]:
                                actual_response[itrial]=1
                            elif response==self.dsgn.key_infreq[0]:
                                actual_response[itrial]=0

                            #see if the response was correct or not
                            print itrial
                            print np.shape(correct_response)
                            print np.shape(actual_response)
                            acc[itrial]=(correct_response[itrial]==actual_response[itrial])*1
                            print response, self.dsgn.key_freq[0], actual_response[itrial], correct_response[itrial],acc[itrial]

                            #redraw the display
                            self.disp.fix.fillColor = (1,1,1)
                            self.disp.fix.draw()
                            if itrial==2 or itrial==8:
                                for i in range(self.dsgn.setsize):
                                    self.disp.square.fillColor = practicetrials_colorrgb[itrial,i]
                                    self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                                    self.disp.square.draw() 
                            else:
                                for i in range(self.dsgn.setsize):
                                    self.disp.circle.fillColor = practicetrials_colorrgb[itrial,i]
                                    self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                                    self.disp.circle.draw()           
                            self.disp.win.flip()

                if itrial==7:
                    self.disp.fix.fillColor = (-1,-1,-1)
                    self.disp.fix.draw() 
                    self.disp.win.flip()
                    core.wait(self.dsgn.t_poststim)               

                    #draw response screen                   
                    self.disp.fix.draw()
                    for i in range(self.dsgn.setsize):
                        for j in range(9):
                            self.disp.wholereportsquare[i,j].draw()
                    self.disp.win.flip()

                    #display mouse at fixation
                    self.disp.mouse.setPos(newPos=(0,0))
                    self.disp.mouse.setVisible(1)

                    count = 0
                    responded = np.zeros(self.dsgn.setsize)
                    respondedind = np.zeros(self.dsgn.setsize)
                    wholereportacc = np.zeros(self.dsgn.setsize)
                    while np.any(responded==0):
                        for i in range(6):
                            if responded[i]==0:
                                for j in range(9):
                                    if self.disp.mouse.isPressedIn(self.disp.wholereportsquare[i,j],buttons=[0]):
                                        #log the time that they responded to that item
                                        responded[i] = 1
                                        respondedind[i] = j
                                        if respondedind[i]==practicetrials_colorind[7][i]:
                                            wholereportacc[i] = 1

                                        #redraw the screen
                                        self.disp.fix.draw()
                                        for i in range(self.dsgn.setsize):
                                            for j in range(9):
                                                self.disp.wholereportsquare[i,j].draw()
                                                #add black lines around the responded item
                                                if responded[i]==1:
                                                    self.disp.square_pos[i].draw()
                                        self.disp.win.flip()

                    #remove memory probe after response
                    self.disp.fix.draw()
                    self.disp.win.flip()
                    self.disp.mouse.setVisible(0)
                    core.wait(self.dsgn.t_iti)

            self.disp.fix.draw()
            self.disp.win.flip() 
            if tk is not None:
                pylink.pumpDelay(100)
                tk.stopRecording()
                pylink.pumpDelay(100)
                tk.setOfflineMode()
            core.wait(.5)   

            print int(np.round(np.mean(acc),decimals=2)*100)
            print int(np.sum(wholereportacc))
            self.disp.text_instructions = visual.TextStim(self.disp.win, 
                text="On the main task, you got " + str(int(np.round(np.mean(acc),decimals=2)*100)) + " percent correct\n\n" + 
                "On the mouse clicking, you got " + str(int(np.sum(wholereportacc))) + " correct.\n\n"+
                "Press r to repeat these trials if you want more practice on both tasks.",
                alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
                fontFiles=[self.dsgn.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
            self.disp.text_instructions.draw()     #instructions text at top of screen
            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw()
            self.disp.text_advance.draw()           #advance text at bottom of screen
            self.disp.win.flip() 
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape+self.dsgn.key_repeat,f=f)  

            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw()
            self.disp.win.flip()
            core.wait(.5)

            if response == 'space':
                break

        #fixation and blinking instruction
        self.disp.text_instructions = visual.TextStim(self.disp.win, 
            text="In general, you should keep your head still and in the chin rest for eye tracking. You can blink when necessary but try to blink infrequently, about once every couple seconds.\n"+
            "One time when it's really important not to blink is when the screen goes blank before Rubik's Cubes appear. "+
            "If you blink during that time, the Rubik's Cube will be aborted and you'll see a message like the one in red below.\n\n"+
            "But, once the Rubik's Cube appears and while you are clicking with your mouse, you can blink as much as you want.",
            alignHoriz='center',wrapWidth=22,height=.75,pos=(0,5),color=(-1,-1,-1),
            fontFiles=[self.disp.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
        self.disp.text_instructions.draw()     #instructions text at top of screen
        self.disp.fix.fillColor = (-1,-1,-1)
        self.disp.fix.draw()
        self.disp.text_blink.draw()
        self.disp.text_advance.draw()           #advance text at bottom of screen
        self.disp.win.flip() 
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape+self.dsgn.key_repeat,f=f)  

        self.disp.fix.fillColor = (-1,-1,-1)
        self.disp.fix.draw()
        self.disp.win.flip() 
        core.wait(.5)


        self.disp.text_instructions = visual.TextStim(self.disp.win, 
            text="Please find the experimenter if you have any questions. " + 
                "Otherwise press spacebar to start! Good luck!",
            alignHoriz='center',wrapWidth=22,height=.75,pos=(0,7),color=(-1,-1,-1),
            fontFiles=[self.disp.project_dir + 'display/fonts/Lato-Reg.ttf'],font=['Lato'])
        self.disp.text_instructions.draw()      #instructions text at top of screen
        self.disp.fix.draw()
        self.disp.text_advance.draw()           #advance text at bottom of screen
        self.disp.win.flip() 
        wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)  
              
        self.disp.fix.draw()
        self.disp.win.flip() 
        core.wait(.5)   

    def startOfBlock(self,iblock,clock,f=None,tk=None):

        #display the start of block screen
        self.disp.text_block_start.setText('Block ' + str(iblock+1) + ' of ' + str(self.dsgn.nblocks+1))   #set the block text to be the right block number
        self.disp.text_block_start.draw()                       #draw the block text
        self.disp.fix.fillColor = (-1,-1,-1)
        self.disp.fix.draw()
        drawpath(self.disp,self.dsgn.nblocks+1,i_fill=iblock-1)   #display progress
        self.disp.text_advance.pos = (0,-4)
        self.disp.text_advance.draw()                           #text advance
        self.disp.win.flip()                                    #flip the window
        
        #wait for the response to start the block
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.dsgn.time_blockstart[iblock]=clock.getTime()


        self.disp.fix.draw()
        self.disp.win.flip()
#        if self.disp.env=='eeg': #Fix that #Start Block Portcode
#            parallel.setData(240+iblock)
        
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            tk.sendMessage(('SYNCTIME %d' %pylink.currentTime()))
            tk.sendCommand(("record_status_message 'Block %d of %d'" %(iblock+1,self.dsgn.nblocks)))
            tk.sendMessage(('Block %d of %d' %(iblock+1,self.dsgn.nblocks)))    
            
        core.wait(self.dsgn.t_poststim)
#        if self.disp.env=='eeg':
#            parallel.setData(0) #reset portcode


    def encodingArray(self,iblock,itrial,clock,f=None,tk=None):

        if tk is not None:
            tk.sendCommand(("record_status_message 'Block %d of %d Trial %d of %d'" %(iblock+1,self.dsgn.nblocks, itrial+1,self.dsgn.ntrials_perblock)))
            tk.sendMessage(('Trial %d of %d' %(itrial+1,self.dsgn.ntrials_perblock)))
        
#        for i in range(2): #present stimtrak with first flip, then turn stimtrak off with second flip
        self.disp.fix.fillColor = (-1,-1,-1)
        self.disp.fix.draw()
#            if i==0:
#                self.disp.stimtrak.draw()
        if self.dsgn.freq_trials[iblock,itrial]==1:
            for i in range(self.dsgn.setsize):
                self.disp.circle.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                self.disp.circle.draw()
        else:
            for i in range(self.dsgn.setsize):
                self.disp.square.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                self.disp.square.draw()         
        self.disp.win.flip()
#        if self.disp.env=='eeg':
#            parallel.setData(10+self.dsgn.freq_trials[iblock,itrial]) #freq/infreq portcode 10=infreq 11=freq
        if tk is not None:
            tk.sendMessage('encArray')
#        else:
#            if self.disp.env=='eeg':
#                parallel.setData(100+itrial%100) #trial portcode
                    
        #save timing information
        self.dsgn.time_encarray_onset[iblock,itrial] = clock.getTime()
        
        response = None
        response_cal = None
        event.clearEvents(eventType='keyboard')

        self.dsgn.eyeL[iblock,itrial],self.dsgn.eyeR[iblock,itrial] = PupilSize(tk)
        self.dsgn.pupil_runningavg[iblock,itrial] = np.mean((np.nanmean(self.dsgn.eyeR[iblock,0:itrial]),np.nanmean(self.dsgn.eyeL[iblock,0:itrial]))) #cumulative running average
        self.dsgn.pupil_runningstd[iblock,itrial] = np.mean((np.nanstd(self.dsgn.eyeR[iblock,0:itrial]),np.nanstd(self.dsgn.eyeL[iblock,0:itrial]))) #cumulative running standard deviation 
        self.dsgn.pupil_runninglowthresh[iblock,itrial] = self.dsgn.pupil_runningavg[iblock,itrial]-self.dsgn.pupil_runningstd[iblock,itrial] #small pupil size threshold #p of norm dist = .159
        self.dsgn.pupil_runninghighthresh[iblock,itrial] = self.dsgn.pupil_runningavg[iblock,itrial]+self.dsgn.pupil_runningstd[iblock,itrial] #large pupil size threshold
        self.dsgn.pupil_trailingavg[iblock,itrial] = np.mean((np.nanmean(self.dsgn.eyeR[iblock,(itrial-2):(itrial+1)]),np.nanmean(self.dsgn.eyeL[iblock,(itrial-2):(itrial+1)])))
        while clock.getTime()<(self.dsgn.time_encarray_onset[iblock,itrial]+self.dsgn.t_stim):
            if itrial>=self.dsgn.n_trials_atstart and itrial<self.dsgn.ntrials_perblock-2:
                if self.dsgn.within_tol[iblock,itrial] == 1 or np.isnan(self.dsgn.within_tol[iblock,itrial]):
                    self.dsgn.within_tol[iblock,itrial] = eyetracker_gazecontrol(self.disp,tk,fix_tol=10,fb=False)
                    
            if response is None:
                response = poll_for_response(self.disp.win,self.dsgn.key_freq+self.dsgn.key_infreq+self.dsgn.key_escape)
                if response is not None:
                    #log time response was made
                    self.dsgn.time_response[iblock,itrial] = clock.getTime()
                    self.dsgn.rts[iblock,itrial] = self.dsgn.time_response[iblock,itrial]-self.dsgn.time_encarray_onset[iblock,itrial]
#                    if self.disp.env=='eeg':
#                        parallel.setData(40) #response portcode
                    #calculate real time details
                    if itrial>=self.dsgn.n_trials_atstart and itrial<self.dsgn.ntrials_perblock-2 and tk is not None:
                        if np.all(self.dsgn.probe_trials[iblock,(itrial-2):(itrial+1)] == 0):
                            if self.dsgn.pupil_trailingavg[iblock,itrial]<self.dsgn.pupil_runninglowthresh[iblock,itrial] and np.nansum(self.dsgn.pupil_triggered_low[iblock,0:itrial])<self.dsgn.n_pupil_triggered_low:
                                self.dsgn.pupil_triggered_low[iblock,itrial] = 1
                            if self.dsgn.pupil_trailingavg[iblock,itrial]>self.dsgn.pupil_runninghighthresh[iblock,itrial] and np.nansum(self.dsgn.pupil_triggered_high[iblock,0:itrial])<self.dsgn.n_pupil_triggered_high:
                                self.dsgn.pupil_triggered_high[iblock,itrial] = 1

                        if self.dsgn.pupil_triggered_high[iblock,itrial] == 1 or self.dsgn.pupil_triggered_low[iblock,itrial]==1:
                            self.dsgn.pupil_triggered[iblock,itrial+1] = 1
                            self.dsgn.probe_trials[iblock,itrial+1] = 1
                            self.dsgn.freq_trials[iblock,itrial+1] = 1
                        
                    #log response that was made
                    self.dsgn.actual_response_string[iblock][itrial] = response[0]
                    
                    if response==self.dsgn.key_freq[0]:
                        self.dsgn.actual_response[iblock,itrial]=1
                    elif response==self.dsgn.key_infreq[0]:
                        self.dsgn.actual_response[iblock,itrial]=0

                    #see if the response was correct or not
                    self.dsgn.acc[iblock,itrial]=(self.dsgn.correct_response[iblock,itrial]==self.dsgn.actual_response[iblock,itrial])*1
                    
                    if np.any(self.dsgn.within_tol[iblock,(itrial-3):(itrial+1)] == 0):
                        if self.dsgn.probe_trials[iblock,itrial]==1 or self.dsgn.pupil_triggered_low[iblock,itrial] == 1 or self.dsgn.pupil_triggered_high[iblock,itrial] == 1:
                            self.dsgn.probe_trials[iblock,(itrial):(itrial+2)] = 0
                            self.dsgn.pupil_triggered[iblock,itrial] = 0
                            self.dsgn.pupil_triggered_high[iblock,(itrial-1):(itrial+1)] = 0
                            self.dsgn.pupil_triggered_low[iblock,(itrial-1):(itrial+1)] = 0

                    # if self.dsgn.acc[iblock,itrial]==0:
                    #     self.dsgn.probe_trials[iblock,itrial] = 0
                    #     self.dsgn.rt_triggered[iblock,itrial] = 0
                    #     self.dsgn.rt_triggered_fast[iblock,itrial] = 0
                    #     self.dsgn.rt_triggered_slow[iblock,itrial] = 0
                        # self.dsgn.rt_triggered_med[iblock,itrial] = 0
                    
                    #redraw the display
                    self.disp.fix.fillColor = (1,1,1)
                    self.disp.fix.draw()
                    if self.dsgn.freq_trials[iblock,itrial]==1:
                        for i in range(self.dsgn.setsize):
                            self.disp.circle.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                            self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                            self.disp.circle.draw()
                    else:
                        for i in range(self.dsgn.setsize):
                            self.disp.square.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                            self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                            self.disp.square.draw()
                    self.disp.win.flip()   

        if response is None:
            if itrial>=self.dsgn.n_trials_atstart:
                self.dsgn.probe_trials[iblock,itrial] = 0
                self.dsgn.pupil_triggered[iblock,itrial] = 0
                self.dsgn.pupil_triggered_low[iblock,itrial-1] = 0
                self.dsgn.pupil_triggered_high[iblock,itrial-1] = 0                 

        if response_cal is None:
            response_cal = poll_for_response(self.disp.win,self.dsgn.key_force_iti)
            if response_cal is not None:
                tk.stopRecording()
                pylink.pumpDelay(100)
                tk.setOfflineMode()
                pylink.pumpDelay(100)
                tk.doTrackerSetup()
                self.disp.fix.draw()
                self.disp.win.flip()
                core.wait(self.dsgn.t_iti)
                tk.startRecording(1,1,1,1) #start 

#
#        if self.disp.env=='eeg':
#                parallel.setData(0) #reset portcode
        
        
    def retentionInterval(self,iblock,itrial,clock,tk=None):
        self.disp.fix.fillColor = (-1,-1,-1)
        self.disp.fix.draw()
        self.disp.win.flip()
#        if self.disp.env=='eeg':
#            parallel.setData(20) #retention interval portcode
        #save timing information
        self.dsgn.time_encarray_offset[iblock,itrial] = clock.getTime()
        
        while clock.getTime()<(self.dsgn.time_encarray_offset[iblock,itrial]+self.dsgn.t_poststim):
            if self.dsgn.within_tol[iblock,itrial]==0:
                if itrial>=80:
                    break
            else:
                self.dsgn.within_tol[iblock,itrial] = eyetracker_gazecontrol(self.disp,tk,fix_tol=10,fb=False)
                
                
        if self.dsgn.within_tol[iblock,itrial]==0:
            print "text blink",iblock,itrial,self.dsgn.within_tol[iblock,itrial]
            self.dsgn.probe_trials[iblock,(itrial-1):(itrial+1)] = 0
            self.dsgn.pupil_triggered[iblock,itrial] = 0
            self.dsgn.pupil_triggered_high[iblock,(itrial-2):(itrial+1)] = 0
            self.dsgn.pupil_triggered_low[iblock,(itrial-2):(itrial+1)] = 0
            self.disp.fix.fillColor = (-1,-1,-1)
            self.disp.fix.draw()
            self.disp.text_blink.draw()
            self.disp.win.flip()
            core.wait(.5)
#
#        core.wait(self.dsgn.t_poststim)
        
#        if self.disp.env=='eeg':
#            parallel.setData(0) #reset portcode
        

    def memoryProbe(self,iblock,itrial,clock,f=None,tk=None):

        #draw response screen                   
        self.disp.fix.draw()
        #self.disp.stimtrak.draw()
        for i in range(self.dsgn.setsize):
            for j in range(9):
                self.disp.wholereportsquare[i,j].draw()
        self.disp.win.flip()
#        if self.disp.env=='eeg':
        if itrial<80:
#                parallel.setData(34) #initial 4 probe trials
            print 34
        elif self.dsgn.pupil_triggered_high[iblock,itrial-1]==1: #rt triggered fast portcode
#                parallel.setData(31)
            print 31
        elif self.dsgn.pupil_triggered_low[iblock,itrial-1]==1: #rt triggered slow portcode
#                parallel.setData(32)
            print 32
            # if self.dsgn.rt_triggered_med[iblock,itrial-1]==1: #rt triggered med portcode
            #     parallel.setData(33) 
        
        if tk is not None:
            if self.dsgn.pupil_triggered_high[iblock,itrial-1]==1:
                tk.sendMessage('largeTrigger')
            if self.dsgn.pupil_triggered_low[iblock,itrial-1]==1:
                tk.sendMessage('smallTrigger')
            tk.sendMessage('memProbe') 
            pylink.pumpDelay(100)
            tk.stopRecording()
            pylink.pumpDelay(100)
            tk.setOfflineMode()

        #save timing information
        self.dsgn.time_memprobe_onset[iblock,itrial] = clock.getTime()

        #display mouse at fixation
        self.disp.mouse.setPos(newPos=(0,0))
        self.disp.mouse.setVisible(1)

        count = 0
        event.clearEvents(eventType='keyboard')
        while np.any(np.isnan(self.dsgn.wholereport_resp[iblock,itrial])):
            response = poll_for_response(self.disp.win,self.dsgn.key_escape)
            for i in range(6):
                if np.isnan(self.dsgn.wholereport_resp[iblock,itrial][i]):
                    for j in range(9):
                        if self.disp.mouse.isPressedIn(self.disp.wholereportsquare[i,j],buttons=[0]):
                            #log the time that they responded to that item
                            self.dsgn.wholereport_rts[iblock,itrial][i] = clock.getTime()-self.dsgn.time_memprobe_onset[iblock,itrial]

                            #which color within that square did they click
                            self.dsgn.wholereport_respcolorind[iblock,itrial][i] = j
                            self.dsgn.wholereport_resprgb[iblock,itrial][i] = self.dsgn.colors_rgb[j]
                            self.dsgn.wholereport_respacc[iblock,itrial][i] = (self.dsgn.wholereport_respcolorind[iblock,itrial][i]==
                                                                            self.dsgn.trials_colorind[iblock,itrial,i])*1
                            #log that they have responded to that item
                            self.dsgn.wholereport_resp[iblock,itrial][i] = 1

                            #log the order that they responded to that item
                            count = count+1
                            self.dsgn.wholereport_resporder[iblock,itrial][i] = count

                            #redraw the screen
                            self.disp.fix.draw()
                            for i in range(self.dsgn.setsize):
                                for j in range(9):
                                    self.disp.wholereportsquare[i,j].draw()
                                    #add black lines around the responded item
                                    if self.dsgn.wholereport_resp[iblock,itrial][i]==1:
                                        self.disp.square_pos[i].draw()
                            self.disp.win.flip()

        #remove memory probe after response
        self.disp.fix.draw()
        self.disp.win.flip()
#        if self.disp.env=='eeg':
#            parallel.setData(0) #reset portcode

        #save timing information
        self.dsgn.time_memprobe_offset[iblock,itrial] = clock.getTime()
        self.disp.mouse.setVisible(0)

    def itiWM(self,iblock,itrial,clock,tk=None):
        self.disp.fix.draw()
        self.disp.win.flip()

#        if self.disp.env=='eeg':
#            parallel.setData(71)

        #save timing information
        self.dsgn.time_iti_onset[iblock,itrial] = clock.getTime()

        #core.wait(self.dsgn.t_iti)
        while clock.getTime()<(self.dsgn.time_iti_onset[iblock,itrial]+self.dsgn.t_iti):
            response = poll_for_response(self.disp.win,keylist=self.dsgn.key_calib[0])
            if response is not None:
                tk.doTrackerSetup()
                self.disp.fix.draw()
                self.disp.win.flip()
                core.wait(self.dsgn.t_iti)

        #start recording
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            tk.sendMessage(('SYNCTIME %d' %pylink.currentTime()))
            tk.sendCommand(("record_status_message 'iti'" ))
            tk.sendMessage(('iti'))
            
#        if self.disp.env=='eeg':
#            parallel.setData(0)

        #save timing information
        self.dsgn.time_iti_offset[iblock,itrial] = clock.getTime()

        
    def endOfBlock(self,iblock,clock,f=None,tk=None):

        if tk is not None: 
            tk.sendMessage('EndBlock') 
            pylink.pumpDelay(100)
            tk.stopRecording()
            pylink.pumpDelay(100)
            tk.setOfflineMode()

        #display the start of block screen
        self.disp.text_block_start.setText('Finished part '+str(iblock+1))   #set the block text to be the right block number
        self.disp.text_block_start.draw()
        #self.disp.text_blockAvgCorrect.setText('You got ' + str(int(np.round(np.nanmean(self.dsgn.wm_acc[iblock])*100))) + ' percent correct in the first part of that bloack')
        #self.disp.text_blockAvgCorrect.draw()
        #self.disp.fix.draw()
        drawpath(self.disp,self.dsgn.nblocks+1,i_fill=iblock-1)   #display progress
        self.disp.text_advance.pos = (0,-4)
        self.disp.text_advance.draw()                           #text advance
        self.disp.win.flip()     #flip the window
        
        #wait for the response to start the block
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)

        response_cal = poll_for_response(self.disp.win,keylist=self.dsgn.key_calib[0])
        if response_cal is not None:
            if tk is not None:
                tk.doTrackerSetup()
        
        self.disp.fix.draw()
        self.disp.win.flip()
        #self.dsgn.time_blockend[iblock] = clock.getTime()
        core.wait(.5)


    def changeDetectionInstructionsIntro(self,f=None):
        
        ##################################################
        ##INSTRUCTIONS SCREEN 0
        #display text
        self.disp.text_CDinstructions0.draw()
        
        #draw_example square
        self.disp.square.pos = (0,0)
        color_ind = np.random.randint(0,self.dsgn.cd_ncolors-3)
        self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind]
        self.disp.square.draw()
        self.disp.text_advance.draw()
        self.disp.win.flip()
        
        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        ##################################################
        ##INSTRUCTIONS SCREEN 1
        #display text
        self.disp.text_CDinstructions1.draw()
        
        #draw example_square
        self.disp.square.draw()
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        ##################################################
        ##INSTRUCTIONS SCREEN 2
        #display text
        self.disp.text_CDinstructions2.draw()
        
        #draw example_square
        new_color_ind = np.random.choice(np.setdiff1d(np.arange(self.dsgn.cd_ncolors-3),color_ind),replace=False)
        self.disp.square.fillColor = self.dsgn.cd_colors_rgb[new_color_ind]
        self.disp.square.draw()
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        ##################################################
        ##INSTRUCTIONS SCREEN 3
        #display text
        self.disp.text_CDinstructions3.draw()
        
        #draw example_square
        self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind]
        self.disp.square.draw()
        
        #draw other things
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)

    def changeDetectionInstructionsPracticeSS1(self,f=None,tk=None):
        #display instructions screen
        self.disp.text_CDinstructions4.draw()
        #drawpath(self.disp,self.dsgn.nblocks)
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        
        #start recording
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            
        core.wait(.5)
        
        acc = 0
        
        while acc == 0:
        
            #encoding array
            self.disp.fix.draw()
            ss = 1
            r =  self.dsgn.cd_rad_min + ( self.dsgn.cd_rad_max - self.dsgn.cd_rad_min)*np.random.random(ss) #radius of eccentricity for square position
            theta = np.random.random(ss)*90+90*np.random.randint(0,4) #[degrees] for square position
            x, y = tools.coordinatetools.pol2cart(theta,r,units='deg')
            self.disp.square.pos = (x[0],y[0])
            color_ind = np.random.randint(0,self.dsgn.cd_ncolors-3)
            self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind]
            self.disp.square.draw()
            self.disp.win.flip()
            core.wait(self.dsgn.cd_t_encarray)
            
            #retention interval
            self.disp.fix.draw()
            self.disp.win.flip()
            
            if tk is not None:
                pylink.pumpDelay(100)
                tk.stopRecording()
                pylink.pumpDelay(100)
                tk.setOfflineMode()
                
            core.wait(self.dsgn.cd_t_retinterval)
            
            #memory probe
            self.disp.fix.draw()
            self.disp.square.draw()
            self.disp.text_diff.draw()  #draw response for different trials
            self.disp.text_same.draw()  #draw response for same trials
            self.disp.win.flip()        #flip window
            
            #wait for responses
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
            
            if [response] == self.dsgn.key_same:
                self.disp.text_CDinstructionsCorrect.draw()
                acc = 1
            else:
                self.disp.text_CDinstructionsIncorrect.draw()
            self.disp.text_advance.draw()
            self.disp.win.flip()
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)

    def changeDetectionInstructionsPracticeSS2Diff(self,f=None,tk=None):
        
        #display instructions screen
        self.disp.text_CDinstructions5.draw()
        #drawpath(self.disp,self.dsgn.nblocks)
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        
        #start recording
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            
        core.wait(.5)
        
        acc = 0
        
        while acc == 0:
        
            #encoding array
            self.disp.fix.draw()
            ss = 2
            current_min_dist = 0
            while current_min_dist<self.dsgn.cd_stim_mindist:
                    
                #polar angle for square position
                r =  self.dsgn.cd_rad_min + ( self.dsgn.cd_rad_max - self.dsgn.cd_rad_min)*np.random.random(ss) #radius of eccentricity for square position
                theta = np.random.random(ss)*90+90*np.random.choice(np.arange(4),size=ss,replace=False)#[degrees] for square position
               
                #cartesian coordinates corresponding to polar position for square position
                x, y = tools.coordinatetools.pol2cart(theta,r,units='deg')
                
                #check that they all satisfy min_dist
                current_min_dist = np.min(pdist(np.transpose(np.vstack((x,y)))))
                
            color_ind = np.random.choice(np.arange(self.dsgn.cd_ncolors),size=ss,replace=False)
            for i in range(ss):
                self.disp.square.pos = (x[i],y[i])
                self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind[i]]
                self.disp.square.draw()
            self.disp.win.flip()
            
            if tk is not None:
                pylink.pumpDelay(100)
                tk.stopRecording()
                pylink.pumpDelay(100)
                tk.setOfflineMode()
                
            core.wait(self.dsgn.cd_t_encarray)
            
            #retention interval
            self.disp.fix.draw()
            self.disp.win.flip()
            core.wait(self.dsgn.cd_t_retinterval)
            
            #memory probe
            self.disp.fix.draw()
            self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind[0]]
            self.disp.square.draw()
            self.disp.text_diff.draw()  #draw response for different trials
            self.disp.text_same.draw()  #draw response for same trials
            self.disp.win.flip()        #flip window
            
            #wait for responses
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
            
            if [response] == self.dsgn.key_diff:
                self.disp.text_CDinstructionsCorrect.draw()
                acc = 1
            else:
                self.disp.text_CDinstructionsIncorrect.draw()
            self.disp.text_advance.draw()
            self.disp.win.flip()
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)

    def changeDetectionInstructionsPracticeSS6Same(self,f=None,tk=None):
            
            #display instructions screen
            self.disp.text_CDinstructions6.draw()
            self.disp.text_advance.draw()
            self.disp.win.flip()

            #wait for response
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            
            #start recording
            if tk is not None:
                tk.startRecording(1,1,1,1) #start 
                pylink.pumpDelay(100)
            
            core.wait(.5)
            
            acc = 0
            
            while acc == 0:
            
                #encoding array
                self.disp.fix.draw()
                ss = 6
                current_min_dist = 0
                while current_min_dist<self.dsgn.cd_stim_mindist:
                        
                    #polar angle for square position
                    r =  self.dsgn.cd_rad_min + ( self.dsgn.cd_rad_max - self.dsgn.cd_rad_min)*np.random.random(ss) #radius of eccentricity for square position
                    theta = np.random.random(ss)*90+90*np.random.choice(np.append(np.arange(4),np.arange(4)),size=ss,replace=False)#[degrees] for square position
                   
                    #cartesian coordinates corresponding to polar position for square position
                    x, y = tools.coordinatetools.pol2cart(theta,r,units='deg')
                    
                    #check that they all satisfy min_dist
                    current_min_dist = np.min(pdist(np.transpose(np.vstack((x,y)))))
                    
                color_ind = np.random.choice(np.arange(self.dsgn.cd_ncolors),size=ss,replace=False)
                for i in range(ss):
                    self.disp.square.pos = (x[i],y[i])
                    self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind[i]]
                    self.disp.square.draw()
                self.disp.win.flip()
                core.wait(self.dsgn.cd_t_encarray)
                
                #retention interval
                self.disp.fix.draw()
                self.disp.win.flip()
                
                if tk is not None:
                    pylink.pumpDelay(100)
                    tk.stopRecording()
                    pylink.pumpDelay(100)
                    tk.setOfflineMode()
                
                core.wait(self.dsgn.cd_t_retinterval)
                
                #memory probe
                self.disp.fix.draw()
                self.disp.square.draw()
                self.disp.text_diff.draw()  #draw response for different trials
                self.disp.text_same.draw()  #draw response for same trials
                self.disp.win.flip()        #flip window
                
                #wait for responses
                response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
                
                if [response] == self.dsgn.key_same:
                    self.disp.text_CDinstructionsCorrect.draw()
                    acc = 1
                else:
                    self.disp.text_CDinstructionsIncorrect.draw()
                self.disp.text_advance.draw()
                self.disp.win.flip()
                response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
                self.disp.win.flip()
                core.wait(.5)
                
            #display instructions screen
            self.disp.text_CDinstructions7.draw()
            self.disp.win.flip()

            #wait for response
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)

    def startOfChangeDetection(self,iblock,clock,f=None,tk=None):
        
        #display the start of block screen
        self.disp.text_block_start.setText('Start of last block')   #set the block text to be the right block number
        self.disp.text_block_start.draw()                       #draw the block text
        drawpath(self.disp,self.dsgn.nblocks+1,i_fill=iblock-1)   #display progress
        self.disp.text_either_key.draw()                        #text advance
        self.disp.win.flip()                                    #flip the window
#        if self.disp.env=='eeg':
#            parallel.setData(250) #cd block start portcode
            
        #start recording
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            tk.sendMessage(('SYNCTIME %d' %pylink.currentTime()))
            tk.sendCommand(("record_status_message 'cd_Block"))
            tk.sendMessage(('cd_Block'))

        #save timing information
        self.dsgn.cd_time_blockstart = clock.getTime()

        #wait for the response to start the block
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)

        core.wait(1)

    def changeDetectionEncodingArray(self,iblock,itrial,clock,f=None,tk=None):

#        if self.disp.env=='eeg':
#            parallel.setData(100+itrial) #cd trial portcode
        
        #draw memory array
        self.disp.fix.draw()    #draw fixation
        for i in range(self.dsgn.cd_setsize):
            self.disp.square.pos = (self.dsgn.cd_encarray_x[itrial][i],self.dsgn.cd_encarray_y[itrial][i])
            self.disp.square.fillColor = self.dsgn.cd_encarray_colorrgb[itrial][i]
            self.disp.square.draw()
        self.disp.win.flip()    #flip the window
#        if self.disp.env=='eeg':
#            parallel.setData(80) #cd trial portcode
            
        if tk is not None:
            tk.sendCommand(("record_status_message 'cd_Block of cd_Trial %d of %d'" %(itrial+1,self.dsgn.ntrials_perblock)))
            tk.sendMessage(('cd_Trial %d of %d' %(itrial+1,self.dsgn.ntrials_perblock)))
        

        #save timing information
        self.dsgn.cd_time_encarray_onset[itrial] = clock.getTime()
        
        #wait for the duration that stimuli should be on the screen
        core.wait(self.dsgn.cd_t_encarray)
#        if self.disp.env=='eeg':
#            parallel.setData(0) #reset portcode

    def changeDetectionRetentionInterval(self,iblock,itrial,clock,f=None):
        
        #draw blank retention interval
        self.disp.fix.draw()
        self.disp.win.flip()
#        if self.disp.env=='eeg':
#            parallel.setData(90) #cd retention interval
        
        #save timing information
        self.dsgn.cd_time_encarray_offset[itrial] = clock.getTime()
        
        #wait for the duration of the retention interval
        core.wait(self.dsgn.cd_t_retinterval)
#        if self.disp.env=='eeg':
#            parallel.setData(0) #reset portcode


    def changeDetectionMemoryProbe(self,iblock,itrial,clock,f=None,tk=None):
        
        #display probe image
        self.disp.fix.draw()        #fixation dot
        self.disp.square.pos = (self.dsgn.cd_memprobe_x[itrial],self.dsgn.cd_memprobe_y[itrial])
        self.disp.square.fillColor = self.dsgn.cd_memprobe_probecolorrgb[itrial]
        self.disp.square.draw()
        self.disp.text_diff.draw()  #draw response for different trials
        self.disp.text_same.draw()  #draw response for same trials
        self.disp.win.flip()        #flip window
#        if self.disp.env=='eeg':
#            parallel.setData(91) #cd probe portcode
            
        if tk is not None:
            tk.sendMessage('cd_memProbe') 
            pylink.pumpDelay(100)
            tk.stopRecording()
            pylink.pumpDelay(100)
            tk.setOfflineMode()
        
        #save timing information
        self.dsgn.cd_time_memprobe_onset[itrial] = clock.getTime()

        #wait for responses
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
        
        #collect responses and RT
        self.dsgn.cd_time_memprobe_offset[itrial] = clock.getTime()
        self.dsgn.cd_rt[itrial]=self.dsgn.cd_time_memprobe_offset[itrial]-self.dsgn.cd_time_memprobe_onset[itrial]
        self.dsgn.cd_actual_response_string[itrial] = response
        if [self.dsgn.cd_actual_response_string[itrial]] == self.dsgn.key_same:
            self.dsgn.cd_actual_response[itrial] = 1
        else:
            self.dsgn.cd_actual_response[itrial] = 0
        
        self.dsgn.cd_acc[itrial] = (self.dsgn.cd_actual_response[itrial]==self.dsgn.cd_correct_response[itrial])*1 

        #remove memory probe after response
        self.disp.fix.draw()
        self.disp.win.flip()
#        if self.disp.env=='eeg':
#            parallel.setData(0) #reset portcode
        
        #save timing information
        self.dsgn.cd_time_memprobe_offset[itrial] = clock.getTime()

    def changeDetectionITI(self,iblock,itrial,clock,f=None,tk=None):
        
        #draw blank iti
        self.disp.fix.draw()
        self.disp.win.flip()
        
        #start recording
        if tk is not None:
            tk.startRecording(1,1,1,1) #start 
            pylink.pumpDelay(100)
            tk.sendMessage(('SYNCTIME %d' %pylink.currentTime()))
            tk.sendCommand(("record_status_message 'cd_iti'" ))
            tk.sendMessage(('cd_iti'))
        
        #log time
        self.dsgn.cd_time_iti_onset[itrial] = clock.getTime()
        
        #wait for the duration of the iti
        core.wait(self.dsgn.cd_t_iti[itrial])
        
        #log time
        self.dsgn.cd_time_iti_offset[itrial] = clock.getTime()

    def endOfExperiment(self,f=None,tk=None):
        self.disp.text_end_expt.draw()
        self.disp.win.flip()
        if tk is not None:
            pylink.pumpDelay(100)
            tk.stopRecording()
            pylink.pumpDelay(100)
            tk.setOfflineMode()
        keys = event.waitKeys(keyList=['space','escape']) 
        