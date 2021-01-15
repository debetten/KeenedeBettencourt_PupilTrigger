1. The Pylink library packaged in Psychopy may not be the one that is compatible with the most recent
Eyelink API (included in the Eyelink Developer's Kit). To check this, go to the "Shell" window of Psychopy
and execute the following two commands

      import pylink
      print pylink.__version__

If the returned value is something 1.0.0.X, please follow the instructions below to replace the Pylink 
library in the Psychopy installation directory. Also, it is a good indication that the Pylink library
in Psychopy is too old if you see a tiny camera image.
      
      a) Find the Pylink library installed with the Developer's Kit. On Windows, go to Start menu -> 
         All Programs -> SR Research -> Eyelink Examples -> Pylink Examples. On Mac, please go to 
         /Applications/Eyelink/pylink. 

      b) On both systems, you will find multiple pylink folders, of which
         the folder name indicates which version of Python it is compatible with. Please rename the folder
         that is compatible with your Python version (should be Python 2.7 in Psychopy) to "pylink" and
         use it to replace the one in your Psychopy installation folder. On Windows, it should be ...\
         Program Files\Psychopy2\Lib\site-packages\. On Mac, it should be under /Applications/Psychopy2/
         Contents/Resources\lib\python2.7.


2. Please note that this library requires a relatively newer version of Psychopy (v 1.84.2 and above).
One of the functionalities implemented is to use Alt + arrow key combinations to adjust the search
limits shown on the host PC display. This requires the detection of key modifiers, which are not 
implemented in earlier versions of Psychopy.


3. The .wav files are for the beeps given during calibration/validation and drift-check/correction.
The .wav files (beeps) won't be played because, for compatibility considerations, the audio playback
routines in the EyeLinkCoreGraphicsPsychoPy.py library have all been commented out. So, it is safe
to simply delete these .wav file.

4. See below for a brief summary of the example scripts posted together with this library.

gc_example.py: This script implements a free viewing task. It shows how to overlay a gaze cursor on
               an image. It also contains pointers to different sections of the user manuals. The
               code is heavily commented; users should be easily figure out how to connect to the 
               tracker, to set up tracking options, to send messages, to set interest areas and trial
               variables, to evoke the calibration/drift-check routines, to retrieve current gaze 
               position, to draw reference visuals on the host PC, to open/close data file, to transfer
	       data file, etc.

animated_target.py: Occasionally, users may want to use animated calibration targets, rather than the
               boring calibration dots. We could tweak the EyelinkCoreGraphicsPsychoPy.py library a 
               bit to use the powerful drawing functions to create animation targets that are native
               to Psychopy.

#### The following two examples may not work in Psychopy v. 1.85.1
video_example.py: This example shows how to send !V VFRAME messages to the tracker, so we can overlay
               gaze on videos in Data Viewer.


stroop_example.py: This script shows the bitmap backdropping function. Backdropping can be slow when 
               the image we want to show on the host PC is large. It is not recommended for timing 
               critical tasks.               