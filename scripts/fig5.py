import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib.font_manager as fm
from matplotlib import rc
import scipy.io as sio
os.chdir('prep')
from req_fns import prettify_plot
os.chdir('../') 

#font defaults
plt.rcParams.update({'font.size': 28})
rc('text', usetex=False)
plt.rcParams['pdf.fonttype'] = 42
if os.path.isfile("/Library/Fonts/HelveticaNeue.ttf"): 
    prop = fm.FontProperties(fname="/Library/Fonts/HelveticaNeue.ttf",size=28)
    prop_light = fm.FontProperties(fname='/Library/Fonts/HelveticaNeue-Light.ttf',size=28) 
else:
    prop = fm.FontProperties(size=28)
    prop_light = fm.FontProperties(size=28)

#color defaults
col_incorr = [30/255.,50/255.,188/255.]
col_corr = [30/255.,188/255.,50/255.]

#project directory
os.chdir('../') 
base_project_dir= os.getcwd()
figure_dir = base_project_dir + '/figures/'

#EXPERIMENT 1
project_name = 'eyeSustAttnWM01'
project_dir = base_project_dir + '/expt1/'
beh_dir = project_dir + '/beh/'
eye_dir = project_dir + '/eye/'
beh_fn = [f for f in sorted(os.listdir(beh_dir )) if f.endswith(project_name + '_explog.csv')]
eye_fn_stim = [f for f in sorted(os.listdir(eye_dir )) if f.endswith(project_name + '_ENCARRAY_EYE_SEG.mat')]

#choose example subject
isubj = 7 

#behavioral data
beh_dat = pd.read_csv(beh_dir+beh_fn[isubj] ,header=12) 
e1_rts_trail = (np.ravel(beh_dat.rt_trailingavg)-np.ravel(beh_dat.rt_runningavg))*1000 #calculate trailing RT 
e1_i_rts_trail = np.logical_and(~np.isnan(e1_rts_trail),e1_rts_trail!=0)
e1_i_low = np.ravel(beh_dat.rt_trailingavg)<(np.ravel(beh_dat.rt_runningavg)-beh_dat.rt_runningstd)
e1_i_high = np.ravel(beh_dat.rt_trailingavg)>(np.ravel(beh_dat.rt_runningavg)+beh_dat.rt_runningstd)
e1_rts_low = (np.ravel(beh_dat.rt_trailingavg)[e1_i_low]-np.ravel(beh_dat.rt_runningavg)[e1_i_low])*1000
e1_rts_high = (np.ravel(beh_dat.rt_trailingavg)[e1_i_high]-np.ravel(beh_dat.rt_runningavg)[e1_i_high])*1000

#pupil data
mat_contents = sio.loadmat(eye_dir+eye_fn_stim[isubj],struct_as_record=False)
tpts=np.arange(0,901)
s_stim = mat_contents['eyeData'][0][0].trial[0][0].pa[:,0,tpts]
e1_pupil_trail = np.zeros(np.shape(np.ravel(beh_dat.rts)))
e1_pupil_trail[:] = np.nan
for itrial in range(len(e1_pupil_trail)): #calculate average pretrial pupil
    if itrial%800>80:
        e1_pupil_trail[itrial] = np.nanmean(s_stim[(itrial-3):itrial])-np.nanmean(s_stim[:itrial])
e1_i_pupil_trail = np.logical_and(~np.isnan(e1_pupil_trail),e1_pupil_trail!=0)
e1_pupil_low = e1_pupil_trail[e1_i_low]
e1_pupil_high = e1_pupil_trail[e1_i_high]


## EXPERIMENT 2
project_name = 'eyeSustAttnWM03'
project_dir = base_project_dir + '/expt2/'
beh_dir = project_dir + 'beh/'
beh_fn = [f for f in os.listdir(beh_dir ) if f.endswith(project_name + '_explog.csv')]

#choose example subject
isubj=0  

#pupil data
beh_dat = pd.read_csv(beh_dir+beh_fn[isubj] ,header=12) 
e2_small_triggered_ind = np.where(np.ravel(beh_dat.pupil_triggered_low)==1)[0]
e2_large_triggered_ind = np.where(np.ravel(beh_dat.pupil_triggered_high)==1)[0]
e2_pupil_trail = np.ravel(beh_dat.pupil_trailingavg)-np.ravel(beh_dat.pupil_runningavg)
e2_i_low = np.ravel(beh_dat.pupil_trailingavg)<(np.ravel(beh_dat.pupil_runningavg)-beh_dat.pupil_runningstd)
e2_pupil_low = np.ravel(beh_dat.pupil_trailingavg)[e2_i_low]-np.ravel(beh_dat.pupil_runningavg)[e2_i_low]
e2_i_high = np.ravel(beh_dat.pupil_trailingavg)>(np.ravel(beh_dat.pupil_runningavg)+beh_dat.pupil_runningstd)
e2_pupil_high = np.ravel(beh_dat.pupil_trailingavg)[e2_i_high]-np.ravel(beh_dat.pupil_runningavg)[e2_i_high]

#behavioral data
e2_rts_trail = np.zeros(np.shape(np.ravel(beh_dat.rts)))
e2_rts_trail[:] = np.nan
for itrial in range(len(e2_rts_trail)):
    if itrial%800>80:
        e2_rts_trail[itrial] = np.nanmean(np.ravel(beh_dat.rts)[(itrial-3):itrial])*1000-np.nanmean(np.ravel(beh_dat.rts)[:itrial])*1000
e2_rts_low = e2_rts_trail[e2_i_low]
e2_rts_high = e2_rts_trail[e2_i_high]


#make figure
fig,ax = plt.subplots(2,2,figsize=(10,8))

#upper left: behavioral triggering, behavior
bw=30
bins = np.arange(np.nanmin(e1_rts_trail[e1_i_rts_trail]),np.nanmax(e1_rts_trail[e1_i_rts_trail]),bw)
y,b = np.histogram(e1_rts_trail[e1_i_rts_trail],bins=bins,density=False)
ax[0, 0].bar(b[:-1],y,np.int(bw*.8),facecolor='gray',alpha=.5,edgecolor='None')
y,b = np.histogram(e1_rts_low,bins=bins,density=False)
ax[0, 0].bar(b[:-1],y,np.int(bw*.8),facecolor=col_incorr,alpha=1,edgecolor='None')
y,b = np.histogram(e1_rts_high,bins=bins,density=False)
ax[0, 0].bar(b[:-1],y,np.int(bw*.8),facecolor=col_corr,alpha=1,edgecolor='None')
prettify_plot(ax[0, 0],xlim=(-500,500),ylim=[0,500],
  xt=([-500,0,500]),xtl=([-500,0,500]),xl="RTs (ms, relative to mean)",
  yt=([0,250,500]),ytl=([0,250,500]),yl="# trials")

#lower left: behavioral triggering, pupil
bw=100
bins = np.arange(np.nanmin(e1_pupil_trail[e1_i_pupil_trail]),np.nanmax(e1_pupil_trail[e1_i_pupil_trail]),bw)
y,b = np.histogram(e1_pupil_trail[~np.isnan(e1_pupil_trail)],bins=bins,density=False)
ax[1,0].bar(b[:-1],y,np.int(bw*.8),facecolor='gray',alpha=.5,edgecolor='None')
y,b = np.histogram(e1_pupil_high,bins=bins,density=False)
ax[1,0].bar(b[:-1],y,np.int(bw*.8),facecolor=col_corr,alpha=1,edgecolor='None')
y,b = np.histogram(e1_pupil_low,bins=bins,density=False)
ax[1,0].bar(b[:-1],y,np.int(bw*.8),facecolor=col_incorr,alpha=1,edgecolor='None')
prettify_plot(ax[1,0],xlim=(-1500,1500),ylim=[0,500],
  xt=([-1500,0,1500]),xtl=([-1500,0,1500]),xl="Pupil size (au, relative to mean)",
  yt=([0,250,500]),ytl=([0,250,500]),yl="# trials")

#upper right: pupil triggering, behavior
bw=30
bins = np.arange(np.nanmin(e2_rts_trail),np.nanmax(e2_rts_trail),bw)
y,b = np.histogram(e2_rts_trail[~np.isnan(e2_rts_trail)],bins=bins,density=False)
ax[0, 1].bar(b[:-1],y,np.int(bw*.8),facecolor='gray',alpha=.5,edgecolor='None')
y,b = np.histogram(e2_rts_low,bins=bins,density=False)
ax[0, 1].bar(b[:-1],y,np.int(bw*.8),facecolor=col_incorr,alpha=1,edgecolor='None')
y,b = np.histogram(e2_rts_high,bins=bins,density=False)
ax[0, 1].bar(b[:-1],y,np.int(bw*.8),facecolor=col_corr,alpha=1,edgecolor='None')
prettify_plot(ax[0, 1],xlim=(-500,500),ylim=[0,500],
  xt=([-500,0,500]),xtl=([-500,0,500]),xl="RTs (ms, relative to mean)",
  yt=([0,250,500]),ytl=([0,250,500]),yl="# trials")

#lower right: pupil triggering, pupil
bw=100
bins = np.arange(np.nanmin(e2_pupil_trail),np.nanmax(e2_pupil_trail),bw)
y,b = np.histogram(e2_pupil_trail[~np.isnan(e2_pupil_trail)],bins=bins,density=False)
ax[1,1].bar(b[:-1],y,np.int(bw*.8),facecolor='gray',alpha=.5,edgecolor='None')
y,b = np.histogram(e2_pupil_low,bins=bins,density=False)
ax[1,1].bar(b[:-1],y,np.int(bw*.8),facecolor=col_incorr,alpha=1,edgecolor='None')
y,b = np.histogram(e2_pupil_high,bins=bins,density=False)
ax[1,1].bar(b[:-1],y,np.int(bw*.8),facecolor=col_corr,alpha=1,edgecolor='None')
prettify_plot(ax[1,1],xlim=(-1500,1500),ylim=[0,500],
  xt=([-1500,0,1500]),xtl=([-1500,0,1500]),xl="Pupil size (au, relative to mean)",
  yt=([0,250,500]),ytl=([0,250,500]),yl="# trials")

#resize and save figure
plt.subplots_adjust(wspace=2,hspace=2)
fig.savefig(figure_dir + 'figure5.pdf', bbox_inches='tight')

