import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import csv
import matplotlib.font_manager as fm
from matplotlib import rc
import scipy.io as sio
import glob
import platform 
import matplotlib.gridspec as gridspec
import matplotlib.pylab as pylab
from scipy.stats import spearmanr
os.chdir('prep')
from req_fns import prettify_plot, scatter_plot_data
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
#project directory
os.chdir('../') 
base_project_dir= os.getcwd()
figure_dir = base_project_dir + '/figures/'
project_dir = base_project_dir + 'expt2/'
beh_dir = project_dir + 'beh/'
eye_dir = project_dir + 'eye/'
results_dir = project_dir + 'results/'

#index subject names
beh_fn = [f for f in os.listdir(beh_dir ) if f.endswith(project_name + '_explog.csv')]
beh_dat = dict.fromkeys(beh_fn)
nsubj = len(beh_fn)

#example data
isubj=5
ifn = beh_fn[isubj]
beh_dat = pd.read_csv(beh_dir+ifn ,header=12) 
ex_fast_triggered_ind = np.ravel(beh_dat.pupil_triggered_low)==1
ex_slow_triggered_ind = np.ravel(beh_dat.pupil_triggered_high)==1
ex_rt_thresh_fast = np.ravel(beh_dat.pupil_runninglowthresh)
ex_rt_thresh_slow = np.ravel(beh_dat.pupil_runninghighthresh)
ex_rt_trail = np.ravel(beh_dat.pupil_trailingavg)
ex_rt_runningavg = np.ravel(beh_dat.pupil_runningavg)
    
wmdat = pd.read_csv(results_dir + 'expt2_attnwmresults.csv')

#pupil size
s_large = np.array([ np.mean(wmdat.s_trail_online[(np.logical_and(
    wmdat.subject_num==isubj,wmdat.triggered_high==1))]) for isubj in range(nsubj)])
s_small = np.array([ np.mean(wmdat.s_trail_online[(np.logical_and(
    wmdat.subject_num==isubj,wmdat.triggered_low==1))]) for isubj in range(nsubj)])
s_mean = (s_large+s_small)/2

#rt
rts_large = np.array([ np.mean(wmdat.rt_trail[(np.logical_and(
    wmdat.subject_num==isubj,wmdat.triggered_high==1))])*1000 for isubj in range(nsubj)])
rts_small = np.array([ np.mean(wmdat.rt_trail[(np.logical_and(
    wmdat.subject_num==isubj,wmdat.triggered_low==1))])*1000 for isubj in range(nsubj)])
rts_mean = (rts_large+rts_small)/2

#working memory
m_large = np.array([ np.mean(wmdat.wm_perf[(np.logical_and(
    wmdat.subject_num==isubj,wmdat.triggered_high==1))]) for isubj in range(nsubj)])
m_small = np.array([ np.mean(wmdat.wm_perf[(np.logical_and(
    wmdat.subject_num==isubj,wmdat.triggered_low==1))]) for isubj in range(nsubj)])
m_mean = (m_large+m_small)/2


#make figure
fig = plt.figure(figsize=(16,12))
gs = fig.add_gridspec(2, 3,height_ratios=[1,1.5])


#figure 1a: triggering example
f3_ax1 = fig.add_subplot(gs[0, :])
xmin = 2750
xmax = 2800
tw = np.arange(xmin,xmax)
f3_ax1.plot(tw,ex_rt_thresh_fast[tw],color=col_incorr,ls='--',lw=2)
f3_ax1.plot(tw,ex_rt_thresh_slow[tw],color=col_corr,ls='--',lw=2)
f3_ax1.plot(tw,ex_rt_trail[tw],color='k',lw=2)
f3_ax1.scatter(ex_fast_triggered_ind[tw]*tw,np.zeros((len(tw)))+850,s=150,color=col_incorr,edgecolor='None')
f3_ax1.scatter(ex_slow_triggered_ind[tw]*tw,np.zeros((len(tw)))+1280,s=150,color=col_corr,edgecolor='None')
prettify_plot(f3_ax1,xlim=(xmin,xmax),ylim=[800,1300],
yt=([800,900,1000,1100,1200,1300]),ytl=([800,900,1000,1100,1200,1300]),yl="Pupil size\n(au)",
xt=([2750,2760,2770,2780,2790,2800]),xtl=([2750,2760,2770,2780,2790,2800]),xl="Trial (#)")

#figure 1b: pupil size
f3_ax4 = fig.add_subplot(gs[1, 0])
f3_ax4.plot([-.5,1.5],[0,0],'--',color='gray',zorder=1)
scatter_plot_data(f3_ax4,s_small-s_mean,x=0,c=col_incorr)
scatter_plot_data(f3_ax4,s_large-s_mean,x=1,c=col_corr)
prettify_plot(f3_ax4,xlim=(-.5,1.5),ylim=[-600,600],
 yt=([-600,-300,0,300,600]),ytl=([-600,-300,0,300,600]),yl="Pupil size\n(au, relative to mean)",
 xt=([0,1]),xtl=(['Small','Large']),xl="Probe")

#figure 1c: RTs 
f3_ax2 = fig.add_subplot(gs[1, 1])
f3_ax2.plot([-.5,1.5],[0,0],'--',color='gray',zorder=1)
scatter_plot_data(f3_ax2,rts_small-rts_mean,x=0,c=col_incorr)
scatter_plot_data(f3_ax2,rts_large-rts_mean,x=1,c=col_corr)
prettify_plot(f3_ax2,xlim=(-.5,1.5),ylim=[-40,40],
 yt=([-40,-20,0,20,40]),ytl=([-40,-20,0,20,40]),yl="Response time\n(ms, relative to mean)",
 xt=([0,1]),xtl=(['Small','Large']),xl="Probe")

#figure 1d: working memory
f3_ax3 = fig.add_subplot(gs[1, 2])
f3_ax3.plot([-.5,1.5],[0,0],'--',color='gray',zorder=1)
scatter_plot_data(f3_ax3,(m_small-m_mean),x=0,c=col_incorr)
scatter_plot_data(f3_ax3,(m_large-m_mean),x=1,c=col_corr)
prettify_plot(f3_ax3,xlim=(-.5,1.5),ylim=[-.4,.4],
 yt=([-.4,-.2,0,.2,.4]),ytl=([-.4,-.2,0,.2,.4]),yl="Working memory\n(# correct, relative to mean)",
 xt=([0,1]),xtl=(['Small','Large']),xl="Probe")


plt.subplots_adjust(wspace=.75)
plt.subplots_adjust(hspace=.75)

fig.savefig(figure_dir + 'figure4.pdf', bbox_inches='tight')

