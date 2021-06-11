import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib.font_manager as fm
from matplotlib import rc
import scipy.io as sio
import platform 
os.chdir('prep')
from req_fns import resampling_statistics1d, prettify_plot, scatter_plot_data, resampling_statistics
os.chdir('../../') 

#font defaults
plt.rcParams.update({'font.size': 32})
rc('text', usetex=False)
plt.rcParams['pdf.fonttype'] = 42
if os.path.isfile("/Library/Fonts/HelveticaNeue.ttf"): 
    prop = fm.FontProperties(fname="/Library/Fonts/HelveticaNeue.ttf",size=24)
    prop_light = fm.FontProperties(fname='/Library/Fonts/HelveticaNeue-Light.ttf',size=28) 
else:
    prop = fm.FontProperties(size=24)

#color defaults
col_incorr = [30/255.,50/255.,188/255.]
col_corr = [30/255.,188/255.,50/255.]
col_lure = [30/255.,119/255.,119/255.]


############## Driver #################

#project directory
base_project_dir=os.getcwd()
project_name = 'eyeSustAttnWM01'
project_dir = base_project_dir + 'expt1/'
beh_dir = project_dir + 'beh/'
eye_dir = project_dir + 'eye/'
figure_dir = base_project_dir + 'figures/'
results_dir = project_dir + 'results/'

#index subject names
beh_fn = [f for f in os.listdir(beh_dir ) if f.endswith(project_name + '_explog.csv')]
eye_fn_stim = [f for f in os.listdir(eye_dir ) if f.endswith(project_name + '_ENCARRAY_EYE_SEG.mat')]
eye_fn_ret = [f for f in os.listdir(eye_dir ) if f.endswith(project_name + '_EYE_SEG.mat')]

#preallocate
nt_stim = 2601 #number of timepoints
tpts_ret = np.arange(0,4500)
nt_ret = np.size(tpts_ret)

nsubj = 30
s_infreq = np.zeros((nsubj,nt_stim))
s_freq = np.zeros((nsubj,nt_stim))
s_ret = np.zeros((nsubj,nt_ret))
## Calculate average pupil timecourse based on sustained attention trial type
for isubj,ifn in enumerate(beh_fn):
    print(isubj,ifn,eye_fn_stim[isubj])

    #load behavioral data
    if ifn=='0905181_eyeSustAttnWM01_explog.csv' or ifn=='0917181_eyeSustAttnWM01_explog.csv':
        beh_dat = pd.read_csv(beh_dir+ifn ,header=11)
    else:
        beh_dat = pd.read_csv(beh_dir+ifn ,header=12)

    #trial conditions
    acc = np.ravel(beh_dat.acc)
    freq_bool = np.ravel(beh_dat.freq_trials==1)
    infreq_bool = np.ravel(beh_dat.freq_trials==0)
    end_block_bool = np.ravel(beh_dat.trial<798) #for choosing eye data for freq trials
    probe_ind = np.where(beh_dat.probe_trials==1)[0] #probes

    #missing eyetracking data 
    if (ifn == '0907182_eyeSustAttnWM01_explog.csv'): 
        j=np.arange(719,721)
        k=np.arange(719,721)
    elif (ifn == '0911184_eyeSustAttnWM01_explog.csv'):
        j=np.arange(3094,3104)
        k=np.arange(3094,3104)
    elif (ifn == '0912183_eyeSustAttnWM01_explog.csv'):
        j=np.append(np.arange(256,261),np.arange(1711,1721))
        k=np.append(np.arange(256,261),np.arange(1711,1721))
    elif (ifn == '0917182_eyeSustAttnWM01_explog.csv'):
        j=np.arange(204,208)
        k=np.arange(203,208)
    else:
        j=[]
        k=[]

    if np.size(j)>0:
        acc = np.delete(acc,j)
        freq_bool = np.delete(freq_bool,j)
        infreq_bool = np.delete(infreq_bool,j)
        end_block_bool = np.delete(end_block_bool,j)
    if np.size(k)>0:
        probe_ind = np.setdiff1d(probe_ind,k)

    #frequent/infrequent trials
    mat_contents = sio.loadmat(eye_dir+eye_fn_stim[isubj],struct_as_record=False)
    eye_data = mat_contents['eyeData'][0][0]
    times_stim = eye_data.trial[0][0].times[0]
    arf_stim = eye_data.arf[0][0].parserBlinks[0]+eye_data.arf[0][0].missingPupil[0]+eye_data.arf[0][0].saccadeX[0]+eye_data.arf[0][0].saccadeY[0]
    s_stim = eye_data.trial[0][0].pa[:,0,:]
    bl = np.nanmean(s_stim[:,0:100],axis=1)
    #calculate mean timecourse for each stimulus array
    if len(freq_bool)==len(arf_stim) and np.shape(s_stim)[1]==2601:
        i = np.logical_and(np.logical_and(infreq_bool==1,arf_stim==0),end_block_bool==1)
        s_infreq[isubj] = np.nanmean(s_stim[i]-np.tile(bl[i],(nt_stim,1)).T,axis=0)

        i = np.logical_and(np.logical_and(infreq_bool==0,arf_stim==0),end_block_bool==1)
        s_freq[isubj] = np.nanmean(s_stim[i]-np.tile(bl[i],(nt_stim,1)).T,axis=0)
    else:
        print('STIM ERROR'+ifn,np.shape(s_stim)[1])
        s_infreq[isubj] = np.nan
        s_freq[isubj] = np.nan

    #probe trials
    mat_contents = sio.loadmat(eye_dir+eye_fn_ret[isubj],struct_as_record=False)
    eye_data = mat_contents['eyeData'][0][0]
    times_ret = eye_data.trial[0][0].times[0,tpts_ret]
    arf_ret = eye_data.arf[0][0].parserBlinks[0]+eye_data.arf[0][0].missingPupil[0]+eye_data.arf[0][0].saccadeX[0]+eye_data.arf[0][0].saccadeY[0]
    s_ret_subj = eye_data.trial[0][0].pa[arf_ret==0,1][:,tpts_ret]
    bl_ret = np.nanmean(s_ret_subj[:,1600:1700],axis=1)     
    if len(probe_ind)==len(arf_ret):#calculate mean timecourse during retention interval
        s_ret[isubj] = np.nanmean(s_ret_subj-np.tile(bl_ret,(nt_ret,1)).T,axis=0)
    else:
        print('RET ERROR',ifn)


#behavioral data

attndat = pd.read_csv(results_dir + 'expt1_attn_results.csv')
wmdat = pd.read_csv(results_dir + 'expt1_wm_results.csv')
acc_nonlure = np.array([ np.mean(attndat.acc[np.logical_and(attndat.freq_trials==1,attndat.subject_num==isubj)]==1) for isubj in range(nsubj)])
acc_lure = np.array([ np.mean(attndat.acc[np.logical_and(attndat.freq_trials==0,attndat.subject_num==isubj)]==1) for isubj in range(nsubj)])
wmacc = np.array([ np.mean(wmdat.wm_perf[wmdat.subject_num==isubj]) for isubj in range(nsubj)])


#make figure
fig, ax = plt.subplots(2,2,figsize=(18,10), gridspec_kw={'width_ratios': [1,4]})


#figure 1b accuracy
ax[0,0].plot([0,1],[hit_rate[~np.isnan(hit_rate)]*100,
    (1-fa_rate[~np.isnan(fa_rate)])*100],color='gray',alpha=.1)
scatter_plot_data(ax[0,0],hit_rate[~np.isnan(hit_rate)]*100,x=0)
scatter_plot_data(ax[0,0],(1-fa_rate[~np.isnan(fa_rate)])*100,x=1,c=col_lure)


#fig 1c working memory performance histograms
x = np.ndarray.flatten(wmacc)
y = np.empty((nsubj,7))
b = np.empty((nsubj,8))
y[:] = np.nan
b[:] = np.nan
for isubj in range(nsubj):
    x = wmacc[isubj]
    if not np.all(np.isnan(x)):
        y[isubj],b[isubj] = np.histogram(x[~np.isnan(x)],bins=np.arange(8),density=True)
        ax[1,0].scatter(b[isubj,:-1],y[isubj],s=10,color='gray',alpha=.1,zorder=20,clip_on='False')
        ax[1,0].plot(b[isubj,:-1],y[isubj],color='gray',alpha=.1,zorder=20,clip_on='False')
ax[1,0].plot(np.arange(7),np.nanmean(y,axis=0),color='k',linewidth=3,zorder=20)
ax[1,0].errorbar(np.arange(7),np.nanmean(y,axis=0),yerr=np.nanstd(y,axis=0)/np.sqrt(np.sum(~np.isnan(y[:,0]))),fmt='none',color='k',linewidth=3,capsize=7,capthick=3,zorder=20)


y1 = s_infreq
y2 = s_freq
avg_low_high = (y1+y2)/2.
norm_infreq = y1+(np.nanmean(avg_low_high,axis=0)-avg_low_high)
norm_freq = y2+(np.nanmean(avg_low_high,axis=0)-avg_low_high)


times_stim = np.arange(0,nt_stim)-100
ax[0,1].plot([0,0],[-60,120],'--',color='gray')
ax[0,1].plot([800,800],[-60,120],'--',color='gray')
ax[0,1].plot([1600,1600],[-60,120],'--',color='gray')
#ax[0,1].plot([2400,2400],[-60,120],'--',color='gray')
ax[0,1].plot(times_stim,np.zeros(np.size(times_stim)),'--',color='gray')
ax[0,1].fill_between(times_stim,
          np.nanmean(s_freq,axis=0)-np.nanstd(s_freq,axis=0)/np.sqrt(nsubj),
          np.nanmean(s_freq,axis=0)+np.nanstd(s_freq,axis=0)/np.sqrt(nsubj),
          facecolor='k',alpha=.25)
ax[0,1].fill_between(times_stim,
          np.nanmean(s_infreq,axis=0)-np.nanstd(s_infreq,axis=0)/np.sqrt(nsubj),
          np.nanmean(s_infreq,axis=0)+np.nanstd(s_infreq,axis=0)/np.sqrt(nsubj),
          facecolor=col_lure,alpha=.5)

#plt.fill_between(np.arange(4411),np.mean(np.nanmean(pa_slow,axis=1),axis=0),lw=3,c=[30/255.,188/255.,50/255.])#,c=[68/255.,12/255.,204/255.])
ax[0,1].plot(times_stim,np.nanmean(s_freq,axis=0),lw=3,c='k')#,c=[68/255.,12/255.,204/255.])
ax[0,1].plot(times_stim,np.nanmean(s_infreq,axis=0),lw=3,c=col_lure)#,c=[68/255.,12/255.,204/255.])

y = s_ret


col_highwm = [128/255.,13/255.,86/255.]
col_lowwm = [234/255.,109/255.,188/255.]

n=np.sum(~np.isnan(s_ret[:,0]))

ax[1,1].plot([0,0],[-60,180],'--',color='gray')
ax[1,1].plot([-800,-800],[-60,180],'--',color='gray')
ax[1,1].plot([-1600,-1600],[-60,180],'--',color='gray')
ax[1,1].plot([-2400,-2400],[-60,180],'--',color='gray')
ax[1,1].fill_between([0,2000],[-60,-60],[180,180],facecolor='gray',alpha=.25,edgecolor='None')
ax[1,1].plot(2000+times_ret,np.zeros(np.size(times_ret)),'--',color='gray')
ax[1,1].fill_between(2000+times_ret,
            np.nanmean(s_ret,axis=0)-np.nanstd(s_ret,axis=0)/np.sqrt(n),
            np.nanmean(s_ret,axis=0)+np.nanstd(s_ret,axis=0)/np.sqrt(n),
            facecolor='k',alpha=.25)
ax[1,1].plot(2000+times_ret,np.nanmean(s_ret,axis=0),lw=3,c='k')#,c=[68/255.,12/255.,204/255.])



#spruce things up
prettify_plot(ax[0,0],xlim=(-.5,1.5),ylim=[0,100],
    yl=('Accuracy (%)'),
    yt=([0,25,50,75,100]),ytl=([0,25,50,75,100]),
    xt=([0,1]))
ax[0,0].set_xticklabels(["Non-lure","Lure"], rotation = 60, ha="right",fontsize=28,fontproperties=prop_light)
prettify_plot(ax[1,0],xlim=(-1,7),ylim=[0,.5],
   yt=([0,.1,.2,.3,.4,.5]),ytl=([0,.1,.2,.3,.4,.5]),
   yl='Proportion of trials',xt=([0,1,2,3,4,5,6]),xtl=([0,1,2,3,4,5,6]),
   xl='Number correct')
prettify_plot(ax[0,1],xlim=[-100,2400],
    xt=[0,800,1600,2400],xtl=[0,800,1600,2400],xl='Time relative to stimulus onset (ms)',
    ylim=[-50,100],yt=[-50,0,50,100],ytl=[-50,0,50,100],yl='Pupil size (a.u.)')
prettify_plot(ax[1,1],xlim=[-900,2000],
    xt=[-2400,-1600,-800,0,1000,2000],
    xtl=[-2400,-1600,-800,0,1000,2000],xl='Time relative to retention interval onset (ms)',
    ylim=[-50,150],yt=[-50,0,50,100,150],ytl=[-50,0,50,100,150],yl='Pupil size (a.u.)')

plt.subplots_adjust(wspace=.5,hspace=1.25)
#plt.show(block=False)

fig.savefig(figure_dir + 'figure1.pdf', bbox_inches='tight')

