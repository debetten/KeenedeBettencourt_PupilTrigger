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
#from scipy.stats import ttest_1samp

#supress scientific notation
#np.set_printoptions(suppress=True) 

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

def prettify_plot(ax,
					 xlim=None,xt=None,xtl=None,xl=None,xaxoffset=None,
					 ylim=None,yt=None,ytl=None,yl=None,ylrot=None,yaxoffset=None,
					 t=None,legend=None,legendloc=None):
    '''
    This is a plotting script that makes the default matplotlib plots a little prettier
    '''

    if os.path.isfile("/Library/Fonts/HelveticaNeue-Light.ttf"): 
        prop_light = fm.FontProperties(fname='/Library/Fonts/HelveticaNeue-Light.ttf')    
    else: 
        prop_light = fm.FontProperties()

    if os.path.isfile("/Library/Fonts/HelveticaNeue.ttf"): 
        prop_reg = fm.FontProperties(fname='/Library/Fonts/HelveticaNeue.ttf')    
    else: 
        prop_reg = fm.FontProperties()

    ax.spines['bottom'].set_linewidth(1)
    ax.spines['bottom'].set_color("gray")
    if xaxoffset is not None: ax.spines['bottom'].set_position(('outward', 10))
    if yaxoffset is not None: ax.spines['left'].set_position(('outward', 10))

    ax.spines['left'].set_linewidth(1)
    ax.spines['left'].set_color("gray")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.yaxis.set_ticks_position("left")   
    ax.tick_params(axis='y',direction='out',length=5,width=1,color='gray')
    ax.xaxis.set_ticks_position("bottom") 
    ax.tick_params(axis='x',direction='out',length=5,width=1,color='gray')
    
    if yt is not None: ax.set_yticks(yt)
    if ytl is not None: ax.set_yticklabels((ytl),fontsize=28,fontproperties=prop_light) 
    if yl is not None: h = ax.set_ylabel(yl,fontsize=36,fontproperties=prop_reg,labelpad=12)
    if ylim is not None: ax.set_ylim(ylim) 
        
    if xt is not None: ax.set_xticks(xt)
    if xtl is not None: ax.set_xticklabels((xtl),fontsize=28,fontproperties=prop_light,horizontalalignment='center')
    if xl is not None: ax.set_xlabel(xl,fontsize=36,fontproperties=prop_reg,labelpad=12)
    if xlim is not None: ax.set_xlim(xlim)
    

    if t is not None: ax.set_title(t,y=1.08,fontsize=36,fontproperties=prop_reg)
    if legend is not None: 
        if legendloc is None: L = ax.legend(loc='center left', bbox_to_anchor=(0,.85))
        else: L = ax.legend(loc='center right', bbox_to_anchor=legendloc)
        plt.setp(L.texts,fontsize='large',fontproperties=prop)
    ax.tick_params(axis='both',pad=10)

    for t in ax.get_xticklabels():    #get_xticklabels will get you the label objects, same for y
        t.set_fontsize(28)
    for t in ax.get_yticklabels():
        t.set_fontsize(28)
    ax.yaxis.label.set_size(32)
    ax.xaxis.label.set_size(32)

    plt.locator_params(nbins=8)
    plt.tight_layout()

def scatter_plot_data(ax,data,e=0,x=0,c='k',w=.4):
    n = np.size(data)
    ax.scatter(np.zeros(n)+x,data,facecolor='k',alpha=.1,clip_on=False,edgecolor='None')#data points
    ax.bar(x,np.mean(data),w,color='None',edgecolor=c,linewidth=3)
    ax.errorbar(x,np.mean(data),yerr=np.std(data)/np.sqrt(n),color=c,linewidth=3,capsize=7,capthick=3)#error bar

def resampling_statistics(data,chance,nsubj=None,nsamples=100000,skipfig=True):
    '''
    Nonparametric resampling statistics
    '''
    if nsubj is None: 
        nsubj = np.shape(data)[0]
    nt=np.shape(data)[1]
    #resample subjects with replacement 100,000 times
    subj_resampled = np.random.randint(0,nsubj,(nsamples,nsubj)) 

    #the empy matrix of resampled data for each of these resampled iterations
    data_resampled = np.empty((nsamples,nt)) 

    #recalculate mean given the resampled subjects
    for i in range(0,nsamples): 
        data_resampled[i] = np.mean(data[subj_resampled[i]],axis=0)

    #calculate p value 
    p = np.sum(np.less(data_resampled,chance),axis=0)/float(nsamples) #count number of resample iterations below chance
    p[np.equal(p,0)]= 1./float(nsamples)

    if skipfig is False:
        plt.figure(figsize=(4,3))
        ax = plt.subplot(111)
        plt.hist(data_resampled,normed=0,facecolor='gray',edgecolor='gray')
        plt.axvline(np.mean(data_resampled),color='b',lw=2,label='resampled mean')
        plt.axvline(np.mean(data),color='m',lw=1,label='original mean')
        plt.axvline(chance,color='c',lw=2,label='chance')
        make_plot_pretty(ax,ylrot=90,yl='Count (#)',legend='1') 
        plt.show()
    
    return p

############## Driver #################

#project directory 
if platform.system()=='Darwin':
    base_project_dir='/Users/megan/Dropbox/writing/articles/2020_pupil/'
else:
    base_project_dir='/Users/Paul Keene/documents/github/pupilometry_paper/'
project_name = 'eyeSustAttnWM03'
project_dir = base_project_dir + 'expt2/'
beh_dir = project_dir + 'beh/'
eye_dir = project_dir + 'eye/'
figure_dir = base_project_dir + 'figures/'
results_dir = project_dir + 'results/'

#index subject names
beh_fn = [f for f in os.listdir(beh_dir ) if f.endswith(project_name + '_explog.csv')]
eye_fn_stim = [f for f in os.listdir(eye_dir ) if f.endswith(project_name + '_ENCARRAY_EYE_SEG.mat')]
eye_fn_ret = [f for f in os.listdir(eye_dir ) if f.endswith(project_name + '_EYE_SEG.mat')]
nsubj = np.size(beh_fn)

#initialize arrays
attndat = pd.read_csv(results_dir + 'expt2_attn_results.csv')
wmdat = pd.read_csv(results_dir + 'expt2_wm_results.csv')
acc_nonlure = np.array([ np.mean(attndat.acc[np.logical_and(attndat.freq_trials==1,attndat.subject_num==isubj)]==1) for isubj in range(nsubj)])
acc_lure = np.array([ np.mean(attndat.acc[np.logical_and(attndat.freq_trials==0,attndat.subject_num==isubj)]==1) for isubj in range(nsubj)])
wmacc = np.array([ np.mean(wmdat.wm_perf[wmdat.subject_num==isubj]) for isubj in range(nsubj)])

nt = 2901
s_freq = np.zeros((nsubj,2501))
s_infreq = np.zeros((nsubj,2501))
s_ret = np.zeros((nsubj,2901))
for isubj,ifn in enumerate(beh_fn):
    print(isubj,ifn,eye_fn_ret[isubj])

    #load behavioral data from CSV 
    beh_dat = pd.read_csv(beh_dir+ifn ,header=12) 
    itrials = np.arange(np.size(beh_dat.acc))

    #load eye-tracking data
    mat_contents = sio.loadmat(eye_dir+eye_fn_stim[isubj],struct_as_record=False)
    eye_data_stim = mat_contents['eyeData'][0][0]
    s_post = eye_data_stim.trial[0][0].pa[:,0,:2501]
    s_post_bl = np.nanmean(eye_data_stim.trial[0][0].pa[:,0,:100],axis=1)
    arf_stim = ((eye_data_stim.arf[0][0].parserBlinks[0]+eye_data_stim.arf[0][0].missingPupil[0]+eye_data_stim.arf[0][0].saccadeX[0]+eye_data_stim.arf[0][0].saccadeY[0])>0)*1

    #delete trials with missing eye data
    if ifn == '0815192_eyeSustAttnWM03_explog.csv': #eyetracker and beh indecies dont match. found discrepency by hand
        j=[646,645,406]
    else:
        j=[]    
    if np.size(j)>0:
        itrials = np.delete(itrials,j) 

    s_freq[isubj] = np.nanmean(s_post[np.logical_and(beh_dat.freq_trials[itrials]==1,arf_stim==0)]-
                    np.tile(s_post_bl[np.logical_and(beh_dat.freq_trials[itrials]==1,arf_stim==0)],(2501,1)).T,axis=0)
    s_infreq[isubj] = np.nanmean(s_post[np.logical_and(beh_dat.freq_trials[itrials]==0,arf_stim==0)]-
                    np.tile(s_post_bl[np.logical_and(beh_dat.freq_trials[itrials]==0,arf_stim==0)],(2501,1)).T,axis=0) 

    #
    mat_contents = sio.loadmat(eye_dir+eye_fn_ret[isubj],struct_as_record=False)
    eye_data_ret = eye_data = mat_contents['eyeData'][0][0]
    s_ret_raw = eye_data.trial[0][0].pa[:,0,2400:5301]
    s_ret_bl = np.nanmean(s_ret_raw[:,:100],axis=1)
    arf_ret = ((eye_data_ret.arf[0][0].parserBlinks[0]+eye_data_ret.arf[0][0].missingPupil[0]+eye_data_ret.arf[0][0].saccadeX[0]+eye_data_ret.arf[0][0].saccadeY[0])>0)*1
    s_ret[isubj] = np.mean(s_ret_raw[arf_ret==0]-np.tile(s_ret_bl[arf_ret==0],(2901,1)).T,axis=0)




#make figure
fig, ax = plt.subplots(2,2,figsize=(18,10), gridspec_kw={'width_ratios': [1,4]})


#figure 1b accuracy
ax[0,0].plot([0,1],[acc_nonlure*100,acc_lure*100],color='gray',alpha=.1)
scatter_plot_data(ax[0,0],acc_nonlure*100,x=0)
scatter_plot_data(ax[0,0],acc_lure*100,x=1,c=col_lure)


#fig 1c working memory performance histograms
x = np.ndarray.flatten(wmacc)
y = np.empty((nsubj,7))
b = np.empty((nsubj,8))
y[:] = np.nan
b[:] = np.nan
for isubj in range(nsubj):
    x = wmdat.wm_perf[wmdat.subject_num==isubj]
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


times_stim = np.arange(0,2501)-100
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

times_ret = np.arange(2901)-900
ax[1,1].plot([0,0],[-60,180],'--',color='gray')
ax[1,1].plot([-800,-800],[-60,180],'--',color='gray')
ax[1,1].fill_between([0,2000],[-60,-60],[180,180],facecolor='gray',alpha=.25,edgecolor='None')
ax[1,1].plot(times_ret,np.zeros(np.size(times_ret)),'--',color='gray')
ax[1,1].fill_between(times_ret,
            np.nanmean(s_ret,axis=0)-np.nanstd(s_ret,axis=0)/np.sqrt(n),
            np.nanmean(s_ret,axis=0)+np.nanstd(s_ret,axis=0)/np.sqrt(n),
            facecolor='k',alpha=.25)
ax[1,1].plot(times_ret,np.nanmean(s_ret,axis=0),lw=3,c='k')#,c=[68/255.,12/255.,204/255.])



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
    xt=[-800,0,1000,2000],
    xtl=[-800,0,1000,2000],xl='Time relative to retention interval onset (ms)',
    ylim=[-50,150],yt=[-50,0,50,100,150],ytl=[-50,0,50,100,150],yl='Pupil size (a.u.)')

plt.subplots_adjust(wspace=.5,hspace=1.25)
#plt.show(block=False)

fig.savefig(figure_dir + 'figure3_expt2_v1.pdf', bbox_inches='tight')

