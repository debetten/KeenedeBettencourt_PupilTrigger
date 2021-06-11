import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, csv
import matplotlib.font_manager as fm
from matplotlib import rc
import scipy.io as sio
import platform 
from scipy.stats import ttest_rel, spearmanr

#supress scientific notation
np.set_printoptions(suppress=True) 

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
col_corr = [0/255.,98/255.,100/255.]
col_incorr = [218/255.,66/255.,36/255.]

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

def scatter_plot_data(ax,data,e=0,x=0,c='k',w=.8):
    n = np.size(data)
    ax.scatter(np.zeros(n)+x,data,facecolor='k',edgecolor='None',alpha=.1,clip_on=False)#data points
    ax.bar(x,np.nanmean(data),width=w,color='None',edgecolor=c,linewidth=3)
    ax.errorbar(x,np.nanmean(data),yerr=np.nanstd(data)/np.sqrt(n),color=c,linewidth=3,capsize=10,capthick=3)#error bar

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
    p[np.equal(p,0)] = 1./float(nsamples)

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
project_name = 'eyeSustAttnWM01'
project_dir = base_project_dir + 'expt1/'
beh_dir = project_dir + 'beh/'
eye_dir = project_dir + 'eye/'
eye_dir_old = base_project_dir + 'expt1_old/' 'eye/'
figure_dir = base_project_dir + 'figures/'

#index subject names
beh_fn = [f for f in os.listdir(beh_dir ) if f.endswith(project_name + '_explog.csv')]
eye_fn_stim = [f for f in os.listdir(eye_dir ) if f.endswith(project_name + '_ENCARRAY_EYE_SEG.mat')]
eye_fn_ret = [f for f in os.listdir(eye_dir ) if f.endswith(project_name + '_EYE_SEG.mat')]
beh_dat = dict.fromkeys(beh_fn)

#preallocate
nsubj = len(beh_fn) #number of subjects
nt_stim = 2601 #number of timepoints
tpts_ret = np.arange(0,4500)
nt_ret = np.size(tpts_ret)

s_infreq = np.zeros((nsubj,nt_stim))
s_infreqinacc = np.zeros((nsubj,nt_stim))
s_infreqacc = np.zeros((nsubj,nt_stim))
s_freq = np.zeros((nsubj,nt_stim))
s_ret_wmacc_high = np.zeros((nsubj,nt_ret))
s_ret_wmacc_low = np.zeros((nsubj,nt_ret))
s_ret_wmacc_high[:] = np.nan
s_ret_wmacc_low[:] = np.nan
r_leftright = np.zeros((nsubj))
## Calculate average pupil timecourse based on sustained attention trial type
for isubj,ifn in enumerate(beh_fn):
    print(isubj,ifn,eye_fn_stim[isubj])

    #load behavioral data
    if ifn=='0905181_eyeSustAttnWM01_explog.csv' or ifn=='0917181_eyeSustAttnWM01_explog.csv':
        beh_dat = pd.read_csv(beh_dir+ifn ,header=11)
    else:
        beh_dat = pd.read_csv(beh_dir+ifn ,header=12)

    #sustained attention behavior
    freq_bool = np.ravel(beh_dat.freq_trials==1)
    infreq_bool = np.ravel(beh_dat.freq_trials==0)
    end_block_bool = np.ravel(beh_dat.trial<798) #for choosing eye data for freq trials

    acc = np.ravel(beh_dat.acc)
    #memory probe behavior
    probe_ind = np.where(beh_dat.probe_trials==1)[0] #probes

    #missing data
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

    wmacc = np.ravel(beh_dat.wholereportrespacc[probe_ind])
    wmacc_high = wmacc>=np.nanmedian(wmacc)
    wmacc_low = wmacc<np.nanmedian(wmacc)

    mat_contents = sio.loadmat(eye_dir+eye_fn_stim[isubj],struct_as_record=False)
    eye_data = mat_contents['eyeData'][0][0]
    times_stim = eye_data.trial[0][0].times[0]
    #arf_stim = eye_data.arf[0][0].combined[0] #find artifacts
    arf_stim = eye_data.arf[0][0].parserBlinks[0]+eye_data.arf[0][0].missingPupil[0]+eye_data.arf[0][0].saccadeX[0]+eye_data.arf[0][0].saccadeY[0]
    s_stim = eye_data.trial[0][0].pa[:,0,:]
    x0=np.nanmean(eye_data.trial[0][0].pa[arf_stim==0,0,:],axis=0)
    x1=np.nanmean(eye_data.trial[0][0].pa[arf_stim==0,1,:],axis=0)
    r_leftright[isubj],_ = spearmanr(x0,x1)
    bl = np.nanmean(s_stim[:,0:100],axis=1)

    mat_contents = sio.loadmat(eye_dir+eye_fn_ret[isubj],struct_as_record=False)
    eye_data = mat_contents['eyeData'][0][0]
    times_ret = eye_data.trial[0][0].times[0,tpts_ret]
    #arf_ret = mat_contents['eyeData'][0][0].arf[0][0].combined[0]
    arf_ret = eye_data.arf[0][0].parserBlinks[0]+eye_data.arf[0][0].missingPupil[0]+eye_data.arf[0][0].saccadeX[0]+eye_data.arf[0][0].saccadeY[0]
    s_ret = eye_data.trial[0][0].pa[:,1,tpts_ret]
    bl_ret = np.nanmean(s_ret[:,0:100],axis=1)    

    #calculate mean timecourse for each stimulus array
    if len(freq_bool)==len(arf_stim) and np.shape(s_stim)[1]==2601:
        i = np.logical_and(np.logical_and(infreq_bool==1,arf_stim==0),end_block_bool==1)
        s_infreq[isubj] = np.nanmean(s_stim[i]-np.tile(bl[i],(nt_stim,1)).T,axis=0)

        i = np.logical_and(np.logical_and(np.logical_and(infreq_bool==1,acc==1),arf_stim==0),end_block_bool==1)
        s_infreqacc[isubj] = np.nanmean(s_stim[i]-np.tile(bl[i],(nt_stim,1)).T,axis=0)

        i = np.logical_and(np.logical_and(np.logical_and(infreq_bool==1,acc==0),arf_stim==0),end_block_bool==1)
        s_infreqinacc[isubj] = np.nanmean(s_stim[i]-np.tile(bl[i],(nt_stim,1)).T,axis=0)

        i = np.logical_and(np.logical_and(np.logical_and(infreq_bool==0,acc==1),arf_stim==0),end_block_bool==1)
        s_freq[isubj] = np.nanmean(s_stim[i]-np.tile(bl[i],(nt_stim,1)).T,axis=0)
    else:
        print('STIM ERROR'+ifn,np.shape(s_stim)[1])
        s_infreqacc[isubj] = np.nan
        s_infreqinacc[isubj] = np.nan
        s_freq[isubj] = np.nan
    
    #calculate mean timecourse during retention interval
    if len(probe_ind)==len(arf_ret):
        s_ret = s_ret-np.tile(bl_ret,(nt_ret,1)).T

        i = np.logical_and(wmacc_high,arf_ret==0)
        s_ret_wmacc_high[isubj] = np.nanmean(s_ret[i],axis=0)

        i = np.logical_and(wmacc_low,arf_ret==0)
        s_ret_wmacc_low[isubj] = np.nanmean(s_ret[i],axis=0)
    else:
        print('RET ERROR',ifn)


#plot figure
fig,ax = plt.subplots(2,2,figsize=(16,12),gridspec_kw={'width_ratios': [4,1]})

y1 = s_infreqacc
y2 = s_infreqinacc
y3 = s_freq
avg_low_high = (y1+y2+y3)/3.
norm_infreqacc = y1+(np.nanmean(avg_low_high,axis=0)-avg_low_high)
norm_infreqinacc = y2+(np.nanmean(avg_low_high,axis=0)-avg_low_high)
norm_freq = y3+(np.nanmean(avg_low_high,axis=0)-avg_low_high)


times_stim = np.arange(0,nt_stim)-100
ax[0,0].plot([0,0],[-60,120],'--',color='gray')
ax[0,0].plot([800,800],[-60,120],'--',color='gray')
ax[0,0].plot([1600,1600],[-60,120],'--',color='gray')
ax[0,0].plot([2400,2400],[-60,120],'--',color='gray')
ax[0,0].plot(times_stim,np.zeros(np.size(times_stim)),'--',color='gray')
ax[0,0].fill_between(times_stim,
          np.nanmean(s_freq,axis=0)-np.nanstd(norm_freq,axis=0)/np.sqrt(nsubj),
          np.nanmean(s_freq,axis=0)+np.nanstd(norm_freq,axis=0)/np.sqrt(nsubj),
          facecolor='gray',alpha=.5)
ax[0,0].fill_between(times_stim,
          np.nanmean(s_infreqinacc,axis=0)-np.nanstd(norm_infreqinacc,axis=0)/np.sqrt(nsubj),
          np.nanmean(s_infreqinacc,axis=0)+np.nanstd(norm_infreqinacc,axis=0)/np.sqrt(nsubj),
          facecolor=[30/255.,188/255.,50/255.],alpha=.5)
ax[0,0].fill_between(times_stim,
           np.nanmean(s_infreqacc,axis=0)-np.nanstd(norm_infreqacc,axis=0)/np.sqrt(nsubj),
           np.nanmean(s_infreqacc,axis=0)+np.nanstd(norm_infreqacc,axis=0)/np.sqrt(nsubj),
           facecolor=[30/255.,50/255.,188/255.],alpha=.5)

#plt.fill_between(np.arange(4411),np.mean(np.nanmean(pa_slow,axis=1),axis=0),lw=3,c=[30/255.,188/255.,50/255.])#,c=[68/255.,12/255.,204/255.])
ax[0,0].plot(times_stim,np.nanmean(s_freq,axis=0),lw=3,c='gray')#,c=[68/255.,12/255.,204/255.])
ax[0,0].plot(times_stim,np.nanmean(s_infreqinacc,axis=0),lw=3,c=[30/255.,188/255.,50/255.])#,c=[68/255.,12/255.,204/255.])
ax[0,0].plot(times_stim,np.nanmean(s_infreqacc,axis=0),lw=3,c=[30/255.,50/255.,188/255.])#,c=[202/255.,52/255.,55/255.])

thr=.005
#tval_stim = np.empty(np.size(times_stim))
#pval_stim = np.empty(np.size(times_stim))
#for i in range(np.size(times_stim)):
#  tval_stim[i],pval_stim[i] = ttest_rel(s_infreqacc[:,i],s_infreqinacc[:,i],nan_policy='omit')
#pval_stim = resampling_statistics(s_infreqinacc-s_infreqacc,0,nsamples=10000)
ax[0,0].scatter((times_stim)[pval_stim<thr],-31+np.ones(np.sum(pval_stim<thr)),marker='s',color='k',zorder=20)
ax[0,0].scatter((times_stim)[pval_stim>(1-thr)],-31+np.ones(np.sum(pval_stim>(1-thr))),marker='s',color='k',zorder=20)

prettify_plot(ax[0,0],xlim=[-100,2500],
    xt=[0,800,1600,2400],xtl=[0,800,1600,2400],xl='Time relative to stimulus onset (ms)',
    ylim=[-60,120],yt=[-60,0,60,120],ytl=[-60,0,60,120],yl='Pupil size')
#ax.legend(["infreq correct","infreq incorrect","frequent"])

scatter_plot_data(ax[0,1],np.nanmean(s_freq[:,(times_stim)>0],axis=1),
    x=0,c='gray')
scatter_plot_data(ax[0,1],np.nanmean(s_infreqacc[:,(times_stim)>0],axis=1),
    x=1,c=[30/255.,50/255.,188/255.])
scatter_plot_data(ax[0,1],np.nanmean(s_infreqinacc[:,(times_stim)>0],axis=1),
    x=2,c=[30/255.,188/255.,50/255.])
prettify_plot(ax[0,1],xlim=[-.5,2.5],
    xt=[0,1,2],
    xtl=['Frequent','Infrequent\naccurate','Infrequent\ninaccurate'],xl='Working memory\nperformance',
    ylim=[-60,120],yt=[-60,0,60,120],ytl=[-60,0,60,120],yl='Pupil size')
ax[0,1].set_xticklabels(['Frequent','Infreq. acc.','Infreq. inacc.'], rotation=90,fontproperties=prop_light)

y1 = s_ret_wmacc_high
y2 = s_ret_wmacc_low
diff_low_high = y1-y2
avg_low_high = (y1+y2)/2.
norm_low = y1+(np.nanmean(avg_low_high,axis=0)-avg_low_high)
norm_high = y2+(np.nanmean(avg_low_high,axis=0)-avg_low_high)

col_highwm = [128/255.,13/255.,86/255.]
col_lowwm = [234/255.,109/255.,188/255.]

n=np.sum(~np.isnan(s_ret_wmacc_low[:,0]))

ax[1,0].plot([0,0],[-60,180],'--',color='gray')
ax[1,0].plot([-800,-800],[-60,180],'--',color='gray')
ax[1,0].plot([-1600,-1600],[-60,180],'--',color='gray')
ax[1,0].plot([-2400,-2400],[-60,180],'--',color='gray')
ax[1,0].fill_between([0,2000],[-60,-60],[180,180],facecolor='gray',alpha=.25,edgecolor='None')
ax[1,0].plot(2000+times_ret,np.zeros(np.size(times_ret)),'--',color='gray')
ax[1,0].fill_between(2000+times_ret,
            np.nanmean(s_ret_wmacc_low,axis=0)-np.nanstd(norm_low,axis=0)/np.sqrt(n),
            np.nanmean(s_ret_wmacc_low,axis=0)+np.nanstd(norm_low,axis=0)/np.sqrt(n),
            facecolor=col_lowwm,alpha=.5)
ax[1,0].fill_between(2000+times_ret,
           np.nanmean(s_ret_wmacc_high,axis=0)-np.nanstd(norm_high,axis=0)/np.sqrt(n),
           np.nanmean(s_ret_wmacc_high,axis=0)+np.nanstd(norm_high,axis=0)/np.sqrt(n),
           facecolor=col_highwm,alpha=.5)
#plt.fill_between(np.arange(4411),np.mean(np.nanmean(pa_slow,axis=1),axis=0),lw=3,c=[30/255.,188/255.,50/255.])#,c=[68/255.,12/255.,204/255.])
ax[1,0].plot(2000+times_ret,np.nanmean(s_ret_wmacc_low,axis=0),lw=3,c=col_lowwm)#,c=[68/255.,12/255.,204/255.])
ax[1,0].plot(2000+times_ret,np.nanmean(s_ret_wmacc_high,axis=0),lw=3,c=col_highwm)#,c=[202/255.,52/255.,55/255.])

# ax[1].plot([-2400,-2400],[-30,150],color='gray',lw=.5,zorder=0)
# ax[1].plot([-1600,-1600],[-30,150],color='gray',lw=.5,zorder=0)
# ax[1].plot([-800,-800],[-30,150],color='gray',lw=.5,zorder=0)
# ax[1].plot([0,0],[-30,150],color='gray',lw=.5,zorder=0)


#tval_ret = np.empty(np.size(times_ret))
#pval_ret = np.empty(np.size(times_ret))
#for i in range(np.size(times_ret)):
#  tval_ret[i],pval_ret[i] = ttest_rel(s_ret_wmacc_high[:,i],s_ret_wmacc_low[:,i],nan_policy='omit')
#pval_ret = resampling_statistics(s_ret_wmacc_high-s_ret_wmacc_low,0,nsamples=10000)

ax[1,0].scatter((2000+times_ret)[pval_ret<thr],-21+np.ones(np.sum(pval_ret<thr)),marker='s',color='k')
ax[1,0].scatter((2000+times_ret)[pval_ret>(1-thr)],-21+np.ones(np.sum(pval_ret>(1-thr))),marker='s',color='k')

prettify_plot(ax[1,0],xlim=[-2500,2000],
    xt=[-2400,-1600,-800,0,1000,2000],
    xtl=[-2400,-1600,-800,0,1000,2000],xl='Time relative to retention interval onset (ms)',
    ylim=[-60,180],yt=[-60,0,60,120,180],ytl=[-60,0,60,120,180],yl='Pupil size')

scatter_plot_data(ax[1,1],np.nanmean(s_ret_wmacc_low[:,(2000+times_ret)>1000],axis=1),
    x=0,c=col_lowwm)
scatter_plot_data(ax[1,1],np.nanmean(s_ret_wmacc_high[:,(2000+times_ret)>1000],axis=1),
    x=1,c=col_highwm)
prettify_plot(ax[1,1],xlim=[-.5,1.5],
    xt=[0,1],
    xtl=['Low','High'],xl='Working memory\nperformance',
    ylim=[-100,400],yt=[-100,0,100,200,300,400],ytl=[-100,0,100,200,300,400],yl='Pupil size')


fig.savefig(figure_dir + 'figure2.pdf', bbox_inches='tight')
plt.show(block=False)
