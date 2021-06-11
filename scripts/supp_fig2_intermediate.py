import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import csv
import matplotlib.font_manager as fm
from matplotlib import rc
import scipy.io as sio
import glob
import scipy
from scipy.stats import ttest_1samp
from scipy.stats import f_oneway


#
plt.rcParams.update({'font.size':200})
font = {'size'   : 32}
rc('font', **font)

#numpy.intersect


#supress scientific notation
np.set_printoptions(suppress=True) 

#font defaults
plt.rcParams.update({'font.size': 32})
rc('text', usetex=False)
plt.rcParams['pdf.fonttype'] = 42
if os.path.isfile("/Library/Fonts/HelveticaNeue.ttf"): 
    prop = fm.FontProperties(fname="/Library/Fonts/HelveticaNeue.ttf",size=24)
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
    if ytl is not None: ax.set_yticklabels((ytl),fontsize=32,fontproperties=prop_light) 
    if yl is not None: h = ax.set_ylabel(yl,fontsize=36,fontproperties=prop_reg,labelpad=12)
    if ylim is not None: ax.set_ylim(ylim) 
        
    if xt is not None: ax.set_xticks(xt)
    if xtl is not None: ax.set_xticklabels((xtl),fontsize=12,fontproperties=prop_light,horizontalalignment='center')
    if xl is not None: ax.set_xlabel(xl,fontsize=36,fontproperties=prop_reg,labelpad=12)
    if xlim is not None: ax.set_xlim(xlim)
    

    if t is not None: ax.set_title(t,y=1.08,fontsize=36,fontproperties=prop_reg)
    if legend is not None: 
        if legendloc is None: L = ax.legend(loc='center left', bbox_to_anchor=(0,.85))
        else: L = ax.legend(loc='center right', bbox_to_anchor=legendloc)
        plt.setp(L.texts,fontsize='large',fontproperties=prop)
    ax.tick_params(axis='both',pad=10)
    plt.locator_params(nbins=8)
    plt.tight_layout()

def scatter_plot_data(ax,data,e=0,x=0,c='k',w=.8):
    n = np.size(data)
    ax.scatter(np.zeros(n)+x,data,facecolor='k',alpha=.1,clip_on=False)#data points
    ax.bar(x,np.nanmean(data),width=w,color='None',edgecolor=c,linewidth=3)
    ax.errorbar(x,np.nanmean(data),yerr=np.nanstd(data)/np.sqrt(n),color=c,linewidth=3,capsize=10,capthick=3)#error bar

@np.vectorize
def calculate_aprime(h,fa):
    '''
    Calculates sensitivity, according to the aprime formula
    '''
    if np.greater_equal(h,fa): a = .5 + (((h-fa) * (1+h-fa)) / (4 * h * (1-fa)))
    else: a = .5 - (((fa-h) * (1+fa-h)) / (4 * fa * (1-h)))
    return a

############## Driver #################

#project directory 
project_name = 'eyeSustAttnWM03'

#where is this project saved on the computer (note you will need to change this if running on a different computer)
base_project_dir='/Users/pkeen/Documents/github/pupilometry_paper/'
project_dir = base_project_dir + 'expt2/'
subjects_dir = project_dir + 'subj_data/'
figure_dir = base_project_dir + 'paper_figs/'

#index subject names
subj_name = [f for f in os.listdir(subjects_dir) if f.endswith(project_name)]
#create dictionary for data
dat = dict.fromkeys(subj_name)
#initialize arrays
nsubj = len(subj_name)
rts_high = np.zeros((nsubj))
rts_low = np.zeros((nsubj))
wmacc_high = np.zeros((nsubj))
wmacc_low = np.zeros((nsubj))
pupil_high = np.zeros((nsubj))
pupil_low = np.zeros((nsubj))
wmacc_midhigh = np.zeros((nsubj))
pupil_midhigh = np.zeros((nsubj))
rts_midhigh = np.zeros((nsubj))
wmacc_midlow = np.zeros((nsubj))
pupil_midlow = np.zeros((nsubj))
rts_midlow = np.zeros((nsubj))
r = np.zeros((nsubj))
sum_high = 0
sum_low = 0

for isubj in range(len(subj_name)):
    #load data
    subj = subj_name[isubj] #string of subjects name
    dat[subj] = pd.read_csv(subjects_dir+subj+'/data/beh/'+subj+'_explog.csv',header=12)
    #index frequent trials
    high_triggered_ind = np.where(dat[subj].pupil_triggered_high==1)[0]+1
    low_triggered_ind = np.where(dat[subj].pupil_triggered_low==1)[0]+1
    triggered_probes = np.concatenate((high_triggered_ind,low_triggered_ind))
    trailing_rts = dat[subj].pupil_trailingavg
    trailing_rts_triggered = trailing_rts[triggered_probes]
    #index probes
    probe_ind = np.where(dat[subj].probe_trials==1)[0]

    #split probes to find intermediates
    [upperlowP,lowerupP,upperupP] = np.percentile(trailing_rts[np.isnan(trailing_rts)==0],[25,50,75])
    lowerRTS = trailing_rts_triggered<upperlowP
    lowerProbes = triggered_probes[lowerRTS]
    lowermidRTS = np.logical_and(upperlowP<trailing_rts_triggered,trailing_rts_triggered<lowerupP)
    lowermidProbes = triggered_probes[lowermidRTS]
    uppermidRTS = np.logical_and(lowerupP<trailing_rts_triggered,trailing_rts_triggered<upperupP)
    uppermidProbes = triggered_probes[uppermidRTS]
    upperRTS = trailing_rts_triggered>upperupP
    upperProbes = triggered_probes[upperRTS]
    #r[isubj] = scipy.stats.spearmanr(trailing_rts_triggered,dat[subj].wholereportrespacc[probe_ind])

    # mat_contents = sio.loadmat(subjects_dir+subj+'/data/eye/'+subj+'_EYE_SEG.mat',struct_as_record=False) ##fix filename
    # #print(mat_contents['eyeData'][0][0])
    # pupil_raw = mat_contents['eyeData'][0][0].trial[0][0].pa 
    # arf = mat_contents['eyeData'][0][0].arf[0][0].combined[0]
    # pupil_size = np.zeros((np.shape(pupil_raw)))

    pupil_high[isubj] = np.nanmean(dat[subj].pupil_trailingavg[upperProbes])
    pupil_low[isubj] = np.nanmean(dat[subj].pupil_trailingavg[lowerProbes])
    pupil_midhigh[isubj] = np.nanmean(dat[subj].pupil_trailingavg[uppermidProbes])
    pupil_midlow[isubj] = np.nanmean(dat[subj].pupil_trailingavg[lowermidProbes])
    wmacc_high[isubj] = np.nanmean(dat[subj].wholereportrespacc[upperProbes])
    wmacc_low[isubj] = np.nanmean(dat[subj].wholereportrespacc[lowerProbes])
    wmacc_midhigh[isubj] = np.nanmean(dat[subj].wholereportrespacc[uppermidProbes])
    wmacc_midlow[isubj] = np.nanmean(dat[subj].wholereportrespacc[lowermidProbes])
    tmp_rts_high = np.zeros(len(upperProbes))
    tmp_rts_low = np.zeros(len(lowerProbes))
    tmp_rts_midhigh = np.zeros(len(uppermidProbes))
    tmp_rts_midlow = np.zeros(len(lowermidProbes))
    for itrial in range(len(upperProbes)):
        tmp_rts_high[itrial] = np.nanmean(dat[subj].rts[(upperProbes[itrial]-3):upperProbes[itrial]])
    for itrial in range(len(lowerProbes)):
        tmp_rts_low[itrial] = np.nanmean(dat[subj].rts[(lowerProbes[itrial]-3):lowerProbes[itrial]])
    for itrial in range(len(lowermidProbes)):
        tmp_rts_midlow[itrial] = np.nanmean(dat[subj].rts[(lowermidProbes[itrial]-3):lowermidProbes[itrial]])
    for itrial in range(len(uppermidProbes)):
        tmp_rts_midhigh[itrial] = np.nanmean(dat[subj].rts[(uppermidProbes[itrial]-3):uppermidProbes[itrial]])
    rts_high[isubj] = np.nanmean(tmp_rts_high)
    rts_low[isubj] = np.nanmean(tmp_rts_low)
    rts_midhigh[isubj] = np.nanmean(tmp_rts_midhigh)
    rts_midlow[isubj] = np.nanmean(tmp_rts_midlow)

#make figure

fig, ax = plt.subplots(1,3,figsize=(24,12))
#fig 1b - attention (A')
scatter_plot_data(ax[1],rts_high,x=0)
scatter_plot_data(ax[1],rts_low,x=3)
scatter_plot_data(ax[1],rts_midhigh,x=1)
scatter_plot_data(ax[1],rts_midlow,x=2)
scatter_plot_data(ax[2],wmacc_high,x=0)
scatter_plot_data(ax[2],wmacc_low,x=3)
scatter_plot_data(ax[2],wmacc_midhigh,x=1)
scatter_plot_data(ax[2],wmacc_midlow,x=2)
scatter_plot_data(ax[0],pupil_high,x=0)
scatter_plot_data(ax[0],pupil_low,x=3)
scatter_plot_data(ax[0],pupil_midhigh,x=1)
scatter_plot_data(ax[0],pupil_midlow,x=2)

# print("high-low ",str(ttest_1samp(wmacc_high-wmacc_low,0))," high-mid ",str(ttest_1samp(wmacc_high-wmacc_intermediate,0))," low-mid ",str(ttest_1samp(wmacc_low-wmacc_intermediate,0)))
# print("high-low ",str(ttest_1samp(rts_high-rts_low,0))," high-mid ",str(ttest_1samp(rts_high-rts_intermediate,0))," low-mid ",str(ttest_1samp(rts_low-rts_intermediate,0)))
# print("wm anova = ",str(f_oneway(wmacc_high,wmacc_midhigh,wmacc_midlow,wmacc_low)))
# print("rt anova = ",str(f_oneway(rts_high,rts_midhigh,rts_midlow,rts_low)))
print("WM: large-midlarge ",str(ttest_1samp(wmacc_high-wmacc_midhigh,0))," small-midsmall ",str(ttest_1samp(wmacc_low-wmacc_midlow,0)))
print("RT: large-midlarge ",str(ttest_1samp(rts_high-rts_midhigh,0))," small-midsmall ",str(ttest_1samp(rts_low-rts_midlow,0)))

#print(np.nanmean(wmacc))

#spruce things up
prettify_plot(ax[0],xlim=(-.5,3.5),ylim=[0,3000],
    yt=([0,500,1000,1500,2000,2500,3000]),ytl=([0,500,1000,1500,2000,2500,3000]),yl="Pretrial Pupil Area (a.u)",
    xt=([0,1,2,3]),xtl=(["Large PA","Mid-Large PA","Mid-Small","Small PA"]),xl="Triggering Condition",t='Pupil Area difference')
prettify_plot(ax[1],xlim=(-.5,3.5),ylim=[0,.8],
    yt=([0,.2,.4,.6,.8]),ytl=([0,200,400,600,800]),yl="Pretrial Response Time (ms)",
    xt=([0,1,2,3]),xtl=(["Large PA","Mid-Large PA","Mid-Small","Small PA"]),xl="Triggering Condition",t='RTs difference')
prettify_plot(ax[2],xlim=(-.5,3.5),ylim=[0,4],
    yt=([0,1,2,3,4,5,6]),ytl=([0,1,2,3,4,5,6]),yl="WM accuracy",
    xt=([0,1,2,3]),xtl=(["Large PA","Mid-Large PA","Mid-Small","Small PA"]),xl="Triggering Condition",t='WM accuaracy')


fig.savefig(figure_dir + 'suppfig2_intermdiate.pdf', bbox_inches='tight')