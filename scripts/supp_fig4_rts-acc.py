import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import os
import csv
import matplotlib.font_manager as fm
from matplotlib import rc
import scipy.io as sio
import glob
import scipy
from scipy.stats import ttest_1samp
#from sklearn.linear_model import LogisticRegression


#
plt.rcParams.update({'font.size':200})
font = {'size'   : 32}
rc('font', **font)


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
    if xtl is not None: ax.set_xticklabels((xtl),fontsize=32,fontproperties=prop_light,horizontalalignment='center')
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
project_name = 'eyeSustAttnWM01'

#where is this project saved on the computer (note you will need to change this if running on a different computer)
base_project_dir='/Users/pkeen/Documents/github/pupilometry_paper/'
project_dir = base_project_dir + 'expt1/'
subjects_dir = project_dir + 'subj_data/'
figure_dir = base_project_dir + 'paper_figs/'

#index subject names
subj_name = [f for f in os.listdir(subjects_dir) if f.endswith(project_name)]
#create dictionary for data
dat = dict.fromkeys(subj_name)
#initialize arrays
nsubj = len(subj_name)
prop_correct = np.zeros((nsubj,10))
p = np.zeros((nsubj,2))
for isubj in range(len(subj_name)):
    #load data
    subj = subj_name[isubj] #string of subjects name
    dat[subj] = pd.read_csv(subjects_dir+subj+'/data/beh/'+subj+'_explog.csv',header=12)
    #index frequent trials
    infreq_trials = np.where(dat[subj].freq_trials==0)[0]
    rts = np.array(dat[subj].rts)
    trailing_rts = np.zeros((len(infreq_trials)))
    for i in range(len(infreq_trials)):
        trailing_rts[i] = np.nanmean(rts[infreq_trials[i]-3:infreq_trials[i]])
    [rt_bins,rt_bin_edges] = np.histogram(trailing_rts[~np.isnan(trailing_rts)],10)
    for i in range(10):
        infreq_in_bin = np.nansum(np.logical_and(rt_bin_edges[i]<trailing_rts,trailing_rts<rt_bin_edges[i+1]))
        infreq_correct_in_bin = np.nansum(np.logical_and(np.logical_and(rt_bin_edges[i]<trailing_rts,trailing_rts<rt_bin_edges[i+1]),dat[subj].acc[infreq_trials]==1))
        prop_correct[isubj,i] = infreq_correct_in_bin/infreq_in_bin
    # clf = LogisticRegression(random_state=0).fit(prop_correct[isubj], np.arange(0,10))
    # p[isubj] = np.nanmean(clf.predict_proba(prop_correct[isubj]))
    p[isubj] = np.polyfit(np.arange(0,10),prop_correct[isubj],1)

print(p)
print(np.nanmean(p,0))

#make figure
fig = plt.figure(figsize=(12,10))
ax = plt.gca()
for i in range(10):
    plt.scatter(np.ones(nsubj)*i,prop_correct[:,i],color='gray',alpha=.3)
plt.plot(np.arange(0,10),np.nanmean(prop_correct,axis=0),color='k')
[t,tp] = ttest_1samp(prop_correct[:,8],prop_correct[:,9],nan_policy='omit')
print(np.nanmean(tp),np.nanmean(prop_correct[:,8])-np.nanmean(prop_correct[:,9]))
prettify_plot(ax,xl='response time bin',yl='prop. of accurate lures (%)',xt=[0,3,6,9],xtl=[1,4,7,10],yt=[0,.2,.4,.6,.8,1],ytl=[0,20,40,60,80,100])

fig.savefig(figure_dir + 'suppfig4_rts-acc-expt1v2.pdf', bbox_inches='tight')

project_name = 'eyeSustAttnWM03'

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
prop_correct = np.zeros((nsubj,10))
p = np.zeros((nsubj,2))
for isubj in range(len(subj_name)):
    #load data
    subj = subj_name[isubj] #string of subjects name
    dat[subj] = pd.read_csv(subjects_dir+subj+'/data/beh/'+subj+'_explog.csv',header=12)
    #index frequent trials
    infreq_trials = np.where(dat[subj].freq_trials==0)[0]
    rts = np.array(dat[subj].rts)
    trailing_rts = np.zeros((len(infreq_trials)))
    for i in range(len(infreq_trials)):
        trailing_rts[i] = np.nanmean(rts[infreq_trials[i]-3:infreq_trials[i]])
    [rt_bins,rt_bin_edges] = np.histogram(trailing_rts[~np.isnan(trailing_rts)],10)
    for i in range(10):
        infreq_in_bin = np.nansum(np.logical_and(rt_bin_edges[i]<trailing_rts,trailing_rts<rt_bin_edges[i+1]))
        infreq_correct_in_bin = np.nansum(np.logical_and(np.logical_and(rt_bin_edges[i]<trailing_rts,trailing_rts<rt_bin_edges[i+1]),dat[subj].acc[infreq_trials]==1))
        prop_correct[isubj,i] = infreq_correct_in_bin/infreq_in_bin
    p[isubj] = np.polyfit(np.arange(0,10),prop_correct[isubj],1)

print(p)
print(np.nanmean(p,0))

#make figure
fig = plt.figure(figsize=(12,10))
ax = plt.gca()
for i in range(10):
    plt.scatter(np.ones(nsubj)*i,prop_correct[:,i],color='gray',alpha=.3)
plt.plot(np.arange(0,10),np.nanmean(prop_correct,axis=0),color='k')
[t,tp] = ttest_1samp(prop_correct[:,8],prop_correct[:,9],nan_policy='omit')
print(np.nanmean(tp),np.nanmean(prop_correct[:,8])-np.nanmean(prop_correct[:,9]))
prettify_plot(ax,xl='response time bin',yl='prop. of accurate lures (%)',xt=[0,3,6,9],xtl=[1,4,7,10],yt=[0,.2,.4,.6,.8,1],ytl=[0,20,40,60,80,100])

fig.savefig(figure_dir + 'suppfig4_rts-acc-expt1v2_expt2.pdf', bbox_inches='tight')