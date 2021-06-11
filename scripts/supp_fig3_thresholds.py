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
fast_thresh_mean = np.zeros((nsubj))
slow_thresh_mean = np.zeros((nsubj))
propSimvAct_fast = np.zeros((nsubj))
propSimvAct_slow = np.zeros((nsubj))
fast_thresh_block = np.zeros((nsubj,5))
slow_thresh_block = np.zeros((nsubj,5))
fast_thresh_chart = np.zeros((nsubj,800))
slow_thresh_chart = np.zeros((nsubj,800))

for isubj in range(len(subj_name)):
    #load data
    subj = subj_name[isubj] #string of subjects name
    dat[subj] = pd.read_csv(subjects_dir+subj+'/data/beh/'+subj+'_explog.csv',header=12)
    #index frequent trials
    fast_thresh_mean[isubj] = np.nanmean(dat[subj].rt_runningfastthresh)
    slow_thresh_mean[isubj] = np.nanmean(dat[subj].rt_runningslowthresh)
    fast_thresh_init = dat[subj].rt_runningfastthresh[80]
    slow_thresh_init = dat[subj].rt_runningslowthresh[80]
    trailingRT = dat[subj].rt_trailingavg
    fast_triggered_ind = np.where(dat[subj].rt_triggered_fast==1)[0]+1
    slow_triggered_ind = np.where(dat[subj].rt_triggered_slow==1)[0]+1
    sum_fast = len(fast_triggered_ind)
    sum_slow = len(slow_triggered_ind)
    sum_fast_sim = sum(trailingRT<fast_thresh_init)
    sum_slow_sim = sum(trailingRT>slow_thresh_init)
    propSimvAct_fast[isubj] = sum_fast_sim/sum_fast
    propSimvAct_slow[isubj] = sum_slow_sim/sum_slow
    for iblock in range(5):
        fast_thresh_block[isubj,iblock] = np.nanmean(dat[subj].rt_runningfastthresh[(800*iblock):(800*(iblock+1))])
        slow_thresh_block[isubj,iblock] = np.nanmean(dat[subj].rt_runningslowthresh[(800*iblock):(800*(iblock+1))])
        fast_thresh_chart[isubj,iblock] = fast_thresh_chart[isubj,iblock]+dat[subj].rt_runningfastthresh[(800*iblock):(800*(iblock+1))]
        slow_thresh_chart[isubj,iblock] = slow_thresh_chart[isubj,iblock]+dat[subj].rt_runningslowthresh[(800*iblock):(800*(iblock+1))]

fast_thresh_diff_subj = np.zeros((nsubj,nsubj))
slow_thresh_diff_subj = np.zeros((nsubj,nsubj))
fast_block_tresh_init = np.zeros((nsubj,5,5))
slow_block_tresh_init = np.zeros((nsubj,5,5))
fast_thresh_chart[isubj,iblock] = fast_thresh_chart[isubj,iblock]/5.0
slow_thresh_chart[isubj,iblock] = slow_thresh_chart[isubj,iblock]/5.0
for isubj in range(len(subj_name)):
    for jsubj in range(nsubj):
        fast_thresh_diff_subj[isubj,jsubj] = abs(fast_thresh_mean[isubj]-fast_thresh_mean[jsubj])
        slow_thresh_diff_subj[isubj,jsubj] = abs(slow_thresh_mean[isubj]-slow_thresh_mean[jsubj])
    for iblock in range(5):
        for jblock in range(5):
            fast_block_tresh_init[isubj,iblock,jblock] = abs(fast_thresh_block[isubj,iblock] - fast_thresh_block[isubj,jblock])
            slow_block_tresh_init[isubj,iblock,jblock] = abs(slow_thresh_block[isubj,iblock] - slow_thresh_block[isubj,jblock])
fast_block_thresh = np.nanmean(fast_block_tresh_init,axis=0)
slow_block_thresh = np.nanmean(slow_block_tresh_init,axis=0)

#make figure
# fig = plt.figure(figsize=(12,10))
# ax = plt.subplot(221)
# ax.imshow(fast_thresh_diff_subj,cmap='hot',interpolation='nearest')
# prettify_plot(ax,xt=[0,nsubj-1],xl="Subject",yt=[0,nsubj-1],yl="Subject",t="fast threshold means")
# ax = plt.subplot(222)
# ax.imshow(slow_thresh_diff_subj,cmap='hot',interpolation='nearest')
# prettify_plot(ax,xt=[0,nsubj-1],xl="Subject",yt=[0,nsubj-1],yl="Subject",t="slow threshold means")
# ax = plt.subplot(223)
# ax.imshow(fast_block_thresh,cmap='hot',interpolation='nearest')
# prettify_plot(ax,xt=[0,4],xl="block",yt=[0,4],yl="block",t="fast threshold block means")
# ax = plt.subplot(224)
# ax.imshow(slow_block_thresh,cmap='hot',interpolation='nearest')
# prettify_plot(ax,xt=[0,4],xl="block",yt=[0,4],yl="block",t="slow threshold block means")
# fig.suptitle("similarities between thresholds between subjects and between blocks (darker is more similar)")

fig = plt.figure(figsize=(12,10))
ax = plt.subplot(211)



print("proportion of fast trials still in with static threshold: ",str(np.nanmean(propSimvAct_fast)))
print("proportion of slow trials still in with static threshold: ",str(np.nanmean(propSimvAct_slow)))

fig.savefig(figure_dir + 'suppfig3_thresholds.pdf', bbox_inches='tight')