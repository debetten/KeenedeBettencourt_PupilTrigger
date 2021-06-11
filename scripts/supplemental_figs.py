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
rts_fast = np.zeros((nsubj))
rts_slow = np.zeros((nsubj))
wmacc_fast = np.zeros((nsubj))
wmacc_slow = np.zeros((nsubj))
pupil_fast = np.zeros((nsubj))
pupil_slow = np.zeros((nsubj))
sum_fast = 0
sum_slow = 0
sum_fast0 = 0
sum_slow0 = 0
sum_fast6 = 0
sum_slow6 = 0
wmacc_fast = np.zeros((nsubj,300))
wmacc_fast[:] = np.nan
wmacc_slow = np.zeros((nsubj,300))
wmacc_slow[:] = np.nan

for isubj in range(len(subj_name)):
    #load data
    subj = subj_name[isubj] #string of subjects name
    dat[subj] = pd.read_csv(subjects_dir+subj+'/data/beh/'+subj+'_explog.csv',header=12)
    #index frequent trials
    fast_triggered_ind = np.where(dat[subj].rt_triggered_fast==1)[0]+1
    slow_triggered_ind = np.where(dat[subj].rt_triggered_slow==1)[0]+1
    sum_fast = sum_fast+ len(fast_triggered_ind)
    sum_slow = sum_slow+ len(slow_triggered_ind)
    #index pr1bes
    probe_ind = np.where(dat[subj].probe_trials==1)[0]
    wmacc_fast[isubj,:np.size(fast_triggered_ind)] = dat[subj].wholereportrespacc[fast_triggered_ind] 
    wmacc_slow[isubj,:np.size(slow_triggered_ind)] = dat[subj].wholereportrespacc[slow_triggered_ind] 
    sum_fast0 = sum_fast0 + np.nansum(wmacc_fast[isubj]==0)
    sum_slow0 = sum_slow0 + np.nansum(wmacc_slow[isubj]==0)
    sum_fast6 = sum_fast6 + np.nansum(wmacc_fast[isubj]==6)
    sum_slow6 = sum_slow6 + np.nansum(wmacc_slow[isubj]==6)

    #rts_fast[isubj] = np.nanmean(dat[subj].rts[fast_triggered_ind])
    #rts_slow[isubj] = np.nanmean(dat[subj].rts[slow_triggered_ind])

print("#fast = ",str(sum_fast),"mean fast = ",str(sum_fast/len(subj_name)))
print("prop fast = ",str(sum_fast/(sum_fast+sum_slow)))
print("#slow = ",str(sum_slow),"mean slow = ",str(sum_slow/len(subj_name)))
print("prop slow = ",str(sum_slow/(sum_fast+sum_slow)))
print("#fast0 = ",str(sum_fast0),"mean fast = ",str(sum_fast0/len(subj_name)))
print("prop fast0 = ",str(sum_fast0/(sum_fast0+sum_slow0)))
print("#slow0 = ",str(sum_slow0),"mean slow = ",str(sum_slow0/len(subj_name)))
print("prop slow0 = ",str(sum_slow0/(sum_fast0+sum_slow0)))
print("#fast6 = ",str(sum_fast6),"mean fast = ",str(sum_fast6/len(subj_name)))
print("prop fast6 = ",str(sum_fast6/(sum_fast6+sum_slow6)))
print("#slow6 = ",str(sum_slow6),"mean slow = ",str(sum_slow6/len(subj_name)))
print("prop slow6 = ",str(sum_slow6/(sum_fast6+sum_slow6)))

#make figure
fig, ax = plt.subplots(1,1,figsize=(12,10))
#fig 1b - attention (A')
x1 = np.ndarray.flatten(wmacc_fast)
ax.hist(x1[~np.isnan(x1)],bins=np.arange(-1,8),align='left',
        color='red',edgecolor='k',linewidth=3,density=True,alpha=.5)
x2 = np.ndarray.flatten(wmacc_slow)
ax.hist(x2[~np.isnan(x2)],bins=np.arange(-1,8),align='left',
        color='blue',edgecolor='k',linewidth=3,density=True,alpha=.5)
t,p = ttest_1samp(np.nanmean(x1),np.nanmean(x2),nan_policy='omit')
print(np.shape(p))
print('t = ',str(t),'p = ',str(p))
for isubj in range(nsubj):
    x1 = wmacc_fast[isubj]
    y,b = np.histogram(x1[~np.isnan(x1)],bins=np.arange(8),density=True)
    ax.scatter(b[:-1],y,color='red',alpha=.33,zorder=20,clip_on='False')
for isubj in range(nsubj):
    x2 = wmacc_slow[isubj]
    y,b = np.histogram(x2[~np.isnan(x2)],bins=np.arange(8),density=True)
    ax.scatter(b[:-1],y,color='blue',alpha=.33,zorder=20,clip_on='False')

prettify_plot(ax,yl="Proportion of probes",xlim=[-1,7],yt=[0,.4,.8],ytl=[0,.4,.8],
    xt=([0,1,2,3,4,5,6]),xtl=([0,1,2,3,4,5,6]),xl="num correct",t='distribution of WM probe accuracy for fast and slow')
fast_handle = mpatches.Patch(color = 'red',label = 'fast')
slow_handle = mpatches.Patch(color = 'blue',label = 'slow')
ax.legend(handles=[fast_handle,slow_handle])

#plt.subplots_adjust(left=.125, right=.9, wspace=1)

#plt.show()

fig.savefig(figure_dir + 'suppfig1_fastvsslowWMdist.pdf', bbox_inches='tight')

###plot with lines instead of bars
###print summary stats and add them in, diff for each #corrent