import numpy as np
import matplotlib.pyplot as plt
import os, copy
import matplotlib.font_manager as fm
from matplotlib import rc
from scipy.stats import binom, beta
from scipy.integrate import quad
from scipy.special import comb, gamma, gammaln
#import scipy
#import scipy.stats as stats
#import pylab as pl
#import matplotlib
#from scipy import stats
#import math

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

@np.vectorize
def calculate_aprime(h,fa):
    '''
    Calculates sensitivity, according to the aprime formula
    '''
    if np.greater_equal(h,fa): a = .5 + (((h-fa) * (1+h-fa)) / (4 * h * (1-fa)))
    else: a = .5 - (((fa-h) * (1+fa-h)) / (4 * fa * (1-h)))
    return a

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

def resampling_statistics1d(data,chance,nsubj=None,nsamples=100000,skipfig=True):
    '''
    Nonparametric resampling statistics
    '''
    if nsubj is None: nsubj = np.size(data)
        
    #resample subjects with replacement 100,000 times
    subj_resampled = np.random.randint(0,nsubj,(nsamples,nsubj)) 

    #the empy matrix of resampled data for each of these resampled iterations
    data_resampled = np.empty(nsamples) 

    #recalculate mean given the resampled subjects
    for i in range(0,nsamples): data_resampled[i] = np.mean(data[subj_resampled[i,:]])

    #calculate p value 
    p = np.sum(np.less(data_resampled,chance))/float(nsamples) #count number of resample iterations below chance
    if np.equal(p,0): p = 1./float(nsamples)

    if skipfig is False:
        plt.figure(figsize=(4,3))
        ax = plt.subplot(111)
        plt.hist(data_resampled,normed=0,facecolor='gray',edgecolor='gray')
        plt.axvline(np.mean(data_resampled),color='b',lw=2,label='resampled mean')
        plt.axvline(np.mean(data),color='m',lw=1,label='original mean')
        plt.axvline(chance,color='c',lw=2,label='chance')
        make_plot_pretty(ax,ylrot=90,yl='Count (#)',legend='1') 
        plt.show()
    
    return p, data_resampled

#computational model
betaPdf = lambda x, a, k: beta(a, 1.).pdf(x/k) / k

binomPMF = lambda i, j, p: comb(j,i) * p**i * (1-p) **(j-i)

# return the probablity of occurence k for a beta binomial n, alpha, beta
def beta_binom_density(k, n, a, b):
    return 1.0*gamma(n+1)*gamma(a+k)*gamma(n+b-k)*gamma(a+b) / (gamma(k+1)*gamma(n-k+1)*gamma(a+b+n)*gamma(a)*gamma(b))

def successPdf_betaBinom(nS, a, k):#PDF for the Beta Binomial
    f = lambda nR: beta_binom_density(nR, k, a,1 ) * binomPMF(nS-nR, 6-nR, 1./(9-nR))**(np.logical_or(nR<6,nR>nS))
    result = np.sum([f(nR) for nR in range(k+1)])
    return result

def likelihood_betaBinom(NS, a, k): #calculates the log likelihood for the Beta Binomial
    #NS is the array of the number correct for each subject (e.g. 1x144 for Exp 1a and 1x180 for Exp 1b)
    f = lambda nS: np.log(successPdf_betaBinom(nS, a, k))
    temp = np.empty(np.size(NS))
    for nS in np.unique(NS):
        temp[nS==NS] = f(nS)
    result = np.sum(temp) #loglikelihood
    return result, temp

def runAnalyses(data, likelihood_model,alpha_vals,kmax_vals):
    #preallocate
    result = np.empty((np.size(alpha_vals),np.size(kmax_vals)))
    alpha_mat = copy.deepcopy(result)
    kmax_mat = copy.deepcopy(result)
    for ialpha, alpha in enumerate(alpha_vals):
        for ikmax, kmax in enumerate(kmax_vals):
            alpha_mat[ialpha,ikmax] = alpha
            kmax_mat[ialpha,ikmax] = kmax
            _,temp = likelihood_model(data.astype(int),alpha,kmax)
            result[ialpha,ikmax] = np.sum(temp)#loglikelihood

    #identify the max log likelihood
    ll_max = np.max(result)
    idx_max = int(np.argmax(result))
    i,j = np.unravel_index(idx_max,np.shape(result))
    alpha= alpha_mat[i,j]
    kmax = kmax_mat[i,j]
    
    return alpha, kmax, alpha_mat, kmax_mat, ll_max, result