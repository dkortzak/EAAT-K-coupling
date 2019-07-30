#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from scipy import stats
import sys, os, time
from random import gauss, randint, sample, choice
#The following functions were taken from the pmx package
#
#-----------------------------------------------------------------------
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2011 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
def data_to_gauss( data ):
    m = np.mean( data )
    dev = np.std( data )
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A
def gauss_func( A, mean, dev, x):
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2)/(2.0*(dev**2))))
    return y
def gauss_intersection( g1, g2 ):
    A1, m1, s1 = g1
    A2, m2, s2 = g2
    
    p1 = m1/s1**2-m2/s2**2
    p2 = np.sqrt(1/(s1**2*s2**2)*(m1-m2)**2+2*(1/s1**2-1/s2**2)*np.log(s2/s1))
    p3 = 1/s1**2-1/s2**2
    x1 = (p1+p2)/p3
    x2 = (p1-p2)/p3
    # determine which solution to take
    if x1 > m1 and x1 < m2 or \
       x1 > m2 and x1 < m1:
        return x1
    elif x2 > m1 and x2 < m2 or \
       x2 > m2 and x2 < m1:
        return x2
    else:
        return False # we do not take the intersection
def smooth(x,window_len=11,window='hanning'):

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')
    y=convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]      
def make_plot(  data1, data2, result, err, nbins, dpi=300 ):

    plt.figure( figsize = (3, 2) )
    mf, devf, Af = data_to_gauss( data1 )
    mb, devb, Ab = data_to_gauss( data2 )
    
    maxi = max( list(data1)+list(data2) )
    mini = min( list(data1)+list(data2) )
    n1, bins1, patches1 = plt.hist(data1, range = (mini,maxi),bins=nbins, facecolor='blue', alpha=0.75, normed=True, label='0->1')
    n2, bins2, patches2 = plt.hist(data2, range = (mini,maxi),bins=nbins, facecolor='red', alpha=0.75, normed=True, label='1->0')
    plt.xlabel('W [kJ/mol]', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.title(r'Work Distribution $\lambda$ 0->1 (blue) $\lambda$ 1->0 (red)')
    plt.grid(lw = 2)
    loc, lab = plt.yticks()
    ll = []
    for i in range(len(lab)):
        ll.append("")
    plt.yticks( loc, ll )
    x = np.arange( mini, maxi, .5 )
    y1 = gauss_func( Af, mf, devf, x )
    y2 = gauss_func( Ab, mb, devb, x )
    
    plt.plot(x, y1, 'b--', linewidth=2)
    plt.plot(x, y2, 'r--', linewidth=2)
    
    size = max( [max(y1), max(y2)] )
    res_x = [result, result ]
    res_y = [0, size*1.2 ]
    plt.plot( res_x, res_y, 'k--', linewidth=2, label = r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    plt.legend(shadow=True, fancybox = True)
    plt.ylim(0, size*1.2 )
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.show()
    
    
def make_W_over_time_plot( fname, data1, data2, result, err, nbins, dpi):
    
    
    sb.set_style("white",{'axes.grid': False, 'grid.linestyle': ':','ytick.major.size': 5,'ytick.direction': 'in','axes.linewidth': .5,'lines.linewidth': .1})
    
    plt.figure( figsize = (2.5, 1.5) )
    x1 = range( len(data1) )
    x2 = range( len(data2) )
    if x1>x2: x = x1
    else: x = x2
    mf, devf, Af = data_to_gauss( data1 )
    mb, devb, Ab = data_to_gauss( data2 )
    
    maxi = max( list(data1)+list(data2) )
    mini = min( list(data1)+list(data2) )
    #print mini,maxi
    sm1 = smooth( array(data1) )
    sm2 = smooth( array(data2) )
    subplot(1,2,1)
    plot(x1,data1,'g-',linewidth=1 ,label="Forward")

    plot(x2,data2,'b-',linewidth=1 ,label="Backward")

    xticks([50,100,150])
    legend(loc='upper center')
    ylabel(r'Work (kJ/mol)')
    xlabel(r'Snapshot')
    

    xlim(0,x[-1]+1)

    sb.set_style("white",{'axes.grid': False, 'grid.linestyle': ':','ytick.major.size': 5,'ytick.direction': 'in','axes.linewidth': .5,'lines.linewidth': .1})
    subplot(1,2,2)
    
    xlabel(r'Density')

    sb.distplot(data2,kde=False,norm_hist=True,vertical=True)
    sb.distplot(data1,kde=False,norm_hist=True,vertical=True)
    


    x = np.arange( mini, maxi, .5 )

    y1 = gauss_func( Af, mf, devf, x )
    y2 = gauss_func( Ab, mb, devb, x )

    plt.plot(y1,x, color='g', linewidth=1)
    plt.plot(y2,x, color='b', linewidth=1)
    size = max( [max(y1), max(y2)] )
    res_x = [result, result ]
    res_y = [0, size*1.2 ]
    plt.plot( res_y, res_x, 'k--', linewidth=1, label = r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    legend(shadow=True, fancybox = True, loc='upper center')
    
    
    

    xticks([0.02,0.04])
    xl = gca()

    
    xl.axes.yaxis.set_ticklabels([])
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace = 0.1)

    plt.show()
    
def cgi_error_from_mean(nruns, mu1, sig1, n1, mu2, sig2, n2):
    iseq = []

    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append( gauss(mu1, sig1))
        for i in range(n2):
            g2.append( gauss(mu2, sig2))
        m1 = average(g1)
        s1 = std(g1)
        m2 = average(g2)
        s2 = std(g2)
        p1 = 1./(s1*np.sqrt(2*np.pi))
        p2 = 1./(s2*np.sqrt(2*np.pi))
        iq = (m1+m2)/2.
        iseq.append(iq)
    mean = average(iseq)
    err = std(iseq)
    return err

def cgi_error(nruns, mu1, sig1, n1, mu2, sig2, n2):
    iseq = []
    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append( gauss(mu1, sig1))
        for i in range(n2):
            g2.append( gauss(mu2, sig2))
        m1 = average(g1)
        s1 = std(g1)
        m2 = average(g2)
        s2 = std(g2)
        p1 = 1./(s1*np.sqrt(2*np.pi))
        p2 = 1./(s2*np.sqrt(2*np.pi))
        iq = gauss_intersection([p1,m1,s1],[p2,m2,s2])
        iseq.append(iq)
    mean = average(iseq)
    err = std(iseq)
    return err

#The following function was modified from the pmx package
#

def calc_err_boot2(wf,wr, nboots):
        '''Calculates the standard error of the Crooks Gaussian Intersection
        via non-parametric bootstrap. The work values are resampled randomly
        with replacement multiple (nboots) times, and the CGI free energy
        recalculated for each bootstrap samples. The standard error of
        the estimate is returned as the standard deviation of the bootstrapped
        free energies.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        nboots: int
            number of bootstrap samples to use for the error estimate.

        Returns
        -------
        err : float
            standard error of the mean.
        '''
        nf = len(wf)
        nr = len(wr)

        dg_boots = []
        
        for k in range(nboots):
            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            mf, devf, Af = data_to_gauss( np.array(bootA) )
	    mb, devb, Ab = data_to_gauss( np.array(bootB) )
	    cgi_result = gauss_intersection( [Af, mf, devf], [Ab, mb, devb ] )
            intersection = True
	    if not cgi_result:
	      cgi_result = (mf+mb)*.5
	     
            dg_boots.append(cgi_result)
            
        print "Mean and SD from bootstrapping"
        print np.mean(dg_boots),np.std(dg_boots)
        return np.mean(dg_boots),np.std(dg_boots)
  
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2011 by Daniel Seeliger
#
#                        All Rights Reserved
def do_cgi(anadir):
  integ_ab=tuple(open(anadir+"/integ0.dat"))
  integ_ba=tuple(open(anadir+"/integ1.dat"))
  res_ab=[]
  res_ba=[]
  for i in range(0,len(integ_ab)):
    res_ab.append(float(integ_ab[i].split()[1]))
  for i in range(0,len(integ_ba)):
    res_ba.append(float(integ_ba[i].split()[1]))
  res_ab=np.array(res_ab)
  res_ba=np.array(res_ba)
  cgi_result,cgi_err = calc_err_boot2(res_ab, res_ba, 1000)
  make_plot(  res_ab, res_ba, cgi_result, cgi_err, 10 )
  
print "EAAT1, Na1 selectivity"
anadir="data/cgi/selectivity/eaat1/Na1"
do_cgi(anadir)
print "EAAT1, Na2 selectivity"
anadir="data/cgi/selectivity/eaat1/Na2"
do_cgi(anadir)
print "EAAT1, Na3 selectivity"
anadir="data/cgi/selectivity/eaat1/Na3"
do_cgi(anadir)
print "EAAT1, Na4 selectivity"
anadir="data/cgi/selectivity/eaat1/Na4"
do_cgi(anadir)
print "GltPh, Na1 selectivity"
anadir="data/cgi/selectivity/gltph/Na1"
do_cgi(anadir)
print "GltPh, Na2 selectivity"
anadir="data/cgi/selectivity/gltph/Na2"
do_cgi(anadir)
print "GltPh, Na3 selectivity"
anadir="data/cgi/selectivity/gltph/Na3"
do_cgi(anadir)