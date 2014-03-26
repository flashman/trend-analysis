import time
import gzip
import cPickle
from utils import xlog2x, p
import matplotlib.pyplot as plt
from multiprocessing import Pool
from myextensions import pxlog2x, information
from math import log

def monthlytotals(data):
    '''Get number of words per month'''
    n = len(data[data.keys()[0]])
    totals=n*[0]
    for key,vals in data.iteritems():
        for i,v in enumerate(vals):
            totals[i]+=v
    return totals

def loadAndNormalizeData(VERSION, limit=256):
    '''Read in data from fileLocation and normalize each month count by the total number of words per month'''
    F = 1000000
    data = cPickle.load(gzip.open('../data/'+VERSION+'.pcl.gz'))
    ndata = dict()
    totals = monthlytotals(data)
    for key, vals in data.iteritems():
        ndata[key] = [ F*v/t for v,t in zip(vals,totals) ]
    return ndata[-limit:]

def loadData(VERSION='mcg2'):
    '''Load data'''
    return cPickle.load(gzip.open('../data/'+VERSION+'.pcl.gz'))

def compute_h_avgs(M_j):
    N=sum(M_j)
    h_avgs = dict()
    for n in range(1,1000):
        h_avgs[n] = -sum( pxlog2x(float(n),float(m),float(M)/N) for M in M_j for m in xrange(1,1+int(min(n,M))))
    return h_avgs

def entropy_transform_old(word, ts,tots):
    """Perform entropy transform of the timeseries ts with respect to totals tots"""
    L=len(tots)
    n=float(sum(ts))
    N=float(sum(tots))
    max_width = L
    r=[]

    for w in xrange(1,max_width,2):
        for t in xrange(0,L-w,2):
            s = sum(ts[t:(t+w)])
            S = sum(tots[t:(t+w)])
            m_j = [ s, n-s ]
            M_j = [ S, N-S ]
            h = -sum( xlog2x(m/n) for m in m_j if m>0 )
            h_avg = -sum( p(n,m,M/N)*xlog2x(m/n)
                                for M in M_j
                                for m in xrange(1,1+int(min(n,M))))
            r.append( (w, t, (n/N) * (h_avg - h) ) )
    return (word,r)

def entropy_transform(word, ts,tots):
    """Perform entropy transform of the timeseries ts with respect to totals tots"""
    L=len(tots)
    m=float(sum(ts))
    M=float(sum(tots))
    max_width = L-1
    r=[]

    for w in xrange(1,max_width,2):
        for t in xrange(0,L-w,2):
            s = sum(ts[t:(t+w)])
            S = sum(tots[t:(t+w)])
            m_j = [ s, (m-s) ]
            M_j = [ S, (M-S) ]
            #m_j = [sum(ts[:t]), sum(ts[t:(t+w)]), sum(ts[(t+w):]) ]
            #M_j = [sum(tots[:t]), sum(tots[t:(t+w)]), sum(tots[(t+w):]) ]
            i =  information(m_j,M_j)*log(L/float(w)+1)
            r.append( (w, t, i ) )
    return (word,r)

def plot_entropy_transform(results,ts=None,word=None):
    '''
    Plot the entropy transform, stored in results, for a specific word.
    results is a list of tuples: (bin-size,offset,information).
    ts is the orignal time series data for the.
    word is the word.
    '''
    #Plot heat map of entropy transform
    x = [r[0] for r in results]
    y = [r[1] for r in results]
    z = [r[2] for r in results]
    m = max(results,key= lambda r: r[2])

    plt.subplot(4, 1, 1)
    plt.scatter( y,x, c=z, marker='s', linewidth='0')
    plt.xlim(min(y),max(y))
    plt.ylim(min(x),max(x))
    plt.ylabel('size')
    plt.xlabel('offset')
    plt.title('%s \n max info: (%d,%d) = %.8f ' % (word, m[1],m[0],m[2]) )

    #plot info as a function of lenght scale
    yl = [-1000000]*(max(x)+1)
    for r in results:
        if yl[r[0]]<r[2]:
            yl[r[0]]=r[2]
    yl = [v for v in yl if v > -1000000]
    plt.subplot(4,1,2)
    plt.plot(range(0,len(yl)*2,2), yl)
    plt.xlim(min(x),max(x))
    plt.ylabel('info')
    plt.xlabel('size')

    #plot info as a function of offset
    yo = [-10000000000]*(max(y)+1)
    for r in results:
        if yo[r[1]]<r[2]:
            yo[r[1]]=r[2]
    yo = [v for v in yo if v > -10000000]
    plt.subplot(4,1,3)
    plt.plot(range(0,len(yo)*2,2), yo)
    plt.xlim(min(y),max(y))
    plt.ylabel('info')
    plt.xlabel('offset')


    plt.subplot(4,1,4)
    plt.plot(ts)
    plt.xlim(0,len(ts))
    plt.savefig('../results/'+word+'.png')
    plt.clf()


if __name__ == '__main__':

    data = loadData()
    tots = monthlytotals(data)
    pool = Pool(processes=4)
    results=dict()

    def _collect_results(result):
        word=result[0]
        m = max(result[1], key=lambda v: v[2])
        results[word]=m
        plot_entropy_transform(result[1], data[word], word)
        print "Finished processing %s" % result[0]

    t = time.time()
    for k,vals in data.items()[:100]:
        print "Processing %s" % k
        pool.apply_async(entropy_transform,[k,vals,tots],callback=_collect_results)
    pool.close()
    pool.join()
    print time.time()-t

    f=gzip.open('../results/mcg2.test.gz','wb')
    cPickle.dump(results,f)
    f.close()

    wm = [(k,vals) for k,vals in results.iteritems()]
    wm = sorted(wm,key=lambda w: w[1][2],reverse=True)
    print 'WORDS SORTED BY MAX INFO:'
    for item in wm:
        print '    %50s m=%10.8f l=%d o=%d' % (item[0], item[1][2], item[1][0], item[1][1])
