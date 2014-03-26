import time
import gzip
import cPickle
import matplotlib.pyplot as plt
from multiprocessing import Pool
from myextensions import information

def monthlytotals(data):
    '''Get number of words per month'''
    n = len(data[data.keys()[0]])
    totals=n*[0]
    for key,vals in data.iteritems():
        for i,v in enumerate(vals):
            totals[i]+=v
    return totals

def loadData(VERSION='mcg2'):
    '''Load data'''
    return cPickle.load(gzip.open('../data/'+VERSION+'.pcl.gz'))

def symmantic_information(word,row,month):
    '''
    Compute symantic information content
    '''
    r=[]
    M=len(month)
    for binsize in xrange(1,M/2+1):
        #if M % binsize != 0:
            #continue 
        io = []
        for offset in xrange(0,binsize):
            L=M-offset
            M_j = [sum(month[i:(i+binsize)]) for i in xrange(offset,offset+L,binsize)]
            m_j = [sum( row[i:(i+binsize)]) for i in xrange(offset,offset+L,binsize)]
            #M_j[-1]+=sum(month[:i])
            #m_j[-1]+=sum(row[:i])
            i=information(m_j,M_j)
            io.append(i)
        r.append( (binsize, max(io)) )
    return (word,r)

def plot_symantic_information(results,ts=None,word=None):
    x = [r[0] for r in results]
    y = [r[1] for r in results]
    m = max(results,key= lambda r: r[1])

    plt.subplot(2, 1, 1)
    plt.plot(x,y)
    plt.xlim(min(x),max(x))
    plt.ylim(min(y),max(y))
    plt.title('%s max: (%d, %f)' % (word, m[0],m[1] ) )

    plt.subplot(2,1,2)
    plt.plot(ts)
    plt.xlim(0,len(ts))
    plt.ylim(0,max(ts))
    plt.savefig('../sym-results/'+word+'.png')
    plt.clf()


if __name__ == '__main__':

    data = loadData()
    tots = monthlytotals(data)
    pool = Pool(processes=4)
    results=dict()

    def _collect_results(result):
        word=result[0]
        m = max(result[1], key=lambda v: v[1])
        results[word]=m
        plot_symantic_information(result[1], data[word], word)
        print "Finished processing %s" % result[0]

    t = time.time()
    for k,vals in data.items()[:100]:
        print "Processing %s" % k
        pool.apply_async(symmantic_information,[k,vals[:256],tots[:256]],callback=_collect_results)
    pool.close()
    pool.join()
    print time.time()-t

    f=gzip.open('../sym-results/mcg4.gz','wb')
    cPickle.dump(results,f)
    f.close()

    wm = [(k,vals) for k,vals in results.iteritems()]
    wm = sorted(wm,key=lambda w: w[1][1],reverse=True)
    print 'WORDS SORTED BY MAX INFO:'
    for item in wm:
        print '    %50s m=%10.5f l=%d ' % (item[0], item[1][1], item[1][0])
