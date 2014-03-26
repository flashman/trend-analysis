from math import log, lgamma, exp

def xlog2x(val):
    """Compute x*log2(x).  Does NOT handle case of x=0"""
    try:
        return val*log(val)/log(2)
    except:
        return 0

def choose(n,m):
    """Probability of m success out of n events, where an individual event succeeds with probability P... Useful for calculating <H^hat(j|w)> in semantic_information fun\
ction below"""
    return exp(lgamma(n+1) - lgamma(m+1) - lgamma(n-m+1))


def p(n,m,p):
    """Probability of m success out of n events, where an individual event succeeds with probability P... Useful for calculating <H^hat(j|w)> in semantic_information fun\
ction below"""
    try:
        return exp(lgamma(n+1)  - lgamma(n-m+1) - lgamma(m+1) +  m*log(p) + (n-m)*log(1.0-p))
    except:
        print "WARNING: domain range errer...returning 0"
        return 0

def centerofmass(d):
    """Compute the index of the center of mass for a dict of values"""
    m = sum(d.values())
    if m != 0:
        cm = 1.0*sum(i*n for i,n in d.iteritems())/m
    else:
        cm=(len(d)-1)/2.0
    return int(round(cm))


def information(row,month,binsize,offset=0,freqthresh=10,exact=False):
    '''
    Compute symantic information content
    '''
    M=len(month)
    L=( (M-offset)/binsize )*binsize
    B=1.0*L/binsize                         #actual number of bins after applying offset                                                                                 

    if B<2:
        return 0.0
    M_j = [sum(month[i:i+binsize]) for i in xrange(offset,offset+L,binsize)]
    N = float(sum(M_j))
    M_min = min(M_j)
    m_j = [sum( row[j] for j in xrange(i,i+binsize) if j in row) for i in xrange(offset,offset+L,binsize)]
    n = float(sum(m_j))
    if n>freqthresh:
        h = -sum( xlog2x(m/n) for m in m_j if m>0)
        if not exact and n*M_min/N>10:
            h_avg = -sum( xlog2x(M/N) for M in M_j)
        else:
            h_avg = -sum( p(n,m,M/N)*xlog2x(m/n)
                                for M in M_j
                                for m in xrange(1,1+int(min(n,M))))
        return (n/N) * (h_avg - h)
    else:
        return 0.0

def offset(d,r,l):
    '''
    Shift values in dictionary so that the central of mass is centered in some bin of size r.
    Clip dictionary so that r divides its length.  Return shifted and clipped list
    '''
    if l/r < 2 or r < 3:
        return 0
    m = float(max(d.values()))
    ilarge = [i for i,j in d.iteritems() if j>= 0.5*m]
    imaxc = centerofmass(d)
    offset = imaxc % r - r/2
    if len(ilarge)>0.7*l:
        offset = 0
    if offset<0:
        offset += r
    if (l-offset)/r<2 or imaxc-offset<0:
        offset= 0
    return offset

def monthlytotals(data,n=258):
    totals=n*[0]
    for key,vals in data.iteritems():
        for i,v in enumerate(vals):
            totals[i]+=v
    return totals
