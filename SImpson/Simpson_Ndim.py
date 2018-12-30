#routine for Simpson integration on a hyper-triangular grid

def f(length, conf, tmax):
    if length == 1:
        return Simpson(conf[0], tmax)
    else:
        for i in xrange(length):
            l = conf[-1]
            k = conf[-2]
            if k<tmax-2:
                return Simpson(l - k, tmax-k)*f(length-1,conf[:-1], tmax)
            else:
                return weight(l, k, tmax) * f(length-1,conf[:-1], tmax)

def weight(l, k, tmax):
    if k == tmax:
        w =0.
    elif k == tmax-1:
        w= 1/2.
    elif k==tmax-2:
        if l==tmax:
            w= 1/3.
        elif l==tmax-1:
            w=4/3.
        else:
            w=1/3.
    return w

def Simpson(l, k):
    if k % 2 == 0:
        if (l == 0 or l == k):
            weight = 1/3.
        elif l % 2 != 0:
            weight = 4/3.
        else:
            weight = 2/3.
    else:
        if l == 0:
            weight = 1/3.
        elif l == k:
            weight = 1/2.
        elif l == k -1:
            weight = 5/6.
        elif l % 2 != 0:
            weight = 4/3.
        else:
            weight = 2/3.
    return weight

