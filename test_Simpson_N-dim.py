from SImpson import Simpson_Ndim
import numpy as np

def function1(x):
    return x**4

def function2(x, y):
    return x**2 * y + x * y**3

def function3(x, y, z):
    return x**2 * y * z**3

def function4(x, y, z):
    return np.exp(-(x+y+z)) * np.cos(x+y+z)

def function5(y, z, omega):
    return (y+z) * np.cos(omega*(y+z))

def function6(x, y, z, omega):
    return (x+y+z) * np.cos(omega*(x+y+z))

epsilon = 1e-2

def test_triangularsimpson():

    #polynomial 1 variables
    length = 1
    dt  = 1e-3
    tmreal = 5
    a = 0
    nt = int(tmreal/dt)
    for i in xrange(nt+1):
        a += function1(i*dt) * Simpson_Ndim.f(length, np.array([i]), nt)
    res1 = a*dt
    res2 = tmreal**5/5.

    assert abs((res1-res2)/res2)<epsilon


    #polynomial 2 variables
    length = 2
    dt  = 0.01
    tmreal = 5
    b = 0
    nt = int(tmreal/dt)
    for i in xrange(nt+1):
        for j in xrange(i, nt+1):
            b += function2(i*dt, j*dt) * Simpson_Ndim.f(length, np.array([i, j]), nt)
    res1 = b*dt**2
    res2 = tmreal**5*(1/6.-1/10.)+tmreal**6*(1/8.-1/24.)
    assert abs((res1-res2)/res2)<epsilon


    #polynomial 3 variables
    length = 3
    dt  = 0.1
    tmreal = 5
    c = 0
    nt = int(tmreal/dt)
    for i in xrange(nt+1):
        for j in xrange(i, nt+1):
            for k in xrange(j, nt+1):
                c += function3(i*dt, j*dt, k*dt) * f(length, [i, j, k], nt)
    res1 = c*dt**3
    res2 = tmreal**9*(1/24. - 1/72. - 1/40. + 1/216.)
    assert abs((res1-res2)/res2)<epsilon

    ###oscillating decaying functions
    length = 3
    dt = 0.05
    tmax = 15
    nt = int(tmax/dt)
    d=0
    for i in xrange(nt+1):
        for j in xrange(i, nt+1):
            for k in xrange(j, nt+1):
                d += function4(i*dt, j*dt, k*dt) * f(length, [i, j, k], nt)
    res1 = d*dt**3
    res2 = 0.25 * (-1/4. * np.exp(-3*tmax) * (np.sin(3*tmax)+np.cos(3*tmax))+1/6.*np.exp(-3*tmax)*(np.sin(3*tmax)+np.cos(3*tmax))+1/4.*np.exp(-tmax)*(np.sin(tmax)+np.cos(tmax))-1/6.)
    assert abs((res1-res2)/res2)<epsilon

    length = 2
    dt = 0.02
    tmax = 5
    omega=5
    nt = int(tmax/dt)
    e=0
    for i in xrange(nt+1):
        for j in xrange(i, nt+1):
            e += function5(i*dt, j*dt, omega) * f(length, [i, j], nt)
    res1 = e*dt**2
    res2 = (np.sin(2*omega*tmax)-tmax*omega*np.cos(2*omega*tmax)-2*np.sin(omega*tmax)+omega*tmax*np.cos(omega*tmax))/omega**3
    assert abs((res1-res2)/res2)<epsilon

    length = 3
    dt = 0.05
    tmax = 10
    nt = int(tmax/dt)
    d=0
    omega=5
    for i in xrange(nt+1):
        for j in xrange(i, nt+1):
            for k in xrange(j, nt+1):
                d += function6(i*dt, j*dt, k*dt, omega) * f(length, [i, j, k], nt)
    res1 = d*dt**3
    res2 = 1/(2.*omega**4)* (-omega*tmax*np.sin(3*omega*tmax)-np.cos(3*omega*tmax)+2*tmax*omega*np.sin(2*omega*tmax)-omega*tmax*np.sin(omega*tmax)+3*np.cos(2*omega*tmax)-3*np.cos(omega*tmax)+1)
    assert abs((res1-res2)/res2)<epsilon


if __name__== "__main__":
    test_triangularsimpson()