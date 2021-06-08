#%%
import numpy as np

def RungeKuttaSolve(f, interval, y0, npoints=500, t0 = False):
    if t0 == False:
        t0 = interval[0]
    xresult, h = np.linspace(interval[0],interval[1], npoints, retstep=True)
    yresult = np.zeros((len(y0), len(xresult)))
    m = 0; k = 0
    yresult[:,0] = y0
    if t0 != interval[0]:
        while m*h <= t0-interval[0] - h:
            for i in range(0, len(y0)):
                a = yresult[:, m]
                a1 = np.copy(a)
                k1 = -f(t0 - k*h,a1)[i]
                a2 = np.copy(a); a2[i] = a2[i] + h*k1/2
                k2 = -f(t0 - k*h - h/2, a2)[i]
                a3 = np.copy(a); a3[i] = a3[i] + h*k2/2
                k3 = -f(t0 - k*h - h/2, a3)[i]
                a4 = np.copy(a); a4[i] = a4[i] + h*k3
                k4 = -f(t0 - k*h - h, a4)[i]
                yresult[i, m+1] = yresult[i, m] + h*(k1 + 2*k2 + 2*k3 + k4)/6
            m = m + 1
        yresult = np.roll(np.flip(yresult), m+1) 
    while k*h < interval[1] - t0 :
        for i in range(0, len(y0)):
            a = yresult[:, k + m]
            a1 = np.copy(a)
            k1 = f(t0 + k*h,a1)[i]
            a2 = np.copy(a); a2[i] = a2[i] + h*k1/2
            k2 = f(t0 + k*h + h/2, a2)[i]
            a3 = np.copy(a); a3[i] = a3[i] + h*k2/2
            k3 = f(t0 + k*h + h/2, a3)[i]
            a4 = np.copy(a); a4[i] = a4[i] + h*k3
            k4 = f(t0 + k*h + h, a4)[i]
            yresult[i, k + m + 1] = yresult[i, k + m] + h*(k1 + 2*k2 + 2*k3 + k4)/6
        k = k + 1
    
    return xresult, yresult

# %%
