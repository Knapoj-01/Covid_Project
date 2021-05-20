import numpy as np
def RungeKuttaSolve(slope, point,ds, ic, npoints):
        if len(point) != 2:
            raise TypeError("Point must be 1x2 array!")
        if  point[0] >= point[1]:
            raise Exception("Invalid Input")
        yresult = np.zeros(npoints + 1)
        h = (point[1] - point[0])/npoints
        m = 0; k = 0
        yresult[0] = ic
        if ds != point[0]:
            while m*h <= ds-point[0] - h:
                k1 = -1*slope(ds - m*h, yresult[m])
                k2 = -1*slope(ds - m*h - h/2, yresult[m] + h*k1/2)
                k3 = -1*slope(ds - m*h - h/2, yresult[m] + h*k2/2)
                k4 = -1*slope(ds - m*h - h, yresult[m] + h*k3)
                yresult[m+1] = yresult[m] + h*(k1 + 2*k2 + 2*k3 + k4)/6
                m = m + 1
            yresult = np.roll(np.flip(yresult), m+1)      
        while k*h < point[1] - ds:
            k1 = slope(ds + k*h,yresult[k + m])
            k2 = slope(ds + k*h + h/2, yresult[k + m] + h*k1/2)
            k3 = slope(ds + k*h + h/2, yresult[k + m] + h*k2/2)
            k4 = slope(ds + k*h + h, yresult[k + m] + h*k3)
            yresult[k + m + 1] = yresult[k + m] + h*(k1 + 2*k2 + 2*k3 + k4)/6
            k = k + 1
        xresult = np.linspace(point[0],point[1], npoints+1)
        return xresult, yresult
