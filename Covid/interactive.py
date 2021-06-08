# Interactive Plot
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import matplotlib.pyplot as plt

def interactPlotting(obj, plotrange, npoints = 300):
    def plotI():
        Ix, Iy = obj.predictI(plotrange)
        fig, ax = plt.subplots()
        ax.plot(Ix, Iy)
        ax.scatter(obj.getData("day")[plotrange[0]:plotrange[1]],\
                    obj.getData("infecteds")[plotrange[0]:plotrange[1]], 0.8)
    #def changeN(N = obj.N):
    #        obj.resetN(N)
    #rangeN = obj.N - round(0.3*obj.N), obj.N + round(0.3*obj.N), 1000
    #interact(changeN, N = rangeN)
    if hasattr(obj, "k"):
        if hasattr(obj, "dt"):
            deltaTime = obj.dt
        else: deltaTime = 1
        def changeInitials(tr = 1/obj.b, k =obj.k, ds=obj.ds, dt=deltaTime):
            obj.computeIC(tr,k,ds,dt)  
            plotI()   
        rangetr = 1/obj.b -5 ,1/obj.b+5 , 0.5
        rangek = obj.k-3,  obj.k+3, 0.05
        rangeds = obj.ds - 10, obj.ds + 10, 1
        rangedt = 1,10,1
        interact(changeInitials, k=rangek, tr = rangetr, ds = rangeds, dt = rangedt)
    else:
        if hasattr(obj, "dt"):
            deltaTime = obj.dt
        else: deltaTime = 1
        def changeInitials(tr = 1/obj.b, ds=obj.ds, dt=deltaTime):
            obj.computeIC(tr,ds,dt)
            plotI()
        rangetr = 1/obj.b -5 ,1/obj.b+5 , 0.5
        rangeds = obj.ds - 10, obj.ds + 10, 1
        rangedt = 1,10,1
        interact(changeInitials, tr = rangetr, ds = rangeds, dt = rangedt)