import time
import numpy as np 
import pandas as pd
from scipy.special import lambertw
from .rungeKutta import RungeKuttaSolve
class SIR_Model:
    #Class Constructor
    def __init__(self, data, S0):
        self.day = np.copy(data.iloc[:,0])
        self.I = np.copy(data.iloc[:,1])
        self.R = np.copy(data.iloc[:,2] + data.iloc[:,3])
        self.S0=S0
    #Initial Conditions
    def computeIC(self, tr, ds, dt = 1):
        self.b = 1/tr
        self.a = (self.I[ds +dt] - self.I[ds])/(self.S0*self.I[ds]*dt) + self.b/self.S0
        self.ds = ds
        self.dt = dt
        return {"a": self.a, "b":self.b}
    def setIC(self,a,b,ds):
        self.a = a
        self.b = b
        self.ds = ds
    #Internal Functions
    def func(self, t,x):
        f = np.array([
            -self.a*x[0]*x[1],
            self.a*x[0]*x[1] - self.b*x[1]
        ])
        return f
    #Predict Data
    def predictI(self,day):
        self.interval = day
        self.x,self.y = RungeKuttaSolve(
            self.func, day, y0 = [self.S0, self.I[self.ds]], t0 = self.ds, 
            npoints=(day[1] -day[0])*5  
            )
        return self.x,self.y[1]

    def getData(self, kwd): 
        data = {
            "infecteds" : self.I,
            "removeds" : self.R,
            "day": self.day
        }
        return data.get(kwd,'Invalid')

    def report(self):
        print("รายงานผลการพยากรณ์โรคโควิด")
        print("วันที่-เวลาพยากรณ์: {} {}"\
            .format(
                time.strftime("%d-%m-%Y", time.localtime()),
                time.strftime("%H:%M:%S", time.localtime()))
        )
        print('a = {:.2e}'.format(self.a))
        print('b = {:.2e}'.format(self.b))
        print('พยากรณ์ตั้งแต่วันที่ {} จนถึงวันที่ {} ของการระบาด'.format(self.interval[0], self.interval[1]))
        print('วันที่ใช้คำนวณเงื่อนไขเริ่มต้น คือวันที่ {} ของการระบาด'.format(self.ds))
    def resetS0(self, S0):
        self.S0 = S0
    
    def printtable(self, dt = 1):
        predict_val = self.y[1][::5*dt].copy()
        data_day = self.day[self.interval[0]:self.interval[1]:dt].copy()
        data_val = self.I[self.interval[0]:self.interval[1]:dt].copy()
        error_val = np.abs(data_val - predict_val).copy()
        #error_percent = error_val/data_val
        data= np.array([data_day, data_val, predict_val, error_val])
        res = pd.DataFrame(data.T, columns=['วันที่', 'ค่าจริง', 'ค่าพยากรณ์', 'ผลต่าง'])
        return res

class SIRPlus_Model(SIR_Model):
    def __init__(self, data, N):
        super().__init__(data, N)

    def computeIC(self, tr,k,  ds, dt = 1):
        self.b = 1/tr
        self.a = (self.I[ds +dt] - self.I[ds])/(self.S0*self.I[ds]*dt) + self.b/self.S0
        self.k = k
        self.ds = ds
        self.dt = dt
        return {"a": self.a, "b":self.b,"k": k}
        #Internal Functions
    def func(self, t,x):
        f = np.array([
            -self.k*self.a*x[0]*x[1],
            self.a*x[0]*x[1]- self.b*x[1]
        ])
        return f
    def predictI(self,day):
        self.interval = day

        try:
            ds = self.initial_value[0]
            I0 = self.initial_value[1]
        except:
            ds = self.ds
            I0 = self.I[self.ds]

        x,y =  RungeKuttaSolve(
            self.func, day, y0 = [self.S0, I0], t0 = ds, 
            npoints=(day[1] -day[0])*5 
            )

        k = len(x) -1
        self.last_value = np.array([x[k], y[1][k]])
        if hasattr(self, 'initial_value'):
            self.initial_value =np.array([x[k], y[1][k]])
        self.x = x; self.y = y

        return x, y[1]
    def beginCont(self):
        self.initial_value = self.last_value
    def endCont(self):
        try:
            del self.last_value
        except: pass
    def report(self):
        print("รายงานผลการพยากรณ์โรคโควิด")
        print("วันที่-เวลาพยากรณ์: {} {}"\
            .format(
                time.strftime("%d-%m-%Y", time.localtime()),
                time.strftime("%H:%M:%S", time.localtime()))
        )
        print('a = {:.2e}'.format(self.a))
        print('beta = {:.2e}'.format(self.b))
        print('พยากรณ์ตั้งแต่วันที่ {} จนถึงวันที่ {} ของการระบาด'.format(self.interval[0], self.interval[1]))
        print('วันที่ใช้คำนวณเงื่อนไขเริ่มต้น คือวันที่ {} ของการระบาด'.format(self.ds))
