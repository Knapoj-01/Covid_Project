import time
import numpy as np 
from scipy.special import lambertw
from .rungeKutta import RungeKuttaSolve
class SIR_Model:
    #Class Constructor
    def __init__(self, data, N):
        self.day = np.copy(data.iloc[:,0])
        self.I = np.copy(data.iloc[:,1])
        self.R = np.copy(data.iloc[:,2] + data.iloc[:,3])
        self.S = np.ones(len(self.I))*N - self.I -self.R
        self.N = N
    #Initial Conditions
    def computeIC(self, tr, ds, dt = 1):
        b = 1/tr
        a = (self.I[ds +dt] - self.I[ds])/(self.S[ds]*self.I[ds]*dt) + b/self.S[ds] 
        self.rho = b/a
        self.a = a
        self.b = b
        self.ic = self.I[ds] + self.S[ds] -self.rho*np.log(self.S[ds])
        self.ds = ds
        return {"a": a, "b":b, "rho": self.rho}
    def setIC(self,a,b,ds):
        self.rho = b/a
        self.a = a
        self.b = b
        self.ds = ds
        self.ic = self.I[ds] + self.S[ds] -self.rho*np.log(self.S[ds])
    #Internal Functions
    def iprime(self, t,I):
        #if I > -self.rho + self.rho*np.log(self.rho) + self.ic: return 0
        S_If = -self.rho*lambertw(-np.exp((I-self.ic)/self.rho)/self.rho,self.br)
        if np.imag(S_If) != 0:
            self.br = 0
        slope = np.real(self.a*S_If*I - self.b*I)
        return slope
    #Predict Data
    def predictI(self, npoints,day):
        self.br = -1
        x,y = RungeKuttaSolve(
            self.iprime, day, self.ds,
            self.I[self.ds], npoints
            )
        return x,y

    def getData(self, kwd): 
        data = {
            "susceptibles" : self.S,
            "infecteds" : self.I,
            "removeds" : self.R,
            "day": self.day
        }
        return data.get(kwd,'Invalid')
    def shiftData(self, kwd, range, amount):
        try:
            if len(range) == 1:
                self.getData(kwd)[range[0]:] =  self.getData(kwd)[range[0]:] + amount
            elif len(range) == 2:
                self.getData(kwd)[range[0]:range[1]] =  self.getData(kwd)[range[0]:range[1]] + amount
        except: 
            pass
        self.S = np.ones(len(self.I))*self.N - self.I -self.R
    def report(self):
        print("รายงานผลการพยากรณ์โรคโควิด")
        print("วันที่-เวลาพยากรณ์: {} {}"\
            .format(
                time.strftime("%d-%m-%Y", time.localtime()),
                time.strftime("%H:%M:%S", time.localtime()))
        )
        print('a = {:.2e}'.format(self.a))
        print('b = {:.2e}'.format(self.rho*self.a))
        print('rho = {:.2e}'.format(self.rho))
        print('ic = {:.2e}'.format(self.ic))
        print('วันที่ใช้คำนวณเงื่อนไขเริ่มต้น คือวันที่ {} ของการระบาด'.format(self.ds))

class SIRD_Model(SIR_Model):
    def __init__(self, data, N):
        super().__init__(data, N)

    def computeIC(self, tr,k,  ds, dt = 1):
        b = 1/tr
        a = (self.I[ds +dt] - self.I[ds])/(self.S[ds]*self.I[ds]*dt) + b/self.S[ds] 
        self.rho = b/a
        self.a = a
        self.b = b
        self.k = k
        self.ic = k *self.I[ds] + self.S[ds] -self.rho*np.log(self.S[ds])
        self.ds = ds
        return {"a": a, "b":b,"k": k ,"rho": self.rho}
    def iprime(self, t,I):
        S_If = -self.rho*lambertw(-np.exp((self.k*I-self.ic)/self.rho)/(self.rho),self.br)
        if np.imag(S_If) != 0:
            self.br = 0
        slope = self.a*S_If*I - self.b*I
        return np.real(slope)
    def predictI(self, npoints,day):
        self.br = -1
        x,y = RungeKuttaSolve(
            self.iprime, day, self.ds,
            self.I[self.ds], npoints
            )
        return x,y
   