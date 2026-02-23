import os
import sys
import numpy

def floatmap(*args):
  return map(float,args)

class PbandPID:
  """PID that uses Proportional Band for scaling

N.B. Derivative action is not yet implemented

"""
  def __call__(self,PV):
    errorPct = 100.0 * (PV - self.SP) / (self.PVhi - self.PVlo)
    if self.reverseActing: errorPct *= -1.0

    if self.auto:
      self.biasPct += errorPct * (100.0/self.Pband) * (self.dTms/1E3) / self.Is
      CVPct = self.biasPct + (errorPct * 100.0/self.Pband)
      self.CV = self.CVlo + ((self.CVhi - self.CVlo) * CVPct/100.0)

    else:
      ### Bumpless transfer from Manual to Auto
      self.CV = self.CV0
      CVPct = 100.0 * (self.CV - self.CVlo) / (self.CVhi - self.CVlo)
      self.biasPct = CVPct - (errorPct * 100.0 / self.Pband)
      self.auto = True

    self.lastPV = PV
    self.CV0 = self.CV
    return self.CV

  def __init__(self, CV0=3.0, CVlo=0.0, CVhi=100.0
              , PVlo=0.0, PVhi=1E3, SP=400.0
              , Pband=90.0    ### Proportional Band, %
              , Is=15.0       ### Integral time, seconds
              , Ds=0.0        ### Derivative time, seconds
              , dTms=200.0    ### Loop update time, milliseconds
              , reverseActing=True
              , **kwargs
              ):
    self.CV0,self.CVlo,self.CVhi = floatmap(CV0,CVlo,CVhi)
    self.PVlo,self.PVhi,self.SP = floatmap(PVlo,PVhi,SP)
    self.Pband,self.Is,self.Ds,self.dTms = floatmap(Pband,Is,Ds,dTms)
    self.reverseActing = reverseActing
    self.auto = False

class SAWTOOTH:
  """Sawtooth process model"""
  def __call__(self,dTms=0.0):
    if self.V < self.LIMlo  : recalc,self.rising = True,True
    elif self.V > self.LIMhi: recalc,self.rising = True,False
    else                    : recalc = self.rate is None

    if recalc:
      if self.rising:
        self.rate = 1e-3 * (self.LIMhi - self.LIMlo) / self.RISEtime
      else:
        self.rate = 1e-3 * (self.LIMlo - self.LIMhi) / self.FALLtime

    self.V += self.rate * dTms

    return self.V
    
  def __init__(self
              ,V0=411.5
              ,LIMlo=388.5
              ,LIMhi=411.5
              ,FALLtime=56.065
              ,RISEtime=19.759
              ,rising=False
              ,**kwargs
              ):
    self.V,self.LIMlo,self.LIMhi = floatmap(V0,LIMlo,LIMhi)
    self.FALLtime,self.RISEtime = floatmap(FALLtime,RISEtime)
    self.rising = rising
    assert self.LIMlo < self.LIMhi
    assert self.FALLtime > 0.0
    assert self.RISEtime > 0.0
    self.rate = None

def run_system(**kwargs):
  pbpid = PbandPID(**kwargs)
  model = SAWTOOTH(**kwargs)
  dTms = pbpid.dTms
  PVs,CVs,Ts = list(),list(),list()
  for i in range(int(1 + (3600.0 / (dTms/1e3)))):
    Ts.append(i * dTms / 1e3)
    PVs.append(model(dTms=i and dTms or 0.0))
    CVs.append(pbpid(PVs[-1]))
  return Ts,PVs,CVs,pbpid.SP

def plot_system(Ts,PVs,CVs,SP):
  import matplotlib.pyplot as plt
  pvmin,pvmax = min(PVs),max(PVs)
  cvmin,cvmax = min(CVs),max(CVs)
  cvs = pvmin + (numpy.array(CVs) * (pvmax - pvmin) / 100.0)
  cvs2 = pvmin + ((numpy.array(CVs) - cvmin) * (pvmax - pvmin) / (cvmax - cvmin))
  ts = numpy.array(Ts) / 60.0
  plt.axhline(y=SP,color='k',ls=':')
  plt.plot(ts,PVs,'b',lw=1.5)
  plt.plot(ts,cvs,'g')
  plt.plot(ts,cvs2,'k:',lw=0.5)
  plt.show()
  

if "__main__" == __name__:
  plot_system(*run_system())
