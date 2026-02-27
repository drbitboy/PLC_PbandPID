import os
os.environ["QT_XCB_GL_INTEGRATION"] = "none"
import sys
import numpy
from inspect import signature

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

  def __init__(self, CV0=2.0, CVlo=0.0, CVhi=100.0
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
    self.auto,self.CV = False,self.CV0

class oldSAWTOOTH:
  """Sawtooth process model"""
  def __call__(self,dTms=0.0,**kwargs):
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

CVnull=2.0 + (2.1 * 19.759 / (19.759 + 56.065))
ratePerCV=(23.0/56.065)/(CVnull-2.0)

class SAWTOOTH:
  """Sawtooth process model - accumulation based on sticky CV"""
  def __call__(self,dTms=0.0,CV=None,**kwargs):
    CVminus,CVplus = CV-self.CVjump, CV+self.CVjump
    while self.CV0 < CVminus: self.CV0 += self.CVjump
    while self.CV0 > CVplus: self.CV0 -= self.CVjump

    self.V += (self.CV0-self.CVnull) * self.ratePerCV * (dTms/1000.0)

    return self.V,self.CV0

  def __init__(self
              ,V0=411.5
              ,LIMlo=388.5
              ,LIMhi=411.5
              ,CV0=2.0
              ,CVjump=2.1
              ,CVnull=round(CVnull,3)
              ,ratePerCV=round(ratePerCV,3)
              ,**kwargs
              ):
    self.V,self.LIMlo,self.LIMhi = floatmap(V0,LIMlo,LIMhi)
    self.CV0,self.CVjump,self.CVnull = floatmap(CV0,CVjump,CVnull)
    assert self.LIMlo < self.LIMhi
    assert self.CVjump > 0.0
    self.ratePerCV,_ = floatmap(ratePerCV,0.0)

def run_system(**kwargs):
  pbpid = PbandPID(**kwargs)
  if kwargs.get("oldSAWTOOTH",False):
    model = oldSAWTOOTH(**kwargs)
  else:
    model = SAWTOOTH(**kwargs)
  dTms = pbpid.dTms
  PVs,CVs,CVstickies,Ts = list(),list(),list(),list()
  for i in range(int(1 + (3600.0 / (dTms/1e3)))):
    Ts.append(i * dTms / 1e3)
    PV = model(dTms=i and dTms or 0.0,CV=pbpid.CV)
    try:
      CVstickies.append(PV[1])
      PV = PV[0]
    except: CVstem = pbpid.CV
    PVs.append(PV)
    CVs.append(pbpid(PVs[-1]))

  return Ts,PVs,CVs,CVstickies,pbpid.SP,model

def plot_system(Ts,PVs,CVs,CVstickies,SP,model):
  import matplotlib.pyplot as plt
  pvmin,pvmax = min(PVs),max(PVs)
  cvs = pvmin + (numpy.array(CVs) * (pvmax - pvmin) / 100.0)
  ts = numpy.array(Ts) / 60.0

  fig,ax1 = plt.subplots()

  ax1.axhline(y=SP,color='k',ls=':',label='SP')
  ax1.plot(ts,PVs,'b',lw=1.5,label='PV')
  ax1.plot(ts,cvs,'g',label='CVorig')
  ax1.legend(loc='upper left')

  ax2 = ax1.twinx()
  ax2.plot(ts,CVs,'k',label='CV>')
  if len(CVstickies) == len(ts):
    ax2.plot(ts,CVstickies,'m',label='CVstem>')

  ax2.legend(loc='lower right')

  plt.title(f"({signature(type(model))}")
  plt.show()


if "__main__" == __name__:
  kwargs = dict()
  for arg in sys.argv[1:]:
    toks = arg.split('=')
    tok0 = toks.pop(0)
    L = len(toks)
    if L == 0: val = tok0[:1] != '-'
    else     : val = '='.join(toks)
    kwargs[tok0.lstrip('-')] = val

  if kwargs.get('help',False):
    print(f"oldSAWTOOTH class:\n  {signature(oldSAWTOOTH)}")
    print(f"SAWTOOTH class:\n  {signature(SAWTOOTH)}")
    print(f"PbandPDI class:\n  {signature(PbandPID)}")
  else:
    plot_system(*run_system(**kwargs))
