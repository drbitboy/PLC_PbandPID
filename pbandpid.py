import os
os.environ["QT_XCB_GL_INTEGRATION"] = "none"
import sys
import numpy
from inspect import signature

def floatmap(*args):
  return map(float,args)

class PbandPID:
  """PID that uses Proportional Band for dependent gains

N.B. Derivative action is not yet implemented

"""
  def __call__(self,PV):
    """Loop update; bumpless switch to Auto if currently in Manual"""
    errorPct = 100.0 * (PV - self.SP) / (self.PVhi - self.PVlo)
    if self.reverseActing: errorPct *= -1.0

    if self.auto:
      ### Loop update in Aut
      self.biasPct += errorPct * (100.0/self.Pband) * (self.dTms/1E3) / self.Is
      CVPct = self.biasPct + (errorPct * 100.0/self.Pband)
      self.CV = self.CVlo + ((self.CVhi - self.CVlo) * CVPct/100.0)

    else:
      ### Bumpless transfer from Manual to Auto
      self.CV = self.CV0
      CVPct = 100.0 * (self.CV - self.CVlo) / (self.CVhi - self.CVlo)
      self.biasPct = CVPct - (errorPct * 100.0 / self.Pband)
      self.auto = True

    self.lastPV = PV    ### Used for D-term (eventually)
    self.CV0 = self.CV
    return self.CV

  def __init__(self, CV0=2.0, SP=400.0  ### Initial state
              , PVlo=0.0, PVhi=1E3      ### PV/SP input scaling
              , CVlo=0.0, CVhi=100.0    ### CV output scaling
              , Pband=90.0              ### Proportional Band, %
              , Is=15.0                 ### Integral time, seconds
              , Ds=0.0                  ### Derivative time, seconds
              , dTms=200.0              ### Loop update time, milliseconds
              , reverseActing=True
              , **kwargs
              ):
    """Construct PID in Manual mode"""
    self.CV0,self.CVlo,self.CVhi = floatmap(CV0,CVlo,CVhi)
    self.PVlo,self.PVhi,self.SP = floatmap(PVlo,PVhi,SP)
    self.Pband,self.Is,self.Ds,self.dTms = floatmap(Pband,Is,Ds,dTms)
    self.reverseActing = reverseActing
    self.auto,self.CV = False,self.CV0

class oldSAWTOOTH:
  """Sawtooth process model, independent of valve position"""
  def __call__(self,dTms=0.0,**kwargs):
    """Calculate one loop update"""
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
              ,V0=411.5         ### Initial model state value
              ,LIMlo=388.5      ### Minimum model value
              ,LIMhi=411.5      ### Maximum model value
              ,FALLtime=56.065  ### Time from max value to min value
              ,RISEtime=19.759  ### Time from min value to max value
              ,rising=False     ### Whether initial state is rising
              ,**kwargs
              ):
    """Construct sawtooth model"""
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
    """Calculate one loop update based on CV position"""
    ### Implement sticky (jumping) CV behavior
    CVminus,CVplus = CV-self.CVjump, CV+self.CVjump
    while self.CV0 < CVminus: self.CV0 += self.CVjump
    while self.CV0 > CVplus: self.CV0 -= self.CVjump

    self.V += (self.CV0-self.CVnull) * self.ratePerCV * (dTms/1000.0)

    if self.PVsmooth: return self.V,self.CV0

    return round(self.V,0),self.CV0

  def __init__(self
              ,PVsmooth=False                ### Return non-rounded values
              ,V0=411.5                      ### Initial value
              ,ratePerCV=round(ratePerCV,3)  ### Rise/fall rate per CV unit
              ,CV0=2.0                       ### Initial CV
              ,CVjump=2.1                    ### Size of CV jumps
              ,CVnull=round(CVnull,3)        ### CV for steady-state PV
              ,**kwargs
              ):
    """Construct sawtooth model"""
    self.PVsmooth = PVsmooth
    self.V,self.ratePerCV = floatmap(V0,ratePerCV)
    self.CV0,self.CVjump,self.CVnull = floatmap(CV0,CVjump,CVnull)
    assert self.CVjump > 0.0

def run_system(**kwargs):
  """Model 1h of data"""

  ### Construct PID and proces model
  pbpid = PbandPID(**kwargs)
  model = (kwargs.get("oldSAWTOOTH") and oldSAWTOOTH or SAWTOOTH)(**kwargs)

  ### Initialize loop parameters
  dTms = pbpid.dTms
  PVs,CVs,CVstickies,Ts = list(),list(),list(),list()

  ### Run loop for 1h
  for i in range(int(1 + (3600.0 / (dTms/1e3)))):
    ### Calculate current time
    Ts.append(i * dTms / 1e3)
    ### Run process model; N.B. CV is ignored in oldSAWTOOTH model
    PV = model(dTms=i and dTms or 0.0,CV=pbpid.CV)
    try:
      ### Assume newer sawtooth model returned tuple (PV,CVsticky,)
      ### - append CVsticky; extract PV
      CVstickies.append(PV[1])
      PV = PV[0]
    except: pass    ### oldSAWTOOTH returns only PV

    ### Append new PV, append CV from PID with new PV
    PVs.append(PV)
    CVs.append(pbpid(PVs[-1]))

  return Ts,PVs,CVs,CVstickies,pbpid.SP,model

def plot_system(Ts,PVs,CVs,CVstickies,SP,model):
  import matplotlib.pyplot as plt
  pvmin,pvmax = min(PVs),max(PVs)
  ### Scale CVs [0,100] range to PV range, to approximate trend from OP
  ### Convert times to minutes
  cvs = pvmin + (numpy.array(CVs) * (pvmax - pvmin) / 100.0)
  ts = numpy.array(Ts) / 60.0

  fig,ax1 = plt.subplots()

  ### ax1 uses left axis; plot PVs and scaled CVs
  ax1.axhline(y=SP,color='k',ls=':',label='SP')
  ax1.plot(ts,PVs,'b',lw=1.5,label='PV')
  ax1.plot(ts,cvs,'g',label='CVorig')
  ax1.legend(loc='upper left')

  ### Create right axis for unscaled CVs
  ax2 = ax1.twinx()

  ### Plot CVs and, if present, the sticky CV values
  ax2.plot(ts,CVs,'k',label='CV>')
  if len(CVstickies) == len(ts):
    ax2.plot(ts,CVstickies,'m',label='CVstem>')

  ### Annotate plot
  ax2.legend(loc='lower right')
  plt.title(f"{type(model).__name__}{signature(type(model))}")
  ax1.set_xlabel('Time, minutes')
  ax1.set_ylabel('PV, PSI')
  ax2.set_ylabel('CV, %')

  plt.show()


if "__main__" == __name__:
  ### Parse command line to dictionary
  kwargs = dict()
  for arg in sys.argv[1:]:
    toks = arg.split('=')
    tok0 = toks.pop(0)
    L = len(toks)
    if L == 0: val = tok0[:1] != '-' ### Bool:  keyword => T; -keyword => F
    else     : val = '='.join(toks)  ### String value
    kwargs[tok0.lstrip('-')] = val

  if kwargs.get('help',False):
    print(f"oldSAWTOOTH class:\n  {signature(oldSAWTOOTH)}")
    print(f"SAWTOOTH class:\n  {signature(SAWTOOTH)}")
    print(f"PbandPDI class:\n  {signature(PbandPID)}")
  else:
    plot_system(*run_system(**kwargs))
