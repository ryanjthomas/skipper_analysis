#!/usr/bin/ipython -i

import sys
import re
#To deep copy our classes
import copy
import numpy as np
pyver=sys.version_info[0]

if pyver==2:
  from ConfigParser import ConfigParser
if pyver==3:
  from configparser import ConfigParser

import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, lfilter_zi

#I like interactive plotting
plt.ion()
#Units
ns=1.0
us=1000*ns

#RC constants for our clocks
RC_vclock = 2.35*us
RC_hclock = 33*ns
RC_sclock = 100*ns  

instr_delay=40*ns

def hex_to_time(hexcode):
  code=format(int(hexcode,16), '0>8b')
  if (code[0]!='1'):
    time=int(code[1:],2)*40./1000 #40 nanoseconds units, in microseconds
  else:
    time=int(code[1:],2)*320./1000 #640 nanosecond units, in microseconds
  return time*us

def clean_line(line):
  #Only take stuff up to the comment, strip leading whitespace
  return line.split(";")[0].lstrip()

def lowpass_filter(x,RC, step):
  '''
  Wrapper to perform low pass filtering on data x. Step is the length of time between samples in the data series.
  Assumes RC and step are in the same time units.
  '''
  b,a=butter(1,1./(np.pi*RC/step))
  return lfilter(b,a,x,zi=lfilter_zi(b,a)*x[0])[0]             

class Delays:
  def __init__(self):
    self.len_strs=[]
    self.delays=[]
    
  def add_delay(self, len_str, delay):
    self.len_strs.append(len_str)
    self.delays.append(delay)

  def parse_delay(self, len_str, sequencer):
    for line in sequencer:
      line=clean_line(line)
      if re.match(len_str,line):
        self.len_strs.append(len_str)
        self.delays.append(hex_to_time(re.search('\$[0-9a-fA-F][0-9a-fA-F]', line).group()[1:]))
        return

    print("Warning, delay " + len_str + " not found")
    return -1

  def parse_line(self, line):
    delay=instr_delay
    line=clean_line(line)
    if len(self.delays)<1:
      print("Warning, no delays defined, returning default...")
      return delay
    
    for i,len_str in enumerate(self.len_strs):
      if len_str in line:
        delay+=self.delays[i]
    return delay
  
class Clock:
  def __init__(self,name,high_val_str, low_val_str, high_val, low_val, config=None,RC=1, initial_val=None, show=True, color=None, linestyle=None):
    if config is not None:
      low_val=config.getfloat("clocks",low_val)
      high_val=config.getfloat("clocks",high_val)
    self.name=name
    self.high_val=high_val
    self.low_val=low_val
    self.high_val_str=high_val_str
    self.low_val_str=low_val_str
    self.show=show
    if initial_val is not None:
      self.curr_val=initial_val
    else:
      self.curr_val=high_val
    self.RC=RC
    if color is not None:
      self.color=color
    else:
      self.color='b'
    if linestyle is not None:
      self.linestyle=linestyle
    else:
      self.linestyle='solid'
      
      
  def parse_line(self, line):
    if self.high_val_str in line:
      self.curr_val=self.high_val
      return self.low_val
    elif self.low_val_str in line:
      self.curr_val=self.low_val
      return self.high_val
    else:
      return self.curr_val

class Sequence:
  def __init__(self, name, start_str, end_str, clocks, delays,time_step=10*ns, valid_cmds=["VIDEO", "CLK2","CLK3"], sequencer=None):
    self.name=name
    self.start_str=start_str
    self.clocks=clocks
    self.end_str=end_str
    self.time_step=time_step
    self.clock_seq=[]
    self.time_seq=[]
    self.delays=delays
    self.valid_cmds=valid_cmds
    if sequencer is not None:
      self.parse_seq(sequencer)
    
  def parse_seq(self, seq):
    if self.start_str==None:
      print("Error, no start string, cannot parse sequencer...")
      return -1

    self.clock_seq=[]
    #Find the starting location
    start_line=-1
    for i,line in enumerate(seq):
      line=clean_line(line)
      if self.start_str==line.strip():
        start_line=i
        break

    if start_line==-1:
      print("Error, start string for sequencer not found...")
      return -1

    for i,line in enumerate(seq[start_line:]):
      if self.end_str==line.strip():
        break
      #Skip whitespace and comments
      line=clean_line(line)
      #Maybe second condition not necessary after above
      if line.lstrip()=='' or line.lstrip()[0]==";" or not any(x in line for x in self.valid_cmds):
        continue
      delay=delays.parse_line(line)
      steps=int(delay/self.time_step)
      for clock in self.clocks:
        clock.parse_line(line)

      for j in range(steps):
        step=[]
        for clock in self.clocks:
          step.append(clock.curr_val)
        self.clock_seq.append(step)
    #Cast to a numpy array to make it easier to work with
    self.clock_seq=np.array(self.clock_seq)
    #Filter our clock sequence through the low-pass filter
    self.filter_clocks()
    self.time_seq=np.arange(len(self.clock_seq))*self.time_step

  def filter_clocks(self):
    fclocks=[]
    for i,clock in enumerate(self.clock_seq.T):
      fclocks.append(lowpass_filter(clock, self.clocks[i].RC, self.time_step))
    fclocks=np.array(fclocks)
    self.clock_seq_filtered=fclocks.T
        
  def show(self, *args, **kwargs):
    self.plot(*args, **kwargs)
    
  def plot(self, filtered=True, unit=us,legend=True):
    '''
    Plots the sequence. By default shows the sequencer after RC filtering, use filtered=False to show the raw sequence. Can specify units of either us or ns. Setting legend true displays the plot label as a legend, otherwise it's shown on the y axis.
    '''
    fig=plt.figure()
    nplots=0
    for clock in self.clocks:
      if clock.show:
        nplots+=1
        
    currplot=1
    for i,clock in enumerate(self.clocks):
      if clock.show:
        plt.subplot(nplots,1,currplot, label=clock.name)
        if filtered:
          data=self.clock_seq_filtered[:,i]
        else:
          data=self.clock_seq[:,i]
        plt.plot(self.time_seq/unit, data, label=clock.name, color=clock.color, linestyle=clock.linestyle)
        plt.tick_params(labelbottom=False)
        if legend:
          plt.legend()
        else:
          plt.ylabel(clock.name)
          
        currplot+=1
    if unit==us:
      plt.xlabel("Time (us)")
    elif unit==ns:
      plt.xlabel("Time (ns)")
    else:
      print("Warning, unit not recognized...")
      plt.xlabel("Time")
    #Enable bottom labels only for our bottom plot
    plt.tick_params(labelbottom=True)
    fig.suptitle(self.name)
    plt.show()
    return fig
      
  def __add__(self, b):
    '''
    We can add two sequences together as long as they share the same clocks. Just appends the filtered and unfiltered clock and time sequence.
    '''
    #First check our clocks are the same
    if self.clocks!=b.clocks:
      raise ValueError("Cannot add sequences with different sets of clocks")
    if self.time_step!=b.time_step:
      #TODO: somehow adapt time steps of different sequencers to work together
      #For now just don't combine sequencers with different time steps (no reason to anyways)
      print("Warning, time steps of sequences are not equal, filtered clocks will not be correct...")
    newseq=copy.copy(self)
    newseq.clock_seq=np.append(self.clock_seq, b.clock_seq,axis=0)
    #This needs to be regenerated to properly stich together clock transitions
    #    newseq.clock_seq_filtered=np.append(self.clock_seq_filtered, b.clock_seq_filtered,axis=0)
    newseq.time_seq=np.append(self.time_seq, b.time_seq+self.time_seq[-1],axis=0)
    newseq.name=self.name+"+"+b.name
    #These no longer make sense
    newseq.start_str=None
    newseq.end_str=None
    newseq.filter_clocks()
    return newseq
    
if __name__=="__main__":
  if len(sys.argv)>1:
    cfname=sys.argv[1]
  else:
    cfname="Config.ini"
  if len(sys.argv)>2:
    sfname=sys.argv[2]
  else:
    sfname="pit_super_sequencer_UW2.waveforms"
    
  if pyver==3:
    config=ConfigParser(inline_comment_prefixes=";")
  else:
    config=ConfigParser()
  config.read(cfname)

  int_delay=config.getfloat("timing","IntegralTime")
  ped_delay=config.getfloat("timing","PedestalIntgWait")
  sig_delay=config.getfloat("timing","SignalIntgWait")
  sw_delay=config.getfloat("timing","SWPulseWidth")
  og_delay=config.getfloat("timing","OGWidth")
  dg_delay=config.getfloat("timing","DGWidth")
  rg_delay=config.getfloat("timing","SkippingRGWidth")
  
  with open(sfname) as f:
    sequencer=f.readlines()
    
  delays=Delays()
  delays.parse_delay("P_DELAY", sequencer)
  delays.parse_delay("R_DELAY", sequencer)
  delays.parse_delay("S_DELAY", sequencer)
  delays.parse_delay("SW_DELAY", sequencer)
  delays.parse_delay("VE_DELAY", sequencer)
  delays.parse_delay("DCRST_DELAY", sequencer)

  delays.add_delay("RDL0",rg_delay*us)
  delays.add_delay("PRD0",ped_delay*us)
  delays.add_delay("DLY0",int_delay*us)
  delays.add_delay("SWD0",sw_delay*us)
  delays.add_delay("POD0",sig_delay*us)
  delays.add_delay("DLY1",int_delay*us)
  delays.add_delay("OGD0",og_delay*us)
  delays.add_delay("DDL0",dg_delay*us)
  
  v12=Clock("V1_2","TwoV1H", "TwoV1L","two_vclock_hi","two_vclock_lo",config=config, RC=RC_vclock)
  v22=Clock("V2_2","TwoV2H", "TwoV2L","two_vclock_hi","two_vclock_lo",config=config, RC=RC_vclock)
  v32=Clock("V3_2","TwoV3H", "TwoV3L","two_vclock_hi","two_vclock_lo",config=config, RC=RC_vclock)
  tg2=Clock("TG_2","TwoTH", "TwoTL","tg_hi","tg_lo",config=config, RC=RC_vclock, color='red')
  
  v11=Clock("V1_1","OneV1H", "OneV1L","one_vclock_hi","one_vclock_lo",config=config, RC=RC_vclock)
  v21=Clock("V2_1","OneV2H", "OneV2L","one_vclock_hi","one_vclock_lo",config=config, RC=RC_vclock)
  v31=Clock("V3_1","OneV3H", "OneV3L","one_vclock_hi","one_vclock_lo",config=config, RC=RC_vclock)
  tg1=Clock("TG_1","OneTH", "OneTL","tg_hi","tg_lo",config=config, RC=RC_vclock, color='red')

  vclocks=[]
  vclocks.append(v11)
  vclocks.append(v21)
  vclocks.append(v31)
  vclocks.append(tg1)
  
  vclocks.append(v12)
  vclocks.append(v22)
  vclocks.append(v32)
  vclocks.append(tg2)

  #Vertical clock sequences
  seq_v1=Sequence("Vertical 1","PARALLEL_1","PARALLEL_2",vclocks, delays, sequencer=sequencer)
  seq_v2=Sequence("Vertical 2","PARALLEL_2","PARALLEL_12",vclocks, delays, sequencer=sequencer)
  seq_v12=Sequence("Vertical 12","PARALLEL_12","END_PARALLEL",vclocks, delays, sequencer=sequencer)


  hu3=Clock("HU3","HU3H", "HU3L","u_hclock_hi","u_hclock_lo",config=config, RC=RC_hclock)
  hu2=Clock("HU2","HU2H", "HU2L","u_hclock_hi","u_hclock_lo",config=config, RC=RC_hclock)
  hu1=Clock("HU1","HU1H", "HU1L","u_hclock_hi","u_hclock_lo",config=config, RC=RC_hclock)

  hl3=Clock("HL3","HL3H", "HL3L","l_hclock_hi","l_hclock_lo",config=config, RC=RC_hclock)
  hl2=Clock("HL2","HL2H", "HL2L","l_hclock_hi","l_hclock_lo",config=config, RC=RC_hclock)
  hl1=Clock("HL1","HL1H", "HL1L","l_hclock_hi","l_hclock_lo",config=config, RC=RC_hclock)

  hclocks=[hu1,hu2,hu3,hl1,hl2,hl3]

  seq_hl=Sequence("Read L", "SERIAL_READ_L_STAGE1","END_SERIAL_READ_L_STAGE1", hclocks, delays, sequencer=sequencer)
  seq_hu=Sequence("Read U", "SERIAL_READ_R_STAGE1","END_SERIAL_READ_R_STAGE1", hclocks, delays, sequencer=sequencer)
  seq_hul=Sequence("Read UL", "SERIAL_READ_LR_STAGE1","END_SERIAL_READ_LR_STAGE1", hclocks, delays, sequencer=sequencer)
  
  rg=Clock("RG","RH", "RL","rg_hi","rg_lo",config=config, RC=RC_sclock,color="blue")
  og=Clock("OG","OUTGH", "OUTGL","og_hi","og_lo",config=config, RC=RC_sclock,color="magenta")
  sw=Clock("SW","WH", "WL","sw_hi","sw_lo",config=config, RC=RC_sclock,color="red")
  dg=Clock("DG","DUMPGH","DUMPGL", "dg_hi","dg_lo", config=config, RC=RC_hclock,color="cyan")

  sclocks=[hu3,rg,og,sw,dg]
  
  seq_skip=Sequence("Skipper Read","PIT_SK_NDCR_SERIAL_READ","END_PIT_SK_NDCR_SERIAL_READ",
                    sclocks, delays, sequencer=sequencer)

  seq_clrchrg=Sequence("Clear Charge","SERIAL_READ_CLRCHG_STAGE_2","END_SERIAL_READ_CLRCHG_STAGE_2",
                       sclocks, delays, sequencer=sequencer)

  allhclocks=[hu1,hu2,hu3,hl1,hl2,hl3,sw,og,rg,dg]

  seq_hor=Sequence("Read UL", "SERIAL_READ_LR_STAGE1","END_SERIAL_READ_LR_STAGE1", allhclocks, delays, sequencer=sequencer)

  seq_read=Sequence("Skipper Read","PIT_SK_NDCR_SERIAL_READ","END_PIT_SK_NDCR_SERIAL_READ",
                    allhclocks, delays, sequencer=sequencer)

  seq_done=Sequence("Clear Charge","SERIAL_READ_CLRCHG_STAGE_2","END_SERIAL_READ_CLRCHG_STAGE_2",
                    allhclocks, delays, sequencer=sequencer)

  seq_full=seq_hor+seq_read+seq_read+seq_done


  seq_full.plot()
