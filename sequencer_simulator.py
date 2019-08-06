#!/usr/bin/python



import sys
import re
import numpy as np
pyver=sys.version_info[0]

if pyver==2:
  from ConfigParser import ConfigParser
if pyver==3:
  from configparser import ConfigParser

import matplotlib.pyplot as plt
  
#Units
ns=1.0
us=1000*ns

instr_length=40*ns

def hex_to_time(hexcode):
  code=format(int(hexcode,16), '0>8b')
  if (code[0]!='1'):
    time=int(code[1:],2)*40./1000 #40 nanoseconds units, in microseconds
  else:
    time=int(code[1:],2)*320./1000 #640 nanosecond units, in microseconds
  return time*us


class Lengths:
  def __init__(self):
    self.len_strs=[]
    self.lengths=[]
    
  def add_length(self, len_str, length):
    self.len_strs.append(len_str)
    self.lengths.append(length)

  def parse_length(self, len_str, sequencer):
    for line in sequencer:
      if re.match(len_str,line.lstrip()):
        self.len_strs.append(len_str)
        self.lengths.append(hex_to_time(re.search('[0-9a-fA-F][0-9a-fA-F]', line).group()))
        return

    print("Warning, delay " + len_str + " not found")
    return -1

  def parse_line(self, line):
    length=instr_length

    if len(self.lengths)<1:
      print("Warning, no lengths defined, returning default...")
      return length
    
    for i,len_str in enumerate(self.len_strs):
      if len_str in line:
        length+=self.lengths[i]
    return length
  
class Clock:
  def __init__(self,name,high_val_str, low_val_str, high_val, low_val, initial_val=None,RC=0, show=True):
    if initial_val is None:
      inital_val=low_val
    self.name=name
    self.high_val=high_val
    self.low_val=low_val
    self.high_val_str=high_val_str
    self.low_val_str=low_val_str
    self.show=show
    self.curr_val=initial_val
    self.RC=RC

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
  def __init__(self, name, start_str, end_str, clocks, lengths,time_step=10*ns, valid_cmds=["VIDEO", "CLK2","CLK3"], sequencer=None):
    self.name=name
    self.start_str=start_str
    self.clocks=clocks
    self.end_str=end_str
    self.time_step=time_step
    self.clock_seq=[]
    self.time_seq=[]
    self.lengths=lengths
    self.valid_cmds=valid_cmds
    if sequencer is not None:
      self.parse_seq(sequencer)
    
  def parse_seq(self, seq):
    self.clock_seq=[]
    #Find the starting location
    start_line=-1
    for i,line in enumerate(seq):
      if self.start_str in line:
        start_line=i
        break

    if start_line==-1:
      print("Error, start string for sequencer not found...")
      return -1

    for i,line in enumerate(seq[start_line:]):
      if self.end_str==line.strip():
        break
      #Skip whitespace and comments
      if line.lstrip()=='' or line.lstrip()[0]==";" or not any(x in line for x in self.valid_cmds):
        continue
      
      length=lengths.parse_line(line)
      steps=int(length/self.time_step)
      for clock in self.clocks:
        clock.parse_line(line)

      for j in range(steps):
        step=[]
        for clock in self.clocks:
          step.append(clock.curr_val)
        self.clock_seq.append(step)
    #Cast to a numpy array to make it easier to work with
    self.clock_seq=np.array(self.clock_seq)
    self.time_seq=np.arange(len(self.clock_seq))*self.time_step/us
  
  def plot(self):
    plt.figure()
    nplots=0
    for clock in self.clocks:
      if clock.show:
        nplots+=1
        
    currplot=1
    for i,clock in enumerate(self.clocks):
      if clock.show:
        plt.subplot(nplots,1,currplot, label=clock.name)
        plt.plot(self.time_seq, self.clock_seq[:,i], label=clock.name)
        plt.ylabel(clock.name)

        currplot+=1
    plt.xlabel("Time (us)")

        
if __name__=="__main__":
  if pyver==3:
    config=ConfigParser(inline_comment_prefixes=";")
  else:
    config=ConfigParser()
  config.read("Config.ini")

  vhi=config.getfloat("clocks","one_vclock_hi")
  vlo=config.getfloat("clocks","one_vclock_lo")

  hhi=config.getfloat("clocks","u_hclock_hi")
  hlo=config.getfloat("clocks","u_hclock_lo")

  
  with open('pit_super_sequencer_UW2_Paolo.waveforms') as f:
    sequencer=f.readlines()

    
  lengths=Lengths()
  lengths.parse_length("P_DELAY", sequencer)
  

  v12=Clock("V1_2","TwoV1H", "TwoV1L",vhi, vlo)
  v22=Clock("V2_2","TwoV2H", "TwoV2L",vhi, vlo)
  v32=Clock("V3_2","TwoV3H", "TwoV3L",vhi, vlo)
  tg2=Clock("TG_2","TwoTH", "TwoTL",vhi, vlo)
  
  v11=Clock("V1_1","OneV1H", "OneV1L",vhi, vlo)
  v21=Clock("V2_1","OneV2H", "OneV2L",vhi, vlo)
  v31=Clock("V3_1","OneV3H", "OneV3L",vhi, vlo)
  tg1=Clock("TG_1","OneTH", "OneTL",vhi, vlo)

  vclocks=[]
  vclocks.append(v11)
  vclocks.append(v21)
  vclocks.append(v31)
  vclocks.append(tg1)
  
  vclocks.append(v12)
  vclocks.append(v22)
  vclocks.append(v32)
  vclocks.append(tg2)
  
  seq_v1=Sequence("Vertical 1","PARALLEL_1","PARALLEL_2",vclocks, lengths, sequencer=sequencer)
  seq_v2=Sequence("Vertical 2","PARALLEL_2","PARALLEL_12",vclocks, lengths, sequencer=sequencer)
  seq_v12=Sequence("Vertical 12","PARALLEL_12","END_PARALLEL",vclocks, lengths, sequencer=sequencer)
  

  

  
