#!/usr/bin/python



import sys
import re
pyver=sys.version_info[0]

if pyver==2:
  from ConfigParser import ConfigParser
if pyver==3:
  from configparser import ConfigParser

import matplotlib.pyplot as plot
  
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
        self.lengths.append(hex_to_time(re.search('\d\d', line).group()))
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
  def __init__(self, name, start_str, end_str, clocks, lengths,time_step=10*ns, valid_cmds=["VIDEO", "CLK2","CLK3"]):
    self.name=name
    self.start_str=start_str
    self.clocks=clocks
    self.end_str=end_str
    self.time_step=time_step
    self.clock_seq=[]
    self.lengths=lengths
    self.valid_cmds=valid_cmds
    
  def parse_seq(self, seq):
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
        for clock in clocks:
          step.append(clock.curr_val)
        self.clock_seq.append(step)        

  def plot(self):
    plt.figure()
    

        
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
  lengths.add_length("P_DELAY", 0.1*us)
  
  clocks=[]
  v1U=Clock("V1","TwoV1H", "TwoV1L",vhi, vlo)
  v2U=Clock("V2","TwoV2H", "TwoV2L",vhi, vlo)
  v3U=Clock("V3","TwoV3H", "TwoV3L",vhi, vlo)

  v1L=Clock("V1","OneV1H", "OneV1L",vhi, vlo)
  v2L=Clock("V2","OneV2H", "OneV2L",vhi, vlo)
  v3L=Clock("V3","OneV3H", "OneV3L",vhi, vlo)
  
  clocks.append(v1U)
  clocks.append(v2U)
  clocks.append(v3U)

  clocks.append(v1L)
  clocks.append(v2L)
  clocks.append(v3L)
  
  seq_v1=Sequence("Vertical","PARALLEL_1","PARALLEL_2",clocks, lengths)

  seq_v1.parse_seq(sequencer) 


  
