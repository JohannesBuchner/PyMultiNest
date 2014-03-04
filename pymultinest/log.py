# Back-up MultNest output at regular intervals.

from __future__ import absolute_import, unicode_literals, print_function
import threading
import os

class backUp(threading.Thread):
  """
  Class for periodically backing up the MultiNest output.
  Creates new folders with date-time name with backups.
  """
  def __init__(self, interval_hr = 1, outputfiles_basename = "chains/1-"):
    """
    Find all data file names. Set running to True.
    """
    threading.Thread.__init__(self)
    self.outputfiles_basename = outputfiles_basename
    self.interval = interval_hr*60*60 # Convert from hours to seconds.
    self.running = True
    
    # List of file suffixes for the data files to be backed up.
    self.suffix = [".txt",
              "post_separate.dat",
              "stats.dat",
              "post_equal_weights.dat",
              "phys_live.points",
              "live.points",
              "ev.dat",
              "resume.dat",
              "summary.txt"
              ]
    
    # Make list of file names.
    self.path, self.base = os.path.split(self.outputfiles_basename)
    self.path = os.path.abspath(self.path)
    self.file_name = [ "%s%s" % (self.base , suffix) for suffix in self.suffix ]  
    self.full_name = [ "%s%s" % (self.outputfiles_basename , suffix) for suffix in self.suffix ] 
  
  def stop(self):
    """
    Stop backing up; set running to False. Set to to True in __init__.
    """
    self.running = False


  def run(self):
    """ 
    Back up the data files every set length of time.
    """
    import shutil
    import time
          
    while self.running:

      if not self.running: break # break before sleeping!
      time.sleep(self.interval) # in seconds, but interval already converted from hours. 
            
      # Make backup directory.
      self.backup = os.path.join(self.path,"%4d-%02d-%02d-%02d-%02d-%02d" % time.localtime()[0:6])
      if not os.path.exists(self.backup):
        os.makedirs(self.backup)
        
      # Make sure that backup was successful.
      while not self.checkBackUp():
        
        print("Attempting to backup files.")
      
        # Copy files to a new backup.
        for name in self.file_name:

          # Join file names with paths.
          orig = os.path.join(self.path,name)
          backup = os.path.join(self.backup,name)
          
          if os.path.exists(orig):
            shutil.copy2(orig,backup)

              
  def checkBackUp(self):
    """
    Check that a backup was successful. i.e. that it the file-lenghts
    are consistent with the [root]resume.dat.
    """
    
    # Relevant files.
    self.resume = os.path.join(self.backup,self.base+"resume.dat")
    self.rejected = os.path.join(self.backup,self.base+"ev.dat")
    self.live = os.path.join(self.backup,self.base+"live.points")
    self.physlive = os.path.join(self.backup,self.base+"phys_live.points")

    # Check that backups exists.
    for name in self.file_name:
      if not os.path.exists(os.path.join(self.backup,name)): return False

    # Pick relevant line from resume file, the second line, and read contents.
    self.nlive = int(open(self.resume).readlines()[1].split()[3])
    self.nresume = int(open(self.resume).readlines()[1].split()[0]) - self.nlive 

    # Check that file lengths are correct.
    if len(file(self.rejected, 'r').readlines()) != self.nresume: return False
    if len(file(self.live, 'r').readlines()) != self.nlive: return False
    if len(file(self.physlive, 'r').readlines()) != self.nlive: return False

    return True