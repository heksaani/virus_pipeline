"""Watcher module for Virus Piper pipeline
This module handles:
1) Detect that a folder appeared

2) Identify which parent platform directory it belongs to

3) Trigger the pipeline associated with that platform
"""

import logging
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from .fastq_files import extract_barcode, get_sample_name
#from .nextseq
#from .minion
#from .miseq


class DirectoryMonitor(FileSystemEventHandler):
  """Handler for file system events"""

  def on_created(self, event):
    """Called when a file or directory is created"""
    if event.is_directory:
      logging.info(f"New directory created: {event.src_path}")
      
    else:
      path_str  = str(event.src_path)
      platform = path_str.split('/')[4] # Assumes platform is the 5th element in the path set in the config
      sample_reads = {}
      if path_str.endswith('.fastq.gz'):
        logging.info(f"Fastq file detected: {event.src_path}")
        if platform not in ['minion', 'nextseq','miseq']:
          logging.error(f"Invalid platform detected {platform}")
        elif platform == 'minion':
          barcode = extract_barcode(path_str)
          logging.info(f"Barcode : {barcode}")
          #combine_minion_reads(path_str)
          #lauch_irma()
        elif platform == 'nextseq':
          # jos data on sekvensoitu NextSeqillä (= kaikki fastq samassa kansiossa) 
          sample_name = get_sample_name(path_str)
          if sample_name not in sample_reads.keys():
             sample_reads[sample_name] = {"R1": None, "R2": None}
          if "R1" in path_str:
            sample_reads[sample_name]["R1"] = path_str
            logging.info(f"Sample name: {sample_name} R1 added to dictionary")
          elif "R2" in path_str:
            sample_reads[sample_name]["R2"] = path_str
            logging.info(f"Sample name: {sample_name} R2 added to dictionary")
          else:
            logging.error("Invalid fastq filename")
          elif platform == 'miseq':
            sample_name = get_sample_name(path_str)
            if sample_name not in sample_reads.keys():
               sample_reads[sample_name] = {"R1": None, "R2": None}
            if "R1" in path_str:
              sample_reads[sample_name]["R1"] = path_str
              logging.info(f"Sample name: {sample_name} R1 added to dictionary")
            elif "R2" in path_str:
              sample_reads[sample_name]["R2"] = path_str
              logging.info(f"Sample name: {sample_name} R2 added to dictionary")
            else:
              logging.error("Invalid fastq filename")
      #lauch_irma()

         


def main(path:str, recursive: bool = True):
  """
  Start monitoring the specified path.
    
  Args:
        path: Directory path to monitor
        recursive: Whether to monitor subdirectories"""
  
  logging.info('Starting watcher script...')
  event_handler = DirectoryMonitor()
  observer = Observer()
  observer.schedule(event_handler, path, recursive=True)
  logging.info(f"Monitoring started on path: {path}")

  observer.start()
  logging.info("Observer started.")
  try:
     while observer.is_alive():
        observer.join(1) # Wait 1 second, check if alive, repeat
  except KeyboardInterrupt:
    logging.info("Stopping observer...")
    observer.stop()
  observer.join()
  logging.info("Observer stopped.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        print("Usage: python watcher.py /path/to/monitor")