"""Watcher module for Virus Piper pipeline
This module handles:
1) Detect that a folder appeared

2) Identify which parent platform directory it belongs to

3) Trigger the pipeline associated with that platform
"""

import logging
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

class DirectoryMonitor(FileSystemEventHandler):
  """Handler for file system events"""

  def on_created(self, event):
    """Called when a file or directory is created"""
    if event.is_directory:
      logging.info(f"New directory created: {event.src_path}")
    else:
      logging.info(f"New file created: {event.src_path}")


def main(path:str, recursive: bool = True):
  """
  Start monitoring the specified path.
    
  Args:
        path: Directory path to monitor
        recursive: Whether to monitor subdirectories"""
  
  logging.info('Starting script...')
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