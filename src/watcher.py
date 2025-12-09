"""Watcher module for Virus Piper pipeline
This module handles directory monitoring and file event handling.
"""
import logging
from watchdog.events import FileSystemEventHandler

class DirectoryMonitor(FileSystemEventHandler):

  def __init__(self, callback):
    super().__init__()
    self.callback = callback

  def on_created(self, event):
    if not event.is_directory:
      path_str = str(event.src_path)
      if path_str.endswith('.fastq.gz'):
          self.callback(path_str)
      else:
        logging.info(f"Ignored non-fastq file: {event.src_path}")
    else:
      logging.info(f"Ignored directory creation: {event.src_path}")