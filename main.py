"""Main script for virus piper."""
import config
import os
import logging
from src.watcher import DirectoryMonitor
from src.sample_handler import SampleHandler
from src.pipeline_launcher import PipelineLauncher
from watchdog.observers import Observer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('/mnt/volume/logs/virus_piper.log'),
        logging.StreamHandler()])

def main():
    logging.info('Starting main script...')

    path = config.DATA_PATH
    if not os.path.exists(path):
        logging.error(f"Path {path} does not exist.")
        return
    if not os.path.isdir(path):
        print(f"Error: Path is not a directory: {path}")
        return
    launcher = PipelineLauncher()
    sample_handler = SampleHandler(launcher)
    event_handler = DirectoryMonitor(callback=sample_handler.process_file)

    observer = Observer()
    observer.schedule(event_handler, path, recursive=True)
    observer.start()
    logging.info(f"Monitoring started on path: {path}")

    try:
        while observer.is_alive():
            observer.join(1)
    except KeyboardInterrupt:
        logging.info("Stopping observer...")
        observer.stop()
    
    observer.join()
    logging.info("Observer stopped.")
if __name__ == "__main__":
    main()
