"""Main script for virus piper."""
import config
import os
#import subprocess
import src.watcher as watcher
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('/mnt/volume/logs/virus_piper.log'),
        logging.StreamHandler()
    ]
    )

def main():
    logging.info('Starting main script...')
    os.umask(0) # THIS PREVENTS FOLDERS TO BE CREATED WITH NOT ALL PERMISSIONS
    # umask(0) sets the file mode creation mask to 0, allowing all permissions for newly created files and directories.
    path = config.DATA_PATH
    if not os.path.exists(path):
        logging.error(f"Path {path} does not exist.")
        return
    if not os.path.isdir(path):
        print(f"Error: Path is not a directory: {path}")
        return
    watcher.main(path)
if __name__ == "__main__":
    main()
