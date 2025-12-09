import logging
import os 
import config 
from datetime import datetime
import subprocess

class PipelineLauncher:
    def __init__(self):
        self.host_data = config.DATA_PATH # /mnt/volume/data
        self.container_data = getattr(config, 'DOCKER_MOUNT_PATH', '/data')
        self.results = getattr(config, 'RESULTS_PATH', '/mnt/volume/results')
        self.docker_command = getattr(config, 'DOCKER_COMMAND', 'docker run --rm -u $(id -u):$(id -g) -v ' + self.host_data + ':/data:ro PROXY_ENV')  
        self.docker_image = getattr(config, 'DOCKER_IMAGE', 'irma')
        self.uid = os.getuid()
        self.gid = os.getgid()
        logging.info(f"PipelineLauncher initialized - UID:GID = {self.uid}:{self.gid}")
        #self.
    def _to_container_path(self, host_path:str) -> str:
        """
        Changing the host path to the container path
        """
        return host_path.replace(self.host_data, self.container_data)
    ## TODO, implement minion launching
    ## TODO generalize so that irma is not the only pipeline so that other species can be used

    def launch(self, platform:str, sample_id:str, R1_path:str, R2_path:str):
        """Launch the analysis pipeline inside a Docker container.
        
        Args:
            platform (str): Sequencing platform ('minion', 'nextseq', 'miseq').
            sample_id (str): Identifier for the sample.
            R1_path (str): Path to the R1 FASTQ file.
            R2_path (str): Path to the R2 FASTQ file.
        """
        logging.info(f"Launching pipeline for {sample_id}")
        if platform != 'minion':
            container_R1 = self._to_container_path(R1_path)
            container_R2 = self._to_container_path(R2_path)
            date = datetime.now().strftime("%Y%m%d")
            results = self._to_container_path(self.results) + f'/{date}_IRMA_RESULTS/{sample_id}'
            command = self.docker_command + f" {self.docker_image} IRMA FLU {container_R1} {container_R2} {results}"
            # turn command into a list for subprocess
            command = command.split()
            logging.info(f"Command to run: {command}")
            subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True
            )
            logging.info(f"IRMA completed for {sample_id}")

        else:
            #command = docker_command + f" {self.docker_image} IRMA-minion {container_R1} {container_R2}"
            logging.error("Minion platform not yet implemented")
            return
        logging.info(f"Running command: {command}")