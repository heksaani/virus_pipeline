import logging
import os 
import config 
from datetime import datetime
import subprocess
    ## TODO, implement minion launching
    ## TODO generalize so that irma is not the only pipeline so that other species can be used
class PipelineLauncher:
    def __init__(self):
        self.host_data = config.DATA_PATH # /mnt/volume/data
        self.container_data = getattr(config, 'DOCKER_MOUNT_PATH', '/data')
        self.uid = os.getuid()
        self.gid = 1003 #config.GID
        self.user_flag = f"{self.uid}:{self.gid}"
        self.docker_command = (config.DOCKER_COMMAND).format(USER_FLAG=self.user_flag) #+ ' ' + config.PROXY_ENV
        self.docker_image = getattr(config, 'DOCKER_IMAGE', 'irma')
    def _to_container_path(self, host_path:str) -> str:
        """
        Changing the host path to the container path

        Args:
            host_path (str): Path on the host machine.
        Returns:
            str: Corresponding path inside the Docker container.
        """
        return host_path.replace(self.host_data, self.container_data)
    
    def influenza_pipeline(self, platform, sample_id:str, reads:dict):
        """Function that launches the IRMA influenza pipeline inside a Docker container.

        Args:
            platform (str): Sequencing platform ('minion', 'nextseq', 'miseq').
            sample_id (str): Identifier for the sample.
            reads (dict): Dictionary containing paths to R1 and R2 reads.
        """

        if platform != 'minion':
            container_R1 = self._to_container_path(reads['R1'])
            container_R2 = self._to_container_path(reads['R2'])
            date = datetime.now().strftime("%Y%m%d")
            container_results = '/output' + f'/{date}_IRMA_RESULTS/{sample_id}'
            command = self.docker_command + f" {self.docker_image} IRMA FLU {container_R1} {container_R2} {container_results}"
            logging.info(f"Running command: {command}")
            command = command.split()
                #subprocess.run(
                #    command,
                #    check=True,
                #    capture_output=True,
                #    text=True
                #)
            logging.info(f"IRMA completed for {sample_id}")

        else:
            logging.error("Minion platform not yet implemented")
            return
    def launch(self, platform:str, species:str, sample_id:str, reads:dict):
        """Launch the analysis pipeline inside a Docker container.
        
        Args:
            platform (str): Sequencing platform ('minion', 'nextseq', 'miseq').
            sample_id (str): Identifier for the sample.
            species (str): Species being analyzed (e.g., 'influenza').
            reads (dict): Dictionary containing paths to R1 and R2 reads.
        """
        
        if species == 'influenza':
            self.influenza_pipeline(platform, sample_id, reads)
        else:
            logging.error(f"Species {species} not yet implemented")
            return