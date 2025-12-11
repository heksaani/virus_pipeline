"""Class for species-specific pipeline launching functions."""
import logging
from datetime import datetime
import os
import config
class SpeciesPipelines:
    def __init__(self, docker_command: str, host_data: str, container_data_path: str) -> None:
        """Initialize the SpeciesPipelines with Docker command and paths.   
        Args:
            docker_command (str): Base Docker command.
            host_data (str): Path on the host machine data directory.
            container_data_path (str): Path inside the container.
        """
        self.docker_command = docker_command # base command for all pipelines
        self.host_data = host_data # path on the host machine data directory
        self.container_data_path = container_data_path # path inside the container
    
    def _to_container_path(self, host_path:str) -> str:
        """
        Changing the host path to the container path

        Args:
            host_path (str): Path on the host machine.
        Returns:§
            str: Corresponding path inside the Docker container.
        """
        return host_path.replace(self.host_data, self.container_data_path)
    
    def influenza_pipeline(self, docker_image, platform:str, sample_id:str, reads:dict) -> list[str]:
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
            inner_command = inner_command = f"umask 0002 && IRMA FLU {container_R1} {container_R2} {container_results}"
            uid = os.getuid()
            gid = getattr(config, 'GID', os.getgid())
            command_list = [
                'docker', 
                'run', 
                '--rm', 
                '-u', f"{uid}:{gid}", 
                '-v', '/mnt/volume/data:/data', 
                '-v', '/mnt/volume/results:/output', 
                'irma_thl', 
                'bash', 
                '-c', 
                inner_command
            ]
            logging.info(f"Running command: {command_list}")
            
            return command_list

        else:
            logging.error("Minion platform not yet implemented")
            return []