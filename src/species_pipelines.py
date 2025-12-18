"""Class for species-specific pipeline launching functions."""
import logging
from datetime import datetime
import os
import config
from .fastq_files import check_platform
class SpeciesPipelines:
    def __init__(self, species:str) -> None:
        """Initialize the SpeciesPipelines with Docker command and paths.   
        Args:
            species (str): Species being analyzed.
            getuid (function, optional): Function to get user ID. Defaults to os.getuid.
            docker_command (str, optional): Docker command template. Defaults to config.DOCKER_COMMAND.
        """
        self.species = species

        self.uid = os.getuid()
        self.gid = config.GID
        self.pipeline_config = config.PIPELINE_CONFIGS.get(species, {})
        self.docker_command = config.DOCKER_COMMAND.format(USER_FLAG=f"{self.uid}:{self.gid}", 
                                                           SPECIES_IMAGE=self.pipeline_config.get("image_name", ""),
                                                           INNER_COMMAND="{INNER_COMMAND}")

    def check_sample(self, sample_name:str, sample_reads:dict):
        """ checks if a sample is ready based on platform and available reads"""
        read_path = sample_reads[sample_name]['R1'] or sample_reads[sample_name]['R2']
        platform = check_platform(read_path)
        required = self.pipeline_config['platforms'][platform]['required']
        for read_type in required:
            if sample_reads[sample_name].get(read_type) is None:
                return False
        
        return True

    def build_docker_command(self) -> list[str]:
        """Function that builds the final Docker command as a list of strings.
        
        Returns:
            list[str]: Base Docker command.
        """
        docker_command = self.docker_command
        if self.species == 'influenza':
            #inner_command = self.influenza_pipeline()
            logging.info('launching')
            inner_command = ""
        else:
            logging.error('not yet')
            inner_command = ""
        full_command = docker_command.replace("{INNER_COMMAND}", inner_command)
        return full_command.split()
    
    def influenza_pipeline(self, sample_id:str, reads:dict) -> list[str]:
        """Function that builds the command to launch the Influenza pipeline inside a Docker container.

        Args:
            platform (str): Sequencing platform ('minion', 'nextseq', 'miseq').
            sample_id (str): Identifier for the sample.
            reads (dict): Dictionary containing paths to R1 and R2 reads.
        Returns:
            list[str]: Command to run as a list of strings.
        """
        platform = check_platform(reads['R1'])
        try:
            date = datetime.now().strftime("%Y%m%d")
            container_results = f'{config.DOCKER_RESULTS}/{date}_IRMA_RESULTS/{sample_id}'

            if platform == 'minion':

                logging.info("Not yet implementeds")
            elif platform  in ['nextseq', 'miseq']:
                logging.info("Building Influenza pipeline command for Illumina data")
                container_r1 = f'{config.DOCKER_DATA}/{reads["R1"]}'
                container_r2 = f'{config.DOCKER_DATA}/{reads["R2"]}'
                irma_command = self.pipeline_config['platforms'][platform]['irma_command']
                inner_command = f'IRMA {irma_command} {container_r1} {container_r2} {container_results}'
                return inner_command
            else:
                logging.error(f"Platform {platform} not recognized for Influenza pipeline")
                return []
        except Exception as e:
            logging.error(f"Error building Influenza pipeline command: {e}")
            return []

        else:
            logging.error("Minion platform not yet implemented")
            return []
        