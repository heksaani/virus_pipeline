"""Class for species-specific pipeline launching functions."""
import logging
import os 
from datetime import datetime
import config
from .fastq_files import check_platform
class SpeciesPipelines:
    def __init__(self, species:str, sample_name:str, sample_reads:dict):
        """Initialize the SpeciesPipelines with Docker command and paths.   
        Args:
            species (str): Species being analyzed.
            getuid (function, optional): Function to get user ID. Defaults to os.getuid.
            docker_command (str, optional): Docker command template. Defaults to config.DOCKER_COMMAND.
        """
        self.species = species

        self.uid = os.getuid()
        self.gid = config.GID
        self.user_flag = f"{self.uid}:{self.gid}"
        self.pipeline_config = config.PIPELINE_CONFIGS.get(species, {})
        self.pipeline_name = self.pipeline_config.get('pipeline_name', '')
        read_path = sample_reads[sample_name]['R1'] or sample_reads[sample_name]['R2']
        self.platform = check_platform(read_path)
        self.image_name = self.pipeline_config.get('image_name', '')
        self.required_reads = self.pipeline_config.get('platforms', {}).get(self.platform, {}).get('required', [])
        self.docker_command = config.DOCKER_COMMAND

    def convert_to_docker_path(self, filepath:str) -> str:
        """Function that converts a host file path to a Docker container file path.
        
        Args:
            filepath (str): Original file path on the host system.
        
        Returns:
            str: Converted file path for Docker container.
        """
        docker_path = filepath.replace(config.DATA_PATH, config.DOCKER_DATA_PATH)
        return docker_path
    def check_sample(self, sample_name:str, sample_reads:dict):
        """ checks if a sample is ready based on platform and available reads"""
        read_path = sample_reads[sample_name]['R1'] or sample_reads[sample_name]['R2']
        platform = check_platform(read_path)
        required = self.pipeline_config['platforms'][platform]['required']
        for read_type in required:
            if sample_reads[sample_name].get(read_type) is None:
                return False
    
        return True

    def build_docker_command(self, sample_name:str, sample_reads:dict) -> list[str]:
        """Function that builds the final Docker command as a list of strings.
        
        Returns:
            list[str]: Base Docker command.
        """
        possible_species = config.PIPELINE_CONFIGS.keys()
        if self.species in possible_species:
            logging.info(f"Building Docker command for species: {self.species}")
            inner_command = getattr(self, f"{self.species}_pipeline")(sample_name, sample_reads)
        else:
            raise ValueError(f"Unsupported species: {self.species}")
        #full_command = docker_command.replace("{INNER_COMMAND}", inner_command)
        #parts = full_command.split()
        #bash_c_index = parts.index('-c')
        # Rejoin everything after '-c' into a single element
        #return parts[:bash_c_index + 1] + [' '.join(parts[bash_c_index + 1:])]
        logging.info(f"Inner command: {inner_command}")
        #command = [part.format(USER_FLAG=self.user_flag, SPECIES_IMAGE=self.image_name, INNER_COMMAND=inner_command) for part in self.docker_command]
        command = [part.format(DATA_PATH=config.DATA_PATH, DOCKER_DATA_PATH=config.DOCKER_DATA_PATH, RESULTS_PATH=config.RESULTS_PATH,
                               DOCKER_RESULTS=config.DOCKER_RESULTS, SPECIES_IMAGE=self.image_name, INNER_COMMAND=inner_command) for part in self.docker_command]
        logging.info(f"command list : {command}")

        return command
    
    def influenza_pipeline(self, sample_id:str, reads:dict):
        """Function that builds the command to launch the Influenza pipeline inside a Docker container.

        Args:
            sample_id (str): Identifier for the sample.
            reads (dict): Dictionary containing paths to R1 and R2 reads.
        Returns:
            list[str]: Command to run as a list of strings.
        """
        platform = self.platform
        try:
            date = datetime.now().strftime("%Y%m%d")
            container_results = f'{config.DOCKER_RESULTS}/{date}_IRMA_RESULTS/{sample_id}'

            if platform == 'minion':

                logging.info("Not yet implementeds")
            elif platform  in ['nextseq', 'miseq']:
                logging.info(f"Building command for platform: {platform}")
                container_r1 = self.convert_to_docker_path(reads['R1'])
                container_r2 = self.convert_to_docker_path(reads['R2'])
                irma_command = self.pipeline_config['platforms'][platform]['irma_command']
                inner_command = f'IRMA {irma_command} {container_r1} {container_r2} {container_results}'
                
                return inner_command
            else:
                logging.error(f"Platform {platform} not recognized for Influenza pipeline")
                raise ValueError(f"Platform {platform} not recognized for Influenza pipeline")
        except Exception as e:
            logging.error(f"Error building Influenza pipeline command: {e}")
            raise
        