import logging
import os 
import config 
from .species_pipelines import SpeciesPipelines
import subprocess
    ## TODO, implement minion launching
    ## TODO generalize so that irma is not the only pipeline so that other species can be used
class PipelineLauncher:
    def __init__(self):
        self.uid = os.getuid()
        self.gid = getattr(config, 'GID', os.getgid())
        self.user_flag = f"{self.uid}:{self.gid}"
        self.docker_command = (config.DOCKER_COMMAND).format(USER_FLAG=self.user_flag)
        self.host_data = config.DATA_PATH
        self.container_data_path = config.DOCKER_MOUNT_PATH

    def _subprocess_run(self, command:list):
        """Function that runs the subprocess command

        Args:
            command (list): Command to run as a list of strings.
        """
        subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True)
        logging.info("Subprocess command completed successfully.")

    def launch(self, platform:str, species:str, sample_id:str, reads:dict):
        """Launch the analysis pipeline inside a Docker container.
        
        Args:
            platform (str): Sequencing platform ('minion', 'nextseq', 'miseq').
            sample_id (str): Identifier for the sample.
            species (str): Species being analyzed (e.g., 'influenza').
            reads (dict): Dictionary containing paths to R1 and R2 reads.
        """
        try:
            docker_image = config.SPECIES[species]
        except KeyError:
            logging.error(f"Species {species} is known, but no Docker image is defined in config.")
            return
        if species == 'influenza':
            #crete the pipeline object
            pipeline_factory = SpeciesPipelines(
                docker_command=self.docker_command,
                host_data=self.host_data,
                container_data_path=self.container_data_path
            )
            command = pipeline_factory.influenza_pipeline(docker_image, platform, sample_id, reads)
            self._subprocess_run(command)
        else:
            logging.error(f"Species {species} not yet implemented")
            return