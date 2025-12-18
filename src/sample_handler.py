import logging
import config
from .fastq_files import extract_barcode, get_sample_name, get_species
from .pipeline_launcher import PipelineLauncher
from .species_pipelines import SpeciesPipelines
class SampleHandler:

    """Class that handles samples using the fastq_files module and species"""

    def __init__(self, pipeline_launcher:PipelineLauncher):
        self.pipeline_launcher = pipeline_launcher
        self.sample_reads = {}
    def _add_read_to_sample(self, sample_name:str, read_path_str:str):
        """
        Function that adds the read paths to the sample dictionary
        Args:
            sample_name (str) : Name of the sample
            path_str (str) : Path to the fastq file
        """
        if sample_name not in self.sample_reads:
            logging.info("Sample not found in dict")
            self.sample_reads[sample_name] = {
                'R1': None, 
                'R2': None, 
                'species': None}
            self.sample_reads[sample_name]['species'] = get_species(read_path_str)
        if 'R1' in read_path_str:
            self.sample_reads[sample_name]['R1'] = read_path_str
            logging.info(f"Added R1 read for sample {sample_name}: {read_path_str}")
        elif 'R2' in read_path_str:
            self.sample_reads[sample_name]['R2'] = read_path_str
            logging.info(f"Added R2 read for sample {sample_name}: {read_path_str}")
        else:
            logging.info(f"Read is not R1 or R2: {read_path_str}")
            self.sample_reads[sample_name]['R1'] = read_path_str

    def process_file(self, filepath):
        """Function that processes a single fastq file

            1) Checks if the file is a fastq.gz file
            2) Get sample name
            3) Add read to sample dict, this also addds the species
            4) Check if the sample is ready, which is dependent on the species and platform

        Args:
            filepath (str): Path to the fastq file.
        """
        if filepath.endswith('.fastq.gz'):
            logging.info(f"Fastq file detected: {filepath}")
            sample_name = get_sample_name(filepath)
            self._add_read_to_sample(sample_name, filepath)
            logging.info(self.sample_reads)
            species_obj = SpeciesPipelines(self.sample_reads[sample_name]['species'])
            if species_obj.check_sample(sample_name, self.sample_reads):
                logging.info(f'sample {sample_name} ready...')
                #docker_cmd = species_obj.build_docker_command(sample_name, self.sample_reads[sample_name])
                #self.pipeline_launcher.launch(docker_cmd)
                del self.sample_reads[sample_name]
            else:
                logging.info(f'sample {sample_name} not ready')
        else:
            logging.info(f"File is not a fastq.gz file: {filepath}")
            return
