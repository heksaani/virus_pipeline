import logging
from .fastq_files import extract_barcode, get_sample_name, check_platform, get_species
from .pipeline_launcher import PipelineLauncher

class SampleHandler:
    def __init__(self, pipeline_launcher:PipelineLauncher):
        self.pipeline_launcher = pipeline_launcher
        self.sample_reads = {}
    def _add_reads_to_sample(self, sample_name:str, path_str:str):
        """Function that adds the read paths to the sample dictionary

        Args:
            sample_name (str) : Name of the sample
            path_str (str) : Path to the fastq file
        """
        if sample_name not in self.sample_reads:
            logging.info("Sample not found in dict adding")
            self.sample_reads[sample_name] = {'R1': None, 'R2': None}
        if 'R1' in path_str:
            self.sample_reads[sample_name]['R1'] = path_str
            logging.info(f"Added R1 for {sample_name}")
        elif 'R2' in path_str:
            self.sample_reads[sample_name]['R2'] = path_str
            logging.info(f"Added R2 for {sample_name}")
        else:
            logging.error("No R1 or R2 in file name")
            return

    def process_file(self, filepath):
        logging.info(f"Processing file: {filepath}")
        if filepath.endswith('.fastq.gz'):
            logging.info(f"Fastq file detected: {filepath}")
            platform = check_platform(filepath)
            if platform is None:
                logging.error("Platform could not be determined.")

            elif platform == 'minion':
                barcode = extract_barcode(filepath)
                logging.info(f"Barcode : {barcode}")
                #combine_minion_reads(path_str)
                # self._add_reads_to_sample(sample_name, path_str)
            elif platform == 'nextseq' or platform == 'miseq':
                sample_name = get_sample_name(filepath)
                logging.info(f"Sample {sample_name} ")
                self._add_reads_to_sample(sample_name, filepath)
                try:
                    reads = self.sample_reads[sample_name] # this is dictionary with R1 and R2 paths
                    if reads['R1'] and reads['R2']:
                        logging.info(f"Both reads found for sample {sample_name}, launching pipeline.")
                        species = get_species(filepath)
                        self.pipeline_launcher.launch(platform, species, sample_name, reads)
                        del self.sample_reads[sample_name]
                        logging.info(f"Sample {sample_name} removed from sample_reads dictionary after launching pipeline.")
                except KeyError:
                    logging.error(f"Sample {sample_name} not found in sample_reads dictionary.")
            else:
                logging.error("No R1 or R2 in file name")
        else:
            logging.info(f"File is not a fastq.gz file: {filepath}")