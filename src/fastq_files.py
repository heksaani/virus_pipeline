""""Script that handle fastq files for virus piper."""
import re
import logging

def extract_barcode(filename:str, pattern:str=r'barcode\d+') -> str | None:
    """Function that extracts the barcode from a given filename. 

    Args:
        filename (str): The name of the file.
        pattern (str): The pattern to look for in the filename.

    Returns:
        str: The extracted barcode or None if not found.
    """
    try:
        match = re.search(pattern, filename)
    except re.error as e:
        logging.error(f"Regex error: {e}")
        return None

    return match.group(0) if match else None


def get_sample_name(input_path:str):
    """Function that extracts sample name from fastq file path

    Args:
        input_path (str) : Path to the fastq file 
    
    Returns:
        sample_name (str) : Extracted sample name
    """
    filename = input_path.split('/')[-1]
    sample_name = filename.split('_')[0]
    return sample_name 

def check_platform(input_path:str) -> str | None:
    """Function that checks the sequencing platform from the file path

    Args:
        input_path (str) : Path to the fastq file
    Returns:
        platform (str) : Sequencing platform (minion, nextseq, miseq) or None if not found
    """
    platform = input_path.split('/')[5]
    if platform not in ['minion', 'nextseq','miseq']:
        logging.error(f"Invalid platform detected {platform}")
        return None
    return platform

def get_species(input_path:str) -> str:
    """Function that checks the species from the file path

    Args:
        input_path (str) : Path to the fastq file
    Returns:
        species (str) : Species extracted from the fastq file path
    """
    try:
        species = input_path.split('/')[4]
    except IndexError as e:
        logging.error(f"Error extracting species from path {input_path}: {e}")
        return ""
    return species

def combine_minion_reads(input_path:str):
    """Function that combines all fastqs into one for one minion barcode

    Args:
        path (str) : Path to the folder with fastq files 
    
    Returns: 
        output_path (str) : Path to the combined fastq file
    """
    #folder = '/'.join(input_path.split('/')[:-1])
    pass

