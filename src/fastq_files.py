""""Script that handle fastq files for virus piper."""
import re
import glob
import logging

def extract_barcode(filename:str, pattern:str=r'barcode\d+') -> str | None:
    """Extracts the barcode from a given filename based on a pattern using regex.

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

def combine_minion_reads(input_path:str):
    """Function that combines all fastqs into one for one minion barcode

    Args:
        path (str) : Path to the folder with fastq files 
    
    """
    pass  # Implementation goes here

