import logging
import subprocess

class PipelineLauncher:

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

    def launch(self, command:list):
        """Function that launches the pipeline using the given command.

        Args:
            command (list): Command to run as a list of strings.
        """
        logging.info(f"Launching pipeline with command: {' '.join(command)}")
        try:
            self._subprocess_run(command)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running pipeline: {e.stderr}")
        else:
            logging.info("Pipeline launched successfully.")
        
