import logging
import subprocess
class PipelineLauncher:

    def launch(self, command: list):
        """Function that launches the pipeline using the given command.

        Args:
            command (list): Command to run as a list of strings.
        """
        logging.info(f"Launching pipeline with command: {''.join(command)}")
        try:
            subprocess.run(command,
                           capture_output=True,
                           text=True,
                           check=True)
            logging.info("Pipeline completed successfully.")
        except subprocess.CalledProcessError as e:
            logging.error("Error running pipeline")
            logging.error(f"Return code: {e.returncode}")
            logging.error(f"Output: {e.output}")
            logging.error(f"Error message: {e.stderr}")
            raise