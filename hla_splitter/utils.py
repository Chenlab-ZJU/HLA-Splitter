import subprocess
import shlex
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)

def run_command(command_str, check=True, cwd=None, shell=False):
    """
    Executes a command string using subprocess.run.

    Args:
        command_str (str): The command string to execute.
        check (bool): If True, raise CalledProcessError on non-zero exit code.
        cwd (str, optional): Directory to run command in. Defaults to None.
        shell (bool): Whether to use the shell. Be cautious if command_str
                      includes untrusted input. Usually safer to set shell=False
                      and pass command as a list.

    Returns:
        subprocess.CompletedProcess: The result object.

    Raises:
        subprocess.CalledProcessError: If check is True and command fails.
    """
    log.info(f"Running command: {command_str}")

    # By default, capture stderr for logging purposes.
    # Capture stdout unless shell=True and redirection (> or >>) is used in the command.
    # If redirection is used, we let the shell handle it instead of Python capturing stdout.
    stdout_target = subprocess.PIPE
    stderr_target = subprocess.PIPE
    text_decode = True

    if shell and ('>' in command_str or '>>' in command_str):
        # When the command string contains shell redirection operators,
        # we instruct subprocess.run not to capture stdout.
        # The command_str will be executed directly by the shell, which will handle the file redirection.
        stdout_target = None
        
    try:
        # If not using shell, split the command safely
        cmd_list = shlex.split(command_str) if not shell else command_str
        result = subprocess.run(
            cmd_list,
            check=check,
            stdout=stdout_target, # Determine whether to capture stdout based on redirection checks
            stderr=stderr_target,
            text=text_decode,
            cwd=cwd,
            shell=shell # Use shell=True cautiously if needed, e.g., for pipes in command_str
        )
        if result.stdout:
            log.info(f"Command STDOUT:\n{result.stdout.strip()}")
        if result.stderr:
            log.warning(f"Command STDERR:\n{result.stderr.strip()}")
        log.info(f"Command finished successfully.")
        return result
        
    except subprocess.CalledProcessError as e:
        log.error(f"Command failed with exit code {e.returncode}")
        log.error(f"Command: {e.cmd}")
        log.error(f"STDERR:\n{e.stderr.strip()}")
        if e.stdout:
            log.error(f"STDOUT:\n{e.stdout.strip()}")
        # Re-raise the error if check=True was intended behaviour
        if check:
            raise e
        # Otherwise return the failed result object
        return e
    except FileNotFoundError:
        log.error(f"Error: Command not found. Is the tool installed and in PATH? Command: {command_str}")
        sys.exit(1)
    except Exception as e:
        log.error(f"An unexpected error occurred running command '{command_str}': {e}")
        sys.exit(1)
