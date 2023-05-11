# standard library imports
import logging
import os

# 3rd party library imports
from rich.logging import RichHandler

# local library imports
from gbrs.emase import emase_utils


class DebugLogFilter(logging.Filter):
    def filter(self, record: logging.LogRecord):
        if record.levelno == logging.DEBUG:
            return True
        return False


class NoDebugLogFilter(logging.Filter):
    def filter(self, record: logging.LogRecord):
        if record.levelno != logging.DEBUG:
            return True
        return False


ch1 = RichHandler(level=logging.NOTSET, show_level=True, show_time=True, show_path=False, omit_repeated_times=False)
ch1.addFilter(NoDebugLogFilter())
ch2 = RichHandler(level=logging.NOTSET, show_level=True, show_time=True, show_path=True, omit_repeated_times=False)
ch2.addFilter(DebugLogFilter())

logging.basicConfig(
    level="NOTSET",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[ch1, ch2]
)


def configure_logging(level: int = 0) -> logging.Logger:
    """
    Configure the logger with the specified `level`. Valid `level` values
    are:

    ======  =========================
    level   logging value
    ======  =========================
    0       WARNING is informational
    1       INFO is user debug
    2+      DEBUG is developer debug
    ======  =========================

    Anything greater than 2 is treated as 2.

    Args:
        level: The logging level; defaults to 0.

    Returns:
        logging.Logger: The logging object.
    """
    log = logging.getLogger('gbrs')

    if level == 0:
        log.setLevel(logging.WARNING)
    elif level == 1:
        log.setLevel(logging.INFO)
    elif level > 1:
        log.setLevel(logging.DEBUG)

    return log


def check_file(file_name: str, mode: str | None = "r") -> str:
    """
    Check if file_name exists and accessible for reading or writing.

    Args:
        file_name: The name of the file.
        mode: "r" for reading, "w" for writing.

    Returns:
        The absolute path of the file.

    Raises:
        FileNotFoundError: When the file doesn't exist or cannot be generated.
        OSError: When the file cannot be generated.
    """
    if mode == 'r':
        if file_name and os.path.exists(file_name):
            return os.path.abspath(file_name)

        raise FileNotFoundError(f"The following file does not exist: {file_name}")
    elif mode == 'r':
        file_dir = '.'

        if file_name:
            file_name = os.path.abspath(file_name)
            file_dir = os.path.dirname(file_name)

            if not os.access(file_dir, os.W_OK | os.X_OK):
                raise OSError(f"Cannot generate file: {file_name}")

            return file_name

    raise OSError(f"Unspecified mode to open file, '{mode}'")


def is_comment(s: str) -> bool:
    """
    Check if the line is a comment (starts with a #)

    Args:
        s: the string to check

    Returns:
        True if the string starts with a #.
    """
    return s.startswith('#')


def get_names(id_file: str) -> list[str]:
    """
    Load a file and get the names.

    Args:
        id_file: the name of the file to open

    Returns:
        A list of the names in the file.
    """
    ids = dict()
    master_id = 0
    with open(id_file) as fh:
        for line in fh:
            item = line.rstrip().split('\t')
            g = item[0]
            if g not in ids:
                ids[g] = master_id
                master_id += 1
    num_ids = len(ids)
    names = {index: name for name, index in ids.items()}
    return [names[k] for k in range(num_ids)]
