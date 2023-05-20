# standard library imports
import logging
import os

# 3rd party library imports
from rich.logging import RichHandler

# local library imports
#


def get_logger(logger_name: str = 'gbrs') -> logging.Logger:
    """
    Get the logger.

    Args:
        logger_name: Name of the logger.

    Returns:
        logging.Logger: The logging object.
    """
    return logging.getLogger(logger_name)


def configure_logging(
    logger_name: str = 'gbrs', level: int = 0
) -> logging.Logger:
    """
    Configure the logger with the specified `level`. Valid `level` values
    are:

    ======  =================================
    level   logging value
    ======  =================================
    0       logging.WARNING is informational
    1       logging.INFO is user debug
    2+      logging.DEBUG is developer debug
    ======  =================================

    Anything greater than 2 is treated as 2.

    if the environment variable ENSIMPL_APP_DEBUG is set to 1,
    there will be more detailed debugging information.

    Args:
        logger_name: The name of the logger.
        level: The logging level; defaults to 0.

    Returns:
        logging.Logger: The logging object.
    """
    ensimpl_app_debug = nvli(os.environ.get('GBRS_APP_DEBUG', '0'), -1)

    rich_handler = RichHandler(
        level=logging.NOTSET,
        show_level=False,
        show_time=True,
        show_path=False,
        omit_repeated_times=False,
    )

    if ensimpl_app_debug == 1:
        rich_handler = RichHandler(
            level=logging.NOTSET,
            show_level=True,
            show_time=True,
            show_path=True,
            omit_repeated_times=False,
        )

    # this is configuring the root logger and below
    # set the level to WARNING, as that is the default
    logging.basicConfig(
        level=logging.WARNING,
        format='%(message)s',
        datefmt=f'{logger_name} [%X]',
        handlers=[rich_handler],
    )

    log = logging.getLogger(logger_name)

    # set gbrs's logging level
    if level == 0:
        log.setLevel(logging.WARNING)
    elif level == 1:
        log.setLevel(19)
    elif level > 1:
        log.setLevel(logging.DEBUG)

    return log


def nvli(value, default) -> int:
    """Returns `value` as an int if `value` can be converted, else `default`.

    Args:
        value: The value to evaluate and convert to an it.
        default: The default value.

    Returns:
        Either `value` or `default`.
    """
    ret = default
    if value:
        try:
            ret = int(value)
        except ValueError:
            pass
    return ret


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
