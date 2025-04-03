import logging
from typing import Optional
from pathlib import Path
import sys

def setup_logger(
    name: str = "pocket_extraction",
    level: int = logging.INFO,
    quiet: bool = False,
    debug: bool = False,
    logfile: Optional[str] = None,
    file_level: int = logging.DEBUG
) -> logging.Logger:
    """Configure and return a standardized logger with file and console output.
    
    Args:
        name: Logger name (package name)
        level: Console logging level
        quiet: If True, sets console level to WARNING
        debug: If True, sets console level to DEBUG
        logfile: Path to log file (None disables file logging)
        file_level: Log level for file output
        
    Returns:
        Configured logger instance
        
    Raises:
        PermissionError: If log file cannot be written
    """
    # Determine console log level
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)  # Set to lowest level, control via handlers
    
    # Clear existing handlers
    logger.handlers = []
    
    # Create formatters
    console_formatter = logging.Formatter(
        '%(levelname)-8s %(message)s'
    )
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler (stderr)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if logfile:
        try:
            log_path = Path(logfile).resolve()
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_path, mode='a')
            file_handler.setLevel(file_level)
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
            logger.debug(f"Logging to file: {log_path}")
        except Exception as e:
            logger.error(f"Failed to initialize log file: {str(e)}")
            raise PermissionError(f"Cannot write to log file: {logfile}") from e
    
    return logger

# Default package logger (console only)
logger = setup_logger()
