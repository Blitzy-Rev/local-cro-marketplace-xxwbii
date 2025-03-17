# Import the logging module
import logging  # logging 2.0+

# Define a logger for the molecular test package
logger = logging.getLogger(__name__)

# Indicate that this is a test package for molecular processing components
__test_package__ = True

# Configure basic log formatting for tests if needed
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    logger.info("Molecular test package initialized")