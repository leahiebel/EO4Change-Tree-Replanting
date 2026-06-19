"""Earth Engine authentication helpers."""

from __future__ import annotations

import logging

import ee

LOGGER = logging.getLogger(__name__)


def initialize_earth_engine(project: str, authenticate: bool = True) -> None:
    """Authenticate and initialize Google Earth Engine.

    Parameters
    ----------
    project:
        Google Cloud project registered for Earth Engine.
    authenticate:
        When true, call ``ee.Authenticate()`` before initialization. No
        credentials or service-account secrets are stored in source code.
    """
    if not project or project == "replace-with-project-id":
        raise ValueError("Set earth_engine_project in config/config.yaml.")

    if authenticate:
        LOGGER.info("Starting Earth Engine user authentication.")
        ee.Authenticate()

    LOGGER.info("Initializing Earth Engine with project %s.", project)
    ee.Initialize(project=project)

