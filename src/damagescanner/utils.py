"""DamageScanner - a directe damage assessment toolkit."""

from pathlib import Path


def _check_output_path(given_args: dict) -> Path | str:
    """Ensures the output directory exists and returns its path.
    
    Returns:
        A Path object for the output directory, or an empty string if no path is provided.
    """
    output_path = given_args.get("output_path", "")

    if output_path == "":
        return output_path

    # Convert to Path if it's a string
    output_path = Path(output_path)

    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)

    return output_path


def _check_scenario_name(given_args: dict) -> str:
    """Validates that a scenario name is provided in the arguments.

    Args:
        given_args: Dictionary of keyword arguments, expected to contain 'scenario_name'.

    Returns:
        The scenario name.

    Raises:
        ValueError: If 'scenario_name' is missing from the arguments.
    """
    scenario_name = given_args.get("scenario_name", False)
    if not scenario_name:
        raise ValueError("Required `scenario_name` not defined.")

    return scenario_name
