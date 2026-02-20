"""DamageScanner - a directe damage assessment toolkit."""

from pathlib import Path


def _check_output_path(output_path: str | Path | None = "") -> Path | str:
    """Ensures the output directory exists and returns its path.

    Returns:
        A Path object for the output directory, or an empty string if no path is provided.
    """
    if output_path == "" or output_path is None:
        return ""

    # Convert to Path if it's a string
    output_path = Path(output_path)

    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)

    return output_path


def _check_scenario_name(scenario_name: str | None = None) -> str:
    """Validates that a scenario name is provided in the arguments.

    Args:
        scenario_name: Name of the scenario.

    Returns:
        The scenario name.

    Raises:
        ValueError: If 'scenario_name' is missing or False.
    """
    if not scenario_name:
        raise ValueError("Required `scenario_name` not defined.")

    return scenario_name
