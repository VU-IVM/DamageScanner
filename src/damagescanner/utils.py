"""DamageScanner - a directe damage assessment toolkit"""


def _check_output_path(given_args):
    """Ensures the output directory exists and returns its path.

    Args:
        given_args (dict): Dictionary of keyword arguments, potentially containing 'output_path'.

    Returns:
        str or Path: The output path. Returns an empty string if not specified.

    Raises:
        OSError: If the directory cannot be created (e.g., due to permissions).
    """

    output_path = given_args.get("output_path", "")

    if output_path != "" and not output_path.exists():
        output_path.mkdir(parents=True)
    return output_path


def _check_scenario_name(given_args):
    """Validates that a scenario name is provided in the arguments.

    Args:
        given_args (dict): Dictionary of keyword arguments, expected to contain 'scenario_name'.

    Returns:
        str: The scenario name.

    Raises:
        ValueError: If 'scenario_name' is missing from the arguments.
    """
    scenario_name = given_args.get("scenario_name", False)
    if not scenario_name:
        raise ValueError("Required `scenario_name` not defined.")

    return scenario_name
