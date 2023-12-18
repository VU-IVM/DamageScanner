"""DamageScanner - a directe damage assessment toolkit
"""

def check_output_path(given_args):
    """Ensures given output path exists.

    Arguments:
        *given_args* : dict, of keyword arguments.

    Returns:
        *str* : output_path, which may be empty string ('')
    """
    output_path = given_args.get('output_path', '')

    if output_path != '' and not output_path.exists():
        output_path.mkdir(parents=True)
    return output_path

def check_scenario_name(given_args):
    """Ensures given output path exists.

    Arguments:
        *given_args* : dict, of keyword arguments.

    Returns:
        *str* : scenario_name
    """
    scenario_name = given_args.get('scenario_name', False)
    if not scenario_name:
        raise ValueError("Required `scenario_name` not defined.")

    return scenario_name