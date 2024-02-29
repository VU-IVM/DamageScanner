"""DamageScanner - a directe damage assessment toolkit
"""

import requests
import tempfile
import time
import os
import shutil
from tqdm import tqdm


def fetch_and_save(
    url, file_path, overwrite=False, max_retries=3, delay=5, chunk_size=16384
):
    """
    Fetches data from a URL and saves it to a temporary file, with a retry mechanism.
    Moves the file to the destination if the download is complete.
    Removes the temporary file if the download is interrupted.
    """
    if not overwrite and file_path.exists():
        return True

    attempts = 0
    temp_file = None

    while attempts < max_retries:
        try:
            print(f"Downloading {url} to {file_path}")
            # Attempt to make the request
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raises HTTPError for bad status codes

            # Create a temporary file
            temp_file = tempfile.NamedTemporaryFile(delete=False)

            # Write to the temporary file
            total_size = int(response.headers.get("content-length", 0))
            progress_bar = tqdm(total=total_size, unit="B", unit_scale=True)
            for data in response.iter_content(chunk_size=chunk_size):
                temp_file.write(data)
                progress_bar.update(len(data))
            progress_bar.close()

            # Close the temporary file
            temp_file.close()

            # Move the temporary file to the destination
            shutil.move(temp_file.name, file_path)

            return True  # Exit the function after successful write

        except requests.RequestException as e:
            # Log the error
            print(f"Request failed: {e}. Attempt {attempts + 1} of {max_retries}")

            # Remove the temporary file if it exists
            if temp_file is not None and os.path.exists(temp_file.name):
                os.remove(temp_file.name)

            # Increment the attempt counter and wait before retrying
            attempts += 1
            time.sleep(delay)

    # If all attempts fail, raise an exception
    raise Exception("All attempts to download the file have failed.")


def _check_output_path(given_args):
    """Ensures given output path exists.

    Arguments:
        *given_args* : dict, of keyword arguments.

    Returns:
        *str* : output_path, which may be empty string ('')
    """
    output_path = given_args.get("output_path", "")

    if output_path != "" and not output_path.exists():
        output_path.mkdir(parents=True)
    return output_path


def _check_scenario_name(given_args):
    """Ensures given output path exists.

    Arguments:
        *given_args* : dict, of keyword arguments.

    Returns:
        *str* : scenario_name
    """
    scenario_name = given_args.get("scenario_name", False)
    if not scenario_name:
        raise ValueError("Required `scenario_name` not defined.")

    return scenario_name
