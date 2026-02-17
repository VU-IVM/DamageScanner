import pathlib

data_path = pathlib.Path(__file__).parent.parent / "data"
output_folder = pathlib.Path("tests/output")
tmp_folder = pathlib.Path("tests/tmp")


def create_test_folders():
    """Create output and tmp folders for tests if they don't exist."""
    output_folder.mkdir(exist_ok=True)
    tmp_folder.mkdir(exist_ok=True)


# Run once on import
create_test_folders()
