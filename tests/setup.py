import pathlib

data_path = pathlib.Path(__file__).parent.parent / "data"

output_folder = pathlib.Path("tests/output")
output_folder.mkdir(exist_ok=True)

tmp_folder = pathlib.Path("tests/tmp")
tmp_folder.mkdir(exist_ok=True)
