import pytest
from pathlib import Path

from .helpers import tmp_folder

from damagescanner.utils import _check_output_path, _check_scenario_name


class TestCheckOutputPath:
    """Tests for _check_output_path function."""

    def test_empty_string_returns_empty(self):
        """Test that empty output_path returns empty string."""
        result = _check_output_path({})
        assert result == ""

    def test_empty_string_explicit(self):
        """Test that explicitly empty output_path returns empty string."""
        result = _check_output_path({"output_path": ""})
        assert result == ""

    def test_returns_path_object(self):
        """Test that result is a Path object when path is provided."""
        result = _check_output_path({"output_path": str(tmp_folder)})
        assert isinstance(result, Path)

    def test_accepts_string_path(self):
        """Test that string paths are accepted."""
        result = _check_output_path({"output_path": str(tmp_folder)})
        assert result == tmp_folder

    def test_accepts_path_object(self):
        """Test that Path objects are accepted."""
        result = _check_output_path({"output_path": tmp_folder})
        assert result == tmp_folder

    def test_creates_directory_if_not_exists(self):
        """Test that directory is created if it doesn't exist."""
        new_dir = tmp_folder / "test_new_directory"

        # Make sure it doesn't exist
        if new_dir.exists():
            new_dir.rmdir()

        result = _check_output_path({"output_path": new_dir})

        assert new_dir.exists()
        assert result == new_dir

        # Cleanup
        new_dir.rmdir()

    def test_creates_nested_directories(self):
        """Test that nested directories are created."""
        nested_dir = tmp_folder / "level1" / "level2" / "level3"

        # Make sure it doesn't exist
        if nested_dir.exists():
            import shutil

            shutil.rmtree(tmp_folder / "level1")

        result = _check_output_path({"output_path": nested_dir})

        assert nested_dir.exists()
        assert result == nested_dir

        # Cleanup
        import shutil

        shutil.rmtree(tmp_folder / "level1")

    def test_existing_directory_unchanged(self):
        """Test that existing directory is returned unchanged."""
        result = _check_output_path({"output_path": tmp_folder})
        assert result == tmp_folder
        assert tmp_folder.exists()


class TestCheckScenarioName:
    """Tests for _check_scenario_name function."""

    def test_returns_scenario_name(self):
        """Test that scenario name is returned correctly."""
        result = _check_scenario_name({"scenario_name": "test_scenario"})
        assert result == "test_scenario"

    def test_raises_error_when_missing(self):
        """Test that ValueError is raised when scenario_name is missing."""
        with pytest.raises(ValueError, match="Required `scenario_name` not defined"):
            _check_scenario_name({})

    def test_raises_error_when_false(self):
        """Test that ValueError is raised when scenario_name is False."""
        with pytest.raises(ValueError, match="Required `scenario_name` not defined"):
            _check_scenario_name({"scenario_name": False})

    def test_raises_error_when_empty_string(self):
        """Test that ValueError is raised when scenario_name is empty string."""
        with pytest.raises(ValueError, match="Required `scenario_name` not defined"):
            _check_scenario_name({"scenario_name": ""})

    def test_accepts_various_string_formats(self):
        """Test that various string formats are accepted."""
        assert _check_scenario_name({"scenario_name": "simple"}) == "simple"
        assert (
            _check_scenario_name({"scenario_name": "with_underscore"})
            == "with_underscore"
        )
        assert _check_scenario_name({"scenario_name": "with-dash"}) == "with-dash"
        assert _check_scenario_name({"scenario_name": "CamelCase"}) == "CamelCase"
        assert (
            _check_scenario_name({"scenario_name": "with123numbers"})
            == "with123numbers"
        )

    def test_preserves_scenario_name_type(self):
        """Test that scenario name type is preserved."""
        result = _check_scenario_name({"scenario_name": "test"})
        assert isinstance(result, str)
