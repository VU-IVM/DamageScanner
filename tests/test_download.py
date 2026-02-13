import pytest
from pathlib import Path
from unittest.mock import patch

from .setup import tmp_folder

from damagescanner.download import (
    _create_gf_download_url,
    _download_file,
    get_country_geofabrik,
    get_region_geofabrik,
    get_planet_file,
)
from damagescanner.config import GEOFABRIK_URL, PLANET_URL


class TestCreateGfDownloadUrl:
    """Tests for _create_gf_download_url function."""

    def test_creates_pbf_url(self):
        """Test creating a PBF download URL."""
        url = _create_gf_download_url("NLD", "pbf")

        assert GEOFABRIK_URL in url
        assert url.endswith(".osm.pbf")
        assert "netherlands" in url.lower()

    def test_creates_shp_url(self):
        """Test creating a SHP download URL."""
        url = _create_gf_download_url("NLD", "shp")

        assert GEOFABRIK_URL in url
        assert url.endswith("-free.shp.zip")
        assert "netherlands" in url.lower()

    def test_invalid_iso3_raises_keyerror(self):
        """Test that invalid ISO3 raises KeyError."""
        with pytest.raises(KeyError):
            _create_gf_download_url("INVALID", "pbf")

    def test_russia_raises_specific_error(self):
        """Test that Russia ISO3 raises specific KeyError with guidance."""
        with pytest.raises(KeyError, match="Russia comes in two files"):
            _create_gf_download_url("RUS", "pbf")

    def test_invalid_format_raises_notimplemented(self):
        """Test that invalid format raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="Invalid file format"):
            _create_gf_download_url("NLD", "geojson")

    def test_various_countries(self):
        """Test URL creation for various countries."""
        test_cases = [
            ("DEU", "germany"),
            ("FRA", "france"),
            ("GBR", "britain"),
            ("USA", "us"),
            ("JPN", "japan"),
        ]

        for iso3, expected_name in test_cases:
            try:
                url = _create_gf_download_url(iso3, "pbf")
                # Just verify it creates a valid-looking URL
                assert url.endswith(".osm.pbf")
            except KeyError:
                # Some countries might not be in DICT_GEOFABRIK
                pass


class TestDownloadFile:
    """Tests for _download_file function."""

    @patch("damagescanner.download.urllib.request.urlretrieve")
    def test_downloads_new_file(self, mock_urlretrieve):
        """Test that new file is downloaded."""
        filepath = tmp_folder / "test_download.pbf"

        # Ensure file doesn't exist
        if filepath.exists():
            filepath.unlink()

        _download_file("http://example.com/file.pbf", filepath, overwrite=True)

        mock_urlretrieve.assert_called_once_with(
            "http://example.com/file.pbf", filepath
        )

    @patch("damagescanner.download.urllib.request.urlretrieve")
    def test_overwrites_existing_file(self, mock_urlretrieve):
        """Test that existing file is overwritten when overwrite=True."""
        filepath = tmp_folder / "test_existing.pbf"

        # Create existing file
        filepath.touch()

        _download_file("http://example.com/file.pbf", filepath, overwrite=True)

        mock_urlretrieve.assert_called_once()

        # Cleanup
        if filepath.exists():
            filepath.unlink()

    @patch("damagescanner.download.urllib.request.urlretrieve")
    def test_skips_existing_file_no_overwrite(self, mock_urlretrieve):
        """Test that existing file is skipped when overwrite=False."""
        filepath = tmp_folder / "test_skip.pbf"

        # Create existing file
        filepath.touch()

        _download_file("http://example.com/file.pbf", filepath, overwrite=False)

        mock_urlretrieve.assert_not_called()

        # Cleanup
        if filepath.exists():
            filepath.unlink()


class TestGetCountryGeofabrik:
    """Tests for get_country_geofabrik function."""

    @patch("damagescanner.download._download_file")
    def test_returns_path(self, mock_download):
        """Test that function returns a Path object."""
        result = get_country_geofabrik("NLD", save_path=tmp_folder)

        assert isinstance(result, Path)

    @patch("damagescanner.download._download_file")
    def test_pbf_format(self, mock_download):
        """Test downloading in PBF format."""
        result = get_country_geofabrik("NLD", file_format="pbf", save_path=tmp_folder)

        assert str(result).endswith(".osm.pbf")
        mock_download.assert_called_once()

    @patch("damagescanner.download._download_file")
    def test_shp_format(self, mock_download):
        """Test downloading in SHP format."""
        result = get_country_geofabrik("NLD", file_format="shp", save_path=tmp_folder)

        assert str(result).endswith("-free.shp.zip")
        mock_download.assert_called_once()

    @patch("damagescanner.download._download_file")
    def test_creates_directory(self, mock_download):
        """Test that save directory is created if needed."""
        new_dir = tmp_folder / "new_osm_dir"

        # Remove if exists
        if new_dir.exists():
            new_dir.rmdir()

        get_country_geofabrik("NLD", save_path=new_dir)

        assert new_dir.exists()

        # Cleanup
        if new_dir.exists():
            new_dir.rmdir()

    def test_invalid_iso3_raises_error(self):
        """Test that invalid ISO3 raises KeyError."""
        with pytest.raises(KeyError):
            get_country_geofabrik("INVALID", save_path=tmp_folder)

    @patch("damagescanner.download._download_file")
    def test_overwrite_parameter_passed(self, mock_download):
        """Test that overwrite parameter is passed to download function."""
        get_country_geofabrik("NLD", save_path=tmp_folder, overwrite=True)

        # Check that overwrite=True was passed
        call_args = mock_download.call_args
        assert (
            call_args[1].get(
                "overwrite", call_args[0][2] if len(call_args[0]) > 2 else None
            )
            == True
        )


class TestGetRegionGeofabrik:
    """Tests for get_region_geofabrik function."""

    @patch("damagescanner.download._download_file")
    def test_returns_path(self, mock_download):
        """Test that function returns a Path object."""
        result = get_region_geofabrik("europe", save_path=tmp_folder)

        assert isinstance(result, Path)

    @patch("damagescanner.download._download_file")
    def test_correct_url_format(self, mock_download):
        """Test that correct URL is constructed."""
        get_region_geofabrik("europe", save_path=tmp_folder)

        call_args = mock_download.call_args[0]
        download_url = call_args[0]

        assert "europe-latest.osm.pbf" in download_url

    @patch("damagescanner.download._download_file")
    def test_various_regions(self, mock_download):
        """Test downloading various regions."""
        regions = ["europe", "africa", "asia", "north-america", "south-america"]

        for region in regions:
            mock_download.reset_mock()
            result = get_region_geofabrik(region, save_path=tmp_folder)

            assert isinstance(result, Path)
            assert region in str(result)

    @patch("damagescanner.download._download_file")
    def test_case_insensitive(self, mock_download):
        """Test that region name is converted to lowercase."""
        get_region_geofabrik("EUROPE", save_path=tmp_folder)

        call_args = mock_download.call_args[0]
        download_url = call_args[0]

        assert "europe-latest" in download_url
        assert "EUROPE" not in download_url


class TestGetPlanetFile:
    """Tests for get_planet_file function."""

    @patch("damagescanner.download._download_file")
    def test_returns_path(self, mock_download):
        """Test that function returns a Path object."""
        save_path = tmp_folder / "planet-latest.osm.pbf"
        result = get_planet_file(save_path=save_path)

        assert isinstance(result, Path)

    @patch("damagescanner.download._download_file")
    def test_uses_planet_url(self, mock_download):
        """Test that Planet URL is used."""
        save_path = tmp_folder / "planet-latest.osm.pbf"
        get_planet_file(save_path=save_path)

        call_args = mock_download.call_args[0]
        download_url = call_args[0]

        assert download_url == PLANET_URL

    @patch("damagescanner.download._download_file")
    def test_overwrite_parameter(self, mock_download):
        """Test that overwrite parameter is passed."""
        save_path = tmp_folder / "planet-latest.osm.pbf"
        get_planet_file(save_path=save_path, overwrite=True)

        call_args = mock_download.call_args
        # Check overwrite is True (either as kwarg or positional arg)
        assert call_args[0][2] == True or call_args[1].get("overwrite") == True


class TestIntegration:
    """Integration tests (these actually create files but don't download)."""

    @patch("damagescanner.download.urllib.request.urlretrieve")
    def test_full_workflow_country(self, mock_urlretrieve):
        """Test full workflow for country download."""

        # Mock urlretrieve to create an empty file
        def create_file(url, filepath):
            Path(filepath).touch()

        mock_urlretrieve.side_effect = create_file

        result = get_country_geofabrik("NLD", save_path=tmp_folder, overwrite=True)

        assert result.exists()

        # Cleanup
        if result.exists():
            result.unlink()

    @patch("damagescanner.download.urllib.request.urlretrieve")
    def test_full_workflow_region(self, mock_urlretrieve):
        """Test full workflow for region download."""

        def create_file(url, filepath):
            Path(filepath).touch()

        mock_urlretrieve.side_effect = create_file

        result = get_region_geofabrik("europe", save_path=tmp_folder, overwrite=True)

        assert result.exists()

        # Cleanup
        if result.exists():
            result.unlink()
