from src.metadata import LIMITATIONS


def test_metadata_limitations_include_thermal_sampling_warning() -> None:
    joined = " ".join(LIMITATIONS)
    assert "30-m grid" in joined
    assert "approximately 100 m" in joined
    assert "not 2-m air temperature" in joined
    assert "not a direct measurement of plant physiological stress" in joined

