from src.quality_mask import (
    QA_CIRRUS_BIT,
    QA_CLOUD_BIT,
    QA_CLOUD_SHADOW_BIT,
    QA_DILATED_CLOUD_BIT,
    QA_FILL_BIT,
    QA_SNOW_BIT,
    QA_WATER_BIT,
    qa_pixel_is_clear,
    qa_radsat_is_clear,
)


def test_qa_pixel_clear_when_no_blocked_bits() -> None:
    assert qa_pixel_is_clear(0, mask_water=True)


def test_qa_pixel_masks_required_bits() -> None:
    blocked_bits = [
        QA_FILL_BIT,
        QA_DILATED_CLOUD_BIT,
        QA_CIRRUS_BIT,
        QA_CLOUD_BIT,
        QA_CLOUD_SHADOW_BIT,
        QA_SNOW_BIT,
        QA_WATER_BIT,
    ]
    for bit in blocked_bits:
        assert not qa_pixel_is_clear(1 << bit, mask_water=True)


def test_water_can_be_kept_when_configured() -> None:
    assert qa_pixel_is_clear(1 << QA_WATER_BIT, mask_water=False)


def test_qa_radsat_requires_zero() -> None:
    assert qa_radsat_is_clear(0)
    assert not qa_radsat_is_clear(1)

