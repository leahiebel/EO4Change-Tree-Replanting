from src.temperature import st_b10_dn_to_celsius_value


def test_st_b10_dn_to_celsius_value() -> None:
    assert st_b10_dn_to_celsius_value(0) == -124.14999999999998
    assert round(st_b10_dn_to_celsius_value(44947), 3) == 29.501

