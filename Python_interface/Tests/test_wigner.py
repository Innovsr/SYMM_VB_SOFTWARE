import pytest
from CISVB_GUI import Output


def test_wigner():

    dummy_root = None
    output = Output(dummy_root)

#    print(f"Test 1 result: {output.wigner(0, 5, 2)}")
#    print(f"Test 2 result: {output.wigner(1, 6, 2)}")

    assert output.wigner(0, 5, 2) == (1, 5, 15),f"expected (1, 5, 15), got {output.wigner(0, 5, 2)}"
    assert output.wigner(1, 6, 2) == (6, 5, 90),f"expected (6, 5, 90), got {output.wigner(1, 6, 2)}"

