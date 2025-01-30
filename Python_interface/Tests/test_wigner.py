#import pytest
#from CISVB_GUI import Output
#
#
#def test_wigner():
#
#    dummy_root = None
#    output = Output(dummy_root)
#
##    print(f"Test 1 result: {output.wigner(0, 5, 2)}")
##    print(f"Test 2 result: {output.wigner(1, 6, 2)}")
#
#    assert output.wigner(0, 5, 2) == (1, 5, 15),f"expected (1, 5, 15), got {output.wigner(0, 5, 2)}"
#    assert output.wigner(1, 6, 2) == (6, 5, 90),f"expected (6, 5, 90), got {output.wigner(1, 6, 2)}"

import pytest
import ctypes
from ctypes.util import find_library
from CISVB_GUI import Output

import os

print("Current Working Directory:", os.getcwd())
print("Directory Contents:", os.listdir(os.getcwd()))
print("Library Exists:", os.path.exists("./symm_str.so"))

#def test_library_loading():
#    lib_name = "symm_str"  # Replace with the actual library name
#    print("Find Library:", find_library(lib_name))
#
#    try:
#        ctypes.CDLL(f"./{lib_name}.so")
#        print("Library loaded successfully.")
#    except OSError as e:
#        pytest.fail(f"Error loading library: {e}")
#
#def test_wigner():
#    dummy_root = None
#    output = Output(dummy_root)
#
#    assert output.wigner(0, 5, 2) == (1, 5, 15), f"expected (1, 5, 15), got {output.wigner(0, 5, 2)}"
#    assert output.wigner(1, 6, 2) == (6, 5, 90), f"expected (6, 5, 90), got {output.wigner(1, 6, 2)}"
#
