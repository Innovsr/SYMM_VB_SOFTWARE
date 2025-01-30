import os
import ctypes
import pytest
from ctypes.util import find_library

def test_library_loading():
    lib_name = "symm_str"  # Replace with actual library name
    lib_file = f"{lib_name}.so"

    # Print working directory
    cwd = os.getcwd()
    print(f"Current Working Directory: {cwd}")

    # Print directory contents
    dir_contents = os.listdir(cwd)
    print(f"Directory Contents: {dir_contents}")

    # Check if the shared library exists in the directory
    lib_path = os.path.abspath(lib_file)
    file_exists = os.path.exists(lib_path)
    print(f"Library Exists: {file_exists}, Expected Path: {lib_path}")

    if not file_exists:
        pytest.fail(f"Shared library {lib_file} not found at {lib_path}")

    # Use ctypes.util.find_library to check if the system can find it
    found_lib = find_library(lib_name)
    print(f"find_library('{lib_name}') returned: {found_lib}")

    if found_lib is None:
        print(f"Warning: {lib_name} not found via find_library().")

    # Try loading the shared library with ctypes
    try:
        ctypes.CDLL(lib_path)
        print("Library loaded successfully.")
    except OSError as e:
        pytest.fail(f"Error loading library: {e}")
