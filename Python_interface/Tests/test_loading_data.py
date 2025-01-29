import pytest
from unittest.mock import mock_open, patch
from CISVB_GUI import Output  # Replace with the actual module name

def test_load_data_from_file():
    mock_file_content = """
    | Data 1
    | Data 2
    Set_number 1
    | Data 3
    | Data 4
    Set_number 2
    """
    output = Output("/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface/C5H5_output")

    with patch("builtins.open", mock_open(read_data=mock_file_content)):
        result = output.load_output_file("dummy_file.txt")
        print('result',result)

    assert len(result) == 2
    assert result[0] == ["| Data 1", "| Data 2"]
    assert result[1] == ["| Data 3", "| Data 4"]
