import pytest
import tkinter as tk
from CISVB_GUI import Output


def test_set_info():
    dummy_root = tk.Tk()
    output = Output(dummy_root)
    various_qualities = ['1 2 3 4', '0 5 5 2', '6 8 9 0']
    overall_qualities = [4, 5, 8]
    rumers = ['R,-,-']
    linear_indset = ['set', 1, 1, 3]

    print(various_qualities, rumers)
    result = output.set_info(
            various_qualities, 
            overall_qualities, 
            rumers, 
            linear_indset
            ) 
    print('result', result)
    assert result == 12, f"expected 12 but got {result}"

