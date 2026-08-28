"""
Pytest configuration for rnadnavar tests.

This conftest.py adds the bin directory to the Python path so that
modules like vcf_stats and vcf_utils can be imported in tests.
"""

import os
import sys

bin_path = os.path.join(os.path.dirname(__file__), "..", "bin")
if bin_path not in sys.path:
    sys.path.insert(0, bin_path)
