"""
The common module is used by other ``fuc`` modules such as `VcfFrame` and
`BedFrame`. It also provides many useful methods.
"""

import pathlib
import re
from difflib import SequenceMatcher

def get_script_name(script):
    """Return the script name."""
    return pathlib.Path(script).stem

def fuc_dir():
    """Return the path to the fuc directory."""
    return pathlib.Path(__file__).parent.parent.parent.absolute()

def get_similarity(a, b):
    """Return a value from 0 to 1 representing how similar two strings are."""
    return SequenceMatcher(None, a, b).ratio()

def is_similar(a, b, threshold=0.9):
    """Return True if the similarity is equal to or greater than threshold."""
    return get_similarity(a, b) >= threshold

def get_most_similar(a, l):
    """Return the most similar string in a list."""
    s = [get_similarity(a, x) for x in l]
    m = max(s)
    i = [i for i, x in enumerate(s) if x == m][0]
    return l[i]
