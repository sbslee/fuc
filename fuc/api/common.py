"""
The ``common`` submodule is used by other ``fuc`` submodules such as
``pyvcf`` and ``pybed``. It also provides many day-to-day actions used in
the field of bioinformatics.
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

def sumstat(tp, fp, fn, tn):
    """Return various important statistics."""
    d = {
        'tpr': tp / (tp + fn), # sensitivity, recall, hit rate
        'tnr': tn / (tn + fp), # specificity, selectivity
        'ppv': tp / (tp + fp), # precision
        'npv': tn / (tn + fn),
        'fnr': fn / (fn + tp), # miss rate
        'fpr': fp / (fp + tn), # fall-out rate
        'fdr': fp / (fp + tp),
        'for': fn / (fn + tn),
        'acc': (tp + tn) / (tp + tn + fp + fn),
    }
    return d
