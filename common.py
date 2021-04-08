import pathlib

def is_numeric(s):
    """Return True if the string is numeric."""
    try:
        int(s)
        return True
    except ValueError:
        try:
            float(s)
            return True
        except ValueError:
            return False

def fuc_dir():
    """Return the path to the fuc directory."""
    return pathlib.Path(__file__).parent.absolute()
