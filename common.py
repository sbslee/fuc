import pathlib
import re

def fuc_dir():
    """Return the path to the fuc directory."""
    return pathlib.Path(__file__).parent.absolute()

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

def parse_where(where):
    """Parse the SQLite WHERE clause."""
    key = where.split(' IN ')[0][1:-1]
    values = re.findall("'([^']*)'", where.split(' IN ')[1])
    return key, values
