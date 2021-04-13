import pathlib
import re

def get_script_name(script):
    """Return the script name."""
    return pathlib.Path(script).stem

def fuc_dir():
    """Return the path to the fuc directory."""
    return pathlib.Path(__file__).parent.parent.parent.absolute()

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

def parse_condition(condition):
    """Parse one condition in the SQLite WHERE clause."""
    k = condition.split(' IN ')[0][1:-1]
    v = re.findall("'([^']*)'", condition.split(' IN ')[1])
    return {k: v}

def parse_where(where):
    """Parse the SQLite WHERE clause."""
    conditions = where.split(' AND ')
    results = {}
    for condition in conditions:
        results = {**results, **parse_condition(condition)}
    return results
