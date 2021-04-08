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
