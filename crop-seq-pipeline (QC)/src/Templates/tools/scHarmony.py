def retrieve_protocol(x):
    """
    retrieves cellline
    """
    try:
        split = x.split("-")
        return "-".join(split[2:])
    except (AttributeError, IndexError):
        return x  