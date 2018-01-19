from tempfile import NamedTemporaryFile


def make_tempfile(contents, suffix):
    """Returns a NamedTemporaryFile with text contents

    Parameters
    ----------
    contents: string
        The text to be written to the new temporary file
    suffix: string
        The suffix of the file path (e.g. ".fasta")
    """
    with NamedTemporaryFile(mode='wt', suffix=suffix, delete=False) as f:
        f.file.write(contents)
    return f
