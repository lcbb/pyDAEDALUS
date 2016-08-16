from StringIO import StringIO


def open_string_as_file(string):
    file = StringIO()
    file.write(string)
    file.seek(0)
    return file
