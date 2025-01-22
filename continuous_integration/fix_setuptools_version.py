import os
import platform

print(os.getcwd())

# temporary fix for https://github.com/pypa/distutils/issues/283
if "PyPy" in platform.python_implementation():
    with open("sdpa-python/pyproject.toml", "w") as file:
        file.write('[build-system]\n')
        file.write('requires = [')
        file.write('    "setuptools<72.2.0"')
        file.write(']')
