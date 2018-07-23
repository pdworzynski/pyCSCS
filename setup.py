import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyCSCS",
    version="0.0.1",
    author="Lasse Buur Rasmussen",
    author_email="lassebuurrasmussen@gmail.com",
    description="A python implementation of the CSCS distance",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/askerdb/pyCSCS",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
