import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyCSCS",
    version="0.0.2",
    author="Asker Brejnrod, Piotr Dworzynski",
    author_email="brejnrod@sund.ku.dk",
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
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
)

