import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BEstimate",
    version="0.0.1",
    author="Cansu Dincer",
    author_email="cd@sanger.ac.uk",
    description="Estimate and Analyse Base Editor Target Site - BEstimate",
    long_description_content_type="text/markdown",
    url="https://https://github.com/CansuDincer/BEstimate/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Independent",
    ],
    python_requires='>=3.8',
)
