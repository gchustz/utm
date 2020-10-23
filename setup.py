import setuptools

with open("README.md","r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="utm",
    version="0.0.1",
    author="George Chustz",
    author_email="gchustz@tamu.edu",
    description="Python library for conversions between Latitude, Longitude, and Height and NED/ENU (in meters).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gchustz/utm",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.6'
)