# -*- coding: utf-8 -*-
"""
Setup file for kartopy package.

Created by: Kat Sejan 19th February 2020.

"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="kryopy", # Replace with your own username
    version="0.0.1",
    author="Katarztna M. Sejan",
    author_email="k.m.sejan@uu.nl",
    description="Processing of CryoSat2 radar altimetry data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CryoKatOnIce/kryopy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)