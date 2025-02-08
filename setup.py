import setuptools
import io

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = io.open("requirements.txt").read().split("\n")

setuptools.setup(
    name="FGOT",
    version="0.0.9",
    author="Yang Chenghui",
    license = 'MIT License',  
    author_email="2022282110116@whu.edu.cn",
    description="FGOT is a tool for uncovering cellular heterogeneity and their associated transcriptional regulatory links.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/YCH-bioinf/FGOT-master",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = install_requires
)

