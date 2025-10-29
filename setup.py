from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="bsmotif",
    version="0.0.1",
    author="Vladimir Raditsa",
    author_email="raditsavv@bionet.nsc.ru", 
    description="A pipeline for hierarchical classification of motifs based on the DSDs of their transcription factors",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/vladiMIR28161/BSMotif",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.22.0",
        "openpyxl>=3.1.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
    'console_scripts': [
        'bsmotif = bsmotif.main:main',
    ],
    },
    python_requires=">=3.8",
)