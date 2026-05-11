from setuptools import setup, find_packages

DISTNAME = 'nanopore_10x_multiome'
VERSION = '0.2.0'
MAINTAINER = 'Chris Jackson'
MAINTAINER_EMAIL = 'cjackson@nygenome.org'
LICENSE = 'MIT'

setup(
    name=DISTNAME,
    version=VERSION,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    license=LICENSE,
    packages=find_packages(include=['nanopore_10x_multiome', "nanopore_10x_multiome.*"]),
    install_requires=[
      'numpy',
      'scipy',
      'pandas',
      'joblib',
      'parasail',
      'regex',
      'tqdm',
      'pysam'
    ],
    zip_safe=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta"
    ]
)
