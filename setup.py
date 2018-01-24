from distutils.core import setup

PACKAGE = "inmembrane"
DESCRIPTION = "A bioinformatic pipeline for proteome annotation \
to predict if a protein is exposed on the surface of a bacteria."
AUTHOR = "Andrew Perry & Bosco Ho"
AUTHOR_EMAIL = "ajperry@pansapiens.com"
URL = "http://github.com/boscoh/inmembrane"
# Must be a semantic version number. Also update inmembrane/__init__.py
VERSION = "0.95.0"  # __import__(PACKAGE).__version__

try:
    extra_requires = []
    from collections import OrderedDict
except:
    extra_requires.append("ordereddict")

setup(
    name=PACKAGE,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    description=DESCRIPTION,
    packages=['inmembrane', 'inmembrane.plugins',
              'inmembrane.protocols', 'inmembrane.tests'],
    # NOTE: some packaging filters are also in MANIFEST.in
    package_data={'inmembrane': ['protocols/*/*',
                                 'tests/*/*',
                                 'plugins/*/*'], },
    scripts=['inmembrane_scan'],
    # README, examples & docs are included via MANIFEST.in
    license='BSD',
    long_description=open('README.rst', 'rt').read(),
    install_requires=["BeautifulSoup >= 3.2.1",
                      "bs4",
                      "cssselect",
                      "lxml",
                      "requests >= 2.0.0",
                      "semantic_version",
                      "suds >= 0.4",
                      "twill == 0.9.1",
                      ] + extra_requires,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=False,
)
