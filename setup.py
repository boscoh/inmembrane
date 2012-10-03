from distutils.core import setup

PACKAGE = "inmembrane"
NAME = ""
DESCRIPTION = ""
AUTHOR = "Andrew Perry & Bosco Ho"
URL = "http://github.com/boscoh/inmembrane"
VERSION = __import__(PACKAGE).__version__

setup(
    name='inmembrane',
    version=VERSION,
    author=AUTHOR,
    url=URL,
    description=DESCRIPTION,
    packages=['plugins','protocols','tests'],
    package_data={'protocols': ['*/*'],
                  'tests': ['*/*'],
                  'plugins': ['*/*'], },
    scripts=['inmembrane.py', 'run_example.py','run_tests.py'] + ['helpers.py'],
    # README, examples & docs are included via MANIFEST.in
    license='BSD',
    long_description=open('README.markdown').read(),
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
