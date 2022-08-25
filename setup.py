from setuptools import setup


with open("README.md","r") as fh:
    long_description = fh.read()

setup(
    name = 'dsolve',
    version = '0.0.3',
    description = 'Solver of dynamic equations with forward looking variables',
    long_description = long_description,
    long_description_content_type='text/markdown',
    py_modules = ["dsolve.atoms", "dsolve.expressions", "dsolve.solvers"],
    package_dir={'':'src'},
    author='Marc de la Barrera i Bardalet',
    author_email='mbarrera@mit.edu',
    install_requires = ["scipy >= 1.9.0"],
    extras_require={"dev":["pytest>=7.1.2",],},
    classifiers =[
        "Programming Language :: Python :: 3.10"
    ]
)