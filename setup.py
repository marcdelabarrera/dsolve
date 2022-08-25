from setuptools import setup


with open("README.md","r") as fh:
    long_description = fh.read()

setup(
    name = 'dsolve',
    version = '0.0.2',
    description = 'Solver of dynamic equations with forward looking variables',
    long_description_content_type='text/markdown',
    long_description='test',
    py_modules = ["dsolve.atoms", "dsolve.expressions", "dsolve.solvers"],
    package_dir={'':'src'},
    install_requires = ["scipy >= 1.9.0"],
    extras_require={"dev":["pytest>=7.1.2",],},
    classifiers =[
        "Programming Language :: Python :: 3.10"
    ]
)