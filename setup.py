import setuptools
  
with open("README.md", "r") as fh:
    description = fh.read()
  
setuptools.setup(
    name="dsolve",
    version="0.0.1",
    author="Marc de la Barrera i Bardalet",
    author_email="mbarrera@mit.edu",
    packages=["test_package"],
    description="A package to solve systems of dynamic equations in python",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/marcdelabarrera/dsolve",
    license='MIT',
    python_requires='>=3.10',
    install_requires=[]
)