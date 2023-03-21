# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
sys.path.append('/bbkinghome/mbarrera/git_supply/dsolve/src/dsolve')
sys.path.append('/bbkinghome/mbarrera/git_supply/dsolve')

project = 'dsolve'
copyright = '2023, Marc de la Barrera i Bardalet'
author = 'Marc de la Barrera i Bardalet'
release = '0.1.'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
    'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']


napoleon_numpy_docstring = True