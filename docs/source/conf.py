# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import mock
sys.path.insert(0, os.path.abspath('../../slime'))


# -- Project information -----------------------------------------------------

project = 'slime'
copyright = '2019, Georgia Tsambos'
author = 'Georgia Tsambos'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# MOCK_MODULES = ['numpy', 'numpy.core', 'numpy.core.multiarray']
# for mod_name in MOCK_MODULES:
#     sys.modules[mod_name] = mock.Mock()
autodoc_mock_imports = ["msprime", "pyslim", "pandas"]

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo'
]

# Modules to be mocked.
# autodoc_mock_imports = ["numpy", "numpy.core.multiarray", "tskit", "pandas"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Master ("home") document
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# Show same sidebar on all pages.
# See https://stackoverflow.com/questions/18969093/how-to-include-the-toctree-in-the-sidebar-of-each-page
html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }

#
