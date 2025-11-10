from importlib import metadata
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# Project information
project = "morfeus"
copyright = "2021, Kjell Jorner"
author = "Kjell Jorner"
release = metadata.version("morfeus-ml")
version = ".".join(metadata.version("morfeus-ml").split(".")[:-1])

# General configuration
extensions = [
    "sphinx_copybutton",
    "sphinx_inline_tabs",
    "sphinx.ext.autodoc",
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
]
master_doc = "index"
source_suffix = ".rst"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_theme_options = {
    "light_logo": "logo-text-light.svg",
    "dark_logo": "logo-text-dark.svg",
    "sidebar_hide_name": True,
}
html_title = "ᴍᴏʀғᴇᴜs"
html_static_path = ["_static"]
html_favicon = "_static/logo-icon-dark.svg"

# Extension configuration

# Autodoc
autodoc_mock_imports = [
    "ase",
    "dftd4",
    "pymeshfix",
    "openbabel",
    "pyvista",
    "pyvistaqt",
    "rdkit",
    "vtk",
    "xtb",
]
autodoc_typehints = "description"

# Bibtex
bibtex_bibfiles = ["refs.bib"]
bibtex_bibliography_header = ".. rubric:: References"
bibtex_footbibliography_header = bibtex_bibliography_header
bibtex_default_style = "unsrt"
