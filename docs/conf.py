import os
import sys

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(".."))

project = "patchsim"
author = ""

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]

# MyST configuration
myst_enable_extensions = [
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "colon_fence",
]

# Register a CSV lexer name for Pygments; use named lookup and fall back to text.
from sphinx.highlighting import lexers as sphinx_lexers
from pygments import lexers as pyg_lexers
try:
    lexer = pyg_lexers.get_lexer_by_name("csv")
except Exception:
    lexer = pyg_lexers.get_lexer_by_name("text")

sphinx_lexers["csv"] = lexer
