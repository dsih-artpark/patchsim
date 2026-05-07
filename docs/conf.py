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

# Register a CSV lexer name for Pygments; fall back to plain text if unavailable.
try:
    from pygments.lexers.data import CsvLexer
    from sphinx.highlighting import lexers

    lexers["csv"] = CsvLexer()
except Exception:
    from pygments.lexers.text import TextLexer
    from sphinx.highlighting import lexers

    lexers["csv"] = TextLexer()
