import os
import sys

from pygments import lexers as pyg_lexers
from sphinx.highlighting import lexers as sphinx_lexers  # type: ignore

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(".."))

project = "PatchSim"
author = ""
html_title = "PatchSim"

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
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/dsih-artpark/patchsim",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        }
    ]
}

# MyST configuration
myst_enable_extensions = [
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "dollarmath",
]
myst_heading_anchors = 3

# Register a CSV lexer name for Pygments; use named lookup and fall back to text.
try:
    lexer = pyg_lexers.get_lexer_by_name("csv")
except Exception:
    lexer = pyg_lexers.get_lexer_by_name("text")

sphinx_lexers["csv"] = lexer
