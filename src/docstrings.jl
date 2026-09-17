# DocStringExtensions templates: set defaults for every documented
# object so individual docstrings stay focused on prose.

@template (FUNCTIONS, METHODS, MACROS) = """
    $(TYPEDSIGNATURES)
$(DOCSTRING)
"""

@template (TYPES) = """
    $(TYPEDEF)
$(DOCSTRING)

---
## Fields
$(TYPEDFIELDS)
"""

@template MODULES = """
$(DOCSTRING)

---
## Exports
$(EXPORTS)
---
## Imports
$(IMPORTS)
"""
