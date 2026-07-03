# GitHub Copilot Instructions for raycloudtools Repository

<!--
Adapted from
https://github.com/csiro-internal/cps/blob/master/.github/copilot-instructions.md.
-->

## Language and Spelling

Always use Australian English spelling in all generated code, comments, documentation, and variable names:

-   Write "colour" not "color"
-   Write "centre" not "center"
-   Write "licence" not "license" (noun)
-   Write "organisation" not "organization"
-   Write "categorise" not "categorize"
-   Write "utilise" not "utilize"
-   Write "artefact" not "artifact"

## Coding Standards

### Apply Language-Specific Standards

Generate code following these standards:

-   **Python**: Follow PEP8 guidelines
-   **C++**: Follow Google C++ style guide

### Python Code Quality Standards

Python code is checked with both **flake8** (code style) and **pydocstyle** (docstring style) using ament configurations. Ensure generated Python code adheres to these standards.

**Flake8 compliance (ament_flake8 configuration):**
The repository uses ament_flake8 with these specific settings:

-   **Line length**: 99 characters maximum (not the PEP 8 default of 79)
-   **Import order style**: Google style
-   **Show source**: Enabled for error context
-   **Statistics**: Enabled for error summaries

**Flake8 ignored rules:**

-   **B902**: Invalid first argument used for instance method
-   **C816**: Missing trailing comma in Python 3.6+
-   **D100-D107**: Missing docstring rules (handled by pydocstyle instead)
-   **D203**: 1 blank line required before class docstring
-   **D212**: Multi-line docstring summary should start at the first line
-   **D404**: First word of the docstring should not be "This"
-   **I202**: Additional newline in a group of imports

**Pydocstyle compliance (ament_pep257 configuration):**
The repository uses the **ament** convention for pydocstyle, which ignores these rules:

-   **D100-D107**: Missing docstring rules (docstrings are optional)
-   **D203**: 1 blank line required before class docstring
-   **D212**: Multi-line docstring summary should start at the first line
-   **D404**: First word of the docstring should not be "This"

**When generating Python code:**

-   Use **99 character line length** (not 79)
-   Follow Google import order style
-   Follow PEP 8 style guidelines for indentation, spacing, and naming
-   Avoid unused imports and variables
-   Use proper naming conventions
-   Maintain consistent code formatting

**When writing docstrings:**

-   Use NumPy style, respecting D213, D406, D407, D413
-   Use Australian English spelling
-   Follow PEP 257 formatting for docstrings that are present
-   Multi-line docstring summaries may start on the second line
-   Docstrings may begin with "This" if appropriate
-   No blank line required before class docstrings
-   Docstrings are optional but recommended for public APIs

### Copyright Headers

Always include copyright headers in generated source files:

**Python files:**

```python
# Copyright 2020-2026 Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# All Rights Reserved
```

**C++ files:**

```cpp
// Copyright 2020-2026 Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// All Rights Reserved
```

## Third-Party Code Integration

When suggesting third-party code:

1. Prioritise system packages (OS packages, PyPI)
2. When copying code, generate appropriate licence documentation including:
    - Original licence information
    - Upstream repository references
    - Commit hash or retrieval date

## Testing Code Generation

Generate test code that includes:

-   Tests appropriate to module's intended use and reuse level
-   Verification that functionality works as intended
-   Testing framework integration
-   Regression tests for bug fixes

## Code Quality Guidelines

Generate code that:

-   Uses descriptive, functional names over branded names
-   Balances documentation thoroughness with practicality
-   Respects third-party code licences
-   Follows adequate testing practices for the use case
-   Avoids unnecessary submodule creation unless genuinely needed
