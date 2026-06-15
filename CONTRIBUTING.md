# Contributing to DamageScanner

Welcome to `DamageScanner`! We appreciate your interest in contributing. By contributing, you help make `DamageScanner` better for everyone. Please review the following guidelines before making your contribution.

## Code of Conduct

This project and everyone participating in it is governed by the
[Code of Conduct](https://github.com/VU-IVM/DamageScanner/blob/main/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. Please report unacceptable behavior
to elco.koks@vu.nl.

## Asking Questions and Reporting Issues

If you encounter any bugs or issues while using `DamageScanner`, please report them by opening an *issue* in the GitHub repository. Be sure to provide detailed information about the problem, such as steps to reproduce it, including operating system and Python version.

If you have any suggestions for improvements, or questions, please also raise an *issue*.

## Contributing to Code

### Setting Up Your Development Environment

1. Clone the repository to your local machine:

    ```bash
    git clone https://github.com/VU-IVM/DamageScanner.git
    ```

2. Set up the development environment using `uv` (recommended):

    ```bash
    uv sync --all-groups
    ```

    Alternatively, using conda:

    ```bash
    conda env create -f environment.yml
    conda activate ds-test
    ```

3. Checkout a new branch for your changes from the main branch:

    ```bash
    git checkout -b feature/your-feature-name
    ```

4. Make your changes.

## Submitting a Pull Request

Once you have made your changes and are ready to contribute, follow these steps to submit a pull request:

1. Ensure your code passes linting and formatting checks:

    ```bash
    ruff check .
    ruff format .
    ```

2. Run the test suite:

    ```bash
    pytest tests/
    ```

3. Push your changes back to origin:

    ```bash
    git push origin feature/your-feature-name
    ```

4. Create a pull request to merge your branch into the main branch.

Provide a clear description of your changes in the pull request.

#### Pull Request Guidelines
- Write clear and concise commit messages.
- Test your changes thoroughly before submitting a pull request.
- If the pull request adds functionality, the docs should also be updated. Improving documentation helps users better understand how to use `DamageScanner`.

## Review Process

All pull requests will undergo a review process before being accepted. Reviewers may provide feedback or request changes to ensure the quality of the codebase.