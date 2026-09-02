# Installation

Get ScanMST installed on your system in just a few minutes.

## Prerequisites

!!! info "Requirements" - **Python**: 3.9 or higher - **Operating System**: Linux, macOS

## Installation Methods

Choose your preferred installation method:

=== "pip"

    The easiest way to install ScanMST:

    ```bash
    pip install scanmst
    ```

    Verify the installation:

    ```bash
    scanmst --version
    ```

=== "conda"

    Install using conda:

    ```bash
    conda install -c conda-forge scanmst
    ```

    Or create a new environment:

    ```bash
    conda create -n scanmst python=3.10
    conda activate scanmst
    pip install scanmst
    ```

=== "uv"

    Using the fast `uv` package manager:

    ```bash
    uv pip install scanmst
    ```

    Or with a virtual environment:

    ```bash
    uv venv
    source .venv/bin/activate
    uv pip install scanmst
    ```

=== "from source"

    For development or the latest features:

    ```bash
    # Clone the repository
    git clone https://github.com/ylab-hi/ScanMST.git
    cd ScanMST

    # Install with Poetry
    poetry install

    # Verify installation
    scanmst --version
    ```

## Verification

Confirm ScanMST is installed correctly:

```bash
# Check version
scanmst --version

# View available commands
scanmst --help
```

## Troubleshooting Installation

### Common Issues

??? question "ImportError: No module named 'scanmst'"

    **Solution**: Ensure you've activated the correct Python environment:

    ```bash
    # Check which Python is being used
    which python

    # Reinstall in the current environment
    pip install --force-reinstall scanmst
    ```

??? question "Permission denied errors"

    **Solution**: Install in user space without sudo:

    ```bash
    pip install --user scanmst
    ```

For more issues, see the [Troubleshooting Guide](troubleshooting.md).

## Next Steps

Now that ScanMST is installed, try the [Quick Start](quick-start.md) tutorial to run your first prediction!
