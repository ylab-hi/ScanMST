# Troubleshooting

Common issues and solutions for ScanMST users.

## Installation Issues

### Python Version Errors

??? question "ModuleNotFoundError or ImportError after installation"

**Symptom**: `ModuleNotFoundError: No module named 'scanmst'`

**Cause**: Wrong Python environment or installation failed

**Solution**:

```bash
# Check Python version (must be 3.9 or higher )
python --version

# Verify pip is using correct Python
which pip
python -m pip --version

# Reinstall in current environment
python -m pip install --force-reinstall scanmst
```

### Dependency Conflicts

??? question "ERROR: pip's dependency resolver does not currently take into account all the packages"

**Symptom**: Pip reports dependency conflicts during installation

**Solution**:

```bash
# Install in a clean environment
python -m venv scanmst_env
source scanmst_env/bin/activate
pip install scanmst
```

## Runtime Issues

### BAM File Issues

??? question "PermissionError: [Errno 13] Permission denied"

    **Symptom**: Cannot read or write BAM files

    **Solution**:
    ```bash
    # Check file permissions
    ls -l input.bam

    # Add read permission
    chmod +r input.bam
    ```

## General Help

### Enable Verbose Logging

For debugging, enable detailed output:

```bash
scanmst ... --log-level trace
```

### Check ScanMST Version

Ensure you're using the latest version:

```bash
# Check current version
scanmst --version

# Update to latest
pip install --upgrade scanmst
```

### System Information

Collect system info for bug reports:

```bash
# Python version
python --version

# ScanMST version
scanmst --version
```

## Getting Further Help

If your issue isn't covered here:

1. **Check existing issues**: [GitHub Issues](https://github.com/ylab-hi/ScanMST/issues)
2. **Search discussions**: [GitHub Discussions](https://github.com/ylab-hi/ScanMST/discussions)
3. **Open a new issue**: Include:
   - ScanMST version (`scanmst --version`)
   - Python version (`python --version`)
   - Operating system
   - Complete error message
   - Minimal reproducible example

!!! tip "Before Opening an Issue"

    - Update to the latest version
    - Try with sample data (`tests/data/isoseq_test.bam`)
    - Include full error traceback
    - Describe what you expected vs. what happened
