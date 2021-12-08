# ScanNLS

<div align="center">

[![PyPI](https://img.shields.io/pypi/v/ScanNLS.svg)](https://pypi.org/project/ScanNLS/ "PyPI")
[![Status](https://img.shields.io/pypi/status/ScanNLS.svg)](https://pypi.org/project/ScanNLS/ "Status")
[![Python Version](https://img.shields.io/pypi/pyversions/ScanNLS)](https://pypi.org/project/ScanNLS "Python Version")
[![License](https://img.shields.io/pypi/l/ScanNLS)](https://opensource.org/licenses/MIT "License")
[![Read the Docs](https://img.shields.io/readthedocs/ScanNLS/latest.svg?label=Read%20the%20Docs)](https://ScanNLS.readthedocs.io/ "Read the documentation at https://ScanNLS.readthedocs.io/")
[![Tests](https://github.com/ylab-hi/ScanNLS/workflows/Tests/badge.svg)](https://github.com/ylab-hi/ScanNLS/actions?workflow=Tests "Tests")
[![Codecov](https://codecov.io/gh/ylab-hi/ScanNLS/branch/main/graph/badge.svg)](https://codecov.io/gh/ylab-hi/ScanNLS "Codecov")
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit "pre-commit")
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black "Black")

</div>

## Features

- TODO

## Requirements

<details>
<summary> Python </summary>
</details>

<details>
<summary> Biopython </summary>
</details>

<details>
<summary> gapmis </summary>
We use the <em>gmapmis</em> package, a tool based on <em>C</em>,  to perform the alignment. you can find the gapmis package <a href="https://github.com/xflouris/gapmis">here</a> Install that:

```bash
$ git clone https://github.com/xflouris/gapmis
$ cd gapmis
$ make
```

You may need to add the path to the gapmis folder to your `$PATH` environment variable.

```bash
$ export PATH=$PATH:$HOME/<gapmis folder>
```

</details>

## Installation

You can install _ScanNLS_ via pip* from PyPI*:

```bash
$ pip install scannls
```

## Usage

Please see the `Command-line Reference <Usage_>`\_ for details.

## Contributing

Contributions are very welcome.
To learn more, see the `Contributor Guide`\_.

## License

Distributed under the terms of the `MIT license`\_,
_ScanNLS_ is free and open source software.

## Issues

If you encounter any problems,
please `file an issue`\_ along with a detailed description.

## Credits

This project was generated from `@cjolowicz`_'s `Hypermodern Python Cookiecutter`_ template.
