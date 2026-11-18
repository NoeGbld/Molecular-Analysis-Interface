# Molecular Analysis Interface

A Multi-database Molecular Analyzer consolidating molecular information from various databases into a single, unified platform.

## Overview

This project provides a Python and HTML-based interface designed to facilitate the analysis and visualization of molecular information collected from multiple databases.

The goal is to provide a centralized and user-friendly platform for exploring molecular data.

## Features

- Collection of molecular information from multiple databases
- Centralized molecular data analysis
- User-friendly web interface
- Python-based data processing
- HTML-based interface

## Installation

> **Important:** the **Quickstart** section of `doc_en` and `doc_fr` is outdated
> (`installer.py` and `requirements.txt` no longer exist).
> Use the procedure below instead. The rest of the documentation is still valid.

### Requirements

- **Python 3.12** : it only works with this version
- [uv](https://docs.astral.sh/uv/): it downloads and uses Python 3.12 automatically,
  so you do not need to install Python yourself

### Steps

1. Install uv:
   ```
   pip install uv
   ```

2. In the project folder, install the dependencies:
   ```
   uv sync
   ```
   On Windows ARM, use this instead,
   because some dependencies have no ARM64 build:
   ```
   uv sync --python cpython-3.12-windows-x86_64
   ```

3. Double-click `Windows Interface Launcher.exe`.

The web interface opens in your browser at `http://127.0.0.1:5000/`.

### First launch

The first launch downloads Chromium (about 330 MB), used for web scraping.
Allow a few minutes: the browser opens automatically once the server is ready.

### Running without the exe (any OS)

```
uv run --no-sync src/launcher.py
```

### Closing the application

Use the quit button of the web interface, or press `Ctrl+C` in the console window.
