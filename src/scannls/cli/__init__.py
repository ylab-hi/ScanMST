"""Initialize the scannls.draft module.

@Filename:    __init__.py
@contact:     yangyang.li@northwestern.edu
@Time:        1/1/22 8:28 PM
"""

from .arg import DefaultOptions
from .nls_inference import infer_nls_from_connected_reads

__all__ = ["infer_nls_from_connected_reads", "DefaultOptions"]
