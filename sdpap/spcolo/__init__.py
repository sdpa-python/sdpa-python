"""__init__.py
Initialization file of spcolo

December 2010, Kenta KATO
"""

from .spcolo import *
from .asputils import *
from .clique import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
