# This little bit of magic fills the __all__ list 
# with every plugin name, and means that calling:
# from plugins import *
# within inmembrane.py will import every plugin
import pkgutil

__all__ = []
for p in pkgutil.iter_modules(__path__):
    __all__.append(p[1])
