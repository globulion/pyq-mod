from .HandlerList import HandlerList # utility

from . import XYZ
from . import CML


format_handlers= HandlerList(XYZ.Handler, CML.Handler)
