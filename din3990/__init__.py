"""Implementation of DIN 3990 - Calculation of load capacity of cylindrical gears.

See https://github.com/hanslhansl/din3990."""

from .din3990_1 import Ritzel, Rad
from . import din3990_1
from . import din3990_5
from . import din3990_11


"""
Funktionen die Werkstoffart benötigen:
Z_E
Z_W
Z_Xdyn
Z_NT
sigma_HGdyn
Y_deltarelTstat
sigma_FGdyn

"""

# enable colored console output on windows
import colorama
colorama.just_fix_windows_console()