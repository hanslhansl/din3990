from dataclasses import dataclass
from enum import IntEnum, auto
from typing import Optional


@dataclass
class Werkstoff:
    class Art(IntEnum):
        # 18 total
        Baustahl = auto()
        Stahlguß = auto()
        Vergütungsstahl = auto()
        Grauguß = auto()
        PerlitischesGußeisenMitKugelgraphit = auto()
        BainitischesGußeisenMitKugelgraphit = auto()
        FerritischesGußeisenMitKugelgraphit = auto()
        SchwarzerTemperguß = auto()
        Einsatzstahl = auto()
        InduktionsgehärteterStahl = auto()
        FlammgehärteterStahl = auto()
        InduktionsgehärtetesGußeisen = auto()
        FlammgehärtetesGußeisen = auto()
        Nitrierstahl = auto()
        NitrierterVergütungsstahl = auto()
        NitrierterEinsatzstahl = auto()
        NitrokarburierterVergütungsstahl = auto()
        NitrokarburierterEinsatzstahl = auto()

    art : Art
    sigma_Hlim : float
    sigma_FE : float
    HB : float
    """Nur das Intervall (130; 470) ist relevant"""
    sigma_02 : Optional[float] = None
    """Notwendig für Baustahl, Vergütungsstahl, PerlitischesGußeisenMitKugelgraphit, BainitischesGußeisenMitKugelgraphit"""

    def __post_init__(self):
        if self.art in (self.Art.Baustahl, self.Art.Vergütungsstahl, self.Art.PerlitischesGußeisenMitKugelgraphit, self.Art.BainitischesGußeisenMitKugelgraphit):
            assert self.sigma_02 is not None
