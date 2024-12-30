from dataclasses import dataclass
from enum import Enum


@dataclass
class Werkstoff:
    class Art(Enum):
        # 17 total
        Baustahl = 0
        Vergütungsstahl = 1
        Grauguß = 2
        PerlitischesGußeisenMitKugelgraphit = 3
        BainitischesGußeisenMitKugelgraphit = 4
        FerritischesGußeisenMitKugelgraphit = 5
        SchwarzerTemperguß = 6
        Einsatzstahl = 7
        InduktionsgehärteterStahl = 8
        FlammgehärteterStahl = 9
        InduktionsgehärtetesGußeisen = 10
        FlammgehärtetesGußeisen = 11
        Nitrierstahl = 12
        NitrierterVergütungsstahl = 13
        NitrierterEinsatzstahl = 14
        NitrokarburierterVergütungsstahl = 15
        NitrokarburierterEinsatzstahl = 16

    art : Art
    sigma_Hlim : float
    sigma_FE : float
    HB : float
    """Nur das Intervall (130; 470) ist relevant"""