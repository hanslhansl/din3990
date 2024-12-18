from dataclasses import dataclass
from enum import Enum


@dataclass
class Werkstoff:
    class Art(Enum):
        Stahl = 0
        VergütungsStahl = 1
        Grauguß = 2
        GußeisenMitKugelgraphit = 3
        SchwarzerTemperguß = 4
        Einsatzstahl = 5
        InduktionsgehärteterStahl = 6
        FlammgehärteterStahl = 7
        InduktionsgehärtetesGußeisen = 8
        FlammgehärtetesGußeisen = 9
        Nitrierstahl = 10
        NitrierterVergütungsstahl = 11
        NitrierterEinsatzstahl = 12
        NitrokarburierterVergütungsstahl = 11
        NitrokarburierterEinsatzstahl = 12

        St = Stahl,
        V = VergütungsStahl,
        GG = Grauguß,
        GGG = GußeisenMitKugelgraphit,
        GTS = SchwarzerTemperguß,
        ES = Einsatzstahl,
        IF = InduktionsgehärteterStahl, FlammgehärteterStahl, InduktionsgehärtetesGußeisen, FlammgehärtetesGußeisen
        NT = Nitrierstahl,
        NVnitr = NitrierterVergütungsstahl, NitrierterEinsatzstahl
        NVnitrocar = NitrokarburierterVergütungsstahl, NitrokarburierterEinsatzstahl
        
    art : Art
    sigma_Hlim : float
    sigma_FE : float
    HB : float
    """Nur das Intervall (130; 470) ist relevant"""