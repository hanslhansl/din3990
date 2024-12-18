from dataclasses import dataclass
from enum import Enum


@dataclass
class Werkstoff:
    class Art(Enum):
        Stahl = 0
        Verg�tungsStahl = 1
        Graugu� = 2
        Gu�eisenMitKugelgraphit = 3
        SchwarzerTempergu� = 4
        Einsatzstahl = 5
        Induktionsgeh�rteterStahl = 6
        Flammgeh�rteterStahl = 7
        Induktionsgeh�rtetesGu�eisen = 8
        Flammgeh�rtetesGu�eisen = 9
        Nitrierstahl = 10
        NitrierterVerg�tungsstahl = 11
        NitrierterEinsatzstahl = 12
        NitrokarburierterVerg�tungsstahl = 11
        NitrokarburierterEinsatzstahl = 12

        St = Stahl,
        V = Verg�tungsStahl,
        GG = Graugu�,
        GGG = Gu�eisenMitKugelgraphit,
        GTS = SchwarzerTempergu�,
        ES = Einsatzstahl,
        IF = Induktionsgeh�rteterStahl, Flammgeh�rteterStahl, Induktionsgeh�rtetesGu�eisen, Flammgeh�rtetesGu�eisen
        NT = Nitrierstahl,
        NVnitr = NitrierterVerg�tungsstahl, NitrierterEinsatzstahl
        NVnitrocar = NitrokarburierterVerg�tungsstahl, NitrokarburierterEinsatzstahl
        
    art : Art
    sigma_Hlim : float
    sigma_FE : float
    HB : float
    """Nur das Intervall (130; 470) ist relevant"""