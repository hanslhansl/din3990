from din3990 import *
from din3962 import din3962_2
import diniso21771


S_Hstatmin = 1.3
S_Hdyn_interval = (1.2, 1.5)
S_Fstatmin = 3.5
S_Fdyn_interval = (1.5, 2)

werkstoff = din3990_5.Werkstoff(din3990_5.Werkstoff.Art.Einsatzstahl, 1500, 860, 220)

geometrie = diniso21771.GearGeometry(m_n = 4,
                    z = (19, 104),
                    x = (0.5, 0.15),
                    bezugsprofil = (din3990_11.Normalprofil1, din3990_11.Normalprofil1),
                    beta = 0,
                    k = 0,
                    b_d_1_verhältnis = 0.64)


getriebe = din3990_11.Calculator(geometrie = geometrie, P = 55,
            n_1 = 980,
            verzahnungsqualität = (din3962_2.GearToothQuality.DIN6, din3962_2.GearToothQuality.DIN7),
            werkstoff = (werkstoff, werkstoff),
            K_A = 1.75,
            K_S = 2.5,
            R_z = (5, 5),

            anpassungmaßnahmeUndFlankenlinienkorrektur = din3990_11.AnpassungmaßnahmeUndFlankenlinienkorrektur.ohne,
            f_ma = 0,  # Annahme, siehe Fußnote 5
            s = 0,
            fertigungsverfahren = (din3990_11.Fertigungsverfahren.geläpptGeschliffenGeschabt, din3990_11.Fertigungsverfahren.geläpptGeschliffenGeschabt),

            S_Hstatmin=S_Hstatmin, S_Hdynmin=S_Hdyn_interval, S_Fstatmin=S_Fstatmin, S_Fdynmin=S_Fdyn_interval,

            _assert = True)


assert geometrie.b / geometrie.m_n <= 30 # Konstruktionsvorgaben Tabelle 4
assert getriebe.R_z100 < 4 # Konstruktionsvorgaben Seite 7
