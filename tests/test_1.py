from din3990 import *
from din3962 import din3962_2
import diniso21771

def test_Norm_Beispiel_1():
    werkstoff = din3990_5.Werkstoff(din3990_5.Werkstoff.Art.Einsatzstahl, 1500, 860, 650), din3990_5.Werkstoff(din3990_5.Werkstoff.Art.Vergütungsstahl, 740, 590, 266, 900)

    geometrie = diniso21771.GearGeometry(m_n = 16,
                        z = (23, 113),
                        x = (0.313, -0.071),
                        bezugsprofil = (din3990_11.Protuberanzprofil, din3990_11.Normalprofil1),
                        beta = 7,
                        k = 0,
                        b = 480)

    getriebe = din3990_11.Calculator(geometrie = geometrie, P = 1500,
                n_1 = 275.2,
                verzahnungsqualität = (din3962_2.GearToothQuality.DIN6, din3962_2.GearToothQuality.DIN6),
                werkstoff = werkstoff,
                K_A = 1.25,
                K_S = 1.25,
                R_z = (6, 12),

                anpassungmaßnahmeUndFlankenlinienkorrektur = din3990_11.AnpassungmaßnahmeUndFlankenlinienkorrektur.ohne,
                f_ma = 10,
                s = 0.151 * 1125,
                l = 1125.,
                d_sh = (370., 370.),
                stützwirkung = False,
                kontakttragbild = (din3990_11.Kontakttragbild.a, din3990_11.Kontakttragbild.a),
                ritzelPosition = din3990_11.RitzelPosition.e,
                fertigungsverfahren = (din3990_11.Fertigungsverfahren.geläpptGeschliffenGeschabt, din3990_11.Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt))
    
    assert abs(getriebe.S_Hdyn[Ritzel] - 2.1) < 0.01
    assert abs(getriebe.S_Hdyn[Rad] - 1.2) < 0.01
    assert abs(getriebe.S_Fdyn[Ritzel] - 4.8) < 0.1
    assert abs(getriebe.S_Fdyn[Rad] - 3.3) < 0.01

def test_MEL_2023_24():
    werkstoff = din3990_5.Werkstoff(din3990_5.Werkstoff.Art.Einsatzstahl, 1300, 860, 650)

    geometrie = diniso21771.GearGeometry(m_n = 3,
                        z = (17, 38),
                        x = (0.2, -0.2),
                        bezugsprofil = (din3990_11.Normalprofil2, din3990_11.Normalprofil2),
                        beta = 15,
                        k = 0,
                        b = 40)

    getriebe = din3990_11.Calculator(geometrie = geometrie, P = 7.3,
                n_1 = 1274.11764705882363,
                verzahnungsqualität = (din3962_2.GearToothQuality.DIN6, din3962_2.GearToothQuality.DIN6),
                werkstoff = (werkstoff, werkstoff),
                K_A = 1.1,
                K_S = 1.1,
                R_z = (6, 12),
                
                K_Halpha = (1, 1),
                K_Hbeta = (1.3, 1.3),

                K_Falpha = (1, 1),
                K_Fbeta = (1, 1),
                fertigungsverfahren = (din3990_11.Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt, din3990_11.Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt))
    
    assert abs(getriebe.S_Hdyn[Rad] - 2.058) < 0.01
