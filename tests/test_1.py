from din3990 import *

def test_Norm_Beispiel_1():
    werkstoff = din3990_5.Werkstoff(din3990_5.Werkstoff.Art.Einsatzstahl, 1500, 860, 650), din3990_5.Werkstoff(din3990_5.Werkstoff.Art.Vergütungsstahl, 740, 590, 266, 900)

    geometrie = din3990_11.DIN_21771(m_n = 16,
                        z = (23, 113),
                        x = (0.313, -0.071),
                        bezugsprofil = (din3990_11.Protuberanzprofil, din3990_11.Normalprofil1),
                        beta = 7,
                        k = 0,
                        b = 480)

    getriebe = din3990_11.Calculator(geometrie = geometrie, P = 1500,
                n_1 = 275.2,
                verzahnungsqualität = (din3990_11.Verzahnungsqualität.DIN6, din3990_11.Verzahnungsqualität.DIN6),
                werkstoff = werkstoff,
                K_A = 1.25,
                K_S = 1.25,
                R_z = (6, 12),

                A = din3990_11.Tabelle3_2.ohneBreitenballigkeltOderEndrücknahme,
                f_ma = (10, 10),  # Annahme, siehe Fußnote 5
                s = (0.151 * 1125, 0.151 * 1125),
                l = (1125., 1125.),
                d_sh = (370., 370.),
                stützwirkung = (False, False),
                bild3_1 = (din3990_11.Bild3_1.a, din3990_11.Bild3_1.a),
                bild3_2 = (din3990_11.Bild3_2.e, din3990_11.Bild3_2.e),
                fertigungsverfahren = (din3990_11.Fertigungsverfahren.geläpptGeschliffenGeschabt, din3990_11.Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt))
    
    assert abs(getriebe.S_Hdyn[Ritzel] - 2.1) < 0.01
    assert abs(getriebe.S_Hdyn[Rad] - 1.2) < 0.01
    assert abs(getriebe.S_Fdyn[Ritzel] - 4.8) < 0.1
    assert abs(getriebe.S_Fdyn[Rad] - 3.3) < 0.01

