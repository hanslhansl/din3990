from dataclasses import dataclass
from re import L
from typing import Optional
from enum import Enum, IntEnum
import math as m
from scipy import optimize
from . import din3990_5


def to_float(val) -> float:
    return val

def involute(alpha):
    return m.tan(m.radians(alpha)) - m.radians(alpha)
def inverse_involute(alpha, anfangswert = 20):
    try:
        return float(optimize.newton(lambda x: involute(x) - alpha, anfangswert))
    except RuntimeError:
        assert(False)
def interpolate(a : float, b : float, t : float):
    return a - t * (a - b)

dauerfest = float("inf")
Ritzel = 0
Rad = 1
_indices = (Ritzel, Rad)

@dataclass
class Profil:
    alpha_n : float
    h_aP_s : float
    h_fP_s : float
    rho_fP_s : float
Normalprofil1 =     Profil(20, 1, 1.25, 0.250)
Normalprofil2 =     Profil(20, 1, 1.25, 0.375)
Protuberanzprofil = Profil(20, 1, 1.40, 0.400)

class Verzahnungsqualität(IntEnum):
    """nach DIN 3962 Teil1 bis Teil 3 bzw. ISO 1328 - 1975"""
    DIN6 = 6
    DIN7 = 7
    DIN8 = 8
    DIN9 = 9
    DIN10 = 10
    DIN11 = 11
    DIN12 = 12
    ISO5 = 105
    ISO6 = 106
    ISO7 = 107
    ISO8 = 108
    ISO9 = 109
    ISO10 = 110
    ISO11 = 111
    ISO12 = 112

class Tabelle3_2(float, Enum):
    ohneBreitenballigkeltOderEndrücknahme = 0.023
    mitSinnvollerBreitenballigkeit = 0.012
    mitSinnvollerEndrücknahme = 0.016

class Bild3_1(IntEnum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4
    f = 5

class Bild3_2(IntEnum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4

class Fertigungsverfahren(IntEnum):
    wälzgefrästWälzgestoßenWälzgehobelt = 0
    geläpptGeschliffenGeschabt = 1

class _WerkstoffKategorie:
    Art = din3990_5.Werkstoff.Art

    St = Art.Baustahl,
    V = Art.Vergütungsstahl,
    GG = Art.Grauguß,
    GGG = Art.PerlitischesGußeisenMitKugelgraphit, Art.BainitischesGußeisenMitKugelgraphit, Art.FerritischesGußeisenMitKugelgraphit
    GTS = Art.SchwarzerTemperguß,
    Eh = Art.Einsatzstahl,
    IF = Art.InduktionsgehärteterStahl, Art.FlammgehärteterStahl, Art.InduktionsgehärtetesGußeisen, Art.FlammgehärtetesGußeisen
    NT = Art.Nitrierstahl,
    NVnitr = Art.NitrierterVergütungsstahl, Art.NitrierterEinsatzstahl
    NVnitrocar = Art.NitrokarburierterVergütungsstahl, Art.NitrokarburierterEinsatzstahl
    NTV = Art.Nitrierstahl, Art.NitrierterVergütungsstahl, Art.NitrokarburierterVergütungsstahl, Art.NitrokarburierterEinsatzstahl

    # für Tabelle 4.1
    Stahl = (Art.Baustahl, Art.Vergütungsstahl, Art.Einsatzstahl,
             Art.InduktionsgehärteterStahl, Art.FlammgehärteterStahl,
             Art.Nitrierstahl, Art.NitrierterVergütungsstahl, Art.NitrierterEinsatzstahl,
             Art.NitrokarburierterVergütungsstahl, Art.NitrokarburierterEinsatzst)
    Stahlguß = ()
    GußeisenMitKugelgraphit = Art.PerlitischesGußeisenMitKugelgraphit, Art.BainitischesGußeisenMitKugelgraphit, Art.FerritischesGußeisenMitKugelgraphit
    GußeisenMitLamellengraphit = Art.Grauguß,

def epsilon_alphan(epsilon_alpha : float, beta_b : float):
    return epsilon_alpha / m.cos(m.radians(beta_b))**2

def v(n : float, d : float):
    return n * m.pi * d / 60000
def F_t(P : float, n : float, d : float):
    """Glg 3.01"""
    return 60000000 * P / n / m.pi / d
def T(P : float, n : float):
    """Glg 3.02"""
    return 30000 * P / n / m.pi

def K_1(verzahnungsqualität : Verzahnungsqualität, geradverzahnt : bool):
    """Tabelle 3.1"""
    if geradverzahnt:
        match verzahnungsqualität:
            case Verzahnungsqualität.DIN6:
                return 9.6
            case Verzahnungsqualität.DIN7:
                return 15.3
            case Verzahnungsqualität.DIN8:
                return 24.5
            case Verzahnungsqualität.DIN9:
                return 34.5
            case Verzahnungsqualität.DIN10:
                return 53.6
            case Verzahnungsqualität.DIN11:
                return 76.6
            case Verzahnungsqualität.DIN12:
                return 122.5
            case Verzahnungsqualität.ISO5:
                return 7.5
            case Verzahnungsqualität.ISO6:
                return 14.9
            case Verzahnungsqualität.ISO7:
                return 26.8
            case Verzahnungsqualität.ISO8:
                return 39.1
            case Verzahnungsqualität.ISO9:
                return 52.8
            case Verzahnungsqualität.ISO10:
                return 76.6
            case Verzahnungsqualität.ISO11:
                return 102.6
            case Verzahnungsqualität.ISO12:
                return 146.3
    else:
        match verzahnungsqualität:
            case Verzahnungsqualität.DIN6:
                return 8.5
            case Verzahnungsqualität.DIN7:
                return 13.6
            case Verzahnungsqualität.DIN8:
                return 21.8
            case Verzahnungsqualität.DIN9:
                return 30.7
            case Verzahnungsqualität.DIN10:
                return 47.7
            case Verzahnungsqualität.DIN11:
                return 68.2
            case Verzahnungsqualität.DIN12:
                return 109.1
            case Verzahnungsqualität.ISO5:
                return 6.7
            case Verzahnungsqualität.ISO6:
                return 13.3
            case Verzahnungsqualität.ISO7:
                return 23.9
            case Verzahnungsqualität.ISO8:
                return 34.8
            case Verzahnungsqualität.ISO9:
                return 47.0
            case Verzahnungsqualität.ISO10:
                return 68.2
            case Verzahnungsqualität.ISO11:
                return 91.4
            case Verzahnungsqualität.ISO12:
                return 130.3
    raise ValueError(f"Falsche Verzahnungsqualität {verzahnungsqualität}")
def K_2(geradverzahnt : bool):
    """Tabelle 3.1"""
    if geradverzahnt:
        return 0.0193
    return 0.0087
def K_V(z_1 : int, v : float, u : float, F_t : float, K_A : float, b : float, epsilon_beta : float, geradverzahnt : bool, verzahnungsqualität : Verzahnungsqualität, _print):
    """Abschnitt 3.3"""
    # Glg 3.04
    temp1 = z_1 * v / 100 * m.sqrt(u**2 / (1 + u**2))
    _print("z_1 * v / 100 * sqrt(u^2 / (1 + u^2)) =", temp1)
    assert temp1 < 10

     # Abschnitt 3.3.1b
    temp2 = max(F_t * K_A / b, 100)
    _print("F_t * K_A / b =", temp2)

    def _K_V(geradverzahnt):
        return 1 + (K_1(verzahnungsqualität, geradverzahnt) / temp2 + K_2(geradverzahnt)) * temp1

    if geradverzahnt:
        return _K_V(True)
    elif epsilon_beta >= 1:
        return _K_V(False)
    else:
        K_Valpha = _K_V(True)
        K_Vbeta = _K_V(False)
        return interpolate(K_Valpha, K_Vbeta, epsilon_beta)

def F_m(F_t : float, K_A : float, K_V : float):
    """Glg 3.07"""
    return F_t * K_A * K_V
def K_s(stützwirkung : bool, bild3_2 : Bild3_2):
    """Bild 3.2"""
    match bild3_2:
        case Bild3_2.a:
            return 0.48 if stützwirkung else 0.8
        case Bild3_2.b:
            return -0.48 if stützwirkung else -0.8
        case Bild3_2.c:
            return 1.33
        case Bild3_2.d:
            return -0.36 if stützwirkung else -0.6
        case Bild3_2.e:
            return -0.6 if stützwirkung else -1.0
def f_sh(F_m : float, d : float, b : float, doppelschrägverzahnt : bool, A : Tabelle3_2, s : float,
         stützwirkung : Optional[bool] = None,
         bild3_2 : Optional[Bild3_2] = None,
         l : Optional[float] = None,
         d_sh : Optional[float] = None) -> float:
    """Glg 3.14, 3.15"""
    temp = 0. if s == 0 else K_s(stützwirkung, bild3_2) * l * s / d**2 * (d / d_sh)**4

    if not doppelschrägverzahnt:
        print(F_m)
        print(b)
        print(A)
        return F_m / b * A * (abs(1 + temp - 0.3) + 0.3) * (b / d)**2
    else:
        b_B = b / 2
        return F_m / b * 2 * A * (abs(1.5 + temp - 0.3) + 0.3) * (b_B / d)**2
def F_betax(d : float, f_sh : float, f_ma : float, doppelschrägverzahnt : bool, s : float,
            bild3_1 : Optional[Bild3_1],
            bild3_2 : Optional[Bild3_2] = None,
            stützwirkung : Optional[bool] = None,
            l : Optional[float] = None,
            d_sh : Optional[float] = None) -> float:
    """Glg 3.09"""
    
    if f_ma == 0:
        return abs(1.33 * f_sh)
    else:
        match bild3_1:
            case Bild3_1.a | Bild3_1.f:
                multi = -1
            case Bild3_1.b | Bild3_1.e:
                multi = 1
            case Bild3_1.c:
                B_s = 1.5 if doppelschrägverzahnt else 1
                multi = 1 if abs(K_s(stützwirkung, bild3_2)) * l * s / d**2 * (d / d_sh)**4 <= B_s else -1
            case Bild3_1.d:
                B_s = 1.5 if doppelschrägverzahnt else 1
                multi = 1 if abs(K_s(stützwirkung, bild3_2)) * l * s / d**2 * (d / d_sh)**4 >= B_s - 0.3 else -1
        return abs(1.33 * f_sh + multi * f_ma)
def _y_beta(werkstoff : din3990_5.Werkstoff, v : float, F_betax : float):
    """Abschnitt 3.4.2.6"""
    match werkstoff.art:
        case din3990_5.Werkstoff.Art.Baustahl | din3990_5.Werkstoff.Art.Vergütungsstahl | din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit | din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit:
            # Glg 3.16
            y_beta = 320 / werkstoff.sigma_Hlim * F_betax
            if v <= 5:
                pass
            elif v <= 10:
                assert y_beta <= 25600 / werkstoff.sigma_Hlim
            else:
                assert y_beta <= 12800 / werkstoff.sigma_Hlim
            return y_beta
        case din3990_5.Werkstoff.Art.Grauguß | din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit:
            # Glg 3.17
            y_beta = 0.55 * F_betax
            if v <= 5:
                pass
            elif v <= 10:
                assert y_beta <= 45
            else:
                assert y_beta <= 22
            return y_beta
        case (din3990_5.Werkstoff.Art.Einsatzstahl |
              din3990_5.Werkstoff.Art.InduktionsgehärteterStahl | din3990_5.Werkstoff.Art.FlammgehärteterStahl |
              din3990_5.Werkstoff.Art.Nitrierstahl | din3990_5.Werkstoff.Art.NitrierterEinsatzstahl |
              din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl |
              din3990_5.Werkstoff.Art.InduktionsgehärtetesGußeisen | din3990_5.Werkstoff.Art.FlammgehärtetesGußeisen
              ):
            y_beta = 0.15 * F_betax
            assert y_beta <= 6
            return y_beta
    raise NotImplementedError
def y_beta(werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff], v : float, F_betax : float):
    return (_y_beta(werkstoff[Ritzel], v, F_betax) + _y_beta(werkstoff[Rad], v, F_betax)) / 2
def F_betay(F_betax : float, y_beta : float):
    return F_betax - y_beta
def K_Hbeta(F_t : float, K_A : float, K_V : float, v : float, d : float, b : float, doppelschrägverzahnt : bool, werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff], A : Tabelle3_2, s : float,
            bild3_1 : Optional[Bild3_2] = None,
            bild3_2 : Optional[Bild3_2] = None,
            stützwirkung : Optional[bool] = None,
            l : Optional[float] = None,
            d_sh : Optional[float] = None,
            _print = print):
    """Abschnitt 3.4"""

    # Abschnitt 3.4.1
    assert(F_t / b * K_A >= 100)

    _F_m = F_m(F_t, K_A, K_V)
    _print("F_m =", _F_m)

    _f_sh = f_sh(_F_m, d, b, doppelschrägverzahnt, A, s, stützwirkung, bild3_2, l, d_sh)
    _print("f_sh =", _f_sh)

    _F_betax = F_betax(d, _f_sh, 0, doppelschrägverzahnt, s, bild3_1, stützwirkung, bild3_2, l, d_sh)
    _print("F_betax =", _F_betax)

    _y_beta = y_beta(werkstoff, v, _F_betax)
    _print("y_beta =", _y_beta)

    _F_betay = F_betay(_F_betax, _y_beta)
    _print("F_betay =", _F_betay)

    # Abschnitt 3.4.3.1
    c_gamma = 20
    _print("c_γ =", c_gamma)

    # Glg 3.20, 3.21
    _K_Hbeta = 1 + c_gamma * _F_betay / (2 * _F_m / b)
    if _K_Hbeta <= 2:
        pass
    else:
        _K_Hbeta = m.sqrt(2 * c_gamma * _F_betay / (_F_m / b))
        assert _K_Hbeta > 2
    return _K_Hbeta

def K_Fbeta(F_t : float, K_A : float, K_Hbeta : float, b : float, h : float):
    """Glg 3.22"""
    assert(F_t / b * K_A >= 100)   # Abschnitt 3.4.1
    h_b = min(h / b, 1. / 3.)
    return m.pow(K_Hbeta, (1 / (1 + h_b + h_b**2)))

def K_H_Falpha(K_A : float, F_t : float, b : float, beta_b : float, epsilon_alpha : float, geradverzahnt : bool, werkstoff : din3990_5.Werkstoff,
               verzahnungsqualität : tuple[Verzahnungsqualität, Verzahnungsqualität], Z_epsilon : float, Y_epsilon : float):
    """
    K_Hα und K_Fα
    Tabelle 3.3
    """
    linienbelastung = F_t / b * K_A
    qualität = Verzahnungsqualität(max(vq + 99 if Verzahnungsqualität.DIN6 <= vq <= Verzahnungsqualität.DIN12 else vq for vq in verzahnungsqualität))

    if werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl,
                         din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                         din3990_5.Werkstoff.Art.InduktionsgehärtetesGußeisen, din3990_5.Werkstoff.Art.FlammgehärtetesGußeisen,
                         din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl,
                         din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):
        if geradverzahnt:
            if linienbelastung > 100:
                match qualität:
                    case Verzahnungsqualität.ISO5 | Verzahnungsqualität.ISO6:
                        return 1., 1.
                    case Verzahnungsqualität.ISO7:
                        return 1.1, 1.1
                    case Verzahnungsqualität.ISO8:
                        return 1.2, 1.2
                    case Verzahnungsqualität.ISO9 | Verzahnungsqualität.ISO10 | Verzahnungsqualität.ISO11 | Verzahnungsqualität.ISO12:
                        K_H = 1 / Z_epsilon**2
                        K_F = 1 / Y_epsilon**2
                        assert K_H >= 1.2
                        assert K_F >= 1.2
                        return K_H, K_F
            else:
                if qualität >= Verzahnungsqualität.ISO5:
                    K_H = 1 / Z_epsilon**2
                    K_F = 1 / Y_epsilon**2
                    assert K_H >= 1.2
                    assert K_F >= 1.2
                    return K_H, K_F
        else:
            if linienbelastung > 100:
                match qualität:
                    case Verzahnungsqualität.ISO5:
                        return 1., 1.
                    case Verzahnungsqualität.ISO6:
                        return 1.1, 1.1
                    case Verzahnungsqualität.ISO7:
                        return 1.2, 1.2
                    case Verzahnungsqualität.ISO8:
                        return 1.4, 1.4
                    case Verzahnungsqualität.ISO9 | Verzahnungsqualität.ISO10 | Verzahnungsqualität.ISO11 | Verzahnungsqualität.ISO12:
                        K = epsilon_alpha / m.cos(m.radians(beta_b))**2
                        assert K >= 1.4
                        return K, K
            else:
                if qualität >= Verzahnungsqualität.ISO5:
                    K = epsilon_alpha / m.cos(m.radians(beta_b))**2
                    assert K >= 1.4
                    return K, K
    else:
        if geradverzahnt:
            if linienbelastung > 100:
                match qualität:
                    case Verzahnungsqualität.ISO5 | Verzahnungsqualität.ISO6 | Verzahnungsqualität.ISO7:
                        return 1., 1.
                    case Verzahnungsqualität.ISO8:
                        return 1.1, 1.1
                    case Verzahnungsqualität.ISO9:
                        return 1.2, 1.2
                    case Verzahnungsqualität.ISO10 | Verzahnungsqualität.ISO11 | Verzahnungsqualität.ISO12:
                        K_H = 1 / Z_epsilon**2
                        K_F = 1 / Y_epsilon**2
                        assert K_H >= 1.2
                        assert K_F >= 1.2
                        return K_H, K_F
            else:
                if qualität >= Verzahnungsqualität.ISO5:
                    K_H = 1 / Z_epsilon**2
                    K_F = 1 / Y_epsilon**2
                    assert K_H >= 1.2
                    assert K_F >= 1.2
                    return K_H, K_F
        else:
            if linienbelastung > 100:
                match qualität:
                    case Verzahnungsqualität.ISO5 | Verzahnungsqualität.ISO6:
                        return 1., 1.
                    case Verzahnungsqualität.ISO7:
                        return 1.1, 1.1
                    case Verzahnungsqualität.ISO8:
                        return 1.2, 1.2
                    case Verzahnungsqualität.ISO9:
                        return 1.4, 1.4
                    case Verzahnungsqualität.ISO10 | Verzahnungsqualität.ISO11 | Verzahnungsqualität.ISO12:
                        K = epsilon_alpha / m.cos(m.radians(beta_b))**2
                        assert K >= 1.4
                        return K, K
            else:
                if qualität >= Verzahnungsqualität.ISO5:
                    K = epsilon_alpha / m.cos(m.radians(beta_b))**2
                    assert K >= 1.4
                    return K, K
    raise ValueError(f"Unerwartete Verzahnungsqualität {qualität}")

def M_1(z : tuple[int, int], d_a : tuple[float, float], d_b : tuple[float, float], alpha_wt : float, epsilon_alpha : float):
    """Glg 4.12"""
    return m.tan(m.radians(alpha_wt)) / m.sqrt(
            (m.sqrt(d_a[Ritzel]**2 / d_b[Ritzel]**2 - 1) - 2 * m.pi / z[Ritzel]) *
            (m.sqrt(d_a[Rad]**2 / d_b[Rad]**2 - 1) - (epsilon_alpha - 1) * 2 * m.pi / z[Rad]))
def Z_B(z : tuple[int, int], d_a : tuple[float, float], d_b : tuple[float, float], alpha_wt : float, epsilon_alpha : float, epsilon_beta : float, geradverzahnt : bool, _print = print):
    """Abschnitt 4.2"""
    
    # Glg 4.12
    _M_1 = M_1(z, d_a, d_b, alpha_wt, epsilon_alpha)
    _print("M_1 =", _M_1)

    if geradverzahnt:
        return max(1., _M_1)
    elif epsilon_beta >= 1:
        return 1.
    else:
        return max(1., _M_1 - epsilon_beta * (_M_1 - 1))

def M_2(z : tuple[int, int], d_a : tuple[float, float], d_b : tuple[float, float], alpha_wt : float, epsilon_alpha : float):
    """Glg 4.13"""
    return m.tan(m.radians(alpha_wt)) / m.sqrt(
            (m.sqrt(d_a[Rad]**2 / d_b[Rad]**2 - 1) - 2 * m.pi / z[Rad]) *
            (m.sqrt(d_a[Ritzel]**2 / d_b[Ritzel]**2 - 1) - (epsilon_alpha - 1) * 2 * m.pi / z[Ritzel]))
def Z_D(z : tuple[int, int], d_a : tuple[float, float], d_b : tuple[float, float], alpha_wt : float, epsilon_alpha : float, epsilon_beta : float, geradverzahnt : bool, innenverzahnt : bool, _print = print):
    """Abschnitt 4.2"""

    if innenverzahnt:
        return 1.

    # Glg 4.12
    _M_2 = M_2(z, d_a, d_b, alpha_wt, epsilon_alpha)
    _print("M_2 =", _M_2)

    if geradverzahnt == 0:
        return max(1., _M_2)
    elif epsilon_beta >= 1:
        return 1.
    else:
        return max(1., _M_2 - epsilon_beta * (_M_2 - 1))

def Z_H(alpha_t : float, alpha_wt : float, beta_b : float):
    """Glg 4.14"""
    return m.sqrt(2 * m.cos(m.radians(beta_b)) * m.cos(m.radians(alpha_wt))
                          / m.cos(m.radians(alpha_t))**2 / m.sin(m.radians(alpha_wt)))

def Z_E(werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff]):
    """Tabelle 4.1"""
    ws1, ws2 = werkstoff
    if ws1 in _WerkstoffKategorie.Stahl:
        if ws2 in _WerkstoffKategorie.Stahl:
            return 189.8
        elif ws2 in _WerkstoffKategorie.Stahlguß:
            return 188.9
        elif ws2 in _WerkstoffKategorie.GußeisenMitKugelgraphit:
            return 181.4
        elif ws2 in _WerkstoffKategorie.GußeisenMitLamellengraphit:
            return 165.4
    elif ws1 in _WerkstoffKategorie.Stahlguß:
        if ws2 in _WerkstoffKategorie.Stahl:
            return 188.9
        elif ws2 in _WerkstoffKategorie.Stahlguß:
            return 188.0
        elif ws2 in _WerkstoffKategorie.GußeisenMitKugelgraphit:
            return 180.5
        elif ws2 in _WerkstoffKategorie.GußeisenMitLamellengraphit:
            return 161.4
    elif ws1 in _WerkstoffKategorie.GußeisenMitKugelgraphit:
        if ws2 in _WerkstoffKategorie.Stahl:
            return 181.4
        elif ws2 in _WerkstoffKategorie.Stahlguß:
            return 180.5
        elif ws2 in _WerkstoffKategorie.GußeisenMitKugelgraphit:
            return 173.9
        elif ws2 in _WerkstoffKategorie.GußeisenMitLamellengraphit:
            return 156.6
    elif ws1 in _WerkstoffKategorie.GußeisenMitLamellengraphit:
        if ws2 in _WerkstoffKategorie.Stahl:
            return 165.4
        elif ws2 in _WerkstoffKategorie.Stahlguß:
            return 161.4
        elif ws2 in _WerkstoffKategorie.GußeisenMitKugelgraphit:
            return 156.6
        elif ws2 in _WerkstoffKategorie.GußeisenMitLamellengraphit:
            return 146.0
    raise NotImplementedError

def Z_epsilon(beta : float, epsilon_alpha : float, epsilon_beta : float):
    """Abschnitt 4.5"""
    if beta == 0:
        return m.sqrt((4. - epsilon_alpha) / 3.)
    elif epsilon_beta >= 1.:
        return m.sqrt(1 / epsilon_alpha)
    else:
        return m.sqrt((4. - epsilon_alpha) / 3. * (1. - epsilon_beta) + epsilon_beta / epsilon_alpha)

def Z_beta(beta : float):
    """ Glg 4.18"""
    return m.sqrt(m.cos(m.radians(beta)))

def Z_LVRstat():
    """Glg 4.22"""
    return 1.
def Z_LVRdyn(fertigungsverfahren : tuple[Fertigungsverfahren, Fertigungsverfahren], R_z : tuple[float, float], a : float, _print = print):
    """Abschnitt 4.8a und c"""
    if fertigungsverfahren[0] == Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt and fertigungsverfahren[1] == Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt:
        return 0.85

    # Glg 4.20
    R_z100 = sum(R_z) / 2 * m.pow(100 / a, 1/3)
    _print("R_z100 =", R_z100)

    if fertigungsverfahren[0] == Fertigungsverfahren.geläpptGeschliffenGeschabt and fertigungsverfahren[1] == Fertigungsverfahren.geläpptGeschliffenGeschabt:
        if R_z100 > 4:
            # Glg 4.21
            return 0.92
        else:
            # Glg 4.22
            return 1.
    else:
        assert R_z100 <= 4
        # Glg 4.22
        return 0.92
_Z_LVRdyn = Z_LVRdyn

def _Z_W(werkstoff : din3990_5.Werkstoff, anderer_werkstoff : din3990_5.Werkstoff, andere_R_z : float):
    if andere_R_z <= 6:
        if (werkstoff.art in _WerkstoffKategorie.St or
            werkstoff.art in _WerkstoffKategorie.V or
            werkstoff.art in _WerkstoffKategorie.GGG):
            if (anderer_werkstoff.art in _WerkstoffKategorie.Eh or
                anderer_werkstoff.art in _WerkstoffKategorie.IF or
                anderer_werkstoff.art in _WerkstoffKategorie.NT or
                anderer_werkstoff.art in _WerkstoffKategorie.NVnitr or
                anderer_werkstoff.art in _WerkstoffKategorie.NVnitrocar):
                HB = werkstoff.HB
                if HB <= 130:
                    return 1.2
                elif HB >= 470:
                    return 1.
                else:
                    return 1.2 - (HB - 130) / 1700
    return 1.
def Z_W(werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff], R_z : tuple[float, float]):
    """Abschnitt 4.9"""
    return _Z_W(werkstoff[0], werkstoff[1], R_z[1]), _Z_W(werkstoff[1], werkstoff[0], R_z[0])

def Z_Xstat():
    return 1.
def Z_Xdyn(m_n : float, werkstoff : din3990_5.Werkstoff):
    """Tabelle 4.2"""
    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl, din3990_5.Werkstoff.Art.Grauguß,
                         din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit,
                         din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.SchwarzerTemperguß):  # 7
        return 1.
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl, din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                           din3990_5.Werkstoff.Art.InduktionsgehärtetesGußeisen, din3990_5.Werkstoff.Art.FlammgehärtetesGußeisen):  # 5
        if m_n <= 10:
            return 1.
        elif m_n < 30:
            return 1.05 - 0.005  * m_n
        else:
            return 0.9
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrierterVergütungsstahl,
                           din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):    # 4
        if m_n <= 7.5:
            return 1.
        elif m_n < 30:
            return 1.08 - 0.011  * m_n
        else:
            return 0.75
    raise ValueError

def Z_NT(werkstoff : din3990_5.Werkstoff):
    """
    Z_NTstat, Z_NTdyn
    Tabelle 4.3
    """
    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl,
                         din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit,
                         din3990_5.Werkstoff.Art.SchwarzerTemperguß, din3990_5.Werkstoff.Art.Einsatzstahl,
                         din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                         din3990_5.Werkstoff.Art.InduktionsgehärtetesGußeisen, din3990_5.Werkstoff.Art.FlammgehärtetesGußeisen):
        return 1.6, 1.0
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Grauguß, din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.Nitrierstahl,
                           din3990_5.Werkstoff.Art.NitrierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl):
        return 1.3, 1.0
    elif werkstoff.art in (din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):
        return 1.1, 1.0
    raise ValueError

def sigma_HGstat(werkstoff : din3990_5.Werkstoff, Z_NTstat : float, Z_LVRstat : float, Z_W : float, Z_Xstat : float):
    """Glg 4.03"""
    return werkstoff.sigma_Hlim * Z_NTstat * Z_LVRstat * Z_W * Z_Xstat
def sigma_HGdyn(sigma_HGstat : float, werkstoff : din3990_5.Werkstoff, Z_NTdyn : float, Z_LVRdyn : float, Z_W : float, Z_Xdyn : float, N_L : float, gewisseGrübchenbildung : bool):
    """Glg 4.03 und Abschnitt 4.1.2b"""
    sigma_HGdauer = werkstoff.sigma_Hlim * Z_NTdyn * Z_LVRdyn * Z_W * Z_Xdyn
    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl,
                         din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit,
                         din3990_5.Werkstoff.Art.SchwarzerTemperguß,
                         din3990_5.Werkstoff.Art.Einsatzstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl,
                         din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                         din3990_5.Werkstoff.Art.InduktionsgehärtetesGußeisen, din3990_5.Werkstoff.Art.FlammgehärtetesGußeisen):    # 12
        if gewisseGrübchenbildung:
            if N_L <= 6 * 10**5:
                return sigma_HGstat
            elif 6 * 10**5 < N_L <= 10**7:
                # Glg 4.05
                exp = 0.3705 * m.log10(sigma_HGstat / sigma_HGdauer)
                # Glg 4.04
                return sigma_HGdauer * m.pow(3 * 10**8 / N_L, exp)
            elif 10**7 < N_L <= 10**9:
                # Glg 4.07
                exp = 0.2791 * m.log10(sigma_HGstat / sigma_HGdauer)
                # Glg 4.06
                return sigma_HGdauer * m.pow(10**9 / N_L, exp)
            else:
                return sigma_HGdauer
        else:
            if N_L <= 10**5:
                return sigma_HGstat
            elif 10**5 < N_L <= 5 * 10**7:
                # Glg 4.05
                exp = 0.3705 * m.log10(sigma_HGstat / sigma_HGdauer)
                # Glg 4.08
                return sigma_HGdauer * m.pow(5 * 10**7 / N_L, exp)
            else:
                return sigma_HGdauer
    elif werkstoff.art in (din3990_5.Werkstoff.Art.NitrierterVergütungsstahl, din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl,
                           din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.Grauguß):   # 5
        if N_L <= 10**5:
            return sigma_HGstat
        elif 10**5 < N_L <= 2 * 10**6:
            # Glg 4.10
            exp = 0.7686 * m.log10(sigma_HGstat / sigma_HGdauer)
            # Glg 4.09
            return sigma_HGdauer * m.pow(2 * 10**6 / N_L, exp)
        else:
            return sigma_HGdauer
    raise NotImplementedError

def Y_epsilon(epsilon_alpha : float, beta_b : float):
    """Abschnitt 5.3"""
    return 0.25 + 0.75 / epsilon_alpha * m.cos(m.radians(beta_b))**2

class DIN_21771:
    def __init__(self,
            m_n : float,
            z: tuple[int, int],
            x: tuple[float, float],
            bezugsprofil : Profil,
            beta : float,
            k : int,
            b : Optional[float] = None,
            b_d_1_verhältnis : Optional[float] = None,
            _print = print):
        """
        Parameters:
        - m_n: Modul
        - z: Zähnezahlen. z1 ist immer positiv. Für Außenradpaare ist z2 positiv, für Innenradpaare ist z2 negativ.
        - x: Profilverschiebungsfaktoren
        """

        self.m_n = m_n
        self.z = z
        self.x = x
        self.beta = beta
        self.k = k

        _print("Getriebegeometrie")
        [_print(key, "=", value) for key, value in vars(self).items()]

        self.alpha_n = bezugsprofil.alpha_n
        self.h_aP = bezugsprofil.h_aP_s * self.m_n
        self.h_fP = bezugsprofil.h_fP_s * self.m_n
        self.rho_fP = bezugsprofil.rho_fP_s * self.m_n
        _print("α_n =", self.alpha_n)
        _print("h_aP =", self.h_aP)
        _print("h_fP =", self.h_fP)
        _print("ρ_fP =", self.rho_fP)

        self.alpha_t = m.degrees(m.atan(m.tan(m.radians(self.alpha_n)) / m.cos(m.radians(self.beta))))
        _print("α_t =", self.alpha_t)

        self.alpha_wt = inverse_involute( involute(self.alpha_t) + 2 * sum(self.x) / sum(self.z) * m.tan(m.radians(self.alpha_n)) )
        _print("α_wt =", self.alpha_wt)

        self.u = self.z[Rad] / self.z[Ritzel]
        _print("u =", self.u)

        def d(idx):
            return self.z[idx] * self.m_n / m.cos(m.radians(self.beta))
        self.d = d(Ritzel), d(Rad)
        _print("d =", self.d)

        def d_b(idx):
            return self.d[idx] * m.cos(m.radians(self.alpha_t))
        self.d_b = d_b(Ritzel), d_b(Rad)
        _print("d_b =", self.d_b)

        def d_a(idx):
            return self.d[idx] + 2 * (self.x[idx] * self.m_n + self.h_aP + self.k * self.m_n)
        self.d_a = d_a(Ritzel), d_a(Rad)
        _print("d_a =", self.d_a)

        def d_f(idx):
            return self.d[idx] - 2 * (self.h_fP - self.x[idx] * self.m_n)
        self.d_f = d_f(Ritzel), d_f(Rad)
        _print("d_f =", self.d_f)

        def d_w(idx):
            return self.d_b[idx] / m.cos(m.radians(self.alpha_wt))
        self.d_w = d_w(Ritzel), d_w(Rad)
        _print("d_w =", self.d_w)

        if b != None:
            self.b = b
        elif b_d_1_verhältnis != None:
            self.b = self.d[Ritzel] * b_d_1_verhältnis
        else:
            raise ValueError("either b or b_d_1_verhältnis must be specified as argument")
        _print("b =", self.b)

        self.m_t = m_n / m.cos(m.radians(self.beta))
        _print("m_t =", self.m_t)

        self.beta_b = m.degrees(m.atan(m.tan(m.radians(self.beta)) * m.cos(m.radians(self.alpha_t))))
        _print("β_b =", self.beta_b)

        self.h = self.h_aP + k * m_n + self.h_fP
        _print("h =", self.h)

        self.a_w = sum(self.d_w) / 2
        _print("a_w =", self.a_w)

        # Profilüberdeckung

        self.p_n = self.m_n * m.pi
        _print("p_n =", self.p_n)

        self.p_t = self.p_n / m.cos(m.radians(self.beta))
        _print("p_t =", self.p_t)

        self.epsilon_alpha = (m.sqrt(self.d_a[Ritzel]**2 - self.d_b[Ritzel]**2) + m.sqrt(self.d_a[Rad]**2 - self.d_b[Rad]**2) - sum(self.d_b) * m.tan(m.radians(self.alpha_wt))) / (2 * self.p_t * m.cos(m.radians(self.alpha_t)))
        _print("ε_α =", self.epsilon_alpha)

        self.epsilon_beta = self.b * m.sin(m.radians(self.beta)) / self.p_n
        _print("ε_β =", self.epsilon_beta)

        self.epsilon_gamma = self.epsilon_alpha + self.epsilon_beta
        _print("ε_γ =", self.epsilon_gamma)
        assert self.epsilon_gamma > 1

        # Unterschnitt

        self.h_aP0 = self.h_fP
        _print("h_aP0 =", self.h_aP0)

        def _z_min(idx):
            return 2 * m.cos(m.radians(self.beta)) * (self.h_aP0 / self.m_n - self.x[idx]) / m.sin(m.radians(self.alpha_t))**2
        z_min = _z_min(Ritzel), _z_min(Rad)
        _print("z_min =", z_min)
        assert self.z[Ritzel] > z_min[Ritzel]
        assert self.z[Rad] > z_min[Rad]

        # Spitzwerden

        def _gamma(idx):
            return inverse_involute( m.pi / 2 / self.z[idx] + 2 * self.x[idx] / self.z[idx] * m.tan(m.radians(self.alpha_n)) + involute(self.alpha_t) )
        gamma = _gamma(Ritzel), _gamma(Rad)
        _print("γ =", gamma)

        def _d_amax(idx):
            return self.m_t * self.z[idx] * m.cos(m.radians(self.alpha_t)) / m.cos(m.radians(gamma[idx]))
        d_amax = _d_amax(Ritzel), _d_amax(Rad)
        _print("d_amax =", d_amax)
        assert self.d_a[Ritzel] <= d_amax[Ritzel]
        assert self.d_a[Rad] <= d_amax[Rad]

        _print()
        return
  
_K_V = K_V
_K_Hbeta = K_Hbeta
_K_Fbeta = K_Fbeta
class DIN_3990_11:
    def __init__(self,
                geometrie : DIN_21771,
                P : float,
                n_1 : float,
                verzahnungsqualität : tuple[Verzahnungsqualität, Verzahnungsqualität],
                werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff],
                K_A : float,
                K_S : float,
                R_z: tuple[float, float],
                N_L : float = dauerfest,
                doppelschrägverzahnt: bool = False,
                innenverzahnt: bool = False,
                gewisseGrübchenbildung : bool = False,
                s_pr: float = 0,
                
                K_V : tuple[Optional[float], Optional[float]] = (None, None),
                K_Hbeta : tuple[Optional[float], Optional[float]] = (None, None),
                A : Optional[Tabelle3_2] = None,
                f_ma : tuple[Optional[float], Optional[float]] = (None, None),
                bild3_1 : tuple[Optional[Bild3_1], Optional[Bild3_1]] =  (None, None),
                s : tuple[Optional[float], Optional[float]] = (None, None),
                stützwirkung : tuple[Optional[bool], Optional[bool]] =  (None, None),
                bild3_2 : tuple[Optional[Bild3_2], Optional[Bild3_2]] =  (None, None),
                l : tuple[Optional[float], Optional[float]] = (None, None),
                d_sh : tuple[Optional[float], Optional[float]] = (None, None),
                Z_LVRdyn : Optional[float] = None,
                fertigungsverfahren : Optional[tuple[Fertigungsverfahren, Fertigungsverfahren]] = None,
                K_Fbeta : tuple[Optional[float], Optional[float]] = (None, None),
                K_Halpha : tuple[Optional[float], Optional[float]] = (None, None),
                K_Falpha : tuple[Optional[float], Optional[float]] = (None, None),
                
                S_Hstatmin : float | tuple[float, float] = 1.,
                S_Hdynmin : float | tuple[float, float] = 1.,
                S_Fstatmin : float | tuple[float, float] = 1.,
                S_Fdynmin : float | tuple[float, float] = 1.,

                _print = print,
                _assert: bool = False):
        """
        Parameters:
        - geometrie
        - P: Leistung [kW]
        - n_1: Antriebsdrehzahl [1/min]
        - verzahnungsqualität
        - werkstoff
        - K_A: siehe Tabelle A.1
        - K_S: ersetz K_A für die statische Berechnung
        - R_z: gemittelte Rauhtiefe
        - N_L: Lastspielzahl
        - doppelschrägverzahnt
        - innenverzahnt
        - gewisseGrübchenbildung: ob gewisse Grübchenbildung zulässig ist
        - s_pr: float, Fussfreischnitt
        - _print: the function used for printing
        - _assert if True, safety is asserted which additionally requires the same arguments as passed to Getriebe.is_safe()

        Optional parameters:
        - K_V
        - K_Hbeta
          or 
            - A: siehe Tabelle 3.2
            - f_ma: siehe Abschnitt 3.4.2.4
              if f_ma != 0:
                - bild3_1: siehe Bild 3.1
            - s: siehe Bild 3.2
              if s != 0:
                - stützwirkung: siehe Bild 3.2
                - bild3_2: siehe Bild 3.2
                - l: siehe Bild 3.2
                - d_sh: Wellendurchmesser
        - Z_LVRdyn: float, siehe Abschnitt 4.8
          or
            - fertigungsverfahren: siehe Abschnitt 4.8
        - K_Fbeta
        - K_Halpha
        - K_Falpha
        """
        
        _print("Parameter")
        [_print(key, "=", value) for key, value in locals().items() if key not in ("self", "_print")]
        _print()

        self.geometrie = geometrie
        self.P = P
        self.n = (n_1, n_1 / self.geometrie.u)
        self.verzahnungsqualität = verzahnungsqualität
        self.werkstoff = werkstoff
        self.K_A = K_A
        self.K_S = K_S
        self.R_z = R_z
        self.N_L = N_L
        self.doppelschrägverzahnt = doppelschrägverzahnt
        self.innenverzahnt = innenverzahnt
        self.gewisseGrübchenbildung = gewisseGrübchenbildung
        self.s_pr = s_pr

        self.K_V = K_V
        self.K_Hbeta = K_Hbeta
        self.A = A
        self.f_ma = f_ma
        self.bild3_1 = bild3_1
        self.s = s
        self.stützwirkung = stützwirkung
        self.bild3_2 = bild3_2
        self.l = l
        self.d_sh = d_sh
        self.Z_LVRdyn = Z_LVRdyn
        self.fertigungsverfahren = fertigungsverfahren
        self.K_Fbeta = K_Fbeta
        self.K_Halpha = K_Halpha
        self.K_Falpha = K_Falpha

        def check_value(idx : int, val : tuple[float, float], min_or_interv : float | tuple[float, float], name : str, full_name : str):
            val = val[idx]
            if isinstance(min_or_interv, tuple):
                _print(name, idx + 1, " = ", val, " in ", min_or_interv, sep="")
                res = min_or_interv[0] <= val <= min_or_interv[1]
            else:
                _print(name, idx + 1, " = ", val, " >= ", min_or_interv, sep="")
                res = min_or_interv <= val
            if res == False:
                _print("\033[31m", full_name, " von Rad ", idx + 1, " ist nicht erfüllt\033[0m", sep="")
            if _assert:
                assert res


        geradverzahnt = self.geometrie.beta == 0

        assert all(n <= 3600 for n in self.n), "Drehzahl darf nicht größer als 3600 sein (siehe Abschnitt 1.2)"
        #assert False, "(siehe Abschnitt 1.3c)"
        assert self.geometrie.beta <= 30, "β darf nicht größer als 30° sein (siehe Abschnitt 1.3d)"
        assert not (self.doppelschrägverzahnt and geradverzahnt), "Ein doppelschrägverzahntes Getriebe kann nicht geradverzahnt sein"

        _print("DIN 3990-11")
        _print("n =", self.n)

        # if self.doppelschrägverzahnt:
        #     self.b_B = self.geometrie.b / 2
        #     _print("b_B =", self.b_B)

        self.v = v(self.n[Ritzel], self.geometrie.d[Ritzel])
        _print("v =", self.v)

        self.F_t = F_t(self.P, self.n[Ritzel], self.geometrie.d[Ritzel])
        _print("F_t =", self.F_t)
  
        self.T = tuple(T(self.P, self.n[idx]) for idx in _indices)
        _print("T =", self.T)

        self.K_V = tuple(_K_V(self.geometrie.z[Ritzel], self.v, self.geometrie.u, self.F_t, self.K_A, self.geometrie.b, self.geometrie.epsilon_beta, geradverzahnt, self.verzahnungsqualität[idx],
                            _print) if K_V is None else K_V for idx, K_V in zip(_indices, self.K_V))
        _print("K_V =", self.K_V)

        self.K_Hbeta = tuple(_K_Hbeta(self.F_t, self.K_A, self.K_V[idx], self.v, self.geometrie.d[idx], self.geometrie.b, self.doppelschrägverzahnt, self.werkstoff, self.A, self.s[idx],
                                    self.bild3_1[idx], self.bild3_2, self.stützwirkung, self.l, self.d_sh, _print) if K_Hbeta is None else K_Hbeta for idx, K_Hbeta in zip(_indices, self.K_Hbeta))
        _print("K_Hβ =", self.K_Hbeta)
        
        self.K_Fbeta = tuple(_K_Fbeta(self.F_t, self.K_A, self.K_Hbeta[idx], self.geometrie.b, self.geometrie.h) if K_Fbeta is None else K_Fbeta for idx, K_Fbeta in zip(_indices, self.K_Fbeta))
        _print("K_Fβ =", self.K_Fbeta)
        
        self.Z_epsilon = Z_epsilon(self.geometrie.beta, self.geometrie.epsilon_alpha, self.geometrie.epsilon_beta)
        _print("Z_ε =", self.Z_epsilon)

        self.Y_epsilon = Y_epsilon(self.geometrie.epsilon_alpha, self.geometrie.beta_b)
        _print("Y_ε =", self.Y_epsilon)

        K_H_Falpha_ritzel, K_H_Falpha_rad = (K_H_Falpha(self.K_A, self.F_t, self.geometrie.b, self.geometrie.beta_b, self.geometrie.epsilon_alpha, geradverzahnt, self.werkstoff[idx],
                                       self.verzahnungsqualität, self.Z_epsilon, self.Y_epsilon) for idx in _indices)
        self.K_Halpha = K_H_Falpha_ritzel[0] if self.K_Halpha[0] is None else self.K_Halpha[0], K_H_Falpha_rad[0] if self.K_Halpha[1] is None else self.K_Halpha[1]
        _print("K_Hα =", self.K_Halpha)
        self.K_Falpha = K_H_Falpha_ritzel[1] if self.K_Falpha[0] is None else self.K_Falpha[0], K_H_Falpha_rad[1] if self.K_Falpha[1] is None else self.K_Falpha[1]
        _print("K_Fα =", self.K_Falpha)

        # Grübchentragfähigkeit
        
        self.Z_B = Z_B(self.geometrie.z, self.geometrie.d_a, self.geometrie.d_b, self.geometrie.alpha_wt, self.geometrie.epsilon_alpha, self.geometrie.epsilon_beta, geradverzahnt, _print)
        _print("Z_B =", self.Z_B)

        self.Z_D = Z_D(self.geometrie.z, self.geometrie.d_a, self.geometrie.d_b, self.geometrie.alpha_wt, self.geometrie.epsilon_alpha, self.geometrie.epsilon_beta, geradverzahnt,
                       self.innenverzahnt, _print)
        _print("Z_D =", self.Z_D)

        self.Z_H = Z_H(self.geometrie.alpha_t, self.geometrie.alpha_wt, self.geometrie.beta_b)
        _print("Z_H =", self.Z_H)

        self.Z_E = Z_E(self.werkstoff)
        _print("Z_E =", self.Z_E)

        self.Z_beta = Z_beta(self.geometrie.beta)
        _print("Z_β =", self.Z_beta)

        # Glg 4.19
        self.Z_LVRstat = Z_LVRstat()
        if self.Z_LVRdyn == None:
            self.Z_LVRdyn = _Z_LVRdyn(self.fertigungsverfahren, self.R_z, self.geometrie.a_w, _print)
        _print("Z_LVRstat =", self.Z_LVRstat)
        _print("Z_LVRdyn =", self.Z_LVRdyn)

        self.Z_W = Z_W(self.werkstoff, self.R_z)
        _print("Z_W =", self.Z_W)

        self.Z_Xstat = 1., 1.
        self.Z_Xdyn = Z_Xdyn(self.geometrie.m_n, self.werkstoff[Ritzel]), Z_Xdyn(self.geometrie.m_n, self.werkstoff[Rad])
        _print("Z_Xstat =", self.Z_Xstat)
        _print("Z_Xdyn =", self.Z_Xdyn)

        Z_NTritzel = Z_NT(self.werkstoff[Ritzel])
        Z_NTrad = Z_NT(self.werkstoff[Rad])
        self.Z_NTstat = Z_NTritzel[0], Z_NTrad[0]
        self.Z_NTdyn = Z_NTritzel[1], Z_NTrad[1]
        _print("Z_NTstat =", self.Z_NTstat)
        _print("Z_NTdyn =", self.Z_NTdyn)

        self.sigma_HGstat = tuple(sigma_HGstat(self.werkstoff[idx], self.Z_NTstat[idx], self.Z_LVRstat, self.Z_W, self.Z_Xstat[idx]) for idx in _indices)
        self.sigma_HGdyn = tuple(sigma_HGdyn(sigma_HGstat[idx], self.werkstoff[idx], self.Z_NTdyn[idx], self.Z_LVRdyn, self.Z_W[idx], self.Z_Xdyn[idx], self.N_L, self.gewisseGrübchenbildung)
                                 for idx in _indices)
        _print("σ_HGstat =", self.sigma_HGstat)
        _print("σ_HGdyn =", self.sigma_HGdyn)

        # Glg 4.02
        self.sigma_H0 = self.Z_H * self.Z_E * self.Z_epsilon * self.Z_beta * m.sqrt(self.F_t / self.geometrie.d[Ritzel] / self.geometrie.b * (self.geometrie.u + 1) / self.geometrie.u)
        _print("σ_H0 =", self.sigma_H0)

        # Glg 4.01
        def sigma_Hstat(idx):
            return (self.Z_B if idx == Ritzel else self.Z_D) * self.sigma_H0 * m.sqrt(self.K_S * self.K_V[idx] * self.K_Hbeta[idx] * self.K_Halpha[idx])
        def sigma_Hdyn(idx):
            return (self.Z_B if idx == Ritzel else self.Z_D) * self.sigma_H0 * m.sqrt(self.K_A * self.K_V[idx] * self.K_Hbeta[idx] * self.K_Halpha[idx])
        self.sigma_Hstat = sigma_Hstat(Ritzel), sigma_Hstat(Rad)
        self.sigma_Hdyn = sigma_Hdyn(Ritzel), sigma_Hdyn(Rad)
        _print("σ_Hstat =", self.sigma_Hstat)
        _print("σ_Hdyn =", self.sigma_Hdyn)

        # Glg 4.11
        def S_Hstat(idx):
            return self.sigma_HGstat[idx] / self.sigma_Hstat[idx]
        def S_Hdyn(idx):
            return self.sigma_HGdyn[idx] / self.sigma_Hdyn[idx]
        self.S_Hstat = S_Hstat(Ritzel), S_Hstat(Rad)
        self.S_Hdyn = S_Hdyn(Ritzel), S_Hdyn(Rad)

        check_value(Ritzel, self.S_Hstat, S_Hstatmin, "S_Hstat", "statische Grübchensicherheit")
        check_value(Ritzel, self.S_Hdyn, S_Hdynmin, "S_Hdyn", "dynamische Grübchensicherheit")
        check_value(Rad, self.S_Hstat, S_Hstatmin, "S_Hstat", "statische Grübchensicherheit")
        check_value(Rad, self.S_Hdyn, S_Hdynmin, "S_Hdyn", "dynamische Grübchensicherheit")
        _print()

        # Zahnfußtragfähigkeit

        # Glg D.1.01
        def z_n(idx):
            return self.geometrie.z[idx] / m.cos(m.radians(self.geometrie.beta_b))**2 / m.cos(m.radians(self.geometrie.beta))
        self.z_n = z_n(Ritzel), z_n(Rad)
        _print("z_n =", self.z_n)

        if not self.innenverzahnt:
            # Glg D.5.01
            self.E = (m.pi / 4 * self.geometrie.m_n - self.geometrie.h_fP * m.tan(m.radians(self.geometrie.alpha_n)) + self.s_pr / m.cos(m.radians(self.geometrie.alpha_n))
                      - (1 - m.sin(m.radians(self.geometrie.alpha_n))) * self.geometrie.rho_fP / m.cos(m.radians(self.geometrie.alpha_n))) 
            _print("E =", self.E)

            # Glg D.5.02
            def G(idx):
                return self.geometrie.rho_fP / self.geometrie.m_n - self.geometrie.h_fP / self.geometrie.m_n + self.geometrie.x[idx]
            self.G = G(Ritzel), G(Rad)
            _print("G =", self.G)

            # Glg D.5.03
            def H(idx):
                return 2 / self.z_n[idx] * (m.pi / 2 - self.E / self.geometrie.m_n) - m.pi / 3
            self.H = H(Ritzel), H(Rad)
            _print("H =", self.H)

            # Glg D.5.04
            def theta(idx):
                theta = m.degrees(m.pi / 6)
                for i in range(5):
                    theta = m.degrees(2 * self.G[idx] / self.z_n[idx] * m.tan(m.radians(theta)) - self.H[idx])
                return theta
            self.theta = theta(Ritzel), theta(Rad)
            _print("ϑ =", self.theta)

            # Glg D.5.05
            def s_Fn(idx):
                return self.geometrie.m_n * (self.z_n[idx] * m.sin(m.pi / 3 - m.radians(self.theta[idx])) + m.sqrt(3) * (self.G[idx] / m.cos(m.radians(self.theta[idx]))
                                                                                                                         - self.geometrie.rho_fP / self.geometrie.m_n))
            self.s_Fn = s_Fn(Ritzel), s_Fn(Rad)
            _print("s_Fn =", self.s_Fn)

            # Glg D.5.06
            def d_n(idx):
                return self.geometrie.m_n * self.z_n[idx]
            self.d_n = d_n(Ritzel), d_n(Rad)
            _print("d_n =", self.d_n)

            # Glg D.5.07
            def d_bn(idx):
                return self.d_n[idx] * m.cos(m.radians(self.geometrie.alpha_n))
            self.d_bn = d_bn(Ritzel), d_bn(Rad)
            _print("d_bn =", self.d_bn)

            # Glg D.5.08
            def d_an(idx):
                return self.d_n[idx] + self.geometrie.d_a[idx] - self.geometrie.d[idx]
            self.d_an = d_an(Ritzel), d_an(Rad)
            _print("d_an =", self.d_an)

            # Glg D.5.09
            def alpha_an(idx):
                return m.degrees(m.acos(self.d_bn[idx] / self.d_an[idx]))
            self.alpha_an = alpha_an(Ritzel), alpha_an(Rad)
            _print("α_an =", self.alpha_an)

            # Glg D.5.10
            def y_a(idx):
                return 1 / self.z_n[idx] * (m.pi / 2 + 2 * self.geometrie.x[idx] * m.tan(m.radians(self.geometrie.alpha_n))) + involute(self.geometrie.alpha_n) - involute(self.alpha_an[idx])
            self.y_a = y_a(Ritzel), y_a(Rad)
            _print("y_a =", self.y_a)

            # Glg D.5.11
            def alpha_Fan(idx):
                return self.alpha_an[idx] - m.degrees(self.y_a[idx])
            self.alpha_Fan = alpha_Fan(Ritzel), alpha_Fan(Rad)
            _print("α_Fan =", self.alpha_Fan)

            # Glg D.5.12
            def h_Fa(idx):
                return self.geometrie.m_n * (0.5 * self.z_n[idx] * (m.cos(m.radians(self.geometrie.alpha_n)) / m.cos(m.radians(self.alpha_Fan[idx])) - m.cos(m.pi / 3 - m.radians(self.theta[idx])))
                                   + 0.5 * (self.geometrie.rho_fP / self.geometrie.m_n - self.G[idx] / m.cos(m.radians(self.theta[idx]))))
            self.h_Fa = h_Fa(Ritzel), h_Fa(Rad)
            _print("h_Fa =", self.h_Fa)

            # Glg D.5.13
            def rho_F(idx):
                return self.geometrie.rho_fP + self.geometrie.m_n * 2 * self.G[idx]**2 / m.cos(m.radians(self.theta[idx])) / (self.z_n[idx] * m.cos(m.radians(self.theta[idx]))**2 - 2 * self.G[idx])
            self.rho_F = rho_F(Ritzel), rho_F(Rad)
            _print("ρ_F =", self.rho_F)
        else:
            raise NotImplementedError("Innenverzahnte Räder sind noch nicht implementiert")

        # Glg D.3.01
        def Y_Fa(idx):
            return (6 * self.h_Fa[idx] / self.geometrie.m_n * m.cos(m.radians(self.alpha_Fan[idx]))) / ((self.s_Fn[idx] / self.geometrie.m_n)**2 * m.cos(m.radians(self.geometrie.alpha_n)))
        self.Y_Fa = Y_Fa(Ritzel), Y_Fa(Rad)
        _print("Y_Fa =", self.Y_Fa)

        # Glg D.4.02
        def L_a(idx):
            return self.s_Fn[idx] / self.h_Fa[idx]
        self.L_a = L_a(Ritzel), L_a(Rad)
        _print("L_a =", self.L_a)

        # Glg D.4.03
        def q_s(idx):
            return self.s_Fn[idx] / 2 / self.rho_F[idx]
        self.q_s = q_s(Ritzel), q_s(Rad)
        _print("q_s =", self.q_s)
        assert all(1 <= q_s < 8 for q_s in self.q_s)

        # Glg D.4.01
        def Y_Sa(idx):
            return (1.2 + 0.13 * self.L_a[idx]) * m.pow(self.q_s[idx], 1 / (1.21 + 2.3 / self.L_a[idx]))
        self.Y_Sa = Y_Sa(Ritzel), Y_Sa(Rad)
        _print("Y_Sa =", self.Y_Sa)

        # Glg 5.08
        def Y_FS(idx):
            return self.Y_Sa[idx] * self.Y_Fa[idx]
        self.Y_FS = Y_FS(Ritzel), Y_FS(Rad)
        _print("Y_FS =", self.Y_FS)

        # Glg 5.10
        self.Y_beta = 1 - min(self.geometrie.epsilon_beta, 1.) * min(self.geometrie.beta, 30) / 120
        _print("Y_β =", self.Y_beta)
 
        # Tabelle 5.8
        epsilon_alphan = epsilon_alphan(self.geometrie.epsilon_alpha, self.geometrie.beta_b)
        _print("ε_αn =", epsilon_alphan)
        def Y_S(idx):
            assert 1 <= self.s_Fn[idx] / self.h_Fa[idx] <= 1.2
            return self.Y_Sa[idx] * (0.6 + 0.4 * epsilon_alphan)
        self.Y_S = Y_S(Ritzel), Y_S(Rad)
        _print("Y_S =", self.Y_S)

        # Abschnitt 5.6
        def Y_deltarelTstat(idx : int):
            wsart = self.werkstoff[idx].art
            if wsart in (din3990_5.Werkstoff.Art.Baustahl, ):
                raise NotImplementedError
            elif wsart in (din3990_5.Werkstoff.Art.vergüteterStahl, ):
                raise NotImplementedError
            elif wsart in (din3990_5.Werkstoff.Art.einsatzgehärteterStahl, ):
                return 0.44 * self.Y_S[idx] + 0.12
            elif wsart in (din3990_5.Werkstoff.Art.nitrierterStahl, din3990_5.Werkstoff.Art.nitrokarburierterStahl ):
                return 0.20 * self.Y_S[idx] + 0.60
            raise NotImplementedError
        # Glg 5.11, 5.12
        def Y_deltarelTdyn(idx):
            if self.q_s[idx] >= 1.5:
                return 1.
            else:
                return 0.95
        self.Y_deltarelTstat = Y_deltarelTstat(Ritzel), Y_deltarelTstat(Rad)
        self.Y_deltarelTdyn = Y_deltarelTdyn(Ritzel), Y_deltarelTdyn(Rad)
        _print("Y_δrelTstat =", self.Y_deltarelTstat)
        _print("Y_δrelTdyn =", self.Y_deltarelTdyn)

        # Abschnitt 5.7
        self.Y_RrelTstat = 1.
        def Y_RrelTdyn(idx):
            if self.R_z[idx] <= 16:
                return 1.
            else:
                return 0.9
        self.Y_RrelTdyn = Y_RrelTdyn(Ritzel), Y_RrelTdyn(Rad)
        _print("Y_RrelTstat =", self.Y_RrelTstat)
        _print("Y_RrelTdyn =", self.Y_RrelTdyn)

        # Tabelle 5.1
        def Y_Xdyn(idx : int):
            wsart = self.werkstoff[idx].art
            if wsart in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.vergüteterStahl):
                if self.geometrie.m_n <= 5:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.03 - 0.006  * self.geometrie.m_n
                else:
                    return 0.85
            elif wsart in (din3990_5.Werkstoff.Art.einsatzgehärteterStahl, din3990_5.Werkstoff.Art.nitrierterStahl, din3990_5.Werkstoff.Art.nitrokarburierterStahl):
                if self.geometrie.m_n <= 5:
                    return 1.
                elif self.geometrie.m_n < 25:
                    return 1.05 - 0.01  * self.geometrie.m_n
                else:
                    return 0.8
            raise NotImplementedError
        self.Y_Xstat = 1., 1.
        self.Y_Xdyn = Y_Xdyn(Ritzel), Y_Xdyn(Rad)
        _print("Y_Xstat =", self.Y_Xstat)
        _print("Y_Xdyn =", self.Y_Xdyn)

        # Tabelle 5.2
        def Y_NT(idx : int):
            wsart = self.werkstoff[idx].art
            if wsart in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.vergüteterStahl):
                return 2.5, 1.0
            elif wsart in (din3990_5.Werkstoff.Art.einsatzgehärteterStahl, ):
                return 2.5, 1.0
            elif wsart in (din3990_5.Werkstoff.Art.nitrierterStahl, ):
                return 1.6, 1.0
            elif wsart in (din3990_5.Werkstoff.Art.nitrokarburierterStahl, ):
                return 1.1, 1.0
            raise NotImplementedError
        Y_NTritzel = Y_NT(Ritzel)
        Y_NTrad = Y_NT(Rad)
        self.Y_NTstat = Y_NTritzel[0], Y_NTrad[0]
        self.Y_NTdyn = Y_NTritzel[1], Y_NTrad[1]
        _print("Y_NTstat =", self.Y_NTstat)
        _print("Y_NTdyn =", self.Y_NTdyn)

        # Glg 5.03
        def sigma_FGstat(idx):
            return self.werkstoff[idx].sigma_FE * self.Y_NTstat[idx] * self.Y_deltarelTstat[idx] * self.Y_RrelTstat * self.Y_Xstat[idx]
        def sigma_FGdyn(idx):
            return self.werkstoff[idx].sigma_FE * self.Y_NTdyn[idx] * self.Y_deltarelTdyn[idx] * self.Y_RrelTdyn[idx] * self.Y_Xdyn[idx]
        self.sigma_FGstat = sigma_FGstat(Ritzel), sigma_FGstat(Rad)
        self.sigma_FGdyn = sigma_FGdyn(Ritzel), sigma_FGdyn(Rad)
        _print("σ_FGstat =", self.sigma_FGstat)
        _print("σ_FGdyn =", self.sigma_FGdyn)

        # Glg 5.02
        def sigma_F0(idx):
            return self.F_t / self.geometrie.b / self.geometrie.m_n * self.Y_FS[idx] * self.Y_epsilon * self.Y_beta
        self.sigma_F0 = sigma_F0(Ritzel), sigma_F0(Rad)
        _print("σ_F0 =", self.sigma_F0)

        # Glg 5.01
        def sigma_Fstat(idx):
            return self.sigma_F0[idx] * self.K_S * self.K_V[idx] * self.K_Fbeta[idx] * self.K_Falpha[idx]
        def sigma_Fdyn(idx):
            return self.sigma_F0[idx] * self.K_A * self.K_V[idx] * self.K_Fbeta[idx] * self.K_Falpha[idx]
        self.sigma_Fstat = sigma_Fstat(Ritzel), sigma_Fstat(Rad)
        self.sigma_Fdyn = sigma_Fdyn(Ritzel), sigma_Fdyn(Rad)
        _print("σ_Fstat =", self.sigma_Fstat)
        _print("σ_Fdyn =", self.sigma_Fdyn)

        # Glg 5.07
        def S_Fstat(idx):
            return self.sigma_FGstat[idx] / self.sigma_Fstat[idx]
        def S_Fdyn(idx):
            return self.sigma_FGdyn[idx] / self.sigma_Fdyn[idx]
        self.S_Fstat = S_Fstat(Ritzel), S_Fstat(Rad)
        self.S_Fdyn = S_Fdyn(Ritzel), S_Fdyn(Rad)

        check_value(Ritzel, self.S_Fstat, S_Fstatmin, "S_Fstat", "statische Zahnbruchsicherheit")
        check_value(Ritzel, self.S_Fdyn, S_Fdynmin, "S_Fdyn", "dynamische Zahnbruchsicherheit")
        check_value(Rad, self.S_Fstat, S_Fstatmin, "S_Fstat", "statische Zahnbruchsicherheit")
        check_value(Rad, self.S_Fdyn, S_Fdynmin, "S_Fdyn", "dynamische Zahnbruchsicherheit")
        _print()
        return


