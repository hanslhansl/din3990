from dataclasses import dataclass
from re import L
from typing import Optional
from enum import Enum, IntEnum, auto
import math as m
from scipy import optimize

from . import din3990_5
from .din3990_1 import Ritzel, Rad

_indices = (Ritzel, Rad)

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


@dataclass
class Profil:
    alpha_n : float
    h_aP_s : float
    h_fP_s : float
    rho_fP_s : float
    s_pr_s : float
Normalprofil1 =     Profil(20, 1, 1.25, 0.250, 0.)
Normalprofil2 =     Profil(20, 1, 1.25, 0.375, 0.)
Protuberanzprofil = Profil(20, 1, 1.4,  0.4,   0.02)

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
    a = auto()
    b = auto()
    c = auto()
    d = auto()
    e = auto()
    f = auto()

class Bild3_2(IntEnum):
    a = auto()
    b = auto()
    c = auto()
    d = auto()
    e = auto()

class Fertigungsverfahren(IntEnum):
    wälzgefrästWälzgestoßenWälzgehobelt = auto()
    geläpptGeschliffenGeschabt = auto()

class _WerkstoffKategorie:
    # für Tabelle 4.1
    Stahl = (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl, din3990_5.Werkstoff.Art.Einsatzstahl,
             din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
             din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl,
             din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl)
    Stahlguß = din3990_5.Werkstoff.Art.Stahlguß,
    GußeisenMitKugelgraphit = (din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit,
                               din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit)
    GußeisenMitLamellengraphit = din3990_5.Werkstoff.Art.Grauguß,

def epsilon_alphan(epsilon_alpha : float, beta_b : float):
    """Bild 5.8"""
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
    _print("\tz_1 * v / 100 * sqrt(u^2 / (1 + u^2)) =", temp1)
    assert temp1 < 10

     # Abschnitt 3.3.1b
    temp2 = max(F_t * K_A / b, 100)
    _print("\tF_t * K_A / b =", temp2)

    def _K_V(geradverzahnt):
        _K_1 = K_1(verzahnungsqualität, geradverzahnt)
        _print("\tK_1 =", _K_1)
        _K_2 = K_2(geradverzahnt)
        _print("\tK_2 =", _K_2)
        return 1 + (_K_1 / temp2 + _K_2) * temp1

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
    raise ValueError(f"{bild3_2}")
def f_sh(F_m : float, d_1 : float, b : float, doppelschrägverzahnt : bool, A : Tabelle3_2, s : float,
         stützwirkung : Optional[bool] = None,
         bild3_2 : Optional[Bild3_2] = None,
         l : Optional[float] = None,
         d_sh : Optional[float] = None) -> float:
    """Glg 3.14, 3.15"""
    if s == 0:
        temp = 0.
    else:
        assert isinstance(stützwirkung, bool)
        assert isinstance(bild3_2, Bild3_2)
        assert isinstance(l, float)
        assert isinstance(d_sh, float)
        temp = K_s(stützwirkung, bild3_2) * l * s / d_1**2 * (d_1 / d_sh)**4

    if not doppelschrägverzahnt:
        return F_m / b * A * (abs(1 + temp - 0.3) + 0.3) * (b / d_1)**2
    else:
        b_B = b / 2
        return F_m / b * 2 * A * (abs(1.5 + temp - 0.3) + 0.3) * (b_B / d_1)**2
def F_betax(d_1 : float, f_sh : float, f_ma : float, doppelschrägverzahnt : bool, s : float,
            bild3_1 : Optional[Bild3_1],
            bild3_2 : Optional[Bild3_2] = None,
            stützwirkung : Optional[bool] = None,
            l : Optional[float] = None,
            d_sh1 : Optional[float] = None) -> float:
    """Glg 3.09"""
    if f_ma == 0:
        return abs(1.33 * f_sh)
    else:
        assert isinstance(bild3_1, Bild3_1)
        B_s = 1.5 if doppelschrägverzahnt else 1.
        match bild3_1:
            case Bild3_1.a | Bild3_1.f:
                multi = -1
            case Bild3_1.b | Bild3_1.e:
                multi = 1
            case Bild3_1.c | Bild3_1.d:
                assert isinstance(stützwirkung, bool)
                assert isinstance(bild3_2, Bild3_2)
                assert isinstance(l, float)
                assert isinstance(d_sh1, float)
                if bild3_1 == Bild3_1.c:
                    multi = 1 if abs(K_s(stützwirkung, bild3_2)) * l * s / d_1**2 * (d_1 / d_sh1)**4 <= B_s else -1
                else:
                    multi = 1 if abs(K_s(stützwirkung, bild3_2)) * l * s / d_1**2 * (d_1 / d_sh1)**4 >= B_s - 0.3 else -1
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
def K_Hbeta(F_t : float, K_A : float, K_V : float, v : float, d_1 : float, b : float, doppelschrägverzahnt : bool, werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff],
            A : Tabelle3_2, f_ma : float, s : float,
            bild3_1 : Optional[Bild3_2] = None,
            bild3_2 : Optional[Bild3_2] = None,
            stützwirkung : Optional[bool] = None,
            l : Optional[float] = None,
            d_sh : Optional[float] = None,
            d_sh1 : Optional[float] = None,
            _print = print):
    """Abschnitt 3.4"""

    # Abschnitt 3.4.1
    assert(F_t / b * K_A >= 100)

    _F_m = F_m(F_t, K_A, K_V)
    _print("\tF_m =", _F_m)
    _print("\tF_m / b =", _F_m / b)

    _f_sh = f_sh(_F_m, d_1, b, doppelschrägverzahnt, A, s, stützwirkung, bild3_2, l, d_sh)
    _print("\tf_sh =", _f_sh)

    _F_betax = F_betax(d_1, _f_sh, f_ma, doppelschrägverzahnt, s, bild3_1, stützwirkung, bild3_2, l, d_sh1)
    _print("\tF_betax =", _F_betax)

    _y_beta = y_beta(werkstoff, v, _F_betax)
    _print("\ty_beta =", _y_beta)

    _F_betay = F_betay(_F_betax, _y_beta)
    _print("\tF_betay =", _F_betay)

    # Abschnitt 3.4.3.1
    c_gamma = 20
    _print("\tc_γ =", c_gamma)

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
    _print("\tM_1 =", _M_1)

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
    _print("\tM_2 =", _M_2)

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
    ws1 = werkstoff[0].art
    ws2 = werkstoff[1].art
    if ws1 in _WerkstoffKategorie.Stahl:
        if ws2 in _WerkstoffKategorie.Stahl:
            return 189.8
        elif ws2 in _WerkstoffKategorie.Stahlguß:
            return 188.9
        elif ws2 in _WerkstoffKategorie.GußeisenMitKugelgraphit:
            return 181.4
        elif ws2 in _WerkstoffKategorie.GußeisenMitLamellengraphit:
            return 165.4
        else:
            ws2
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

def R_z100(R_z : tuple[float, float], a : float):
    """Glg 4.20"""
    return sum(R_z) / 2 * m.pow(100 / a, 1/3)

def Z_LVRstat():
    """Glg 4.22"""
    return 1.
def Z_LVRdyn(fertigungsverfahren : tuple[Fertigungsverfahren, Fertigungsverfahren], R_z100 : float):
    """Abschnitt 4.8a und c"""
    if fertigungsverfahren[0] == Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt and fertigungsverfahren[1] == Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt:
        return 0.85

    if fertigungsverfahren[0] == Fertigungsverfahren.geläpptGeschliffenGeschabt and fertigungsverfahren[1] == Fertigungsverfahren.geläpptGeschliffenGeschabt:
        if R_z100 > 4:
            # Glg 4.21
            return 0.92
        else:
            # Glg 4.22
            return 1.
    else:
        # Glg 4.22
        return 0.92

def _Z_W(werkstoff : din3990_5.Werkstoff, anderer_werkstoff : din3990_5.Werkstoff, andere_R_z : float):
    if andere_R_z <= 6:
        if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl, din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit,
                             din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit):
            if anderer_werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl, 
                                         din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                                         din3990_5.Werkstoff.Art.InduktionsgehärtetesGußeisen, din3990_5.Werkstoff.Art.FlammgehärtetesGußeisen,
                                         din3990_5.Werkstoff.Art.Nitrierstahl,
                                         din3990_5.Werkstoff.Art.NitrierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl,
                                         din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):
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
            elif N_L <= 10**7:
                # Glg 4.05
                exp = 0.3705 * m.log10(sigma_HGstat / sigma_HGdauer)
                # Glg 4.04
                return sigma_HGdauer * m.pow(3 * 10**8 / N_L, exp)
            elif N_L <= 10**9:
                # Glg 4.07
                exp = 0.2791 * m.log10(sigma_HGstat / sigma_HGdauer)
                # Glg 4.06
                return sigma_HGdauer * m.pow(10**9 / N_L, exp)
            else:
                return sigma_HGdauer
        else:
            if N_L <= 10**5:
                return sigma_HGstat
            elif N_L <= 5 * 10**7:
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
        elif N_L <= 2 * 10**6:
            # Glg 4.10
            exp = 0.7686 * m.log10(sigma_HGstat / sigma_HGdauer)
            # Glg 4.09
            return sigma_HGdauer * m.pow(2 * 10**6 / N_L, exp)
        else:
            return sigma_HGdauer
    raise NotImplementedError

def sigma_H0(F_t : float, d_1 : float, b : float, u : float, Z_H : float, Z_E : float, Z_epsilon : float, Z_beta : float):
    """Glg 4.02"""
    return Z_H * Z_E * Z_epsilon * Z_beta * m.sqrt(F_t / d_1 / b * (u + 1) / u)

def sigma_Hstat(K_S : float, K_V : float, K_Halpha : float, K_Hbeta : float, sigma_H0 : float, Z_B_D : float):
    """
    Z_B gilt für das Ritzel, Z_D für das Rad.
    Glg 4.01
    """
    return Z_B_D * sigma_H0 * m.sqrt(K_S * K_V * K_Hbeta * K_Halpha)
def sigma_Hdyn(K_A : float, K_V : float, K_Halpha : float, K_Hbeta : float, sigma_H0 : float, Z_B_D : float):
    """
    Z_B gilt für das Ritzel, Z_D für das Rad.
    Glg 4.01
    """
    return Z_B_D * sigma_H0 * m.sqrt(K_A * K_V * K_Hbeta * K_Halpha)

def S_Hstat(sigma_HGstat : float, sigma_Hstat : float):
    """Glg 4.11"""
    return sigma_HGstat / sigma_Hstat
def S_Hdyn(sigma_HGdyn : float, sigma_Hdyn : float):
    """Glg 4.11"""
    return sigma_HGdyn / sigma_Hdyn

def z_n(z : int, beta : float, beta_b : float):
    """Glg D.1.01"""
    return z / m.cos(m.radians(beta_b))**2 / m.cos(m.radians(beta))

def Y_epsilon(epsilon_alpha : float, beta_b : float):
    """Abschnitt 5.3"""
    return 0.25 + 0.75 / epsilon_alpha * m.cos(m.radians(beta_b))**2

def Y_beta(beta : float, epsilon_beta : float):
    """Glg 5.10"""
    return 1 - min(epsilon_beta, 1.) * min(beta, 30) / 120

def Y_deltarelTstat(werkstoff : din3990_5.Werkstoff, Y_S : float):
    """Abschnitt 5.6"""
    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, ):
        return (1 + 0.93 * (Y_S - 1) * m.pow(200 / werkstoff.sigma_02, 1/4)) / (1 + 0.93 * m.pow(200 / werkstoff.sigma_02, 1/4))
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Vergütungsstahl, din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit):
        return (1 + 0.82 * (Y_S - 1) * m.pow(300 / werkstoff.sigma_02, 1/4)) / (1 + 0.82 * m.pow(300 / werkstoff.sigma_02, 1/4))
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl,
                           din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl):
        return 0.44 * Y_S + 0.12
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrierterVergütungsstahl,
                           din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):
        return 0.20 * Y_S + 0.60
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Stahlguß, ):
        return 0.07 * Y_S + 0.86
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Grauguß, din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit):
        return 1.
    raise NotImplementedError(werkstoff.art)
def Y_deltarelTdyn(q_s : float):
    """Glg 5.11, 5.12"""
    if q_s >= 1.5:
        return 1.
    else:
        return 0.95
    
def Y_RrelTstat():
    """Glg 5.18"""
    return 1.
def Y_RrelTdyn(R_z : float):
    """Glg 5.18 & 5.19"""
    if R_z <= 16:
        return 1.
    else:
        return 0.9

def Y_Xstat():
    """Tabelle 5.1"""
    return 1.
def Y_Xdyn(werkstoff : din3990_5.Werkstoff, m_n : float):
    """Tabelle 5.1"""
    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl,
                         din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.SchwarzerTemperguß):
        if m_n <= 5:
            return 1.
        elif m_n < 30:
            return 1.03 - 0.006  * m_n
        else:
            return 0.85
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl, din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                           din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrierterVergütungsstahl,
                           din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):
        if m_n <= 5:
            return 1.
        elif m_n < 25:
            return 1.05 - 0.01  * m_n
        else:
            return 0.8
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Grauguß, din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit):
        if m_n <= 5:
            return 1.
        elif m_n < 25:
            return 1.075 - 0.015  * m_n
        else:
            return 0.85
    raise NotImplementedError

def Y_S(Y_Sa : float, s_Fn : float, h_Fa : float, epsilon_alphan : float):
    """Bild 5.8"""
    assert 1 <= s_Fn / h_Fa <= 1.2
    return Y_Sa * (0.6 + 0.4 * epsilon_alphan)

def Y_NT(werkstoff : din3990_5.Werkstoff):
    """Tabelle 5.2"""
    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl,
                         din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.SchwarzerTemperguß):
        return 2.5, 1.0
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl,
                           din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl):
        return 2.5, 1.0
    elif werkstoff.art in (din3990_5.Werkstoff.Art.Nitrierstahl, din3990_5.Werkstoff.Art.NitrierterVergütungsstahl,
                           din3990_5.Werkstoff.Art.Grauguß, din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit):
        return 1.6, 1.0
    elif werkstoff.art in (din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl):
        return 1.1, 1.0
    raise NotImplementedError

def sigma_FGstat(werkstoff : din3990_5.Werkstoff, Y_NTstat : float, Y_deltarelTstat : float, Y_RrelTstat : float, Y_Xstat : float):
    """Glg 5.03"""
    return werkstoff.sigma_FE * Y_NTstat * Y_deltarelTstat * Y_RrelTstat * Y_Xstat
def sigma_FGdyn(sigma_FGstat : float, werkstoff : din3990_5.Werkstoff, Y_NTdyn : float, Y_deltarelTdyn : float, Y_RrelTdyn : float, Y_Xdyn : float, N_L : float):
    """Glg 5.03 und Abschnitt 5.1.3b"""
    sigma_FGdauer = werkstoff.sigma_FE * Y_NTdyn * Y_deltarelTdyn * Y_RrelTdyn * Y_Xdyn

    if werkstoff.art in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.Vergütungsstahl,
                         din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit,
                         din3990_5.Werkstoff.Art.SchwarzerTemperguß):

        if N_L <= 10**4:
            return sigma_FGstat
        elif N_L <= 3 * 10**6:
            # Glg 5.05
            exp = 0.4037 * m.log10(sigma_FGstat / sigma_FGdauer)
            # Glg 5.04
            return sigma_FGdauer * m.pow(3 * 10**6 / N_L, exp)
        else:
            return sigma_FGdauer


    elif werkstoff.art in (din3990_5.Werkstoff.Art.Einsatzstahl, din3990_5.Werkstoff.Art.NitrierterEinsatzstahl,
                           din3990_5.Werkstoff.Art.InduktionsgehärteterStahl, din3990_5.Werkstoff.Art.FlammgehärteterStahl,
                           din3990_5.Werkstoff.Art.NitrierterVergütungsstahl, din3990_5.Werkstoff.Art.Nitrierstahl,
                           din3990_5.Werkstoff.Art.NitrokarburierterVergütungsstahl, din3990_5.Werkstoff.Art.NitrokarburierterEinsatzstahl,
                           din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit, din3990_5.Werkstoff.Art.Grauguß):
        if N_L <= 10**3:
            return sigma_FGstat
        elif N_L <= 3 * 10**6:
            # Glg 5.06
            exp = 0.2876 * m.log10(sigma_FGstat / sigma_FGdauer)
            # Glg 5.04
            return sigma_FGdauer * m.pow(3 * 10**6 / N_L, exp)
        else:
            return sigma_FGdauer

    raise NotImplementedError

def sigma_F0(F_t : float, b : float, m_n : float, Y_FS : float, Y_epsilon : float, Y_beta : float):
    """Glg 5.02"""
    return F_t / b / m_n * Y_FS * Y_epsilon * Y_beta

def sigma_Fstat(sigma_F0 : float, K_S : float, K_V : float, K_Falpha : float, K_Fbeta : float):
    """Glg 5.01"""
    return sigma_F0 * K_S * K_V * K_Fbeta * K_Falpha
def sigma_Fdyn(sigma_F0 : float, K_A : float, K_V : float, K_Falpha : float, K_Fbeta : float):
    """Glg 5.01"""
    return sigma_F0 * K_A * K_V * K_Fbeta * K_Falpha

def S_Fstat(sigma_FGstat : float, sigma_Fstat : float):
    """Glg 5.07"""
    return sigma_FGstat / sigma_Fstat
def S_Fdyn(sigma_FGdyn : float, sigma_Fdyn : float):
    """Glg 5.07"""
    return sigma_FGdyn / sigma_Fdyn

def E(m_n : float, h_fP : float, rho_fP : float, alpha_n : float, s_pr : float):
    """Glg D.5.01"""
    return m.pi / 4 * m_n - h_fP * m.tan(m.radians(alpha_n)) + s_pr / m.cos(m.radians(alpha_n)) - (1 - m.sin(m.radians(alpha_n))) * rho_fP / m.cos(m.radians(alpha_n))
def G(x : float, m_n : float, h_fP : float, rho_fP : float):
    """Glg D.5.02"""
    return rho_fP / m_n - h_fP / m_n + x
def H(m_n : float, z_n : float, E : float):
    """Glg D.5.03"""
    return 2 / z_n * (m.pi / 2 - E / m_n) - m.pi / 3
def theta(z_n : float, G : float, H : float):
    """Glg D.5.04"""
    theta = m.degrees(m.pi / 6)
    for _ in range(5):
        theta = m.degrees(2 * G / z_n * m.tan(m.radians(theta)) - H)
    return theta
def s_Fn(m_n : float, h_fP : float, rho_fP : float, alpha_n : float, s_pr : float, z_n : float, G : float, theta : float, innenverzahnt : bool):
    """Glg D.5.05"""
    if not innenverzahnt:
        return m_n * (z_n * m.sin(m.pi / 3 - m.radians(theta)) + m.sqrt(3) * (G / m.cos(m.radians(theta)) - rho_fP / m_n))
    else:
        return m_n * 2 * (m.pi / 4 + m.tan(m.radians(alpha_n)) * (h_fP - rho_fP) / m_n + (rho_fP - s_pr) / m_n / m.cos(m.radians(alpha_n)) - rho_fP / m_n * m.cos(m.pi / 6))
def d_n(m_n : float, z_n : float):
    """Glg D.5.06"""
    return m_n * z_n
def d_bn(alpha_n : float, d_n : float):
    """Glg D.5.07"""
    return d_n * m.cos(m.radians(alpha_n))
def d_an(d : float, d_a : float, d_n : float):
    """Glg D.5.08"""
    return d_n + d_a - d
def d_fn(d : float, d_f : float, d_n : float):
    return d_n + d_f - d
def alpha_an(d_bn : float, d_an : float):
    """Glg D.5.09"""
    return m.degrees(m.acos(d_bn / d_an))
def y_a(x : float, alpha_n : float, z_n : float, alpha_an : float):
    """Glg D.5.10"""
    return 1 / z_n * (m.pi / 2 + 2 * x * m.tan(m.radians(alpha_n))) + involute(alpha_n) - involute(alpha_an)
def alpha_Fan(x : float, alpha_n : float, z_n : float, d_n : float, d_an : float, innenverzahnt : bool, _print = print):
    """Glg D.5.11"""
    if not innenverzahnt:
        _d_bn = d_bn(alpha_n, d_n)
        _print("d_bn =", _d_bn)

        _alpha_an = alpha_an(_d_bn, d_an)
        _print("α_an =", _alpha_an)

        _y_a = y_a(x, alpha_n, z_n, _alpha_an)
        _print("y_a =", _y_a)

        return _alpha_an - m.degrees(_y_a)
    else:
        return alpha_n
def h_Fa(m_n : float, h_fP : float, rho_fP : float, alpha_n : float, z_n : float, G : float, theta : float, d_an : float, d_fn : float, alpha_Fan : float, innenverzahnt : bool):
    """Glg D.5.12 & D.5.15"""
    if not innenverzahnt:
        return m_n * (0.5 * z_n * (m.cos(m.radians(alpha_n)) / m.cos(m.radians(alpha_Fan)) - m.cos(m.pi / 3 - m.radians(theta))) + 0.5 * (rho_fP / m_n - G / m.cos(m.radians(theta))))
    else:
        return m_n * ((d_an - d_fn) / 2 / m_n - (m.pi / 4 + (h_fP / m_n - (d_an - d_fn) / 2 / m_n) * m.tan(m.radians(alpha_n))) * m.tan(m.radians(alpha_n)) - rho_fP / m_n * (1 - m.sin(m.pi / 6)))
def rho_F(m_n : float, rho_fP : float, z_n : float, G : float, theta : float, innenverzahnt : bool):
    """Glg D.5.13 & D.5.16"""
    if not innenverzahnt:
        return rho_fP + m_n * 2 * G**2 / m.cos(m.radians(theta)) / (z_n * m.cos(m.radians(theta))**2 - 2 * G)
    else:
        return rho_fP / 2
def AbschnittD_5(z_n : float, x : float, m_n : float, d : float, d_a : float, d_f : float, h_fP : float, rho_fP : float, s_pr : float, alpha_n : float, innenverzahnt : bool, _print = print):
    """
    s_Fn, h_Fa, rho_F, alpha_Fan
    Abschnitt D.5a und b
    """
    _E = E(m_n, h_fP, rho_fP, alpha_n, s_pr)
    _print("E =", _E)

    _G = G(x, m_n, h_fP, rho_fP)
    _print("G =", _G)

    _H = H(m_n, z_n, _E)
    _print("H =", _H)

    _theta = theta(z_n, _G, _H)
    _print("ϑ =", _theta)

    _d_n = d_n(m_n, z_n)
    _print("d_n =", _d_n)

    _d_an = d_an(d, d_a, _d_n)
    _print("d_an =", _d_an)

    _d_fn = d_fn(d, d_f, _d_n)
    _print("d_fn =", _d_fn)

    _alpha_Fan = alpha_Fan(x, alpha_n, z_n, _d_n, _d_an, innenverzahnt, _print)

    _s_Fn = s_Fn(m_n, h_fP, rho_fP, alpha_n, s_pr, z_n, _G, _theta, innenverzahnt)

    _h_Fa = h_Fa(m_n, h_fP, rho_fP, alpha_n, z_n, _G, _theta, _d_an, _d_fn, _alpha_Fan, innenverzahnt)

    _rho_F = rho_F(m_n, rho_fP, z_n, _G, _theta, innenverzahnt)

    return _s_Fn, _h_Fa, _rho_F, _alpha_Fan

def Y_Fa(m_n : float, alpha_n : float, s_Fn : float, h_Fa : float, alpha_Fan : float):
    """D.3.01"""
    return (6 * h_Fa / m_n * m.cos(m.radians(alpha_Fan))) / ((s_Fn / m_n)**2 * m.cos(m.radians(alpha_n)))

def L_a(s_Fn : float, h_Fa : float):
    """Glg D.4.02"""
    return s_Fn / h_Fa
def q_s(s_Fn : float, rho_F : float):
    """Glg D.4.03"""
    return s_Fn / 2 / rho_F
def Y_Sa(q_s : float, L_a : float):
    """Glg D.4.01"""
    assert 1 <= q_s < 8
    return (1.2 + 0.13 * L_a) * m.pow(q_s, 1 / (1.21 + 2.3 / L_a))

def Y_FS(Y_Sa : float, Y_Fa : float):
    """Glg 5.08"""
    return Y_Sa * Y_Fa

def d(z : int, m_t : float):
    return abs(z) * m_t
def d_b(d : float, alpha_t : float):
    return d * m.cos(m.radians(alpha_t))
def d_a(m_n : float, z : int, x : float, h_aP : float, k : float, d : float):
    return d + 2 * z / abs(z) * (x * m_n + h_aP + k * m_n)
def d_f(m_n : float, z : int, x : float, h_fP : float, d : float):
    return d - 2 * z / abs(z) * (h_fP - x * m_n)
def d_w(d_b : float, alpha_wt : float):
    return d_b / m.cos(m.radians(alpha_wt))
def h(m_n : float, k : float, h_aP : float, h_fP : float):
    return h_aP + k * m_n + h_fP
def z_min(m_n : float, x : float, h_aP0, alpha_t : float, beta : float):
    return 2 * m.cos(m.radians(beta)) * (h_aP0 / m_n - x) / m.sin(m.radians(alpha_t))**2
def gamma(z : int, x : float, alpha_n : float, alpha_t : float):
    return inverse_involute( m.pi / 2 / z + 2 * x / z * m.tan(m.radians(alpha_n)) + involute(alpha_t))
def d_amax(z : int, m_t : float, alpha_t : float, gamma : float):
    return m_t * z * m.cos(m.radians(alpha_t)) / m.cos(m.radians(gamma))
def alpha_t(alpha_n : float, beta : float):
    return m.degrees(m.atan(m.tan(m.radians(alpha_n)) / m.cos(m.radians(beta))))
def alpha_wt(z : tuple[int, int], x : tuple[float, float], alpha_n : float, alpha_t : float):
    return inverse_involute(involute(alpha_t) + 2 * sum(x) / sum(z) * m.tan(m.radians(alpha_n)) )
def u(z : tuple[int, int]):
    return z[Rad] / z[Ritzel]
def m_t(m_n : float, beta : float):
    return m_n / m.cos(m.radians(beta))
def beta_b(beta : float, alpha_t : float):
    return m.degrees(m.atan(m.tan(m.radians(beta)) * m.cos(m.radians(alpha_t))))
def a_w(z_2 : int, d_w : float):
    return (d_w[Rad] + z_2 / abs(z_2) * d_w[Ritzel]) / 2
def p_n(m_n : float):
    return m_n * m.pi
def p_t(p_n : float, beta : float):
    return p_n / m.cos(m.radians(beta))
def epsilon_alpha(z : tuple[int, int], d_a : tuple[float, float], d_b : tuple[float, float], p_t : float, alpha_t : float, alpha_wt : float):
    return (m.sqrt(d_a[Ritzel]**2 - d_b[Ritzel]**2) + z[Rad] / abs(z[Rad]) * m.sqrt(d_a[Rad]**2 - d_b[Rad]**2)
                              - (d_b[Ritzel] + z[Rad] / abs(z[Rad]) * d_b[Rad]) * m.tan(m.radians(alpha_wt))) / (2 * p_t * m.cos(m.radians(alpha_t)))
def epsilon_beta(b : float, p_n : float, beta : float):
    return b * m.sin(m.radians(beta)) / p_n
def epsilon_gamma(epsilon_alpha : float, epsilon_beta : float):
    return epsilon_alpha + epsilon_beta

class DIN_21771:
    def __init__(self,
            m_n : float,
            z: tuple[int, int],
            x: tuple[float, float],
            bezugsprofil : tuple[Profil, Profil],
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

        assert all(z > 0 for z in z), "Innenverzahnung nicht implementiert"

        self.m_n = m_n
        self.z = z
        self.x = x
        self.beta = beta
        self.k = k

        _print("Getriebegeometrie")
        [_print(key, "=", value) for key, value in vars(self).items()]

        assert all(bezugsprofil[idx].alpha_n == bezugsprofil[Ritzel].alpha_n for idx in _indices)
        self.alpha_n = bezugsprofil[Ritzel].alpha_n
        self.h_aP = tuple(bezugsprofil[idx].h_aP_s * self.m_n for idx in _indices)
        self.h_fP = tuple(bezugsprofil[idx].h_fP_s * self.m_n for idx in _indices)
        self.rho_fP = tuple(bezugsprofil[idx].rho_fP_s * self.m_n for idx in _indices)
        self.s_pr = tuple(bezugsprofil[idx].s_pr_s * self.m_n for idx in _indices)
        _print("α_n =", self.alpha_n)
        _print("h_aP =", self.h_aP)
        _print("h_fP =", self.h_fP)
        _print("ρ_fP =", self.rho_fP)

        self.alpha_t = alpha_t(self.alpha_n, self.beta)
        _print("α_t =", self.alpha_t)

        self.alpha_wt = alpha_wt(self.z, self.x, self.alpha_n, self.alpha_t)
        _print("α_wt =", self.alpha_wt)

        self.u = u(self.z)
        _print("u =", self.u)

        self.m_t = m_t(self.m_n, self.beta)
        _print("m_t =", self.m_t)

        self.d = tuple(d(self.z[idx], self.m_t) for idx in _indices)
        _print("d =", self.d)

        self.d_b = tuple(d_b(self.d[idx], self.alpha_t) for idx in _indices)
        _print("d_b =", self.d_b)

        self.d_a = tuple(d_a(self.m_n, self.z[idx], self.x[idx], self.h_aP[idx], self.k, self.d[idx]) for idx in _indices)
        _print("d_a =", self.d_a)

        self.d_f = tuple(d_f(self.m_n, self.z[idx], self.x[idx], self.h_fP[idx], self.d[idx]) for idx in _indices)
        _print("d_f =", self.d_f)

        self.d_w = tuple(d_w(self.d_b[idx], self.alpha_wt) for idx in _indices)
        _print("d_w =", self.d_w)

        if b != None:
            self.b = b
        elif b_d_1_verhältnis != None:
            self.b = self.d[Ritzel] * b_d_1_verhältnis
        else:
            raise ValueError("either b or b_d_1_verhältnis must be specified as argument")
        _print("b =", self.b)

        self.beta_b = beta_b(self.beta, self.alpha_t)
        _print("β_b =", self.beta_b)
        
        self.h = tuple(h(self.m_n, self.k, self.h_aP[idx], self.h_fP[idx]) for idx in _indices)
        _print("h =", self.h)

        self.a_w = a_w(self.z[Rad], self.d_w)
        _print("a_w =", self.a_w)

        # Profilüberdeckung

        self.p_n = p_n(self.m_n)
        _print("p_n =", self.p_n)

        self.p_t = p_t(self.p_n, self.beta)
        _print("p_t =", self.p_t)

        self.epsilon_alpha = epsilon_alpha(self.z, self.d_a, self.d_b, self.p_t, self.alpha_t, self.alpha_wt)
        _print("ε_α =", self.epsilon_alpha)

        self.epsilon_beta = epsilon_beta(self.b, self.p_n, self.beta)
        _print("ε_β =", self.epsilon_beta)

        self.epsilon_gamma = epsilon_gamma(self.epsilon_alpha, self.epsilon_beta)
        _print("ε_γ =", self.epsilon_gamma)
        assert self.epsilon_gamma > 1

        # Unterschnitt

        self.h_aP0 = self.h_fP
        _print("h_aP0 =", self.h_aP0)

        _z_min = tuple(z_min(self.m_n, self.x[idx], self.h_aP0[idx], self.alpha_t, self.beta) for idx in _indices)
        _print("z_min =", _z_min)
        assert self.z[Ritzel] > _z_min[Ritzel]
        assert self.z[Rad] > _z_min[Rad]

        # Spitzwerden

        _gamma = tuple(gamma(self.z[idx], self.x[idx], self.alpha_n, self.alpha_t) for idx in _indices)
        _print("γ =", _gamma)

        _d_amax = tuple(d_amax(self.z[idx], self.m_t, self.alpha_t, _gamma[idx]) for idx in _indices)
        _print("d_amax =", _d_amax)
        assert self.d_a[Ritzel] <= _d_amax[Ritzel]
        assert self.d_a[Rad] <= _d_amax[Rad]

        _print()
        return
  
_K_V = K_V
_K_Hbeta = K_Hbeta
_K_Fbeta = K_Fbeta
_Z_LVRdyn = Z_LVRdyn
class Calculator:
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
        
        assert not innenverzahnt, "Innenverzahnte Getriebe sind nicht implementiert"

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

        self.v = v(self.n[Ritzel], self.geometrie.d[Ritzel])
        _print("v =", self.v)

        self.F_t = F_t(self.P, self.n[Ritzel], self.geometrie.d[Ritzel])
        _print("F_t =", self.F_t)
  
        self.T = tuple(T(self.P, self.n[idx]) for idx in _indices)
        _print("T =", self.T)

        self.K_V = tuple(_K_V(self.geometrie.z[Ritzel], self.v, self.geometrie.u, self.F_t, self.K_A, self.geometrie.b, self.geometrie.epsilon_beta, geradverzahnt, self.verzahnungsqualität[idx],
                            _print) if K_V is None else K_V for idx, K_V in zip(_indices, self.K_V))
        _print("K_V =", self.K_V)

        self.K_Hbeta = tuple(_K_Hbeta(self.F_t, self.K_A, self.K_V[idx], self.v, self.geometrie.d[Ritzel], self.geometrie.b, self.doppelschrägverzahnt, self.werkstoff,
                                      self.A, self.f_ma[idx], self.s[idx], self.bild3_1[idx], self.bild3_2[idx], self.stützwirkung[idx], self.l[idx], self.d_sh[idx], self.d_sh[Ritzel], _print)
                             if K_Hbeta is None else K_Hbeta for idx, K_Hbeta in zip(_indices, self.K_Hbeta))
        _print("K_Hβ =", self.K_Hbeta)
        
        self.K_Fbeta = tuple(_K_Fbeta(self.F_t, self.K_A, self.K_Hbeta[idx], self.geometrie.b, self.geometrie.h[idx]) if K_Fbeta is None else K_Fbeta for idx, K_Fbeta in zip(_indices, self.K_Fbeta))
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

        self.R_z100 = R_z100(self.R_z, self.geometrie.a_w)
        _print("R_z100 =", self.R_z100)

        self.Z_LVRstat = Z_LVRstat()
        if self.Z_LVRdyn == None:
            self.Z_LVRdyn = _Z_LVRdyn(self.fertigungsverfahren, self.R_z100)
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

        self.sigma_HGstat = tuple(sigma_HGstat(self.werkstoff[idx], self.Z_NTstat[idx], self.Z_LVRstat, self.Z_W[idx], self.Z_Xstat[idx]) for idx in _indices)
        self.sigma_HGdyn = tuple(sigma_HGdyn(self.sigma_HGstat[idx], self.werkstoff[idx], self.Z_NTdyn[idx], self.Z_LVRdyn, self.Z_W[idx], self.Z_Xdyn[idx], self.N_L, self.gewisseGrübchenbildung)
                                 for idx in _indices)
        _print("σ_HGstat =", self.sigma_HGstat)
        _print("σ_HGdyn =", self.sigma_HGdyn)

        self.sigma_H0 = sigma_H0(self.F_t, self.geometrie.d[Ritzel], self.geometrie.b, self.geometrie.u, self.Z_H, self.Z_E, self.Z_epsilon, self.Z_beta)
        _print("σ_H0 =", self.sigma_H0)

        self.sigma_Hstat = tuple(sigma_Hstat(self.K_S, self.K_V[idx], self.K_Halpha[idx], self.K_Hbeta[idx], self.sigma_H0, self.Z_B if idx == Ritzel else self.Z_D) for idx in _indices)
        self.sigma_Hdyn = tuple(sigma_Hdyn(self.K_A, self.K_V[idx], self.K_Halpha[idx], self.K_Hbeta[idx], self.sigma_H0, self.Z_B if idx == Ritzel else self.Z_D) for idx in _indices)
        _print("σ_Hstat =", self.sigma_Hstat)
        _print("σ_Hdyn =", self.sigma_Hdyn)

        self.S_Hstat = tuple(S_Hstat(self.sigma_HGstat[idx], self.sigma_Hstat[idx]) for idx in _indices)
        self.S_Hdyn = tuple(S_Hdyn(self.sigma_HGdyn[idx], self.sigma_Hdyn[idx]) for idx in _indices)


        check_value(Ritzel, self.S_Hstat, S_Hstatmin, "S_Hstat", "statische Grübchensicherheit")
        check_value(Ritzel, self.S_Hdyn, S_Hdynmin, "S_Hdyn", "dynamische Grübchensicherheit")
        check_value(Rad, self.S_Hstat, S_Hstatmin, "S_Hstat", "statische Grübchensicherheit")
        check_value(Rad, self.S_Hdyn, S_Hdynmin, "S_Hdyn", "dynamische Grübchensicherheit")
        _print()
        
        # Zahnfußtragfähigkeit

        self.z_n = tuple(z_n(self.geometrie.z[idx], self.geometrie.beta, self.geometrie.beta_b) for idx in _indices)
        _print("z_n =", self.z_n)

        _AbschnittD_5 = tuple(AbschnittD_5(self.z_n[idx], self.geometrie.x[idx], self.geometrie.m_n, self.geometrie.d[idx], self.geometrie.d_a[idx], self.geometrie.d_f[idx], self.geometrie.h_fP[idx],
                                           self.geometrie.rho_fP[idx], self.geometrie.s_pr[idx], self.geometrie.alpha_n, self.innenverzahnt if idx == Rad else False, _print) for idx in _indices)
        self.s_Fn, self.h_Fa, self.rho_F, self.alpha_Fan  = zip(_AbschnittD_5[0], _AbschnittD_5[1])
        _print("s_Fn =", self.s_Fn)
        _print("h_Fa =", self.h_Fa)
        _print("ρ_F =", self.rho_F)
        _print("α_Fan =", self.alpha_Fan)

        self.Y_Fa = tuple(Y_Fa(self.geometrie.m_n, self.geometrie.alpha_n, self.s_Fn[idx], self.h_Fa[idx], self.alpha_Fan[idx]) for idx in _indices)
        _print("Y_Fa =", self.Y_Fa)

        self.L_a = tuple(L_a(self.s_Fn[idx], self.h_Fa[idx]) for idx in _indices)
        _print("L_a =", self.L_a)

        self.q_s = tuple(q_s(self.s_Fn[idx], self.rho_F[idx]) for idx in _indices)
        _print("q_s =", self.q_s)

        self.Y_Sa = tuple(Y_Sa(self.q_s[idx], self.L_a[idx]) for idx in _indices)
        _print("Y_Sa =", self.Y_Sa)

        self.Y_FS = tuple(Y_FS(self.Y_Sa[idx], self.Y_Fa[idx]) for idx in _indices)
        _print("Y_FS =", self.Y_FS)

        self.Y_beta = Y_beta(self.geometrie.beta, self.geometrie.epsilon_beta)
        _print("Y_β =", self.Y_beta)
 
        _epsilon_alphan = epsilon_alphan(self.geometrie.epsilon_alpha, self.geometrie.beta_b)
        _print("ε_αn =", _epsilon_alphan)

        self.Y_S = tuple(Y_S(self.Y_Sa[idx], self.s_Fn[idx], self.h_Fa[idx], _epsilon_alphan) for idx in _indices)
        _print("Y_S =", self.Y_S)

        self.Y_deltarelTstat = tuple(Y_deltarelTstat(self.werkstoff[idx], self.Y_S[idx]) for idx in _indices)
        self.Y_deltarelTdyn = tuple(Y_deltarelTdyn(self.q_s[idx]) for idx in _indices)
        _print("Y_δrelTstat =", self.Y_deltarelTstat)
        _print("Y_δrelTdyn =", self.Y_deltarelTdyn)

        self.Y_RrelTstat = tuple(Y_RrelTstat() for _ in _indices)
        self.Y_RrelTdyn = tuple(Y_RrelTdyn(self.R_z[idx]) for idx in _indices)
        _print("Y_RrelTstat =", self.Y_RrelTstat)
        _print("Y_RrelTdyn =", self.Y_RrelTdyn)

        self.Y_Xstat = tuple(Y_Xstat() for _ in _indices)
        self.Y_Xdyn = tuple(Y_Xdyn(self.werkstoff[idx], self.geometrie.m_n) for idx in _indices)
        _print("Y_Xstat =", self.Y_Xstat)
        _print("Y_Xdyn =", self.Y_Xdyn)

        Y_NTritzel = Y_NT(self.werkstoff[Ritzel])
        Y_NTrad = Y_NT(self.werkstoff[Rad])
        self.Y_NTstat = Y_NTritzel[0], Y_NTrad[0]
        self.Y_NTdyn = Y_NTritzel[1], Y_NTrad[1]
        _print("Y_NTstat =", self.Y_NTstat)
        _print("Y_NTdyn =", self.Y_NTdyn)

        self.sigma_FGstat = tuple(sigma_FGstat(self.werkstoff[idx], self.Y_NTstat[idx], self.Y_deltarelTstat[idx], self.Y_RrelTstat[idx], self.Y_Xstat[idx]) for idx in _indices)
        self.sigma_FGdyn = tuple(sigma_FGdyn(self.sigma_FGstat[idx], self.werkstoff[idx], self.Y_NTdyn[idx], self.Y_deltarelTdyn[idx], self.Y_RrelTdyn[idx], self.Y_Xdyn[idx], self.N_L)
                                 for idx in _indices)
        _print("σ_FGstat =", self.sigma_FGstat)
        _print("σ_FGdyn =", self.sigma_FGdyn)

        self.sigma_F0 = tuple(sigma_F0(self.F_t, self.geometrie.b, self.geometrie.m_n, self.Y_FS[idx], self.Y_epsilon, self.Y_beta) for idx in _indices)
        _print("σ_F0 =", self.sigma_F0)

        self.sigma_Fstat = tuple(sigma_Fstat(self.sigma_F0[idx], self.K_S, self.K_V[idx], self.K_Falpha[idx], self.K_Fbeta[idx]) for idx in _indices)
        self.sigma_Fdyn = tuple(sigma_Fdyn(self.sigma_F0[idx], self.K_A, self.K_V[idx], self.K_Falpha[idx], self.K_Fbeta[idx]) for idx in _indices)
        _print("σ_Fstat =", self.sigma_Fstat)
        _print("σ_Fdyn =", self.sigma_Fdyn)

        self.S_Fstat = tuple(S_Fstat(self.sigma_FGstat[idx], self.sigma_Fstat[idx]) for idx in _indices)
        self.S_Fdyn = tuple(S_Fdyn(self.sigma_FGdyn[idx], self.sigma_Fdyn[idx]) for idx in _indices)

        check_value(Ritzel, self.S_Fstat, S_Fstatmin, "S_Fstat", "statische Zahnbruchsicherheit")
        check_value(Ritzel, self.S_Fdyn, S_Fdynmin, "S_Fdyn", "dynamische Zahnbruchsicherheit")
        check_value(Rad, self.S_Fstat, S_Fstatmin, "S_Fstat", "statische Zahnbruchsicherheit")
        check_value(Rad, self.S_Fdyn, S_Fdynmin, "S_Fdyn", "dynamische Zahnbruchsicherheit")
        _print()
        return


