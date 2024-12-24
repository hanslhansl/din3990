from dataclasses import dataclass
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

Ritzel = 0
Rad = 1

@dataclass
class Profil:
    alpha_n : float
    h_aP_s : float
    h_fP_s : float
    rho_fP_s : float
Normalprofil1 =     Profil(20, 1, 1.25, 0.250)
Normalprofil2 =     Profil(20, 1, 1.25, 0.375)
Protuberanzprofil = Profil(20, 1, 1.40, 0.400)


class Tabelle_3_2(Enum):
    ohneBreitenballigkeltOderEndrücknahme = 0.023
    mitSinnvollerBreitenballigkeit = 0.012
    mitSinnvollerEndrücknahme = 0.016

class Bild_3_1(IntEnum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4
    f = 5

class Bild_3_2(IntEnum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4

class Fertigungsverfahren(IntEnum):
    wälzgefrästWälzgestoßenWälzgehobelt = 0
    geläpptGeschliffenGeschabt = 1

def v(n : float, d : float):
    return n * m.pi * d / 60000
def F_t(P : float, n : float, d : float):
    """Glg 3.01"""
    return 60000000 * P / n / m.pi / d
def T(P : float, n : float):
    """Glg 3.02"""
    return 30000 * P / n / m.pi

def K_1(verzahnungsqualität : int, geradverzahnt : bool):
    if geradverzahnt:
        match verzahnungsqualität:
            case 6:
                return 9.6
            case 7:
                return 15.3
            case 8:
                return 24.5
            case 9:
                return 34.5
            case 10:
                return 53.6
            case 11:
                return 76.6
            case 12:
                return 122.5
    else:
        match verzahnungsqualität:
            case 6:
                return 8.5
            case 7:
                return 13.6
            case 8:
                return 21.8
            case 9:
                return 30.7
            case 10:
                return 47.7
            case 11:
                return 68.2
            case 12:
                return 109.1
    raise ValueError(f"Falsche Verzahnungsqualität {verzahnungsqualität}")
def K_2(geradverzahnt : bool):
    if geradverzahnt:
        return 0.0193
    return 0.0087
def K_V(z_1 : int, v : float, u : float, F_t : float, K_A : float, b : float, epsilon_beta : float, geradverzahnt : bool, verzahnungsqualität : tuple[int, int], _print):
    """Abschnitt 3.3"""
    # Glg 3.04
    temp1 = z_1 * v / 100 * m.sqrt(u**2 / (1 + u**2))
    _print("z_1 * v / 100 * sqrt(u^2 / (1 + u^2)) =", temp1)
    assert temp1 < 10

     # Abschnitt 3.3.1 b
    temp2 = max(F_t * K_A / b, 100)
    _print("F_t * K_A / b =", temp2)

    def _K_V(verzahnungsqualität : int, geradverzahnt : bool):
        return 1 + (K_1(verzahnungsqualität, geradverzahnt) / temp2 + K_2(geradverzahnt)) * temp1

    if geradverzahnt:
        return _K_V(verzahnungsqualität[Ritzel], True), _K_V(verzahnungsqualität[Rad], True)
    elif epsilon_beta >= 1:
        return _K_V(verzahnungsqualität[Ritzel], False), _K_V(verzahnungsqualität[Rad], False)
    else:
        K_Valpha = _K_V(verzahnungsqualität[Ritzel], True), _K_V(verzahnungsqualität[Rad], True)
        K_Vbeta = _K_V(verzahnungsqualität[Ritzel], False), _K_V(verzahnungsqualität[Rad], False)
        return interpolate(K_Valpha[Ritzel], K_Vbeta[Ritzel], epsilon_beta), interpolate(K_Valpha[Rad], K_Vbeta[Rad], epsilon_beta)

def F_m(F_t : float, K_A : float, K_V : float):
    """Glg 3.07"""
    return F_t * K_A * K_V
def K_s(stützwirkung : bool, bild_3_2 : Bild_3_2):
    """Bild 3.2"""
    match bild_3_2:
        case Bild_3_2.a:
            return 0.48 if stützwirkung else 0.8
        case Bild_3_2.b:
            return -0.48 if stützwirkung else -0.8
        case Bild_3_2.c:
            return 1.33
        case Bild_3_2.d:
            return -0.36 if stützwirkung else -0.6
        case Bild_3_2.e:
            return -0.6 if stützwirkung else -1.0
def K_Hbeta(F_t : float, K_A : float, b : float, K_V : tuple[float, float], _print):
    """Abschnitt 3.4"""

    # Abschnitt 3.4.1
    assert(F_t / b * K_A >= 100)

    _F_m = F_m(F_t, K_A, K_V[Ritzel]), F_m(F_t, K_A, K_V[Rad])
    _print("F_m =", _F_m)