from dataclasses import dataclass
from typing import Optional
from enum import Enum
import math as m
from scipy import optimize


def to_float(val) -> float:
    return val

def involute(alpha):
    return m.tan(m.radians(alpha)) - m.radians(alpha)
def inverse_involute(alpha, anfangswert = 20):
    try:
        return float(optimize.newton(lambda x: involute(x) - alpha, anfangswert))
    except RuntimeError:
        assert(False)

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

class Bild_3_1(Enum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4
    f = 5

class Bild_3_2(Enum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4

class Fertigungsverfahren(Enum):
    wälzgefrästWälzgestoßenWälzgehobelt = 0
    geläpptGeschliffenGeschabt = 1

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

        self.u = self.z[Rad2] / self.z[Rad1]
        _print("u =", self.u)

        def d(idx):
            return self.z[idx] * self.m_n / m.cos(m.radians(self.beta))
        self.d = d(Rad1), d(Rad2)
        _print("d =", self.d)

        def d_b(idx):
            return self.d[idx] * m.cos(m.radians(self.alpha_t))
        self.d_b = d_b(Rad1), d_b(Rad2)
        _print("d_b =", self.d_b)

        def d_a(idx):
            return self.d[idx] + 2 * (self.x[idx] * self.m_n + self.h_aP + self.k * self.m_n)
        self.d_a = d_a(Rad1), d_a(Rad2)
        _print("d_a =", self.d_a)

        def d_f(idx):
            return self.d[idx] - 2 * (self.h_fP - self.x[idx] * self.m_n)
        self.d_f = d_f(Rad1), d_f(Rad2)
        _print("d_f =", self.d_f)

        def d_w(idx):
            return self.d_b[idx] / m.cos(m.radians(self.alpha_wt))
        self.d_w = d_w(Rad1), d_w(Rad2)
        _print("d_w =", self.d_w)

        if b != None:
            self.b = b
        elif b_d_1_verhältnis != None:
            self.b = self.d[Rad1] * b_d_1_verhältnis
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

        self.epsilon_alpha = (m.sqrt(self.d_a[Rad1]**2 - self.d_b[Rad1]**2) + m.sqrt(self.d_a[Rad2]**2 - self.d_b[Rad2]**2) - sum(self.d_b) * m.tan(m.radians(self.alpha_wt))) / (2 * self.p_t * m.cos(m.radians(self.alpha_t)))
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
        z_min = _z_min(Rad1), _z_min(Rad2)
        _print("z_min =", z_min)
        assert self.z[Rad1] > z_min[Rad1]
        assert self.z[Rad2] > z_min[Rad2]

        # Spitzwerden

        def _gamma(idx):
            return inverse_involute( m.pi / 2 / self.z[idx] + 2 * self.x[idx] / self.z[idx] * m.tan(m.radians(self.alpha_n)) + involute(self.alpha_t) )
        gamma = _gamma(Rad1), _gamma(Rad2)
        _print("γ =", gamma)

        def _d_amax(idx):
            return self.m_t * self.z[idx] * m.cos(m.radians(self.alpha_t)) / m.cos(m.radians(gamma[idx]))
        d_amax = _d_amax(Rad1), _d_amax(Rad2)
        _print("d_amax =", d_amax)
        assert self.d_a[Rad1] <= d_amax[Rad1]
        assert self.d_a[Rad2] <= d_amax[Rad2]

        _print()
        return
  
class DIN_3990_11:
    def __init__(self,
                geometrie : DIN_21771,
                P : float,
                n_1 : float,
                verzahnungsqualität : tuple[int, int],
                werkstoff : tuple[Werkstoff, Werkstoff],
                K_A : float,
                K_S : float,
                R_z: tuple[float, float],
                doppelschrägverzahnt: bool = False,
                innenverzahnt: bool = False,
                s_pr: float = 0,
                
                K_V : Optional[tuple[float, float]] = None,
                K_Hbeta : Optional[tuple[float, float]] = None,
                A : Optional[Tabelle_3_2] = None,
                f_ma : Optional[tuple[float, float]] = None,
                bild_3_1 : Optional[tuple[Bild_3_1, Bild_3_1]] = None,
                s : Optional[tuple[float, float]] = None,
                stützwirkung : Optional[tuple[bool, bool]] = None,
                bild_3_2 : Optional[tuple[Bild_3_2, Bild_3_2]] = None,
                l : Optional[tuple[float, float]] = None,
                d_sh : Optional[tuple[float, float]] = None,
                Z_LVRdyn : Optional[float] = None,
                fertigungsverfahren : Optional[tuple[Fertigungsverfahren, Fertigungsverfahren]] = None,
                K_Fbeta : Optional[tuple[float, float]] = None,
                K_Halpha : Optional[tuple[float, float]] = None,
                K_Falpha : Optional[tuple[float, float]] = None,
                
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
        - verzahnungsqualität: siehe Tabelle 3.1
        - werkstoff
        - K_A: siehe Tabelle A.1
        - K_S: ersetz K_A für die statische Berechnung
        - R_z: gemittelte Rauhtiefe
        - doppelschrägverzahnt
        - innenverzahnt
        - s_pr: float, Fussfreischnitt
        - _print: the function used for printing
        - _assert if True, safety is asserted which additionally requires the same arguments as passed to Getriebe.is_save()

        Optional parameters:
        - K_V
        - K_Hbeta
          or 
            - A: Tabelle_3_2, siehe Tabelle 3.2
            - f_ma: tuple[float, float], siehe Abschnitt 3.4.2.4
              if f_ma != 0:
                - bild_3_1: tuple[Bild_3_1, Bild_3_1], siehe Bild 3.1
            - s: tuple[float, float], siehe Bild 3.2
              if s != 0:
                - stützwirkung: tuple[bool, bool], siehe Bild 3.2
                - bild_3_2: tuple[Bild_3_2, Bild_3_2], siehe Bild 3.2
                - l: tuple[float, float], siehe Bild 3.2
                - d_sh: tuple[float, float], Wellendurchmesser
        - Z_LVRdyn: float, siehe Abschnitt 4.8
          or
            - fertigungsverfahren: tuple[Fertigungsverfahren, Fertigungsverfahren], siehe Abschnitt 4.8
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
        self.doppelschrägverzahnt = doppelschrägverzahnt
        self.innenverzahnt = innenverzahnt
        self.s_pr = s_pr

        self.K_V = K_V
        self.K_Hbeta = K_Hbeta
        self.A = A
        self.f_ma = f_ma
        self.bild_3_1 = bild_3_1
        self.s = s
        self.stützwirkung = stützwirkung
        self.bild_3_2 = bild_3_2
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

        assert all(vq in range(6, 13) for vq in self.verzahnungsqualität)
        
        _print("DIN 3990-11")
        _print("n =", self.n)

        if self.doppelschrägverzahnt:
            self.b_B = self.geometrie.b / 2
            _print("b_B =", self.b_B)

        self.v = self.n[Rad1] * m.pi * self.geometrie.d[Rad1] / 60000
        _print("v =", self.v)

        self.F_t = 1000 * self.P / self.v
        _print("F_t =", self.F_t)
  
        def T(idx):
            return self.F_t * self.geometrie.d[idx] / 2000
        self.T = T(Rad1), T(Rad2)
        _print("T =", self.T)

        if self.K_V == None:
            # Glg 3.04
            temp1 = self.geometrie.z[Rad1] * self.v / 100 * m.sqrt(self.geometrie.u**2 / (1 + self.geometrie.u**2))
            print("z_1 * v / 100 * sqrt(u^2 / (1 + u^2)) =", temp1)
            assert temp1 < 10  

            # Abschnitt 3.3.1 b
            temp2 = max(self.F_t * self.K_A / self.geometrie.b, 100)
            print("F_t * K_A / b =", temp2)

            def K_V(idx, geradverzahnt) -> float:
                if geradverzahnt:
                    match self.verzahnungsqualität[idx]:
                        case 6:
                            K_1 = 9.6
                        case 7:
                            K_1 = 15.3
                        case 8:
                            K_1 = 24.5
                        case 9:
                            K_1 = 34.5
                        case 10:
                            K_1 = 53.6
                        case 11:
                            K_1 = 76.6
                        case 12:
                            K_1 = 122.5
                    K_2 = 0.0193
                else:
                    match self.verzahnungsqualität[idx]:
                        case 6:
                            K_1 = 8.5
                        case 7:
                            K_1 = 13.6
                        case 8:
                            K_1 = 21.8
                        case 9:
                            K_1 = 30.7
                        case 10:
                            K_1 = 47.7
                        case 11:
                            K_1 = 68.2
                        case 12:
                            K_1 = 109.1
                    K_2 = 0.0087
                print("K_1 =", K_1)
                print("K_2 =", K_2)
                return 1 + (K_1 / temp2 + K_2) * temp1
            if self.geometrie.beta == 0:
                self.K_V = K_V(Rad1, True), K_V(Rad2, True)
            elif self.geometrie.epsilon_beta >= 1:
                self.K_V = K_V(Rad1, False), K_V(Rad2, False)
            else:
                K_Valpha = K_V(Rad1, True), K_V(Rad2, True)
                K_Vbeta = K_V(Rad1, False), K_V(Rad2, False)
                _print("K_Vα =", K_Valpha)
                _print("K_Vβ =", K_Vbeta)
                def interp(idx) -> float:
                    return K_Valpha[idx] - self.geometrie.epsilon_beta * (K_Valpha[idx] - K_Vbeta[idx])
                self.K_V = interp(Rad1), interp(Rad2)
        _print("K_V =", self.K_V)

        if self.K_Hbeta == None:
            # Abschnitt 3.4.1
            assert(self.F_t / self.geometrie.b * self.K_A >= 100)

            # Glg 3.07
            def F_m(idx) -> float:
                return self.F_t * self.K_A * self.K_V[idx]
            self.F_m = F_m(Rad1), F_m(Rad2)
            _print("F_m =", self.F_m)

            def K_s(idx):
                stütz = self.stützwirkung[idx]
                match self.bild_3_2[idx]:
                    case Bild_3_2.a:
                        return 0.48 if stütz else 0.8
                    case Bild_3_2.b:
                        return -0.48 if stütz else -0.8
                    case Bild_3_2.c:
                        return 1.33
                    case Bild_3_2.d:
                        return -0.36 if stütz else -0.6
                    case Bild_3_2.e:
                        return -0.6 if stütz else -1.0
                raise ValueError("Unknown bild_3_2 option ", self.bild_3_2[idx])

            # Glg 3.14, 3.15
            def f_sh(idx) -> float:
                def temp(idx) -> float:
                    return 0. if self.s[idx] == 0 else K_s(idx) * self.l[idx] * self.s[idx] / self.geometrie.d[Rad1]**2 * (self.geometrie.d[Rad1] / self.d_sh[idx])**4
                if self.A == Tabelle_3_2.ohneBreitenballigkeltOderEndrücknahme:
                    A = 0.023
                elif self.A == Tabelle_3_2.mitSinnvollerBreitenballigkeit:
                    A = 0.012
                elif self.A == Tabelle_3_2.mitSinnvollerEndrücknahme:
                    A = 0.016

                if not self.doppelschrägverzahnt:
                    return self.F_m[idx] / self.geometrie.b * A * (abs(1 + temp(idx) - 0.3) + 0.3) * (self.geometrie.b / self.geometrie.d[Rad1])**2
                else:
                    return self.F_m[idx] / self.geometrie.b * 2 * A * (abs(1.5 + temp(idx) - 0.3) + 0.3) * (self.b_B / self.geometrie.d[Rad1])**2
            self.f_sh = f_sh(Rad1), f_sh(Rad2)
            _print("f_sh =", self.f_sh)

            # Glg 3.09
            def F_betax(idx) -> float:
                def multi(idx):
                    match self.bild_3_1[idx]:
                        case Bild_3_1.a | Bild_3_1.f:
                            return -1
                        case Bild_3_1.b | Bild_3_1.e:
                            return 1
                        case Bild_3_1.c:
                            B_s = 1.5 if self.doppelschrägverzahnt else 1
                            return 1 if abs(self.K_s) * self.l * self.s / self.d[Rad1]**2 * (self.d[Rad1] / self.d_sh[Rad1])**4 <= B_s else -1
                        case Bild_3_1.d:
                            B_s = 1.5 if self.doppelschrägverzahnt else 1
                            return 1 if abs(self.K_s) * self.l * self.s / self.d[Rad1]**2 * (self.d[Rad1] / self.d_sh[Rad1])**4 >= B_s - 0.3 else -1
                        case _:
                            raise ValueError("Unknown bild_3_1 option ", self.bild_3_1[idx])
                if self.f_ma[idx] == 0:
                    return abs(1.33 * self.f_sh[idx])
                else:
                    return abs(1.33 * self.f_sh[idx] + multi(idx) * self.f_ma[idx])
            self.F_betax = F_betax(Rad1), F_betax(Rad2)
            _print("F_betax =", self.F_betax)

            # Abschnitt 3.4.2.6
            def y_beta(idx : int):
                werkstoff = self.werkstoff[idx]
                match werkstoff.art:
                    case Werkstoff.Art.Baustahl | Werkstoff.Art.vergüteterStahl:
                        y_beta = 320 / werkstoff.sigma_Hlim * self.F_betax[idx]
                        if self.v <= 5:
                            pass
                        elif self.v <= 10:
                            assert y_beta <= 25600 / werkstoff.sigma_Hlim
                        else:
                            assert y_beta <= 12800 / werkstoff.sigma_Hlim
                        return y_beta
                    case Werkstoff.Art.einsatzgehärteterStahl | Werkstoff.Art.nitrierterStahl | Werkstoff.Art.nitrokarburierterStahl:
                        y_beta = 0.15 * self.F_betax[idx]
                        assert y_beta <= 6
                        return y_beta
                raise NotImplementedError
            self.y_beta = (y_beta(Rad1) + y_beta(Rad2)) / 2
            _print("y_beta =", self.y_beta)

            # Glg 3.08
            def F_betay(idx):
                return self.F_betax[idx] - self.y_beta
            self.F_betay = F_betay(Rad1), F_betay(Rad2)
            _print("F_betay =", self.F_betay)

            # Abschnitt 3.4.3.1
            self.c_gamma = 20
            _print("c_gamma =", self.c_gamma)

            # Glg 3.20, 3.21
            def K_Hbeta(idx):
                val = 1 + self.c_gamma * self.F_betay[idx] / (2 * self.F_m[idx] / self.geometrie.b)
                if val <= 2:
                    return val
                return m.sqrt(2 * self.c_gamma * self.F_betay[idx] / (self.F_m[idx] / self.geometrie.b))
            self.K_Hbeta = K_Hbeta(Rad1), K_Hbeta(Rad2)
        _print("K_Hβ =", self.K_Hbeta)
        
        # Glg 3.22
        if self.K_Fbeta == None:
            def K_Fbeta(idx):
                h_b = min(self.geometrie.h / self.geometrie.b, 1. / 3.)
                return m.pow(self.K_Hbeta[idx], (1 / (1 + h_b + h_b**2)))

            assert(self.F_t / self.geometrie.b * self.K_A >= 100)   # Abschnitt 3.4.1
            self.K_Fbeta = K_Fbeta(Rad1), K_Fbeta(Rad2)
        _print("K_Fβ =", self.K_Fbeta)
        
        # Abschnitt 4.5
        if self.geometrie.beta == 0:
            self.Z_epsilon = m.sqrt((4. - self.geometrie.epsilon_alpha) / 3.)
        elif self.geometrie.epsilon_beta >= 1.:
            self.Z_epsilon = m.sqrt(1 / self.geometrie.epsilon_alpha)
        else:
            self.Z_epsilon = m.sqrt((4. - self.geometrie.epsilon_alpha) / 3. * (1. - self.geometrie.epsilon_beta) + self.geometrie.epsilon_beta / self.geometrie.epsilon_alpha)
        _print("Z_ε =", self.Z_epsilon)

        # Glg 5.09
        self.Y_epsilon = 0.25 + 0.75 / self.geometrie.epsilon_alpha * m.cos(m.radians(self.geometrie.beta_b))**2
        _print("Y_ε =", self.Y_epsilon)

        # Tabelle 3.3
        def K_Halpha_und_Falpha(idx : int):
            linienbelastung = self.F_t / self.geometrie.b * self.K_A
            art = self.werkstoff[idx].art
            qualität = min(self.verzahnungsqualität)
            Art = Werkstoff.Art

            if art in (Art.einsatzgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl):
                if self.geometrie.beta == 0:
                    if linienbelastung > 100:
                        match qualität:
                            case 6 | 7:
                                return 1., 1.
                            case 8:
                                return 1.1, 1.1
                            case 9:
                                return 1.2, 1.2
                            case 10 | 11 | 12:
                                K_H = 1 / self.Z_epsilon**2
                                K_F = 1 / self.Y_epsilon**2
                                assert K_H >= 1.2
                                assert K_F >= 1.2
                                return K_H, K_F
                    else:
                        if qualität <= 6:
                            K_H = 1 / self.Z_epsilon**2
                            K_F = 1 / self.Y_epsilon**2
                            assert K_H >= 1.2
                            assert K_F >= 1.2
                            return K_H, K_F
                else:
                    if linienbelastung > 100:
                        match qualität:
                            case 6:
                                return 1., 1.
                            case 7:
                                return 1.1, 1.1
                            case 8:
                                return 1.2, 1.2
                            case 9:
                                return 1.4, 1.4
                            case 10 | 11 | 12:
                                K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                                assert K >= 1.4
                                return K, K
                    else:
                        if qualität <= 6:
                            K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                            assert K >= 1.4
                            return K, K
            else:
                if self.geometrie.beta == 0:
                    if linienbelastung > 100:
                        match qualität:
                            case 6 | 7 | 8:
                                return 1., 1.
                            case 9:
                                return 1.1, 1.1
                            case 10:
                                return 1.2, 1.2
                            case 11 | 12:
                                K_H = 1 / self.Z_epsilon**2
                                K_F = 1 / self.Y_epsilon**2
                                assert K_H >= 1.2
                                assert K_F >= 1.2
                                return K_H, K_F
                    else:
                        if qualität <= 6:
                            K_H = 1 / self.Z_epsilon**2
                            K_F = 1 / self.Y_epsilon**2
                            assert K_H >= 1.2
                            assert K_F >= 1.2
                            return K_H, K_F
                else:
                    if linienbelastung > 100:
                        match qualität:
                            case 6 | 7:
                                return 1., 1.
                            case 8:
                                return 1.1, 1.1
                            case 9:
                                return 1.2, 1.2
                            case 10:
                                return 1.4, 1.4
                            case 11 | 12:
                                K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                                assert K >= 1.4
                                return K, K
                    else:
                        if qualität <= 6:
                            K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                            assert K >= 1.4
                            return K, K

            raise ValueError
        K_ritzel = K_Halpha_und_Falpha(Rad1)
        K_rad = K_Halpha_und_Falpha(Rad2)
        if self.K_Halpha == None:
            self.K_Halpha = K_ritzel[0], K_rad[0]
        _print("K_Hα =", self.K_Halpha)
        if self.K_Falpha == None:
            self.K_Falpha = K_ritzel[1], K_rad[1]
        _print("K_Fα =", self.K_Falpha)

        # Grübchentragfähigkeit
        
        # Glg 4.12
        self.M_1 = m.tan(m.radians(self.geometrie.alpha_wt)) / m.sqrt(
            (m.sqrt(self.geometrie.d_a[Rad1]**2 / self.geometrie.d_b[Rad1]**2 - 1) - 2 * m.pi / self.geometrie.z[Rad1]) *
            (m.sqrt(self.geometrie.d_a[Rad2]**2 / self.geometrie.d_b[Rad2]**2 - 1) - (self.geometrie.epsilon_alpha - 1) * 2 * m.pi / self.geometrie.z[Rad2]))
        _print("M_1 =", self.M_1)

        # Glg 4.13
        self.M_2 = m.tan(m.radians(self.geometrie.alpha_wt)) / m.sqrt(
            (m.sqrt(self.geometrie.d_a[Rad2]**2 / self.geometrie.d_b[Rad2]**2 - 1) - 2 * m.pi / self.geometrie.z[Rad2]) *
            (m.sqrt(self.geometrie.d_a[Rad1]**2 / self.geometrie.d_b[Rad1]**2 - 1) - (self.geometrie.epsilon_alpha - 1) * 2 * m.pi / self.geometrie.z[Rad1]))
        _print("M_2 =", self.M_2)

        # Abschnitt 4.2
        if self.geometrie.beta == 0:
            self.Z_B = max(1., self.M_1)
            self.Z_D = max(1., self.M_2)
        elif self.geometrie.epsilon_beta >= 1:
            self.Z_B = 1.
            self.Z_D = 1.
        else:
            self.Z_B = max(1., self.M_1 - self.geometrie.epsilon_beta * (self.M_1 - 1))
            self.Z_D = max(1., self.M_2 - self.geometrie.epsilon_beta * (self.M_2 - 1))
        if self.innenverzahnt:
            self.Z_D = 1.
        _print("Z_B =", self.Z_B)
        _print("Z_D =", self.Z_D)

        # Glg. 4.14
        self.Z_H = m.sqrt(2 * m.cos(m.radians(self.geometrie.beta_b)) * m.cos(m.radians(self.geometrie.alpha_wt))
                          / m.cos(m.radians(self.geometrie.alpha_t))**2 / m.sin(m.radians(self.geometrie.alpha_wt)))
        _print("Z_H =", self.Z_H)

        # Tabelle 4.1
        ws1, ws2 = self.werkstoff
        Art = Werkstoff.Art
        stahl = (Art.Baustahl, Art.vergüteterStahl, Art.einsatzgehärteterStahl, Art.randschichtgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl)
        if ws1.art in stahl and ws2.art in stahl:
            self.Z_E = 189.8
        else:
            raise NotImplementedError
        _print("Z_E =", self.Z_E)

        # Glg 4.18
        self.Z_beta = m.sqrt(m.cos(m.radians(self.geometrie.beta)))
        _print("Z_β =", self.Z_beta)

        # Glg 4.20
        self.R_z100 = sum(self.R_z) / 2 * m.pow(100 / self.geometrie.a_w, 1/3)
        _print("R_z100 =", self.R_z100)

        # Glg 4.19
        self.Z_LVRstat = 1.
        if self.Z_LVRdyn == None:
            if self.fertigungsverfahren == (Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt, Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt):
                self.Z_LVRdyn = 0.85
                self.Z_LVRdyn
            elif self.fertigungsverfahren == (Fertigungsverfahren.geläpptGeschliffenGeschabt, Fertigungsverfahren.geläpptGeschliffenGeschabt):
                if self.R_z100 > 4:
                    self.Z_LVRdyn = 0.92
                else:
                    self.Z_LVRdyn = 1.
            elif Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt in self.fertigungsverfahren and Fertigungsverfahren.geläpptGeschliffenGeschabt in self.fertigungsverfahren:
                assert self.R_z100 <= 4
                self.Z_LVRdyn = 0.92
            else:
                raise ValueError(f"Unknown fertigungsverfahren {self.fertigungsverfahren}")
        _print("Z_LVRstat =", self.Z_LVRstat)
        _print("Z_LVRdyn =", self.Z_LVRdyn)

        # Glg 4.23
        wsart1 = self.werkstoff[Rad1].art
        wsart2 = self.werkstoff[Rad2].art
        Art = Werkstoff.Art
        weich = (Art.Baustahl, Art.vergüteterStahl)
        hart = (Art.einsatzgehärteterStahl, Art.randschichtgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl)
        if (wsart1 in weich or wsart2 in weich) and (wsart1 in hart or wsart2 in hart) and self.R_z <= 6:
            HB = min(self.werkstoff[Rad1].HB, self.werkstoff[Rad2].HB)
            if HB < 130:
                self.Z_W = 1.2
            elif HB > 470:
                self.Z_W = 1.
            else:
                self.Z_W = 1.2 - (HB - 130) / 1700
        else:
            self.Z_W = 1.
        _print("Z_W =", self.Z_W)

        # Tabelle 4.2
        def Z_Xdyn(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, Art.vergüteterStahl):
                return 1.
            elif wsart in (Art.einsatzgehärteterStahl, ):
                if self.geometrie.m_n <= 10:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.05 - 0.005  * self.geometrie.m_n
                else:
                    return 0.9
            elif wsart in (Art.nitrierterStahl, Art.nitrokarburierterStahl):
                if self.geometrie.m_n <= 7.5:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.08 - 0.011  * self.geometrie.m_n
                else:
                    return 0.75
            raise NotImplementedError
        self.Z_Xstat = 1., 1.
        self.Z_Xdyn = Z_Xdyn(Rad1), Z_Xdyn(Rad2)
        _print("Z_Xstat =", self.Z_Xstat)
        _print("Z_Xdyn =", self.Z_Xdyn)

        # Tabelle 4.3
        def Z_NT(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, ):
                return 1.6, 1.0
            elif wsart in (Art.vergüteterStahl, Art.einsatzgehärteterStahl):
                return 1.6, 1.0
            elif wsart in (Art.nitrierterStahl, ):
                return 1.3, 1.0
            elif wsart in (Art.nitrokarburierterStahl, ):
                return 1.1, 1.0
            raise NotImplementedError
        Z_NTritzel = Z_NT(Rad1)
        Z_NTrad = Z_NT(Rad2)
        self.Z_NTstat = Z_NTritzel[0], Z_NTrad[0]
        self.Z_NTdyn = Z_NTritzel[1], Z_NTrad[1]
        _print("Z_NTstat =", self.Z_NTstat)
        _print("Z_NTdyn =", self.Z_NTdyn)

        # Glg 4.03
        def sigma_HGstat(idx):
            return self.werkstoff[idx].sigma_Hlim * self.Z_NTstat[idx] * self.Z_LVRstat * self.Z_W * self.Z_Xstat[idx]
        def sigma_HGdyn(idx) -> float:
            return self.werkstoff[idx].sigma_Hlim * self.Z_NTdyn[idx] * self.Z_LVRdyn * self.Z_W * self.Z_Xdyn[idx]
        self.sigma_HGstat = sigma_HGstat(Rad1), sigma_HGstat(Rad2)
        self.sigma_HGdyn = sigma_HGdyn(Rad1), sigma_HGdyn(Rad2)
        _print("σ_HGstat =", self.sigma_HGstat)
        _print("σ_HGdyn =", self.sigma_HGdyn)

        # Glg 4.02
        self.sigma_H0 = self.Z_H * self.Z_E * self.Z_epsilon * self.Z_beta * m.sqrt(self.F_t / self.geometrie.d[Rad1] / self.geometrie.b * (self.geometrie.u + 1) / self.geometrie.u)
        _print("σ_H0 =", self.sigma_H0)

        # Glg 4.01
        def sigma_Hstat(idx):
            return (self.Z_B if idx == Rad1 else self.Z_D) * self.sigma_H0 * m.sqrt(self.K_S * self.K_V[idx] * self.K_Hbeta[idx] * self.K_Halpha[idx])
        def sigma_Hdyn(idx):
            return (self.Z_B if idx == Rad1 else self.Z_D) * self.sigma_H0 * m.sqrt(self.K_A * self.K_V[idx] * self.K_Hbeta[idx] * self.K_Halpha[idx])
        self.sigma_Hstat = sigma_Hstat(Rad1), sigma_Hstat(Rad2)
        self.sigma_Hdyn = sigma_Hdyn(Rad1), sigma_Hdyn(Rad2)
        _print("σ_Hstat =", self.sigma_Hstat)
        _print("σ_Hdyn =", self.sigma_Hdyn)

        # Glg 4.11
        def S_Hstat(idx):
            return self.sigma_HGstat[idx] / self.sigma_Hstat[idx]
        def S_Hdyn(idx):
            return self.sigma_HGdyn[idx] / self.sigma_Hdyn[idx]
        self.S_Hstat = S_Hstat(Rad1), S_Hstat(Rad2)
        self.S_Hdyn = S_Hdyn(Rad1), S_Hdyn(Rad2)

        check_value(Rad1, self.S_Hstat, S_Hstatmin, "S_Hstat", "statische Grübchensicherheit")
        check_value(Rad1, self.S_Hdyn, S_Hdynmin, "S_Hdyn", "dynamische Grübchensicherheit")
        check_value(Rad2, self.S_Hstat, S_Hstatmin, "S_Hstat", "statische Grübchensicherheit")
        check_value(Rad2, self.S_Hdyn, S_Hdynmin, "S_Hdyn", "dynamische Grübchensicherheit")
        _print()

        # Zahnfußtragfähigkeit

        # Glg D.1.01
        def z_n(idx):
            return self.geometrie.z[idx] / m.cos(m.radians(self.geometrie.beta_b))**2 / m.cos(m.radians(self.geometrie.beta))
        self.z_n = z_n(Rad1), z_n(Rad2)
        _print("z_n =", self.z_n)

        if not self.innenverzahnt:
            # Glg D.5.01
            self.E = (m.pi / 4 * self.geometrie.m_n - self.geometrie.h_fP * m.tan(m.radians(self.geometrie.alpha_n)) + self.s_pr / m.cos(m.radians(self.geometrie.alpha_n))
                      - (1 - m.sin(m.radians(self.geometrie.alpha_n))) * self.geometrie.rho_fP / m.cos(m.radians(self.geometrie.alpha_n))) 
            _print("E =", self.E)

            # Glg D.5.02
            def G(idx):
                return self.geometrie.rho_fP / self.geometrie.m_n - self.geometrie.h_fP / self.geometrie.m_n + self.geometrie.x[idx]
            self.G = G(Rad1), G(Rad2)
            _print("G =", self.G)

            # Glg D.5.03
            def H(idx):
                return 2 / self.z_n[idx] * (m.pi / 2 - self.E / self.geometrie.m_n) - m.pi / 3
            self.H = H(Rad1), H(Rad2)
            _print("H =", self.H)

            # Glg D.5.04
            def theta(idx):
                theta = m.degrees(m.pi / 6)
                for i in range(5):
                    theta = m.degrees(2 * self.G[idx] / self.z_n[idx] * m.tan(m.radians(theta)) - self.H[idx])
                return theta
            self.theta = theta(Rad1), theta(Rad2)
            _print("ϑ =", self.theta)

            # Glg D.5.05
            def s_Fn(idx):
                return self.geometrie.m_n * (self.z_n[idx] * m.sin(m.pi / 3 - m.radians(self.theta[idx])) + m.sqrt(3) * (self.G[idx] / m.cos(m.radians(self.theta[idx]))
                                                                                                                         - self.geometrie.rho_fP / self.geometrie.m_n))
            self.s_Fn = s_Fn(Rad1), s_Fn(Rad2)
            _print("s_Fn =", self.s_Fn)

            # Glg D.5.06
            def d_n(idx):
                return self.geometrie.m_n * self.z_n[idx]
            self.d_n = d_n(Rad1), d_n(Rad2)
            _print("d_n =", self.d_n)

            # Glg D.5.07
            def d_bn(idx):
                return self.d_n[idx] * m.cos(m.radians(self.geometrie.alpha_n))
            self.d_bn = d_bn(Rad1), d_bn(Rad2)
            _print("d_bn =", self.d_bn)

            # Glg D.5.08
            def d_an(idx):
                return self.d_n[idx] + self.geometrie.d_a[idx] - self.geometrie.d[idx]
            self.d_an = d_an(Rad1), d_an(Rad2)
            _print("d_an =", self.d_an)

            # Glg D.5.09
            def alpha_an(idx):
                return m.degrees(m.acos(self.d_bn[idx] / self.d_an[idx]))
            self.alpha_an = alpha_an(Rad1), alpha_an(Rad2)
            _print("α_an =", self.alpha_an)

            # Glg D.5.10
            def y_a(idx):
                return 1 / self.z_n[idx] * (m.pi / 2 + 2 * self.geometrie.x[idx] * m.tan(m.radians(self.geometrie.alpha_n))) + involute(self.geometrie.alpha_n) - involute(self.alpha_an[idx])
            self.y_a = y_a(Rad1), y_a(Rad2)
            _print("y_a =", self.y_a)

            # Glg D.5.11
            def alpha_Fan(idx):
                return self.alpha_an[idx] - m.degrees(self.y_a[idx])
            self.alpha_Fan = alpha_Fan(Rad1), alpha_Fan(Rad2)
            _print("α_Fan =", self.alpha_Fan)

            # Glg D.5.12
            def h_Fa(idx):
                return self.geometrie.m_n * (0.5 * self.z_n[idx] * (m.cos(m.radians(self.geometrie.alpha_n)) / m.cos(m.radians(self.alpha_Fan[idx])) - m.cos(m.pi / 3 - m.radians(self.theta[idx])))
                                   + 0.5 * (self.geometrie.rho_fP / self.geometrie.m_n - self.G[idx] / m.cos(m.radians(self.theta[idx]))))
            self.h_Fa = h_Fa(Rad1), h_Fa(Rad2)
            _print("h_Fa =", self.h_Fa)

            # Glg D.5.13
            def rho_F(idx):
                return self.geometrie.rho_fP + self.geometrie.m_n * 2 * self.G[idx]**2 / m.cos(m.radians(self.theta[idx])) / (self.z_n[idx] * m.cos(m.radians(self.theta[idx]))**2 - 2 * self.G[idx])
            self.rho_F = rho_F(Rad1), rho_F(Rad2)
            _print("ρ_F =", self.rho_F)
        else:
            raise NotImplementedError("Innenverzahnte Räder sind noch nicht implementiert")

        # Glg D.3.01
        def Y_Fa(idx):
            return (6 * self.h_Fa[idx] / self.geometrie.m_n * m.cos(m.radians(self.alpha_Fan[idx]))) / ((self.s_Fn[idx] / self.geometrie.m_n)**2 * m.cos(m.radians(self.geometrie.alpha_n)))
        self.Y_Fa = Y_Fa(Rad1), Y_Fa(Rad2)
        _print("Y_Fa =", self.Y_Fa)

        # Glg D.4.02
        def L_a(idx):
            return self.s_Fn[idx] / self.h_Fa[idx]
        self.L_a = L_a(Rad1), L_a(Rad2)
        _print("L_a =", self.L_a)

        # Glg D.4.03
        def q_s(idx):
            return self.s_Fn[idx] / 2 / self.rho_F[idx]
        self.q_s = q_s(Rad1), q_s(Rad2)
        _print("q_s =", self.q_s)
        assert all(1 <= q_s < 8 for q_s in self.q_s)

        # Glg D.4.01
        def Y_Sa(idx):
            return (1.2 + 0.13 * self.L_a[idx]) * m.pow(self.q_s[idx], 1 / (1.21 + 2.3 / self.L_a[idx]))
        self.Y_Sa = Y_Sa(Rad1), Y_Sa(Rad2)
        _print("Y_Sa =", self.Y_Sa)

        # Glg 5.08
        def Y_FS(idx):
            return self.Y_Sa[idx] * self.Y_Fa[idx]
        self.Y_FS = Y_FS(Rad1), Y_FS(Rad2)
        _print("Y_FS =", self.Y_FS)

        # Glg 5.10
        self.Y_beta = 1 - min(self.geometrie.epsilon_beta, 1.) * min(self.geometrie.beta, 30) / 120
        _print("Y_β =", self.Y_beta)
 
        # Tabelle 5.8
        epsilon_alphan = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
        _print("ε_αn =", epsilon_alphan)
        def Y_S(idx):
            assert 1 <= self.s_Fn[idx] / self.h_Fa[idx] <= 1.2
            return self.Y_Sa[idx] * (0.6 + 0.4 * epsilon_alphan)
        self.Y_S = Y_S(Rad1), Y_S(Rad2)
        _print("Y_S =", self.Y_S)

        # Abschnitt 5.6
        def Y_deltarelTstat(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, ):
                raise NotImplementedError
            elif wsart in (Art.vergüteterStahl, ):
                raise NotImplementedError
            elif wsart in (Art.einsatzgehärteterStahl, ):
                return 0.44 * self.Y_S[idx] + 0.12
            elif wsart in (Art.nitrierterStahl, Art.nitrokarburierterStahl ):
                return 0.20 * self.Y_S[idx] + 0.60
            raise NotImplementedError
        # Glg 5.11, 5.12
        def Y_deltarelTdyn(idx):
            if self.q_s[idx] >= 1.5:
                return 1.
            else:
                return 0.95
        self.Y_deltarelTstat = Y_deltarelTstat(Rad1), Y_deltarelTstat(Rad2)
        self.Y_deltarelTdyn = Y_deltarelTdyn(Rad1), Y_deltarelTdyn(Rad2)
        _print("Y_δrelTstat =", self.Y_deltarelTstat)
        _print("Y_δrelTdyn =", self.Y_deltarelTdyn)

        # Abschnitt 5.7
        self.Y_RrelTstat = 1.
        def Y_RrelTdyn(idx):
            if self.R_z[idx] <= 16:
                return 1.
            else:
                return 0.9
        self.Y_RrelTdyn = Y_RrelTdyn(Rad1), Y_RrelTdyn(Rad2)
        _print("Y_RrelTstat =", self.Y_RrelTstat)
        _print("Y_RrelTdyn =", self.Y_RrelTdyn)

        # Tabelle 5.1
        def Y_Xdyn(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, Art.vergüteterStahl):
                if self.geometrie.m_n <= 5:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.03 - 0.006  * self.geometrie.m_n
                else:
                    return 0.85
            elif wsart in (Art.einsatzgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl):
                if self.geometrie.m_n <= 5:
                    return 1.
                elif self.geometrie.m_n < 25:
                    return 1.05 - 0.01  * self.geometrie.m_n
                else:
                    return 0.8
            raise NotImplementedError
        self.Y_Xstat = 1., 1.
        self.Y_Xdyn = Y_Xdyn(Rad1), Y_Xdyn(Rad2)
        _print("Y_Xstat =", self.Y_Xstat)
        _print("Y_Xdyn =", self.Y_Xdyn)

        # Tabelle 5.2
        def Y_NT(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, Art.vergüteterStahl):
                return 2.5, 1.0
            elif wsart in (Art.einsatzgehärteterStahl, ):
                return 2.5, 1.0
            elif wsart in (Art.nitrierterStahl, ):
                return 1.6, 1.0
            elif wsart in (Art.nitrokarburierterStahl, ):
                return 1.1, 1.0
            raise NotImplementedError
        Y_NTritzel = Y_NT(Rad1)
        Y_NTrad = Y_NT(Rad2)
        self.Y_NTstat = Y_NTritzel[0], Y_NTrad[0]
        self.Y_NTdyn = Y_NTritzel[1], Y_NTrad[1]
        _print("Y_NTstat =", self.Y_NTstat)
        _print("Y_NTdyn =", self.Y_NTdyn)

        # Glg 5.03
        def sigma_FGstat(idx):
            return self.werkstoff[idx].sigma_FE * self.Y_NTstat[idx] * self.Y_deltarelTstat[idx] * self.Y_RrelTstat * self.Y_Xstat[idx]
        def sigma_FGdyn(idx):
            return self.werkstoff[idx].sigma_FE * self.Y_NTdyn[idx] * self.Y_deltarelTdyn[idx] * self.Y_RrelTdyn[idx] * self.Y_Xdyn[idx]
        self.sigma_FGstat = sigma_FGstat(Rad1), sigma_FGstat(Rad2)
        self.sigma_FGdyn = sigma_FGdyn(Rad1), sigma_FGdyn(Rad2)
        _print("σ_FGstat =", self.sigma_FGstat)
        _print("σ_FGdyn =", self.sigma_FGdyn)

        # Glg 5.02
        def sigma_F0(idx):
            return self.F_t / self.geometrie.b / self.geometrie.m_n * self.Y_FS[idx] * self.Y_epsilon * self.Y_beta
        self.sigma_F0 = sigma_F0(Rad1), sigma_F0(Rad2)
        _print("σ_F0 =", self.sigma_F0)

        # Glg 5.01
        def sigma_Fstat(idx):
            return self.sigma_F0[idx] * self.K_S * self.K_V[idx] * self.K_Fbeta[idx] * self.K_Falpha[idx]
        def sigma_Fdyn(idx):
            return self.sigma_F0[idx] * self.K_A * self.K_V[idx] * self.K_Fbeta[idx] * self.K_Falpha[idx]
        self.sigma_Fstat = sigma_Fstat(Rad1), sigma_Fstat(Rad2)
        self.sigma_Fdyn = sigma_Fdyn(Rad1), sigma_Fdyn(Rad2)
        _print("σ_Fstat =", self.sigma_Fstat)
        _print("σ_Fdyn =", self.sigma_Fdyn)

        # Glg 5.07
        def S_Fstat(idx):
            return self.sigma_FGstat[idx] / self.sigma_Fstat[idx]
        def S_Fdyn(idx):
            return self.sigma_FGdyn[idx] / self.sigma_Fdyn[idx]
        self.S_Fstat = S_Fstat(Rad1), S_Fstat(Rad2)
        self.S_Fdyn = S_Fdyn(Rad1), S_Fdyn(Rad2)

        check_value(Rad1, self.S_Fstat, S_Fstatmin, "S_Fstat", "statische Zahnbruchsicherheit")
        check_value(Rad1, self.S_Fdyn, S_Fdynmin, "S_Fdyn", "dynamische Zahnbruchsicherheit")
        check_value(Rad2, self.S_Fstat, S_Fstatmin, "S_Fstat", "statische Zahnbruchsicherheit")
        check_value(Rad2, self.S_Fdyn, S_Fdynmin, "S_Fdyn", "dynamische Zahnbruchsicherheit")
        _print()
        return


