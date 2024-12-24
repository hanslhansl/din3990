from ._din3990_11 import *
from . import _din3990_11


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
  
class DIN_3990_11:
    def __init__(self,
                geometrie : DIN_21771,
                P : float,
                n_1 : float,
                verzahnungsqualität : tuple[int, int],
                werkstoff : tuple[din3990_5.Werkstoff, din3990_5.Werkstoff],
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
        - _assert if True, safety is asserted which additionally requires the same arguments as passed to Getriebe.is_safe()

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

        self.v = v(self.n[Ritzel], self.geometrie.d[Ritzel])
        _print("v =", self.v)

        self.F_t = F_t(self.P, self.n[Ritzel], self.geometrie.d[Ritzel])
        _print("F_t =", self.F_t)
  
        self.T = T(self.P, self.n[Ritzel]), T(self.P, self.n[Rad])
        _print("T =", self.T)

        if self.K_V == None:
            self.K_V = _din3990_11.K_V(self.z[Ritzel], self.v, self.geometrie.u, self.F_t, self.K_A, self.geometrie.b, self.geometrie.epsilon_beta,
                                       self.geometrie.beta == 0, self.verzahnungsqualität, _print)
        _print("K_V =", self.K_V)

        if self.K_Hbeta == None:
            # Abschnitt 3.4.1
            assert(self.F_t / self.geometrie.b * self.K_A >= 100)

            # Glg 3.07
            def F_m(idx) -> float:
                return self.F_t * self.K_A * self.K_V[idx]
            self.F_m = F_m(Ritzel), F_m(Rad)
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
                    return 0. if self.s[idx] == 0 else K_s(idx) * self.l[idx] * self.s[idx] / self.geometrie.d[Ritzel]**2 * (self.geometrie.d[Ritzel] / self.d_sh[idx])**4
                if self.A == Tabelle_3_2.ohneBreitenballigkeltOderEndrücknahme:
                    A = 0.023
                elif self.A == Tabelle_3_2.mitSinnvollerBreitenballigkeit:
                    A = 0.012
                elif self.A == Tabelle_3_2.mitSinnvollerEndrücknahme:
                    A = 0.016

                if not self.doppelschrägverzahnt:
                    return self.F_m[idx] / self.geometrie.b * A * (abs(1 + temp(idx) - 0.3) + 0.3) * (self.geometrie.b / self.geometrie.d[Ritzel])**2
                else:
                    return self.F_m[idx] / self.geometrie.b * 2 * A * (abs(1.5 + temp(idx) - 0.3) + 0.3) * (self.b_B / self.geometrie.d[Ritzel])**2
            self.f_sh = f_sh(Ritzel), f_sh(Rad)
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
                            return 1 if abs(self.K_s) * self.l * self.s / self.d[Ritzel]**2 * (self.d[Ritzel] / self.d_sh[Ritzel])**4 <= B_s else -1
                        case Bild_3_1.d:
                            B_s = 1.5 if self.doppelschrägverzahnt else 1
                            return 1 if abs(self.K_s) * self.l * self.s / self.d[Ritzel]**2 * (self.d[Ritzel] / self.d_sh[Ritzel])**4 >= B_s - 0.3 else -1
                        case _:
                            raise ValueError("Unknown bild_3_1 option ", self.bild_3_1[idx])
                if self.f_ma[idx] == 0:
                    return abs(1.33 * self.f_sh[idx])
                else:
                    return abs(1.33 * self.f_sh[idx] + multi(idx) * self.f_ma[idx])
            self.F_betax = F_betax(Ritzel), F_betax(Rad)
            _print("F_betax =", self.F_betax)

            # Abschnitt 3.4.2.6
            def y_beta(idx : int):
                werkstoff = self.werkstoff[idx]
                match werkstoff.art:
                    case din3990_5.Werkstoff.Art.Stahl | din3990_5.Werkstoff.Art.Vergütungsstahl | din3990_5.Werkstoff.Art.PerlitischesGußeisenMitKugelgraphit | din3990_5.Werkstoff.Art.BainitischesGußeisenMitKugelgraphit:
                        y_beta = 320 / werkstoff.sigma_Hlim * self.F_betax[idx]
                        if self.v <= 5:
                            pass
                        elif self.v <= 10:
                            assert y_beta <= 25600 / werkstoff.sigma_Hlim
                        else:
                            assert y_beta <= 12800 / werkstoff.sigma_Hlim
                        return y_beta
                    case din3990_5.Werkstoff.Art.Grauguß | din3990_5.Werkstoff.Art.FerritischesGußeisenMitKugelgraphit:
                        y_beta = 0.55 * self.F_betax[idx]
                        if self.v <= 5:
                            pass
                        elif self.v <= 10:
                            assert y_beta <= 45
                        else:
                            assert y_beta <= 22
                        return y_beta
                    case din3990_5.Werkstoff.Art.Einsatzstahl | din3990_5.Werkstoff.Art.nitrierterStahl | din3990_5.Werkstoff.Art.nitrokarburierterStahl:
                        y_beta = 0.15 * self.F_betax[idx]
                        assert y_beta <= 6
                        return y_beta
                raise NotImplementedError
            self.y_beta = (y_beta(Ritzel) + y_beta(Rad)) / 2
            _print("y_beta =", self.y_beta)

            # Glg 3.08
            def F_betay(idx):
                return self.F_betax[idx] - self.y_beta
            self.F_betay = F_betay(Ritzel), F_betay(Rad)
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
            self.K_Hbeta = K_Hbeta(Ritzel), K_Hbeta(Rad)
        _print("K_Hβ =", self.K_Hbeta)
        
        # Glg 3.22
        if self.K_Fbeta == None:
            def K_Fbeta(idx):
                h_b = min(self.geometrie.h / self.geometrie.b, 1. / 3.)
                return m.pow(self.K_Hbeta[idx], (1 / (1 + h_b + h_b**2)))

            assert(self.F_t / self.geometrie.b * self.K_A >= 100)   # Abschnitt 3.4.1
            self.K_Fbeta = K_Fbeta(Ritzel), K_Fbeta(Rad)
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

            if art in (din3990_5.Werkstoff.Art.einsatzgehärteterStahl, din3990_5.Werkstoff.Art.nitrierterStahl, din3990_5.Werkstoff.Art.nitrokarburierterStahl):
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
        K_ritzel = K_Halpha_und_Falpha(Ritzel)
        K_rad = K_Halpha_und_Falpha(Rad)
        if self.K_Halpha == None:
            self.K_Halpha = K_ritzel[0], K_rad[0]
        _print("K_Hα =", self.K_Halpha)
        if self.K_Falpha == None:
            self.K_Falpha = K_ritzel[1], K_rad[1]
        _print("K_Fα =", self.K_Falpha)

        # Grübchentragfähigkeit
        
        # Glg 4.12
        self.M_1 = m.tan(m.radians(self.geometrie.alpha_wt)) / m.sqrt(
            (m.sqrt(self.geometrie.d_a[Ritzel]**2 / self.geometrie.d_b[Ritzel]**2 - 1) - 2 * m.pi / self.geometrie.z[Ritzel]) *
            (m.sqrt(self.geometrie.d_a[Rad]**2 / self.geometrie.d_b[Rad]**2 - 1) - (self.geometrie.epsilon_alpha - 1) * 2 * m.pi / self.geometrie.z[Rad]))
        _print("M_1 =", self.M_1)

        # Glg 4.13
        self.M_2 = m.tan(m.radians(self.geometrie.alpha_wt)) / m.sqrt(
            (m.sqrt(self.geometrie.d_a[Rad]**2 / self.geometrie.d_b[Rad]**2 - 1) - 2 * m.pi / self.geometrie.z[Rad]) *
            (m.sqrt(self.geometrie.d_a[Ritzel]**2 / self.geometrie.d_b[Ritzel]**2 - 1) - (self.geometrie.epsilon_alpha - 1) * 2 * m.pi / self.geometrie.z[Ritzel]))
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
        stahl = (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.vergüteterStahl, din3990_5.Werkstoff.Art.einsatzgehärteterStahl, din3990_5.Werkstoff.Art.randschichtgehärteterStahl, din3990_5.Werkstoff.Art.nitrierterStahl, din3990_5.Werkstoff.Art.nitrokarburierterStahl)
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
        wsart1 = self.werkstoff[Ritzel].art
        wsart2 = self.werkstoff[Rad].art
        weich = (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.vergüteterStahl)
        hart = (din3990_5.Werkstoff.Art.einsatzgehärteterStahl, din3990_5.Werkstoff.Art.randschichtgehärteterStahl, din3990_5.Werkstoff.Art.nitrierterStahl, din3990_5.Werkstoff.Art.nitrokarburierterStahl)
        if (wsart1 in weich or wsart2 in weich) and (wsart1 in hart or wsart2 in hart) and self.R_z <= 6:
            HB = min(self.werkstoff[Ritzel].HB, self.werkstoff[Rad].HB)
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
            if wsart in (din3990_5.Werkstoff.Art.Baustahl, din3990_5.Werkstoff.Art.vergüteterStahl):
                return 1.
            elif wsart in (din3990_5.Werkstoff.Art.einsatzgehärteterStahl, ):
                if self.geometrie.m_n <= 10:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.05 - 0.005  * self.geometrie.m_n
                else:
                    return 0.9
            elif wsart in (din3990_5.Werkstoff.Art.nitrierterStahl, din3990_5.Werkstoff.Art.nitrokarburierterStahl):
                if self.geometrie.m_n <= 7.5:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.08 - 0.011  * self.geometrie.m_n
                else:
                    return 0.75
            raise NotImplementedError
        self.Z_Xstat = 1., 1.
        self.Z_Xdyn = Z_Xdyn(Ritzel), Z_Xdyn(Rad)
        _print("Z_Xstat =", self.Z_Xstat)
        _print("Z_Xdyn =", self.Z_Xdyn)

        # Tabelle 4.3
        def Z_NT(idx : int):
            wsart = self.werkstoff[idx].art
            if wsart in (din3990_5.Werkstoff.Art.Baustahl, ):
                return 1.6, 1.0
            elif wsart in (din3990_5.Werkstoff.Art.vergüteterStahl, din3990_5.Werkstoff.Art.einsatzgehärteterStahl):
                return 1.6, 1.0
            elif wsart in (din3990_5.Werkstoff.Art.nitrierterStahl, ):
                return 1.3, 1.0
            elif wsart in (din3990_5.Werkstoff.Art.nitrokarburierterStahl, ):
                return 1.1, 1.0
            raise NotImplementedError
        Z_NTritzel = Z_NT(Ritzel)
        Z_NTrad = Z_NT(Rad)
        self.Z_NTstat = Z_NTritzel[0], Z_NTrad[0]
        self.Z_NTdyn = Z_NTritzel[1], Z_NTrad[1]
        _print("Z_NTstat =", self.Z_NTstat)
        _print("Z_NTdyn =", self.Z_NTdyn)

        # Glg 4.03
        def sigma_HGstat(idx):
            return self.werkstoff[idx].sigma_Hlim * self.Z_NTstat[idx] * self.Z_LVRstat * self.Z_W * self.Z_Xstat[idx]
        def sigma_HGdyn(idx) -> float:
            return self.werkstoff[idx].sigma_Hlim * self.Z_NTdyn[idx] * self.Z_LVRdyn * self.Z_W * self.Z_Xdyn[idx]
        self.sigma_HGstat = sigma_HGstat(Ritzel), sigma_HGstat(Rad)
        self.sigma_HGdyn = sigma_HGdyn(Ritzel), sigma_HGdyn(Rad)
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
        epsilon_alphan = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
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


