from enum import Enum
from math import sqrt
from dataclasses import dataclass, field
import math

class AnchorInstalledType(Enum):
    Cast_in = 1
    Post_in = 2

class CastInAnchorType(Enum):
    # ACI318-19 Fig.R2.1
    HexHeadBoltWithWasher = 1
    L_bolt = 2
    J_bolt = 3
    WeldedHeadedStud = 4

class PostInAnchorType(Enum):
    # ACI318-19 Fig.R2.1
    Adhesive = 1
    Andercut = 2
    SleeveType = 3 #TorqueControlledExpansion
    StudType = 4   #TorqueControlledExpansion
    DisplacementControlledExpansion = 5
    Screw = 6

@dataclass
class Anchor:
    pass

class BasePlate:
    pass



@dataclass
class AnchorConnection:
    """
    d       : Kolon kesit boyu
    bf      : Kolon kesit eni
    B       : Taban plakasi genisliği
    N       : Taba plakasi yuksekligi
    tf      : Taban plakasi kalinligi
    P_axial : Gelen eksenel yük
    M       : Gelen moment
    fck     : Beton karakteristik dayanimi

    """
    d       : float = field(default_factory=float)
    bf      : float = field(default_factory=float)
    B       : int = field(default_factory=int)
    N       : int = field(default_factory=int)
    tf      : int = field(default_factory=int)
    P_axial : float = field(default_factory=float)
    M       : float = field(default_factory=float)
    fck     : float = field(default_factory=float)

@dataclass
class BasePlate:
    P_u : float
    M_u : float
    V_u : float
    f_c : float
    d   : float
    b_f : float
    f_c : float
    F_y : float
    Case: int

    """
    Taban plakasının genişlik ve yüksekliğinin aynı olması uygulama, malzeme kesim ve temini açısından büyük avantajdır bu nedenle aynı olduğu kabul edilecek. Rijitleştirme levhaları ilk aşamada hesaplarda göz önüne alınmayacak.

    TASARIM ADIMLARI
    1 - Analiz sonucundan LRFD kombinasyonları ile arttırılmış yüklerden gelen Pu,Mu,Vu değerleri
    2 - Eksenel basınç kuvvetinde minimum gerekli taban plakası alanının bulunması
    3 - Taban plakası başlangıç boyutlarının tespiti B==N ==> max(bf+80 , d+80)
    4 - Taban plakası alanı ile minimum gerekli alanın karşılaştırılması.
    5 - Yük dış merkezliği e=Mu/Pu; f_pmax = fi * 0.85 * f_c * (A2/A1)**0.5, q_max= f_pmax * B ve e_crit= N/2 - (P_u / (2*q_max)) bulunması
    6 - m=(N-0.95*d)/2, n=(B-0.8*b_f)/2, lamda=(2*X**0.5) / (1 + (1-X)**0.5) ve bunlara bağlı olarak konsol plaka boyu l=max(m,n,lamb*n') değerinin hesaplanması.
    7 - Çekme için rod gerekiyorsa, ankraj merkezinden plaka orta noktasına olan mesafenin hesabı f,
    8 - Basınç alanı uzunluğu Y 
    9 - Taban plakası altında oluşan basınç kuvveti q nun bulunması,
    11- Taşıma gücü kontrolü q_max = f_p(max) * B(or N) > q = Pu/Y
    12- Taban plakası altında oluşan basınç gerilmesinin bulunması f_p=P_u/(B*Y),
    13- Plaka kalınlığının bulunması tp,
    """
    import math

    def ApproximateBasePlateArea(P_u : float, f_c : float, Case : int = 3, fi : float = 0.65)-> float:
        """
            A1 : Taban plakasi alani
            A2 : Beton yüzey alani

            Case 1 : A1 = A2
            Case 2 : A2 >= 4A1
            Case 3 : A1 < A2 < 4A1

        """
        if Case == 1:
            A1_req = P_u / (fi * 0.85 * f_c)

        if Case == 2 or Case == 3:
            A1_req = P_u / (2 * fi * 0.85 * f_c)
        return round(A1_req,0)

    def FindPlateDimensions(d : float, b_f : float, A1_req : float)-> int:
        delta = (0.95 * d - 0.8 * b_f) / 2
        
        N_1 = A1_req**0.5 + delta
        N_2 = d + 2 * 80#inç - 3inç yerine 80mm olmalı
        N = max(N_1,N_2)

        B_1 = A1_req / N
        B_2 = b_f + 2 * 80#inç - 3inç yerine 80mm olmalı
        B = max(B_1,B_2)

        N = math.ceil(max(B,N))
        B = N
        A_baseplate = B * N
        if A1_req > A_baseplate:
            raise f"A1_req = {A1_req} > {A_baseplate} = A_baseplate"
        else:
            print("Uygulama açısından kare plaka tercih edilmiştir taban plakasının B ve N boyutları eşittir.")
        
        # while N%5 != 0:
        #         N = N+1
        #         B = B + 1
        return N

    def e_Get(M_u : float, P_u : float)-> float:
        return round(M_u/P_u,2)

    def f_pmax_Get(f_c : float, A1 : float, A2 : float, fi : float = 0.65)-> float:
        if A2 < A1:
            raise ValueError("A2, A1'den küçük olamaz!!!")
        if A2 >= A1:
            f_pmax = fi * 0.85 * f_c * (A2/A1)**0.5

        if (A2/A1)**0.5 >= 1 and (A2/A1)**0.5 <= 2:
            print("Beton sargılama katkısı kullanılabilir.")

        return round(f_pmax,2)

    # Türkiye çelik yapılar yönetmeliği denk. 13.22-23
    def Get_P_p(f_pmax : float, A1 : float, f_c : float)-> float:
        P_p = f_pmax * A1
        P_p = min(P_p, 1.7 * f_c * A1)
        return round(P_p,2)

    def q_max_Get(f_pmax : float, B : float)-> float:
        q_max = f_pmax * B
        return round(q_max,2)

    def e_crit_Get(q_max : float, P_u : float, N : float)-> float:
        e_crit = N/2 - (P_u / (2*q_max))
        return round(e_crit,2)

    def Get_m(N : float, d : float)->int:
        return round((N-0.95*d)/2,2)

    def Get_n(B : float, b_f : float)->int:
        return round((B-0.8*b_f)/2,2)

    def Get_X(d : float, b_f : float, P_u : float, P_p : float, fi : float = 0.9)->float:
        X = (4 * d * b_f * P_u) / ((d + b_f)**2 *  P_p)
        return round(X,2)

    def Get_lambda(X : float)->float:
        lambda_x = 1.0
        if X <1:
            lambda_x = (2*X**0.5) / (1 + (1-X)**0.5)
        lamb_ = min(1.0,lambda_x)
        return round(lamb_,2)

    def Get_l(d : float, b_f : float, m : float, n : float, lamb : float)-> float:
        n_lamb = lamb * (d * b_f)**0.5 / 4
        return max(m,n,n_lamb)

    def BasePlateThickness(P_u : float, l : float, B : float, N : float, F_y : float, fi : float = 0.9)->int:
        t_min = l * ((2 * P_u) / (fi * F_y * B * N))**0.5
        return round(t_min,2)

    def Get_Y(e : float, e_crit : float, P_u : float, N : float, f : float, q_max : float) -> float:
        if e <= e_crit:
            Y = N - 2*e
        
        if e > e_crit:
            limit = (f + N/2)**2
            check = (2 * P_u * (e + f)) / q_max
            if check > limit:
                print("Denge denklemi çözümsüz plaka büyütülmeli.")
                return 0.0
            Y1 = (f + N*0.5)  + ( (f + N * 0.5)**2 - (2*P_u * (e + f))/q_max )**0.5
            Y2 = (f + N*0.5)  - ( (f + N * 0.5)**2 - (2*P_u * (e + f))/q_max )**0.5
            Y = min(Y1,Y2)
        return Y

    def Get_f_p(P_u : float, B : float, Y : float) -> float:
        return round(P_u/(B*Y),2)

    def ForMomentsBasePlateThickness(m : float, n : float, f_p : float, Y : float, F_y : float)-> int:
        if Y >= m:
            t_p_req1 = 1.5 * m * (f_p / F_y)**0.5
        if Y < m :
            t_p_req1 = 2.11 * ( (f_p * Y * (m - Y*0.5) ) / F_y)
        
        if Y >= n:
            t_p_req2 = 1.5 * n * (f_p / F_y)**0.5
        if Y < n :
            t_p_req2 = 2.11 * ( (f_p * Y * (n - Y*0.5) ) / F_y)
        
        t_p_req = max(t_p_req2,t_p_req1)
        return round(t_p_req,2)

    def GetBasePlateArea(B : float, N : float) -> float:
        A1 = B * N
        return round(A1,2)

    def Get_t(N : int, x : float) -> float:
        """ankraj merkezinden plaka orta noktasına olan mesafeyi hesaplar

        Args:
            N (int): İlgili doğrultuda plaka boyutu
            x (float): Ankraj rodun merkezinin plaka kenarına olan mesafesi

        Returns:
            float: t
        """
        t = N/2 - x #ankraj merkezinden plaka orta noktasına olan mesafe
        return round(t,2)

    def Get_q(P_u : float, Y : float) -> float:
        return round(P_u/Y,2)






# Betona ankrajlanmış birleşimlerde limit dayanımlar

"""
Abrg        : Net bearing area of the head of stud, anchor bolt, or headed deformed bar,
ANa         : Projected influence area of a single adhesive anchor or group of adhesive anchors, for calculation of bond strength in tension, in.2
ANao        : Projected influence area of a single adhesive anchor, for calculation of bond strength in tension if not limited by edge distance or spacing, in.2
ANc         : Projected concrete failure area of a single anchor or group of anchors, for calculation of strength intension, in.2
ANco        : Projected concrete failure area of a single anchor, for calculation of strength in tension if not limited by edge distance or spacing,
h_ef'       : Limiting value of hef where anchors are located less than 1.5hef from three or more edges refer to ACI318-19 Fig. R17.6.2.1.2
h_ef_sl     : Effective embedment depth of shear lug
h_sl        : Embedment depth of shear lug
h_a         : Thickness of member in which an anchor is located, measured parallel to anchor axis
da          : outside diameter of anchor or shaft diameter of headed stud, headed bolt, or hooked bolt
da'         : Value substituted for d_a if an oversized anchor is used,


adhesive    : chemical components formulated organic polymers, or a combination of polymers and inorganic materials that cure if blended together.
anchor      : a steel element either cast into concrete or post-installed into a hardened concrete member and used to transmit applied loads to the concrete.
anchor, post-installed : anchor installed in hardened concrete; adhesive, expansion, screw, and undercut anchors are examples of post-installed anchors.
anchor, horizontal or upwardly inclined : Anchor installed in a hole drilled horizontally or in a hole drilled at any orientation above horizontal.
"""

# Tension

class SteelStrengthOfAnchorInTension:

    def NominalSteelStrengthOfAnchorInTension(self, A_se: float, f_ya : float,f_uta : float)-> float:
        """_summary_

        Args:
            A_se (float): is the effective cross-sectional area of an anchor in tension mm^2
            f_ya (float): specified yield strength of anchor steel, shall not exceed 1.9f_ya or 862MPa
            f_uta (float): specified tensile strength of anchor steel, shall not exceed 1.9f_ya or 862MPa

        Returns:
            float: N_sa
        """
        f_uta = min(f_uta,1.9*f_ya,862)
        return round(A_se * f_uta,3)

class ConcreteBreakoutStrengthOfAnchorInTension:
    
    def Psi_cpN_Get(self,Cac : float, Ca_min : float, h_ef: float, InstalledType : AnchorInstalledType) -> float:
        """ACI 318-19 17.6.2.6 Breakout splitting factor

        Args:
            Cac (float): critical edge distance required to develop the basic strength as controlled by concrete breakout or bond of a post-installed anchor in tension in uncracked concrete without supplementary reinforcement to control splitting,
            Ca_min (float): minimum distance from center of an anchor shaft to the edge of concrete
            h_ef (float): Effective embedment depth of anchor
            InstalledType (AnchorInstalledType): _description_

        Returns:
            float: _description_
        """
        # 
        Psi_cpN = 1.0
        if InstalledType.value != 1:
            if Ca_min >= Cac:
                Psi_cpN = 1.0
            if Ca_min < Cac:
                Psi_cpN = max((Ca_min/Cac),(1.5*h_ef/Cac))
        
        return Psi_cpN
        
    def Psi_cN_Get(self, InstalledType : AnchorInstalledType, IsConcCracked : bool) -> float:
        # ACI 318-19 17.6.2.5 Breakout cracking factor, 1.25 for cast-in anchors, 1.4 for post-installed anchors,
        Psi_cN = 1.25
        if InstalledType.value != 1:
            Psi_cN = 1.4
        if IsConcCracked:
            Psi_cN = 1.0
        return Psi_cN

    def Psi_edN_Get(self,Ca_min : float, h_ef : float) -> float:
        """ACI 318-19 17.6.2.4 breakout edge effect factor

        Args:
            Ca_min (float): minimum distance from center of an anchor shaft to the edge of concrete
            h_ef (float): Effective embedment depth of anchor

        Returns:
            float: Psi_edN
        """
        Psi_edN = 1.0
        if Ca_min < 1.5*h_ef:
            Psi_edN = 0.7 + 0.3 * (Ca_min/(1.5*h_ef))
        return Psi_edN
    
    def eN_Get(self)-> float:
        pass

    def Psi_ecN_Get(eN : float, h_ef:float)-> float:
        """ACI 318-19 17.6.2.3 Breakout eccentricity factor Psi_ec,N, shall be calculated by Eq.(17.6.2.3.1). Yükün iki dik eksene göre eksantrik olması durumunda Psi_ec,N her eksen için ayrı ayrı hesaplanacak ve bu faktörlerin çarpımı ACI 318-19 Eq.17.6.2.1(b)'de Psi_ec,N olarak kullanılacaktır.

        Args:
            eN (float)  : Ankraj grubuna etkiyen bileşke kuvvetin yeri ile çekme ankrajların yükleme(geometrik merkezi) yeri arasındaki mesafe, Eğer bu eksantriden kaynaklı belirli ankrajlar çekmede oluyorsa eksantrisite hesaplanırken sadece çekmedeki ankrajlar kullanılmalı.
            h_ef (float): Effective embedment depth of anchor

        Returns:
            float: _description_
        """
        Psi_ecN = 1 / (1 + (eN/(1.5*h_ef)))
        return min(1.0,Psi_ecN)

    def BasicSingleAnchorBreakoutStrength(self,kc:float,lambda_a : float, f_c:float, h_ef:float)-> float:
        """ACI318-19 Eq.17.6.2.2.1 kc = 24 for cast-in anchors and 17 for post-installed anchors. 17.6.2.2.1 Basic concrete breakout strength of a single anchor in tension in cracked concrete, Nb, shall be calculated by Eq. (17.6.2.2.1), except as permitted in 17.6.2.2.3

        Args:
            kc (float):   ;SI units kc = 10 or 7 Imperial units kc = 24 or 17
            lambda_a (float): _description_
            f_c (float): Specified compressive strength of concrete
            h_ef (float): _description_

        Returns:
            float: _description_
        """
        Nb = kc * lambda_a * f_c**0.5 * h_ef**1.5
        return round(Nb,3)
    
    def BasicSingleAnchorBreakoutStrength2(self,kc:float,lambda_a : float, f_c:float, h_ef:float)-> float:
        # ACI318-19 17.6.2.2.3 For single cast-in headed studs and headed bolts with 11in(28cm) <= hef <= 25in(63.5cm), Nb shall be calculated by:
        Nb = 3.9 * lambda_a * f_c**0.5 * h_ef**(5/3) # Imperial birimde 3.9 olan katsayı 16 oluyor.
        return round(Nb,3)
    
    def A_Nc(self)-> float:
        pass

    def A_Nco(self,h_ef: float)-> float:
        """Ankraj tüm kenarlara mesafesi en az 1.5 h_ef ise ACI318-19 Eq.17.6.2.1.4

        Args:
            h_ef (float): _description_

        Returns:
            float: _description_
        """
        return 9*h_ef**2

    def ForSingleAnchor(self,A_Nc: float, A_Nco : float, Psi_edN : float, Psi_cN : float, Psi_cpN : float, Nb : float)-> float:
        """_summary_

        Args:
            A_Nc (float): projected concrete failure area of a single anchor or group of anchors, for calculation of strength in tension
            A_Nco (float): projected concrete failure area of a single anchor, for calculation of strength in tension if not limited by edge distance or spacing
            Psi_edN (float): Breakout edge effect factor used to modify tensile strength of anchors based on proximity to edges of concrete member
            Psi_cN (float): breakout cracking factor used to modify tensile strength of anchors based on the influence of cracks in concrete
            Psi_cpN (float): breakout splitting factor used to modify tensile strength of post-installed anchors intended use in uncracked concrete without reinforcement to account for splitting tensile stresses
            Nb (float): basic concrete breakout strength in tension of a single anchor in cracked concrete

        Returns:
            float: _description_
        """
        N_cb = (A_Nc/A_Nco) * Psi_edN * Psi_cN * Psi_cpN * Nb
        return round(N_cb,3)
    
    def ForGroupAnchor(self,A_Nc: float, A_Nco : float,Psi_ecN : float, Psi_edN : float, Psi_cN : float, Psi_cpN : float, Nb : float)-> float:
        """_summary_

        Args:
            A_Nc (float)    : projected concrete failure area of a single anchor or group of anchors, for calculation of strength in tension
            A_Nco (float)   : projected concrete failure area of a single anchor, for calculation of strength in tension if not limited by edge distance or spacing
            Psi_ecN (float) : Breakout edge effect factor used to modify tensile strength of anchors based on proximity to edges of concrete member
            Psi_edN (float) : Breakout edge effect factor used to modify tensile strength of anchors based on proximity to edges of concrete member
            Psi_cN (float)  : breakout cracking factor used to modify tensile strength of anchors based on the influence of cracks in concrete
            Psi_cpN (float) : breakout splitting factor used to modify tensile strength of post-installed anchors intended use in uncracked concrete without reinforcement to account for splitting tensile stresses
            Nb (float)      : basic concrete breakout strength in tension of a single anchor in cracked concrete

        Returns:
            float: _description_
        """
        N_cb = (A_Nc/A_Nco) * Psi_ecN * Psi_edN * Psi_cN * Psi_cpN * Nb
        return round(N_cb,3)

class PulloutStrengthInTension:

    def BasicSingleAnchorPulloutStrength(self,A_brg : float, f_c : float, e_h : float, d_a : float, AnchorType : int, InstalledType : AnchorInstalledType)-> float:
        """Basic Single Anchor Pullout Strength, Np

        Args:
            A_brg (float)   : Net bearing area of the head of stud, anchor bolt, or headed deformed bar,
            f_c (float)     : Specified compressive strength of concrete
            e_h (float)     : distance from the inner surface of the shaft of a J or L-bolt to the outer tip of the J- or L-bolt
            d_a (float)     : outside diameter of anchor or shaft diameter of headed stud, headed bolt, or hooked bolt
            AnchorType (int): 
                                1 : Cast-in headed studs and headed bolts 
                                2 : J or L bolts
            InstalledType (AnchorInstalledType): _description_

        Returns:
            float: _description_
        """
        if InstalledType.value == 2:
            print("ACI318-19 17.6.3.2.1 For post-installed expansion, screw, and undercut anchors, the values of Np shall be based on the 5 percent fractile of results of tests performed and evaluated according to ACI 355.2. It is not permissible to calculate the pullout strength in tension for such anchors.")
        Np = 8 * A_brg * f_c
        if AnchorType == 2 and 3*d_a <= e_h and e_h <= 4.5*d_a:
            Np = 0.9 * f_c * e_h * d_a

        return Np
    
    def Psi_cP_Get(self, IsConcCracked : bool)-> float:
        """ACI318-19 17.6.3.3 Pullout cracking factor,

        Args:
            IsConcCracked (bool): Beton çatlamış(True) mı çatlamamış mı(False)

        Returns:
            float: Psi_cP
        """
        
        Psi_cP = 1.4
        if IsConcCracked:
            Psi_cP = 1.0
        return Psi_cP
    
    def SingleAnchorNominalPulloutStrength(self, Psi_cP: float, Np : float)->float:
        """Single Anchor Nominal Pullout Strength, Npn

        Args:
            Psi_cP (float): ACI318-19 17.6.3.3 Pullout cracking factor
            Np (float): Basic Single Anchor Pullout Strength

        Returns:
            float: Npn
        """
        Npn = Psi_cP * Np
        return round(Npn,3)
    
class ConcreteSideFaceBlowoutStrengthOfHeadedAnchorInTension:
    # Bu kontrol cast-in ankrajlarda genelde görülür post-in ankrajlarda kurulum sırasında ayrılma durumu genelde govern eder(ACI355.2 ye göre gereksinimler karşılanır.) bu nedenle bu kontrol cast-in ankrajlarda yapılır.
    def SingleHeadedAnchorSideFaceBlowoutStrength(self, h_ef : float, C_a1 : float, C_a2 : float, A_brg : float, lambda_a : float, f_c : float)-> float:
        """Single Headed Anchor Side Face Blowout Strength

        Args:
            h_ef (float)    : Effective embedment depth of anchor
            C_a1 (float)    : distance from the center of an anchor shaft to the edge of concrete in one direction, in. If shear is applied to anchor, ca1 is taken in the direction of the applied shear. If tension is applied to the anchor, ca1 is the minimum edge distance. Where anchors subject to shear are located in narrow sections of limited thickness, see R17.7.2.1.2
            C_a2 (float)    : distance from center of an anchor shaft to the edge of concrete in the direction perpendicular to ca1,
            A_brg (float)   : Net bearing area of the head of stud, anchor bolt, or headed deformed bar,
            lambda_a (float): Modification factor to reflect the reduced mechanical properties of lightweight concrete concrete anchorage applications
            f_c (float)     : Specified compressive strength of concrete

        Returns:
            float: N_sb
        """
        N_sb = 0.0
        if h_ef > 2.5*C_a1 :
            N_sb = 13 * C_a1 * A_brg**0.5 * lambda_a * f_c**0.5 # ACI318-19-17.6.4.1 Imperial birimde 13 olan katsayı 160 oluyor.
            if C_a2 < 3*C_a1 and 1.0 <= C_a2/C_a1 and C_a2/C_a1 <= 3.0:
                N_sb = N_sb * (1 + C_a2/C_a1)
        return N_sb
    
    def MultipleHeadedAnchorSideFaceBlowoutStrength(self, s: float, h_ef : float, C_a1 : float, N_sb : float)-> float:
        """ACI318-19 - R17.6.4.2 To calculate nominal side-face blowout strength for multiple headed anchors, only those anchors close to an edge (ca1 < 0.4hef ) that are loaded in tension should be considered. Their strength is compared to the portion of the tensile load applied to those anchors.

        Args:
            s (float)       : is the distance between the outer anchors along the edge
            h_ef (float)    : Effective embedment depth of anchor
            C_a1 (float)    : distance from the center of an anchor shaft to the edge of concrete in one direction, in. If shear is applied to anchor, ca1 is taken in the direction of the applied shear. If tension is applied to the anchor, ca1 is the minimum edge distance. Where anchors subject to shear are located in narrow sections of limited thickness, see R17.7.2.1.2
            N_sb (float)    : Single Headed Anchor Side Face Blowout Strength

        Returns:
            float: N_sbg
        """
        N_sbg = 0.0
        if h_ef > 2.5*C_a1 and s < 6*C_a1:
            N_sbg = (1 + s / (6*C_a1)) * N_sb
        return N_sbg

class BoundStrengthOfAdhesiveAnchorInTension:
    # post-in olan adhesive ankrajlar için yapılmakta. group ankraj ise ankraj aralıkları bond strength için 2*C_Na dan kısa olmalı.
    def MinimumCharacteristicBondStresses(self,ServiceEnviroment : bool, Cracked : bool, DesignIncludes : bool)-> float:
        """Minimum Characteristic Bond Stresses

        Args:
            ServiceEnviroment (bool): True ise içerde uygulanıyor dahil False ise dışarıda uygulanıyor.
            Cracked (bool): True ise beton çatlamamışsa dahil False ise beton çatlamışsa
            DesignIncludes (bool): True ise deprem kuvvetleri dahil False ise sürekli çekme durumu dahil

        Returns:
            float
        """
        if ServiceEnviroment:
            to = 1000/145 #MPa
            if Cracked:
                to = 300/145 #MPa
        if ServiceEnviroment != True:
            to = 650
            if Cracked:
                to = 200

        if DesignIncludes:
            if Cracked:
                to = to * 0.8
            if Cracked != True:
                to = to * 0.8

        if DesignIncludes != True:
            to = to * 0.4

    def C_Na_Get(self, d_a : float, To_uncr : float)-> float:
        """_summary_

        Args:
            d_a (float): _description_
            To_uncr (float): Minimum Characteristic Bond Stresses

        Returns:
            float: _description_
        """
        C_Na = 10 * d_a * (To_uncr / 76)**0.5  # MPa Imperial biriminde bölüm durumundaki 76 katsayısı 1100 dür.
        return round(C_Na,3)
    
    def Psi_edNa_Get(self, C_a_min : float, C_Na : float)->float:
        """Bond edge effect factor

        Args:
            C_a_min (float): _description_
            C_Na (float): _description_

        Returns:
            float: _description_
        """
        Psi_edNa = 1.0
        if C_a_min < C_Na:
            Psi_edNa = 0.7 + 0.3 * C_a_min / C_Na
        
        return round(Psi_edNa ,3)

    def Psi_cpNa_Get(self,  C_a_min : float, C_ac : float, C_Na : float)->float:
        """Bond splitting factor,

        Args:
            C_a_min (float): _description_
            C_ac (float): _description_
            C_Na (float): _description_

        Returns:
            float: _description_
        """
        Psi_cpNa = 1.0
        if C_a_min < C_ac:
            Psi_cpNa = max(C_a_min/C_ac, C_Na/C_ac)
        return round(Psi_cpNa ,3)

    def Psi_ecNa_Get(self,e_N : float, C_Na : float)->float:
        """_summary_

        Args:
            e_N (float): _description_
            C_Na (float): _description_

        Returns:
            float: _description_
        """
        # 17.6.5.3.1 Modification factor for adhesive anchor groups loaded eccentrically in tension, Psi_ec,Na, shall be calculated by Eq (17.6.5.3.1).

        Psi_ecNa = min(1.0 , 1/(1 + e_N/C_Na))
        return round(Psi_ecNa,3)
    
    def A_Na_Get(self, C_Na : float, s1 : float, s2 : float, C_a1 : float, C_a2 : float)-> float:
        """# ACI 318-19 - 17.6.5.1.1 ANa is projected influence area of a single adhesive anchor or an adhesive anchor group that is approximated as a rectilinear area that projects outward a distance cNa from the centerline of the adhesive anchor, or in the case of an adhesive anchor group, from a line through a row of adjacent adhesive anchors. ANa shall not exceed nANao, where n is the number of adhesive anchors in the group that resist tension.
        # ACI 318-19 - 17.6.5.1.1 A_Na, tek bir yapışkan ankrajın veya bir yapışkan ankraj (Yapışkan ankraj dediği epoksi ile sonradan ekilmiş ankraj grubu post-in ankraj grubu yani) grubunun, yapışkan ankrajın merkez hattından veya bir yapışkan ankraj grubu durumunda bir çizgiden dışarıya doğru bir cNa mesafesi kadar dışarı doğru çıkıntı yapan doğrusal bir alan olarak yaklaşık olarak tahmin edilen, öngörülen etki alanıdır. bir sıra bitişik yapışkan ankraj.A_Na, n*A_Nao'yu aşmayacaktır; burada n, gruptaki gerilime direnen yapışkan ankrajların sayısıdır.
        # Fig R17.6.5.1 Calculation of influence areas ANao and ANa.

        Args:
            C_Na (float): projected distance from center of an anchor shaft on one side of the anchor required to develop the full bond strength of a single adhesive anchor
            s1 (float): X doğrultusunda ankraj aralığı
            s2 (float): Y doğrultusunda ankraj aralığı
            C_a1 (float): X doğrultusunda kenara en yakın ankrajın merkezinden kenar ile arasındaki mesafe
            C_a2 (float): Y doğrultusunda kenara en yakın ankrajın merkezinden kenar ile arasındaki mesafe

        Returns:
            float: _description_
        """
        
        A_Na = (C_Na + s1 + C_a1) * (C_Na + s2 + C_a2)
        return round(A_Na,3)

    def A_Nao_Get(self, C_Na : float)-> float:
        """_summary_

        Args:
            C_Na (float): projected distance from center of an anchor shaft on one side of the anchor required to develop the full bond strength of a single adhesive anchor

        Returns:
            float: A_Nao
        """
        # 17.6.5.1.2 ANao is the projected influence area of a single adhesive anchor with an edge distance of at least cNa:
        A_Nao = (2*C_Na)**2
        return round(A_Nao,3)
        
    def N_ba_Get(self, d_a:float, h_ef:float, To_cr:float, lambda_a: float)-> float:
        """_summary_

        Args:
            d_a (float): _description_
            h_ef (float): _description_
            To_cr (float): Minimum Characteristic Bond Stresses
            lambda_a (float): _description_

        Returns:
            float: _description_
        """
        # ACI 318-19 - 17.6.5.2.1 Basic bond strength of a single adhesive anchor in tension in cracked concrete, Nba,
        N_ba = lambda_a * To_cr * 3.14 * d_a * h_ef
        return round(N_ba,3)
    
    def SingleAnchorNominalBondStrengthInTension(self, A_Na : float, A_Nao : float, Psi_edNa : float, Psi_cpNa : float, N_ba : float)-> float:
        Na = (A_Na/A_Nao) * Psi_edNa * Psi_cpNa * N_ba
        return round(Na,3)
    
    def GroupAnchorNominalBondStrengthInTension(self, A_Na : float, A_Nao : float, Psi_edNa : float, Psi_cpNa : float, Psi_ecNa : float, N_ba : float)-> float:
        Nag = (A_Na/A_Nao) * Psi_edNa * Psi_cpNa * Psi_ecNa * N_ba
        return round(Nag,3)

# Shear

class SteelStrengthOfAnchorInShear:
    """ACI318-19 - 17.7.1.2.1 If anchors are used with built-up grout pads, nominal strength Vsa calculated in accordance with 17.7.1.2 shall be multiplied by 0.80.
    Ankrajların yerleşik harç pedleriyle birlikte kullanılması durumunda, 17.7.1.2'ye göre hesaplanan nominal dayanım Vsa, 0,80 ile çarpılacaktır.
    """
    def A_seV(self,d_a : float, n_t : float)-> float:
        """the effective cross-sectional area of an anchor in shear.For threaded rods and headed bolts, ASME B1.1 defines Ase,V as

        Args:
            d_a (float): outside diameter of anchor or shaft diameter of headed stud, headed bolt, or hooked bolt
            n_t (float): is the number of threads per inch or mm

        Returns:
            float: _description_
        """

        A_seV = 3.14/4 * (d_a - (0.9743/n_t))
        return round(A_seV,2)
    
    def NominalSteelStrengthOfCastInHeadedStudAnchorInShear(self, A_seV: float, f_ya : float,f_uta : float)-> float:
        """ACI318-19 17.7.1.2(a) - Nominal Steel Strength Of Cast-In Headed Stud Anchor In Shear

        Args:
            A_seV (float): is the effective cross-sectional area of an anchor in shear mm^2
            f_ya (float): specified yield strength of anchor steel, shall not exceed 1.9f_ya or 862MPa
            f_uta (float): specified tensile strength of anchor steel, shall not exceed 1.9f_ya or 862MPa

        Returns:
            float: V_sa
        """
        f_uta = min(f_uta,1.9*f_ya,862)
        return round(A_seV * f_uta,3)
    
    def NominalSteelStrengthCastInHeadedAndHookedBoltAnchorInShear(self, A_seV: float, f_ya : float,f_uta : float)-> float:
        """ACI318-19 17.7.1.2(b) - For cast-in headed bolt and hooked bolt anchors and for post-installed anchors where sleeves do not extend through the shear plane
            ACI318-19 17.7.1.2(c) - For post-installed anchors where sleeves extend through the shear plane, Vsa shall be based on the 5 percent fractile of results of tests performed and evaluated in accordance with ACI 355.2. Alternatively, Eq. (17.7.1.2b) shall be permitted to be used.

        Args:
            A_seV (float): is the effective cross-sectional area of an anchor in shear mm^2
            f_ya (float): specified yield strength of anchor steel, shall not exceed 1.9f_ya or 862MPa
            f_uta (float): specified tensile strength of anchor steel, shall not exceed 1.9f_ya or 862MPa

        Returns:
            float: V_sa
        """
        f_uta = min(f_uta,1.9*f_ya,862)
        return round(0.6* A_seV * f_uta,3)

class ConcreteBreakoutStrengthOfAnchorInShear:
    def A_Vco_Get(self, C_a1 : float)-> float:
        """ACI318-19 - Fig.R17.7.2.1a

        Args:
            C_a1 (float): _description_

        Returns:
            float: _description_
        """
        A_Vco = 4.5*C_a1**2
        return round(A_Vco,2)
    
    # TODO Hazır değil Fig. R17.7.2.1.2
    def A_Vc_Get(self, C_a1 : float, h_a : float, C_a2 : float, s : float)-> float:
        # ACI318-19 - Fig. R17.7.2.1b—Calculation of Avc for single anchors and anchor groups.

        if h_a < 1.5*C_a1: #For single anchor ortada olan ankraj
            A_Vc = 3 * C_a1 * h_a

        if C_a2 < 1.5*C_a1: #For single anchor köşede olan ankraj
            A_Vc = 1.5 * C_a1 * (1.5*C_a1 + C_a2)
    # TODO Hazır değil Fig. R17.7.2.1.2
    def A_Vcg_Get(self, C_a11 : float, C_a12 : float, h_a : float, C_a2 : float, s : float, V_dir : bool , Case : int = 0)-> float:
        """Fig. R17.7.2.1b—Calculation of Avc for single anchors and anchor groups.

        Args:
            C_a11 (float)       : Kenara en yakın ankrajın merkezinden, kesme kuvvetine dik doğrultudaki beton kenarına olan mesafe
            C_a12 (float)       : Kenara en yakın ikinci ankrajın merkezinden, kesme kuvvetine dik doğrultudaki beton kenarına olan mesafe
            h_a (float)         : Beton kalınlığı
            C_a2 (float)        : Kenara en yakın ankrajın merkezinden, kesme kuvveti doğrultusundaki beton kenarına olan mesafe
            s (float)           : Ankraj aralığı ilgili doğrultuda
            V_dir (bool)        : Kuvvet doğrultusu True paralel False dik
            Case (int, optional): 0,1,2,3. Defaults to 0.
                                    
                                    Case 1: One assumption of the distribution of forces indicates that half of the shear force would be critical on the front anchor and the projected area. For the calculation of concrete breakout, ca1 is taken as ca1,1.
                                    
                                    Case 2: Another assumption of the distribution of forces indicates that the total shear force would be critical on the rear anchor and its projected area. Only this assumption needs to be considered when anchors are welded to a common plate independent of s. For the calculation of concrete breakout, ca1 is taken as ca1,2.

                                    Case 3: Where s < ca1,1, apply the entire shear load V to the front anchor. This case does not apply for anchors welded to a common plate. For the calculation of concrete breakout, ca1 is taken as ca1,1.
        Returns:
            float: _description_
        """

        if h_a < 1.5*C_a11 and s < 3*C_a11 and V_dir != True: #For group anchor V kuvveti ankraj hesap doğrultusuna dik
            A_Vc = (3 * C_a11 + s) * h_a


        if h_a < 1.5*C_a11 and s >= C_a11 and Case == 1 and V_dir: #For group anchor V kuvveti ankraj hesap doğrultusuna paralel Case1
            A_Vc = 3 * C_a11 * h_a

        if h_a < 1.5*C_a12 and s > C_a11 and Case == 2 and V_dir: #For group anchor V kuvveti ankraj hesap doğrultusuna paralel Case2
            A_Vc = 3 * C_a12 * h_a
        
        if h_a < 1.5*C_a11 and s < C_a11 and Case == 3 and V_dir: #For group anchor V kuvveti ankraj hesap doğrultusuna paralel Case3
            A_Vc = 3 * C_a11 * h_a

    def ForSingleAnchor(self,A_Vc: float, A_Vco : float, Psi_edV : float, Psi_cV : float, Psi_hV : float, V_b : float)-> float:
        """ACI318-19 - Eq.17.7.2.1(a)

        Args:
            A_Vc (float)    : 
            A_Vco (float)   : 
            Psi_edV (float) : 
            Psi_cV (float)  : 
            Psi_hV (float)  : 
            V_b (float)     : basic concrete breakout strength in shear of a single anchor in cracked concrete

        Returns:
            float: V_cb
        """
        N_cb = (A_Vc/A_Vco) * Psi_edV * Psi_cV * Psi_hV * V_b
        return round(N_cb,3)
    
    def ForGroupAnchor(self,A_Vc: float, A_Vco : float, Psi_ecV : float, Psi_edV : float, Psi_cV : float, Psi_hV : float, V_b : float)-> float:
        """ACI318-19 - Eq.17.7.2.1(b)

        Args:
            A_Vc (float)    : 
            A_Vco (float)   : 
            Psi_edV (float) : 
            Psi_cV (float)  : 
            Psi_hV (float)  : 
            V_b (float)     : basic concrete breakout strength in shear of a single anchor in cracked concrete

        Returns:
            float: V_cb
        """
        N_cb = (A_Vc/A_Vco) * Psi_ecV * Psi_edV * Psi_cV * Psi_hV * V_b
        return round(N_cb,3)

        
class ConcretePryoutStrengthOfAnchorInShear:
    
    def SingleAnchorConcPryoutStrengthInShear(self, h_ef : float, N_cp : float)-> float:
        k_cp = 1.0
        if h_ef >= 2.5*25.4:
            k_cp = 2.0
        
        V_cp = k_cp * N_cp
        round(V_cp,2)
        
    def GroupAnchorConcPryoutStrengthInShear(self, h_ef : float, N_cpg : float)-> float:
        k_cp = 1.0
        if h_ef >= 2.5*25.4:
            k_cp = 2.0
        V_cpg = k_cp * N_cpg

        return round(V_cpg,2)


# Tension-Shear Interaction




# def main() -> None:
#     # conn = AnchorConnection(d=12.7, bf= 12.2, B= 240, N= 240, tf=10, P_axial=700, M=50, fck= 3)
#     # BasePlat = CheckBasePlate(Contn=conn, ReductionFactor=0.65)
#     # print(f' A1 ={BasePlat.A1}, N ={BasePlat.N}, ')
#     pt = CastInAnchorType(1)
#     print(pt.__class__.__name__.split("In")[0])

# if __name__ == "__main__":
#     main()