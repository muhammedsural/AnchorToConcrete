"""
    h_a     : thickness of member in which an anchor is located,measured parallel to anchor axis
    h_ef    : effective embedment depth of anchor
    h'_ef   : limiting value of hef where anchors are located less than 1.5hef from three or more edges, in.; refer to Fig. R17.6.2.1.2
    h_efsl  : effective embedment depth of shear lug
    h_sl    : embedment depth of shear lug

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


adhesive                                : chemical components formulated organic polymers, or a combination of polymers and inorganic materials that cure if blended together.
anchor                                  : a steel element either cast into concrete or post-installed into a hardened concrete member and used to transmit applied loads to the concrete.
anchor, post-installed                  : anchor installed in hardened concrete; adhesive, expansion, screw, and undercut anchors are examples of post-installed anchors.
anchor, horizontal or upwardly inclined : Anchor installed in a hole drilled horizontally or in a hole drilled at any orientation above horizontal.
"""
from dataclasses import dataclass
from enum import Enum
import math

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
    
class AnchorMaterial:
    fy : float
    fu : float
    fya : float
    futa : float
    
@dataclass
class Anchor:
    h1 : float
    h2 : float
    h3 : float
    h4 : float
    h5 : float #h_eff
    D  : float
    t  : float
    AnchorType : CastInAnchorType | PostInAnchorType
    Material : AnchorMaterial

class ConcArea:
    Bx : float
    By : float
    Ha : float

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

    def ApproximateBasePlateArea(self,P_u : float, f_c : float, Case : int = 3, fi : float = 0.65)-> float:
        """ A1 : Taban plakasi alani
            A2 : Beton yüzey alani

            

        Args:
            P_u (float): Gelen eksenel basınç yükü
            f_c (float): Beton basınç dayanımı
            Case (int, optional): _description_. Defaults to 3.
                                    Case 1 : A1 = A2
                                    Case 2 : A2 >= 4A1
                                    Case 3 : A1 < A2 < 4A1
            fi (float, optional): _description_. Defaults to 0.65.

        Returns:
            float: _description_
        """
        if Case == 1:
            A1_req = P_u / (fi * 0.85 * f_c)

        if Case == 2 or Case == 3:
            A1_req = P_u / (2 * fi * 0.85 * f_c)
        return round(A1_req,0)

    def FindPlateDimensions(self,d : float, b_f : float, A1_req : float)-> int:
        delta = (0.95 * d - 0.8 * b_f) / 2
        
        N_1 = A1_req**0.5 + delta
        # The base plate dimension N × B should be large enough for the installation of four anchor rods, as required by OSHA.
        # N_2 = d + 2 * 80 #80mm OSHA için, Abi O. Aghayere çelik yapı tasarımı kitabı 100mm öneriyor.
        N_2 = 0 
        N = max(N_1,N_2)

        B_1 = A1_req / N
        # B_2 = b_f + 2 * 80 #80mm OSHA için, Abi O. Aghayere çelik yapı tasarımı kitabı 100mm öneriyor.
        B_2 = 0
        B = max(B_1,B_2)

        N = math.ceil(max(B/10,N/10))*10
        B = N*10
        A_baseplate = B * N
        if A1_req > A_baseplate:
            raise ValueError(f"A1_req = {A1_req} > {A_baseplate} = A_baseplate")        
        return N

    def e_Get(self,M_u : float, P_u : float)-> float:
        return round(M_u/P_u,2)

    def f_pmax_Get(self,f_c : float, A1 : float, A2 : float)-> float:
        """Bearing stress max

        Args:
            f_c (float): specified compressive strength of concrete MPa
            A1 (float): area of steel concentrically bearing on a concrete support,(mm2)
            A2 (float): maximum area of the portion of the supporting surface that is geometrically similar to and concentric with the loaded area,(mm2)

        Raises:
            ValueError: A2 notless than A1

        Returns:
            float: _description_
        """
        if A2 < A1:
            raise ValueError("A2, A1'den küçük olamaz!!!")
        if A2 >= A1:
            f_pmax = 0.85 * f_c * (A2/A1)**0.5

        if (A2/A1)**0.5 >= 1 and (A2/A1)**0.5 <= 2:
            print("Beton sargılama katkısı kullanılabilir.")

        return round(f_pmax,2)

    # Türkiye çelik yapılar yönetmeliği denk. 13.22-23
    def NominalBearingStrengthConc(self,f_pmax : float, A1 : float, f_c : float)-> float:
        """The nominal bearing strength of conc

        Args:
            f_pmax (float): Max bearing stress
            A1 (float): area of steel concentrically bearing on a concrete support,(mm2)

        Returns:
            float: _description_
        """
        P_p = f_pmax * A1
        P_p = min(P_p, 1.7 * f_c * A1)
        return round(P_p,2)

    def DesignBearingStrengthConc(self,Pp : float, fi : float = 0.65) -> float:
        return round(fi*Pp,2)




    def q_max_Get(self,f_pmax : float, B : float)-> float:
        q_max = f_pmax * B
        return round(q_max,2)

    def e_crit_Get(self,q_max : float, P_u : float, N : float)-> float:
        e_crit = N/2 - (P_u / (2*q_max))
        return round(e_crit,2)

    def Get_m(self,N : float, d : float)->int:
        return round((N-0.95*d)/2,2)

    def Get_n(self,B : float, b_f : float)->int:
        return round((B-0.8*b_f)/2,2)

    def Get_X(self,d : float, b_f : float, P_u : float, P_p : float, fi : float = 0.9)->float:
        X = (4 * d * b_f * P_u) / ((d + b_f)**2 *  P_p)
        return round(X,2)

    def Get_lambda(self,X : float)->float:
        lambda_x = 1.0
        if X <1:
            lambda_x = (2*X**0.5) / (1 + (1-X)**0.5)
        lamb_ = min(1.0,lambda_x)
        return round(lamb_,2)

    def Get_l(self,d : float, b_f : float, m : float, n : float, lamb : float)-> float:
        n_lamb = lamb * (d * b_f)**0.5 / 4
        return max(m,n,n_lamb)

    def BasePlateThickness(self,P_u : float, l : float, B : float, N : float, F_y : float, fi : float = 0.9)->int:
        t_min = l * ((2 * P_u) / (fi * F_y * B * N))**0.5
        t_min = math.ceil(t_min)
        while t_min%5 != 0:
             t_min += 1
        return round(t_min,2)

    def Get_Y(self,e : float, e_crit : float, P_u : float, N : float, f : float, q_max : float) -> float:
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

    def Get_f_p(self,P_u : float, B : float, Y : float) -> float:
        return round(P_u/(B*Y),2)

    def ForMomentsBasePlateThickness(self,m : float, n : float, f_p : float, Y : float, F_y : float)-> int:
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

    def GetBasePlateArea(self,B : float, N : float) -> float:
        A1 = B * N
        return round(A1,2)

    def Get_t(self,N : int, x : float) -> float:
        """ankraj merkezinden plaka orta noktasına olan mesafeyi hesaplar

        Args:
            N (int): İlgili doğrultuda plaka boyutu
            x (float): Ankraj rodun merkezinin plaka kenarına olan mesafesi

        Returns:
            float: t
        """
        t = N/2 - x #ankraj merkezinden plaka orta noktasına olan mesafe
        return round(t,2)

    def Get_q(self,P_u : float, Y : float) -> float:
        return round(P_u/Y,2)

@dataclass
class AnchorConnection:
    basePlate : BasePlate
    anchor    : Anchor
    concArea  : ConcArea


def Table1753a(SteelType : int, WhichEffect : int) -> float:
    """_summary_

    Args:
        SteelType (int): 0 is ductile others brittle
        WhichEffect (int): 0 is Tension others is Shear

    Returns:
        float: _description_
    """
    if SteelType == 0 :
        if WhichEffect == 0:
            fi = 0.75
        else:
            fi = 0.65
    else:
        if WhichEffect == 0:
            fi = 0.65
        else:
            fi = 0.6
    return fi

def Table1753b(SupReinforcement : bool, AncTypeInstall : bool, WhichEffect : int, AnchorCategory : int) -> float:
    """Anchor strength governed by concrete breakout, bond, and side-face blowout

    Args:
        SupReinforcement (bool): True is yes, False is no
        AncTypeInstall (bool): True is Cast-in. False is Post-Installed
        WhichEffect (int): 0 is Tension others is Shear
        AnchorCategory (int): 1 or 2 or 3 for post-installed

    Returns:
        float: _description_
    """
    
    if SupReinforcement :
        if AncTypeInstall:
            fi = 0.75
        else:
            if WhichEffect != 0:
                fi = 0.75
            else:
                if AnchorCategory == 1:
                    fi = 0.75
                if AnchorCategory == 2:
                    fi = 0.65
                if AnchorCategory == 3:
                    fi = 0.55
                else:
                    fi == 0.0 #Non Defined
    else:
        if AncTypeInstall:
            fi = 0.70 
        else:
            if WhichEffect != 0:
                fi = 0.75
            else:
                if AnchorCategory == 1:
                    fi = 0.75
                if AnchorCategory == 2:
                    fi = 0.65
                if AnchorCategory == 3:
                    fi = 0.55
                else:
                    fi == 0.0 #Non Defined
    return fi

def Table1753c(AncTypeInstall : bool, WhichEffect : int, AnchorCategory : int) -> float:
    """Anchor strength governed by concrete pullout, or pryout strength

    Args:
        AncTypeInstall (bool): True is Cast-in. False is Post-Installed
        WhichEffect (int): 0 is Tension others is Shear
        AnchorCategory (int): 1 or 2 or 3 for post-installed

    Returns:
        float: _description_
    """
    if AncTypeInstall:
        fi = 0.7
    else:
        if WhichEffect != 0:
                fi = 0.70
        else:
            if AnchorCategory == 1:
                fi = 0.65
            if AnchorCategory == 2:
                fi = 0.55
            if AnchorCategory == 3:
                fi = 0.45
            else:
                fi == 0.0 #Non Defined
            
def h_ef_Check(h_eff : float, Ca_x1 : float, Ca_x2 : float, Ca_y1 : float, Ca_y2 : float, sx : float, sy : float) -> float:
    """

    Conditions
                1- Eğer ankraj grubu 3 kenarada 1.5*h_ef ten daha yakınsa h_ef aşağıdaki şarta göre belirlenir. Normal şartlarda gömme derinliğidir.
                        h_ef = max(C_amax/1.5, s/3)

    Returns:
        float: _description_
    """
    edgeLimitCounter = 0
    if Ca_x1 < 1.5*h_eff:
        edgeLimitCounter += 1

    if Ca_x2 < 1.5*h_eff:
        edgeLimitCounter += 1

    if Ca_y1 < 1.5*h_eff:
        edgeLimitCounter += 1

    if Ca_y2 < 1.5*h_eff:
        edgeLimitCounter += 1
    
    if edgeLimitCounter >= 3:
        smax = max(sx, sy)
        Camax = max(Ca_x1,Ca_x2,Ca_y1,Ca_y2)
        h_eff = max(Camax/1.5, smax/3)
        print(f"Ankraj grubu en az 3 beton kenarına 1.5*h_ef ten daha yakın h_ef:{h_eff} olarak düzeltildi...")
        
    return round(h_eff,2)

def A_Nco(h_ef : float) -> float:
    """ACI 318, Eq. 17.6.2.1.4

    Args:
        h_ef (float): _description_

    Returns:
        float: _description_
    """
    return 9 * h_ef**2

#INFO Çekme etkisindeki ankraj sayısına göre alan değişebilmekte bunun için çekme bölgesi doğrultusu belirtilmeli.
#INFO Bu fonksiyonun geliştirilmesi gerekebilir.
def A_Nc(Anc_x : float, Anc_y : float, h_ef : float, Ca_x1 : float, Ca_x2 : float, Ca_y1 : float, Ca_y2 : float, n : int, TensionRegion : int) -> float:
    """_summary_

    Args:
        Anc_x (float): X doğrultusunda dizilmiş iki uçtaki ankrajların merkezleri arasındaki mesafe
        Anc_y (float): Y doğrultusunda dizilmiş iki uçtaki ankrajların merkezleri arasındaki mesafe
        h_ef (float): Ankrajın betona gömülme derinliği. Taban plakasının alt yüzeyinden ankraj alt başlığı veya levhasının üst yüzeyine olan mesafe.
        Ca_x1 (float): _description_
        Ca_x2 (float): _description_
        Ca_y1 (float): _description_
        Ca_y2 (float): _description_
        n (int): number of anchors in the group that resist tension.
        TensionRegion (int): which direction resist tension.
                            1 -> X dir
                            2 -> Y dir
                            3 -> All in tension

    Returns:
        float: _description_
    """
    h_ef_üs = h_ef_Check(h_ef, Ca_x1, Ca_x2, Ca_y1, Ca_y2, Anc_x, Anc_y)
    x1 = Ca_x1
    x2 = Ca_x2
    y1 = Ca_y1
    y2 = Ca_y2

    if 1.5*h_ef_üs < Ca_x1:
        x1 = 1.5*h_ef_üs

    if 1.5*h_ef_üs < Ca_x2:
        x2 = 1.5*h_ef_üs

    if 1.5*h_ef_üs < Ca_y1:
        y1 = 1.5*h_ef_üs

    if 1.5*h_ef_üs < Ca_y1:
        y2 = 1.5*h_ef_üs

    
    X_edge = x1 + x2
    if TensionRegion ==1 or TensionRegion == 3:
        X_edge += Anc_x
    
    Y_edge = y1 + y2
    if TensionRegion == 2 or TensionRegion == 3:
        Y_edge += Anc_y
            
    A_Nc = X_edge * Y_edge

    # n*A_Nco < A_Nc olamaz
    if n*9*h_ef**2 < A_Nc:
        A_Nc = n*9*h_ef**2
    return A_Nc

# TENSION LIMIT STATE
def NominalSteelStrengthOfAnchorInTension(A_se: float, f_ya : float,f_uta : float)-> float:
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
    """ACI 318-19 17.6.2     """
    
    def Psi_cpN_Get(self,Cac : float, Ca_min : float, h_ef: float, InstalledType : CastInAnchorType | PostInAnchorType) -> float:
        """ACI 318-19 17.6.2.6 Breakout splitting factor (For post installed anchors)

        Args:
            Cac (float): critical edge distance required to develop the basic strength as controlled by concrete breakout or bond of a post-installed anchor in tension in uncracked concrete without supplementary reinforcement to control splitting,
            Ca_min (float): minimum distance from center of an anchor shaft to the edge of concrete
            h_ef (float): Effective embedment depth of anchor
            InstalledType (CastInAnchorType | PostInAnchorType): _description_

        Returns:
            float: _description_
        """
        # 
        Psi_cpN = 1.0
        
        if InstalledType.__class__.__name__ != "CastInAnchorType":
            if Ca_min >= Cac:
                Psi_cpN = 1.0
            if Ca_min < Cac:
                Psi_cpN = max((Ca_min/Cac),(1.5*h_ef/Cac))
        
        return Psi_cpN
        
    def Psi_cN_Get(self, InstalledType : CastInAnchorType | PostInAnchorType, IsConcCracked : bool) -> float:
        # ACI 318-19 17.6.2.5 Breakout cracking factor, 1.25 for cast-in anchors, 1.4 for post-installed anchors,
        Psi_cN = 1.0
        if InstalledType.__class__.__name__ == "CastInAnchorType" and IsConcCracked == False:
            Psi_cN = 1.4
        if InstalledType.__class__.__name__ == "PostInAnchorType" and IsConcCracked == False:
            Psi_cN = 1.25
            
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

    def Psi_ecN_Get(self,eN : float, h_ef:float)-> float:
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

        if h_ef >= 11*25.4 and h_ef <= 25*25.4:
            Nb = 3.9 * lambda_a * f_c**0.5 * h_ef**(5/3) # ACI318-19 17.6.2.2.3 For single cast-in headed studs and headed bolts with 11in(28cm) <= hef <= 25in(63.5cm), coeff for Imperial units 16 not 3.9 
        return round(Nb,3)
        
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
        modificationFactor = (A_Nc/A_Nco)
        N_cb = modificationFactor * Psi_edN * Psi_cN * Psi_cpN * Nb
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
        modificationFactor = (A_Nc/A_Nco)
        N_cb = modificationFactor * Psi_ecN * Psi_edN * Psi_cN * Psi_cpN * Nb
        return round(N_cb,3)

class PulloutStrengthInTension:

    def BasicSingleAnchorPulloutStrength(self,A_brg : float, f_c : float, e_h : float, d_a : float, AnchorType : int)-> float:
        """Basic Single Anchor Pullout Strength, Np

        Args:
            A_brg (float)   : Net bearing area of the head of stud, anchor bolt, or headed deformed bar,
            f_c (float)     : Specified compressive strength of concrete
            e_h (float)     : distance from the inner surface of the shaft of a J or L-bolt to the outer tip of the J- or L-bolt
            d_a (float)     : outside diameter of anchor or shaft diameter of headed stud, headed bolt, or hooked bolt
            AnchorType (int): 
                                1 : Cast-in headed studs and headed bolts 
                                2 : J or L bolts

        Returns:
            float: _description_
        """
        if AnchorType.__class__.__name__ == "PostInAnchorType":
            print("ACI318-19 17.6.3.2.1 For post-installed expansion, screw, and undercut anchors, the values of Np shall be based on the 5 percent fractile of results of tests performed and evaluated according to ACI 355.2. It is not permissible to calculate the pullout strength in tension for such anchors.")
        Np = 8 * A_brg * f_c
        if AnchorType == 2 and 3*d_a <= e_h and e_h <= 4.5*d_a:
            Np = 0.9 * f_c * e_h * d_a

        return round(Np,3)
    
    def Psi_cP_Get(self, IsConcCracked : bool)-> float:
        """ACI318-19 17.6.3.3 Pullout cracking factor,

        Args:
            IsConcCracked (bool): Beton çatlamış(True) mı çatlamamış mı(False)

        Returns:
            float: Psi_cP
        """
        Psi_cP = 1.0
        if IsConcCracked == False:
            Psi_cP = 1.4
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
        N_sb = 10**50
        if h_ef > 2.5*C_a1 :
            N_sb = 13 * C_a1 * A_brg**0.5 * lambda_a * f_c**0.5 # ACI318-19-17.6.4.1 Imperial birimde 13 olan katsayı 160 oluyor.
            if C_a2 < 3*C_a1 and 1.0 <= C_a2/C_a1 and C_a2/C_a1 <= 3.0:
                N_sb = N_sb * (1 + C_a2/C_a1)
        return round(N_sb,3)
    
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
        N_sbg = 10**50
        if h_ef > 2.5*C_a1 and s < 6*C_a1:
            N_sbg = (1 + s / (6*C_a1)) * N_sb
        return round(N_sbg,3)

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


# SHEAR LIMIT STATE
# Ca1 = the edge distance in the direction of load

class SteelStrengthOfAnchorInShear:
    """ACI318-19 - 17.7.1.2.1 If anchors are used with built-up grout pads, nominal strength Vsa calculated in accordance with 17.7.1.2 shall be multiplied by 0.80.
    Ankrajların yerleşik harç pedleriyle birlikte kullanılması durumunda, 17.7.1.2'ye göre hesaplanan nominal dayanım Vsa, 0,80 ile çarpılacaktır.
    """
    def A_seV(self,d_a : float, n_t : float)-> float:
        """the effective cross-sectional area of an anchor in shear.

        Args:
            d_a (float): outside diameter of anchor or shaft diameter of headed stud, headed bolt, or hooked bolt
            n_t (float): is the number of threads per inch or mm

        Returns:
            float: _description_
        """

        A_seV = 3.14/4 * (d_a - (0.9743/n_t))**2
        return round(A_seV,2)
    
    def NominalSteelStrengthOfCastInHeadedStudAnchorInShear(self, A_seV: float, f_ya : float,f_uta : float)-> float:
        """ACI318-19 17.7.1.2(a) - Nominal Steel Strength Of Cast-In Headed Stud Anchor In Shear
        Grout harcı standart hale geldiği için 0.8 lik azaltma direkt uygulanmıştır. Design guide yazarlarına göre böyle bir azaltmaya gerek yoktur. ACI318-19 17.7.1.2.1

        Args:
            A_seV (float): is the effective cross-sectional area of an anchor in shear mm^2
            f_ya (float): specified yield strength of anchor steel, shall not exceed 1.9f_ya or 862MPa
            f_uta (float): specified tensile strength of anchor steel, shall not exceed 1.9f_ya or 862MPa

        Returns:
            float: V_sa
        """
        f_uta = min(f_uta,1.9*f_ya,862)
        return round(0.8* A_seV * f_uta,3)
    
    def NominalSteelStrengthCastInHeadedAndHookedBoltAnchorInShear(self, A_seV: float, f_ya : float,f_uta : float)-> float:
        """ACI318-19 17.7.1.2(b) - For cast-in headed bolt and hooked bolt anchors and for post-installed anchors where sleeves do not extend through the shear plane
            ACI318-19 17.7.1.2(c) - For post-installed anchors where sleeves extend through the shear plane, Vsa shall be based on the 5 percent fractile of results of tests performed and evaluated in accordance with ACI 355.2. Alternatively, Eq. (17.7.1.2b) shall be permitted to be used.
            Grout harcı standart hale geldiği için 0.8 lik azaltma direkt uygulanmıştır. Design guide yazarlarına göre böyle bir azaltmaya gerek yoktur. ACI318-19 17.7.1.2.1

        Args:
            A_seV (float): is the effective cross-sectional area of an anchor in shear mm^2
            f_ya (float): specified yield strength of anchor steel, shall not exceed 1.9f_ya or 862MPa
            f_uta (float): specified tensile strength of anchor steel, shall not exceed 1.9f_ya or 862MPa

        Returns:
            float: V_sa
        """
        f_uta = min(f_uta,1.9*f_ya,862)
        return round(0.8* 0.6* A_seV * f_uta,3)

class ConcreteBreakoutStrengthOfAnchorInShear:

    def Get_Ca1(self,
                Anc_x : float, 
                Anc_y : float, 
                ha    : float, 
                Ca_x1 : float, 
                Ca_x2 : float, 
                Ca_y1 : float, 
                Ca_y2 : float, 
                ShearDirection : int, 
                IsWeldCommonPlate : bool) -> float:
        """_summary_

        Args:
                Anc_x (float)       : X doğrultusunda dizilmiş iki uçtaki ankrajların merkezleri arasındaki mesafe
                Anc_y (float)       : Y doğrultusunda dizilmiş iki uçtaki ankrajların merkezleri arasındaki mesafe
                ha (float)          : Beton kalınlığı
                Ca_x1 (float)       : Y doğrultusunda sol kenara olan mesafe
                Ca_x2 (float)       : Y doğrultusunda sağ kenara olan mesafe
                Ca_y1 (float)       : Y doğrultusunda sol kenara olan mesafe
                Ca_y2 (float)       : Y doğrultusunda sağ kenara olan mesafe
                ShearDirection (int): Kesme kuvveti yönü
                                        1 -> +X dir
                                        2 -> -X dir
                                        3 -> +Y dir
                                        4 -> -Y dir
                IsWeldCommonPlate   : Ankrajlar ortak bir plakaya kaynaklı mı? 

        Returns:
                float: _description_
        """
        if ShearDirection == 1 or ShearDirection == 2:
            if ShearDirection == 1:
                    Ca1 = Ca_y1 #Kesme kuvvetine dik olan kenara, kesme kuvvetini aktaran bulon grubunun merkezinin uzaklığı
                    if IsWeldCommonPlate:
                        Ca1 = Ca_y1 + Anc_y #Bulonlar genel sabitleme plakasına kaynaklıysa

            if ShearDirection == 2:
                    Ca1 = Ca_y2 #Kesme kuvvetine dik olan kenara, kesme kuvvetini aktaran bulon grubunun merkezinin uzaklığı
                    if IsWeldCommonPlate:
                        Ca1 = Ca_y2 + Anc_y #Bulonlar genel sabitleme plakasına kaynaklıys

            # ha ve Ca2(Kesme kuvvetini aktaran bulon grubunun uç ankrajlarının, kesme kuvvetine paralel olan kenarlara olan mesafelerinden maximumu) 1.5 Ca1 den küçükse Ref --> ACI 318-19 - 17.7.2.1.2     
            if ha < 1.5*Ca1 and max(Ca_x1,Ca_x2) < 1.5*Ca1:
                    Ca1 = max(max(Ca_x1,Ca_x2)/1.5, ha/1.5, Anc_x/3) 

        if ShearDirection == 3 or ShearDirection == 4:
            if ShearDirection == 3:
                    Ca1 = Ca_x1
                    if IsWeldCommonPlate:
                        Ca1 = Ca_x1 + Anc_y
                    
            if ShearDirection ==4:
                    Ca1 = Ca_x2
                    if IsWeldCommonPlate:
                        Ca1 = Ca_x2 + Anc_y

            if ha < 1.5*Ca1 and max(Ca_y1,Ca_y2) < 1.5*Ca1:
                    Ca1 = max(max(Ca_y1,Ca_y2)/1.5, ha/1.5, Anc_y/3)

        return round(Ca1,3)
    
    def A_Vco(self, C_a1 : float)-> float:
        """ACI318-19 - Fig.R17.7.2.1a

        Args:
            C_a1 (float): _description_

        Returns:
            float: _description_
        """
        A_Vco = 4.5*C_a1**2
        return round(A_Vco,2)

    def A_Vc(self, 
             Anc_x : float, 
             Anc_y : float, 
             ha    : float, 
             Ca_x1 : float, 
             Ca_x2 : float, 
             Ca_y1 : float, 
             Ca_y2 : float, 
             ShearDirection : int, 
             IsWeldCommonPlate : bool) -> float:

        """ 
        Kesme kuvvetine dik kenarda oluşan efektif kesme alanını hesaplar.

        
                    |
                    |
                    |
                    |                    
                 ___|___
                 \  |  /  V
                  \ | /
                   \|/ 
                   
        |    |               |     |  
         Cax1      Ancx       Cax2
         ___________________________  ____
        |\                        /|
        | \                      / |
        |  \                    /  |  Cay2
        |   \                  /   |
        |    \                /    |  
        |     O------O-------O     |  ----               |-X                    
        |     |              |     |                     |              
        |     |              |     |                     |              
        |     O              O     |  Ancy        +Y_____|______-Y              
        |     |              |     |                     |                  
        |     |              |     |                     |                  
        |     O------O-------O     |  ----               |                      
        |    /                \    |                     |+X
        |   /                  \   |  Cay1
        |  /                    \  |    
        | /                      \ |
        |/________________________\|  ____ 
                    

            Args:
                Anc_x (float)       : X doğrultusunda dizilmiş iki uçtaki ankrajların merkezleri arasındaki mesafe
                Anc_y (float)       : Y doğrultusunda dizilmiş iki uçtaki ankrajların merkezleri arasındaki mesafe
                ha (float)          : Beton kalınlığı
                Ca_x1 (float)       : Y doğrultusunda sol kenara olan mesafe
                Ca_x2 (float)       : Y doğrultusunda sağ kenara olan mesafe
                Ca_y1 (float)       : Y doğrultusunda sol kenara olan mesafe
                Ca_y2 (float)       : Y doğrultusunda sağ kenara olan mesafe
                ShearDirection (int): Kesme kuvveti yönü
                                    1 -> +X dir
                                    2 -> -X dir
                                    3 -> +Y dir
                                    4 -> -Y dir
                IsWeldCommonPlate   : Ankrajlar ortak bir plakaya kaynaklı mı? 

            Returns:
                float: _description_
        """
        Ca1 = self.Get_Ca1(Anc_x,Anc_y,ha,Ca_x1,Ca_x2,Ca_y1,Ca_y2,ShearDirection,IsWeldCommonPlate)
        
        if ShearDirection == 1 or ShearDirection == 2:
            
            if ha > 1.5*Ca1:
                ha = 1.5*Ca1 # ha 1.5 Ca1 den büyükse sınırlandırılır.             
            
            if Ca_x2 > 1.5*Ca1:
                Ca_x2 = 1.5*Ca1
                
            if Ca_x1 > 1.5*Ca1:
                Ca_x1 = 1.5*Ca1

            # Avc = (Ankraj grubunun sol ucundan kuvvete paralel sol kenara olan mesafesi + Ankraj grubunun sol ucundan sağ ucuna olan mesafe + Ankraj grubunun sağ ucundan kuvvete paralel sağ kenara olan mesafesi) * ha
            Avc = (Ca_x1 + Anc_x + Ca_x2) * ha

        elif ShearDirection == 3 or ShearDirection == 4:
            if ha > 1.5*Ca1:
                ha = 1.5*Ca1 

            if Ca_y2 > 1.5*Ca1:
                Ca_y2 = 1.5*Ca1
                
            if Ca_y1 > 1.5*Ca1:
                Ca_y1 = 1.5*Ca1
                
            if ha > 1.5*Ca1:
                ha = 1.5*Ca1

            # Avc = (Ankraj grubunun sol ucundan kuvvete paralel sol kenara olan mesafesi + Ankraj grubunun sol ucundan sağ ucuna olan mesafe + Ankraj grubunun sağ ucundan kuvvete paralel sağ kenara olan mesafesi) * ha
            Avc = (Ca_y1 + Ca_y2 + Anc_y) * ha
                
        else:
            return Exception("ShearDirection unknow")
        
        return Avc

    def Get_le(self, AnchorType : int, h_ef : float, d_a : float) -> float:
        if AnchorType.__class__.__name__ == CastInAnchorType.__name__:
            l_e = h_ef
        if AnchorType == PostInAnchorType.SleeveType or AnchorType == PostInAnchorType.StudType:
            l_e = 2*d_a
        else:
            l_e = 8*d_a

        return l_e
    
    def V_b(self,
            AnchorType      : int, 
            l_e             : float, 
            d_a             : float, 
            f_c             : float, 
            Ca1             : float,
            t_attachment    : float,
            s_min           : float,
            Ca2             : float,
            h_ef            : float,
            eV              : float = 0.0,
            lambda_a        : float = 1.0) -> float:

        Vb1 = (7 * (l_e / d_a)**0.2 * d_a**0.5) * lambda_a * f_c**0.5 * Ca1**1.5 #ACI318-19 Eq. 17.7.2.2.1a
        Vb2 = 9 * lambda_a * f_c**0.5 * Ca1**1.5                                 #ACI318-19 Eq. 17.7.2.2.1b
        Vb3 = (8 * (l_e / d_a)**0.2 * d_a**0.5) * lambda_a * f_c**0.5 * Ca1**1.5 #ACI318-19 Eq. 17.7.2.2.1c

        Vb = min(Vb1,Vb2)

        if AnchorType.__class__.__name__ == CastInAnchorType.__name__:
            #For cast-in headed studs, headed bolts, or hooked bolts that are *continuously welded to steel attachments*, basic concrete breakout strength of a single anchor in shear in cracked concrete
            provides = True
            Vb = min(Vb2,Vb3)
            #provided that a through d are satisfied
            # (a) Steel attachment thickness is the greater of 0.5da and 3/8 in.
            # (b) Anchor spacing s is at least 2.5 in.
            # Reinforcement is provided at the corners if ca2 <= 1.5*hef
            #(d) For anchor groups, the strength is calculated based on the strength of the row of anchors farthest from the edge.
            if t_attachment < max(0.5*d_a, (3/8)*25.4):
                print(f"Steel attachments thickness must greater than {round(max(0.5*d_a, (3/8)*25.4), 2)}mm")
                provides = False
            if s_min < 2.5*25.4:
                print(f"Anchor spacing is too small, please make it larger than {round(2.5*25.4,2)}mm")
                provides = False
            if Ca2 <= 1.5*h_ef:
                print("Please provide the reinforcement at the corners" )

            if provides != True:
                Vb = 0.0
                
        if eV > 0.0:
            Vb = Vb3
            
        return Vb

    def Psi_ecV(self, eV : float, Ca1 : float) -> float:
        """Modification factor for anchor groups loaded eccentrically in shear

        Args:
            eV (float): _description_
            Ca1 (float): _description_

        Returns:
            float: _description_
        """
        Psi_ecV = 1 / (1 + (eV / 1.5*Ca1))
        if Psi_ecV > 1:
            Psi_ecV = 1.0
        return Psi_ecV

    def Psi_edV(self, Ca1 : float, Ca2 : float) -> float:
        if Ca2 >= 1.5*Ca1:
            Psi_edV = 1.0
        else:
            Psi_edV = 0.7 + (0.3 * (Ca2/(1.5*Ca1)))
        return Psi_edV

    def Psi_cV(self,ConcCrack : bool, Reinforcement : bool, ReinforcementBarNo : int, StirrupsSpaced : float) -> float:
        """Modification factor for the influence of cracking in anchor regions at service load levels and presence or absence of supplementary reinforcement

        Args:
            ConcCrack (bool): analysis indicates no cracking at service load levels True -> Yes, False -> No
            Reinforcement (bool): supplementary reinforcement use it True -> Yes, False -> No
            ReinforcementBarNo (int): _description_
            StirrupsSpaced (float): _description_

        Returns:
            float: _description_
        """
        if ConcCrack:
            Psi_cV = 1.4
        if ConcCrack != True:
            if Reinforcement != True or ReinforcementBarNo < 4:
                Psi_cV = 1.0
            if Reinforcement and ReinforcementBarNo >= 4:
                print("No 4 bardan büyük takviye donatıları ankrajlar arasındaysada bu katsayı kullanılabilir...")
                Psi_cV = 1.2
            if Reinforcement and ReinforcementBarNo >= 4 and StirrupsSpaced <= 4*25.4:
                print("No 4 bardan büyük takviye donatıları ankrajlar arasındaysada bu katsayı kullanılabilir...")
                Psi_cV = 1.4
        return Psi_cV

    def Psi_hV(self, Ca1 : float, h_a : float) -> float:
        """Modification factor for anchors located in a concrete member

        Args:
            Ca1 (float): _description_
            h_a (float): _description_

        Returns:
            float: _description_
        """
        Psi_hV = 1.0
        if h_a < 1.5*Ca1:
            Psi_hV = ((1.5*Ca1) / h_a)**0.5
        return Psi_hV
    
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
        round(V_cp,3)
        
    def GroupAnchorConcPryoutStrengthInShear(self, h_ef : float, N_cpg : float)-> float:
        k_cp = 1.0
        if h_ef >= 2.5*25.4:
            k_cp = 2.0
        V_cpg = k_cp * N_cpg

        return round(V_cpg,3)

#TENSION-SHEAR INTERACTION
def TensionShearInteractionStrength(Nua : float, Vua : float, Design_Nn : float, Design_Vn : float) -> float:
    """Tension-Shear Interaction capacity ratio
        ratioTensionShearCapacity <= 1.2 must.
    Args:
        Nua (float): Tension force single anchor in tension
        Vua (float): Shear force single anchor in shear
        Design_Nn (float): Design tension capacity conc
        Design_Vn (float): Design shear capacity conc

    Returns:
        float: ratio tension shear interaction capacity
    """
    ratioTensionCapacity = Nua / (Design_Nn)
    ratioShearCapacity = Nua / (Design_Vn)
    ratioTensionShearCapacity = ratioTensionCapacity + ratioShearCapacity

    return round(ratioTensionShearCapacity,3)
        



    
# def main() -> None:
#     # conn = AnchorConnection(d=12.7, bf= 12.2, B= 240, N= 240, tf=10, P_axial=700, M=50, fck= 3)
#     # BasePlat = CheckBasePlate(Contn=conn, ReductionFactor=0.65)
#     # print(f' A1 ={BasePlat.A1}, N ={BasePlat.N}, ')
#     # pt = CastInAnchorType(1)
#     # print(pt.__class__.__name__.split("In")[0])

#     #AISC Design Guide-01 Third Edition Example 4.7.7

#     #INFO ACI318-19 17.5.1.3.1 ankraj aralığı 3*h_ef ten küçük olduğu için tüm ankrajlar grup olarak hareket edileceği düşünülüyor ve n=8 kabul ediliyor.
#     # Inputs 
#     """
#         1 in   = 25.4 mm
#         1 psi  = 0.00689476 N/mm2
#         1 kips = 4448.22 N
#         1 ksi  = 6.8947 MPa
#     """
#     # BASE PLATE
#     Fpy = 50 * 6.8947
#     Fpu = 60 * 6.8947

#     # ANCHOR RODS
#     Fay = 36 * 6.8947
#     Fau = 58 * 6.8947
#     d_a = 1.5 * 25.4

#     # OTHER VARIABLES
#     Anc_x = 6.5 *25.4
#     Anc_y = 28 *25.4
#     Cax1  = 100 *25.4
#     Cax2  = 100 *25.4
#     Cay1  = 100 *25.4
#     Cay2  = 100 *25.4
#     h_ef  = 2*12 *25.4
#     f_c   = 5000 *0.00689476 
#     n = 8
#     nt = 6 #Design Guide Table-4.1
#     fi_breakout  = 0.7  #Because no supplementary reinforcement was specified, ACI 318, Table 17.5.3(b)
#     fi_pullout   = 0.7  #ACI 318, Table 17.5.3(c)
#     fi_steelshear = 0.65 # ACI 318, Table 17.5.3(a)
#     fi_pryoutshear = 0.7 # ACI 318, Table 17.5.3(c)
    
#     Abrg  = 3.12 * 25.4**2 #Design Guide Table-4.2

#     # FORCES
#     Vu = 116 * 4448.22
#     Nu = 138 * 4448.22
#     # Standart delikler ve kesmede etkin olabilecek şekilde ankrajlar beton kenarına yakın olmadığından kesme kuvveti ankrajlar arasında eşit dağıldığı kabul edilir. Aynı şekilde çekme durumunda da bir eksantrisite olmadığı için çekme kuvvetlerinin ankrajlara eşit dağıldığı kabul edilir.
#     Vua_i = Vu/n
#     Nua_i = Nu/n
    
#     print("CALCULATE BASIC PARAMETERS")
#     print(80*"=")
#     h_ef = h_ef_Check(h_ef,Cax1,Cax2,Cay1,Cay2,Anc_x,Anc_y)
#     print(f"h_ef = {round(h_ef/10**3,3)}m - ({round(h_ef/(25.4),3)}in)")
    
#     ANco = A_Nco(h_ef)
#     print(f"ANco = {round(ANco/10**6,3)}m2 - ({round(ANco/25.4**2,3)}in^2)")
    
#     ANc = A_Nc(Anc_x, Anc_y, h_ef, Cax1, Cax2, Cay1, Cay2,n)
#     print(f"ANc = {round(ANc/10**6,3)}m2 - ({round(ANc/25.4**2,3)}in^2)")
#     print(80*"=")

#     print(80*"*")
#     print("TENSION CHECKS")
#     print(80*"=")
    
#     print("CONCRETE BREAKOUT STRENGTH IN TENSION")
#     print(80*"=")
    
#     BreakoutTension = ConcreteBreakoutStrengthOfAnchorInTension()
    
#     Nb = BreakoutTension.BasicSingleAnchorBreakoutStrength(kc=16, lambda_a=1, f_c=f_c, h_ef=h_ef)
#     print(f"Nb = {round(Nb/10**3,3)}kN - ({round(Nb/4448.22,3)}kips)")
    
#     psi_ecN = BreakoutTension.Psi_ecN_Get(0.0, h_ef)
#     print(f"psi_ecN = {psi_ecN}")
    
#     psi_edN = BreakoutTension.Psi_edN_Get(Ca_min=Cax1, h_ef=h_ef)
#     print(f"psi_edN = {psi_edN}")
    
#     psi_cN = BreakoutTension.Psi_cN_Get(InstalledType=AnchorInstalledType.Cast_in, IsConcCracked = None)
#     print(f"psi_cN = {psi_cN}")
    
#     psi_cpN = BreakoutTension.Psi_cpN_Get(Cac=Cax1 , Ca_min=Cax1 , h_ef=h_ef, InstalledType=AnchorInstalledType.Cast_in)
#     print(f"psi_cpN = {psi_cpN}")

#     Ncbg = BreakoutTension.ForGroupAnchor(A_Nc=ANc, A_Nco=ANco, Psi_ecN=psi_ecN, Psi_edN=psi_edN, Psi_cN=psi_cN, Psi_cpN=psi_cpN, Nb=Nb)
#     print(f"Ncbg = {round(Ncbg/10**3,3)}kN - ({round(Ncbg/4448.2,3)}kips)")

#     Design_Ncbg = fi_breakout * Ncbg
#     print(f"fi_d*Ncbg = {round(Design_Ncbg/10**3,3)}kN - ({round(Design_Ncbg/4448.2,3)}kips)")

#     print(f"Is concrete breakout strength in tension capacity enough? ==> {Nu < Design_Ncbg} - {round(Nu/1000,3)}kN < {round(Design_Ncbg/1000,3)}kN - ({round(Nu/4448.2,3)}kips < {round(Design_Ncbg/4448.2,3)}kips)")

#     print(80*"=")

#     print("CONCRETE PULLOUT STRENGTH IN TENSION")
#     print(80*"=")
#     PulloutTension = PulloutStrengthInTension()
    
#     psi_cp = PulloutTension.Psi_cP_Get(IsConcCracked=None)
#     print(f"psi_cp = {psi_cp}")
    
#     Np     = PulloutTension.BasicSingleAnchorPulloutStrength(A_brg=Abrg, 
#                                                              f_c=f_c, 
#                                                              e_h=0.0, 
#                                                              d_a=d_a, 
#                                                              AnchorType= CastInAnchorType.HexHeadBoltWithWasher, 
#                                                              InstalledType= AnchorInstalledType.Cast_in)
#     print(f"Np = {round(Np/10**3,3)}kN - ({round(Np/4448.2,3)}kips)")
    
#     Npn = PulloutTension.SingleAnchorNominalPulloutStrength(Psi_cP=psi_cp, Np=Np)
#     print(f"Npn = {round(Npn/10**3,3)}kN - ({round(Npn/4448.2,3)}kips)")

#     Design_Npn = fi_pullout * Npn
#     print(f"fi_d*Ncbg = {round(Design_Npn/10**3,3)}kN - ({round(Design_Npn/4448.2,3)}kips)")

#     print(f"Is concrete pullout strength in tension capacity enough? ==> {Nua_i < Design_Npn} - {Nua_i/1000}kN < {Design_Npn/1000}kN - ({round(Nua_i/4448.2,3)}kips < {round(Design_Npn/4448.2,3)}kips)")
#     print(80*"=")

#     print("CONCRETE SIDE-FACE BLOWOUT STRENGTH IN TENSION")
#     print(80*"=")
#     print("Determine the concrete side-face blowout strength in tension of the anchor group The available concrete side-face blowout strength in tension of the anchor group is determined according to ACI 318, Section 17.6.4. Because there are no cases with anchor rods close to an edge (hef > 2.5ca1), side-face blowout is not applicable.")
#     print(80*"=")

#     print(80*"*")
#     print("SHEAR CHECKS")
#     print(80*"=")
    
#     print("AVAILABLE STEEL STRENGTH IN SHEAR")
#     print(80*"=")

#     steelStrengthShear = SteelStrengthOfAnchorInShear()
    
#     Asev = steelStrengthShear.A_seV(d_a=d_a, n_t=nt)
#     print(f"A_seV = {round(Asev,3)}mm2 - ({round(Asev/25.4**2,3)}in^2)")

#     # Design guide örneğinde efektif ankrajın çekme alanı alınmış(DG01-Table-4.1 A_seN) ACI da önerilen formül nominal alana yakın sonuç veriyor.(ACI318-19 - R17.7.1.2)
#     Vsa = steelStrengthShear.NominalSteelStrengthCastInHeadedAndHookedBoltAnchorInShear(A_seV=1.41*25.4**2, f_ya=Fay, f_uta= Fau)
#     print(f"Vsa = {round(Vsa/10**3,3)}kN - ({round(Vsa/4448.2,3)}kips)")

#     Design_Vsa = fi_steelshear * Vsa
#     print(f"fi*Vsa = {round(Design_Vsa/10**3,3)}kN - ({round(Design_Vsa/4448.2,3)}kips)")
    
#     print(f"Is available steel strength in shear capacity enough? ==> {Vua_i < Design_Vsa} - {Vua_i/1000}kN < {Design_Vsa/1000}kN - ({round(Vua_i/4448.2,3)}kips < {round(Design_Vsa/4448.2,3)}kips)")
#     print(80*"=")

#     print("CONCRETE BREAKOUT STRENGTH IN SHEAR")
#     print(80*"=")
#     BreakoutShear = ConcreteBreakoutStrengthOfAnchorInShear()
#     print("Because there are no edges adjacent to the base connection, the available concrete breakout strength in shear is not applicable.")

#     print("CONCRETE PRYOUT STRENGTH IN SHEAR")
#     print(80*"=")
#     PryShear = ConcretePryoutStrengthOfAnchorInShear()

#     Ncpg = Ncbg
    
#     Vcpg = PryShear.GroupAnchorConcPryoutStrengthInShear(h_ef=h_ef, Ncpg=Ncpg)
#     print(f"Vcpg = {round(Vcpg/10**3,3)}kN - ({round(Vcpg/4448.2,3)}kips)")

#     Design_Vcpg = fi_pryoutshear * Vcpg
#     print(f"fi*Vsa = {round(Design_Vcpg/10**3,3)}kN - ({round(Design_Vcpg/4448.2,3)}kips)")
    
#     print(f"Is available pryout strength in shear capacity enough? ==> {Vu < Design_Vcpg} - {Vu/1000}kN < {Design_Vcpg/1000}kN - ({round(Vu/4448.2,3)}kips < {round(Design_Vcpg/4448.2,3)}kips)")
#     print(80*"=")
    
# if __name__ == "__main__":
#     # main()
#     cbot = ConcreteBreakoutStrengthOfAnchorInShear()
#     Avc = cbot.A_Vc(Anc_x = 12, Anc_y = 12, ha = 15000, Ca_x1 = 12, Ca_x2 = 12, Ca_y1 = 12, Ca_y2 = 12, ShearDirection = 1, IsWeldCommonPlate = False)
#     print(Avc)