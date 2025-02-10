import pytest
from AnchorToConcrete import ConcreteBreakoutStrengthOfAnchorInShear


@pytest.fixture(scope="class")
def CBS():
    a = ConcreteBreakoutStrengthOfAnchorInShear()
    return a

# def test_Ca1(CBS) -> bool:
#     assert CBS.Get_Ca1(Anc_x = 6, 
#                        Anc_y = 24, 
#                        ha = ha, 
#                        Ca_x1=Cax1, 
#                        Ca_x2=Cax2, 
#                        Ca_y1=Cay1, 
#                        Ca_y2=Cay2, 
#                        ShearDirection=ShearDirection, 
#                        IsWeldCommonPlate=IsWeldCommonPlate) == 12.0

def test_Avco(CBS) -> bool:
    assert CBS.A_Vco(C_a1 = 12) == 648.0
    
def test_Avc(CBS) -> bool:
    assert CBS.A_Vc(4, 4, 10000, 50, 50, 12, 12, 1, False) == 720.0

