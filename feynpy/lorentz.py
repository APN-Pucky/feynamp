from sympy.physics.hep.gamma_matrices import LorentzIndex
from sympy.tensor.tensor import TensorHead, TensorIndexType

# LorentzIndex = TensorIndexType('Lorentz', dummy_name='L')

DiracSpinIndex = TensorIndexType("DiracSpin", dim=4, dummy_name="DS")
SpinIndex = TensorIndexType("Spin", dummy_name="S")

Gamma = TensorHead("gamma", [LorentzIndex, SpinIndex, SpinIndex])
ProjP = TensorHead("ProjP", [SpinIndex, SpinIndex])
ProjM = TensorHead("ProjM", [SpinIndex, SpinIndex])

Eps = TensorHead("Eps", [LorentzIndex, SpinIndex])
Eps_star = TensorHead("Eps_star", [LorentzIndex, SpinIndex])
u = TensorHead("u", [SpinIndex])
u_bar = TensorHead("u_bar", [SpinIndex])
v = TensorHead("v", [SpinIndex])
v_bar = TensorHead("v_bar", [SpinIndex])
