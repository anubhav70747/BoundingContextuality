This folder contains methods for evaluating bounds on tripartite parity oblivious random access codes (PORAC/MPORAC) (Octave/Matlab + Yalmip + sdpt3/sedumi/mosek):

-**parityConSpekkens.m**: Spekkens' parity calculations
-**parityConEven.m**: even parity calculations
- **nto1dpoRACBE_OC.m**: LP for evaluating bounds pertaining to contextual polytope for oblivious communication games. 
- **nto1dpoRACBE_Q1.m**: level = 1 NPA-like SDP for oblivious communication games. 
- **nto1dpoRACBE_Q1BE.m**: level = 1 + BE NPA-like SDP for oblivious communication games.
- **nto1dpoRACBE_Q.m**: See-Saw SDP for oblivious communication games.
- **SuccessBE.m**: An auxiliary function for nto1dpoRACBE_Q.m for calculating success probability given tasks parameters, states and measurements.
- **BinEntropy.m**: Calculate binary entropy.
- **keyRate.m**: Calculate lower bound on the key-rate.










