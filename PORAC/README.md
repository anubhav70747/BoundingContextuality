This folder contains methods for evaluating bounds on parity oblivious random access codes with measurement equivalences (PORAC/MPORAC) (Octave/Matlab + Yalmip + sdpt3/sedumi/mosek):

-**parityConSpekkens.m**: Spekkens' parity calculations
-**parityConEven.m**: even parity calculations
- **nto1dpoRAC_OC.m**: LP for evaluating bounds pertaining to contextual polytope for oblivious communication games. 
- **nto1dpoRAC_Q1.m**: level = 1 NPA-like SDP for oblivious communication games with unitaries. 
- **nto1dpoRAC_UQ1.m**: level = 1 NPA-like SDP for oblivious communication games with projectors. 
- **nto1dpoRAC_Q2.m**: level = 2 NPA-like SDP for oblivious communication games with projectors.
- **nto1dpoRAC_UQ2_beta.m**: level = 2 NPA-like SDP for oblivious communication games with unitaries. 
- **nto1dpoRAC_Q.m**: See-Saw SDP for oblivious communication games.
- **Success.m**: An auxiliary function for nto1dpoRAC_Q.m for calculating success probability given tasks parameters, states and measurements.
- **nto1dpoRAC_NC.m**: See-Saw LP for calculation NC success probability of oblivious communication games.
- **SuccC.m**: An auxiliary function for nto1dpoRAC_NC.m for calculating success probability.










