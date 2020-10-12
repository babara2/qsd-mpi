// Two level atom coupled with vaccum, decaying

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "ACG.h"
#include "Traject.h"
#include "State.h"
#include "Operator.h"
#include "SpinOp.h"
#include "Complex.h"

long double HALF_PI = 1.5707963267948966192313216916397514420985846996875529104874722961;

double gate1(double t)
{
  if ((t<HALF_PI) && (t>0)) return 1.0;
  return 0.0;
}

double gate2(double t)
{
  if ((t<2*HALF_PI) && (t>=HALF_PI)) return 1.0;
  return 0.0;
}

double gate3(double t)
{
  if ((t<3*HALF_PI) && (t>=2*HALF_PI)) return 1.0;
  return 0.0;
}

int main()
{
  RealFunction gate_1 = gate1;
  RealFunction gate_2 = gate2;
  RealFunction gate_3 = gate3;

  // Basic operators
  // idenity operators for DOFs 0,1,2
  IdentityOperator id_0(0);
  IdentityOperator id_1(1);
  IdentityOperator id_2(2);
  // DOF 0
  SigmaX sigma_x_0(0);
  SigmaY sigma_y_0(0);
  SigmaZ sigma_z_0(0);
  SigmaPlus sigma_plus_0(0);
  Operator sigma_minus_0 = sigma_plus_0.hc(); // hermitian conjugate is sigma_minus_0
  Operator hadamard_dof0 = (1/sqrt(2))*(sigma_z_0 + sigma_x_0);
  Operator s_gate_dof0 = sigma_minus_0*sigma_plus_0;
  // DOF 1
  SigmaX sigma_x_1(1);
  SigmaZ sigma_z_1(1);
  SigmaPlus sigma_plus_1(1);
  Operator sigma_minus_1 = sigma_plus_1.hc();
  // DOF 2
  SigmaX sigma_x_2(2);
  SigmaZ sigma_z_2(2);
  SigmaPlus sigma_plus_2(2);
  Operator sigma_minus_2 = sigma_plus_2.hc();

  // Hamiltonians
  // TODO: document and verify energy scale is correct
  double OMEGA_eg = 1.0;
  Operator Ham_id_0 = 0.0*(sigma_x_0);
  // DOF 0
  // X gate for 3 consecutive gates
  Operator Ham_X_0 = OMEGA_eg*(sigma_x_0);
  // Hadamard for 3 consecutive gates
  Operator Ham_H_0 = OMEGA_eg*(hadamard_dof0);
  // Hadamard -> X -> Hadamard
  Operator Ham_HXH_0 = OMEGA_eg*(gate3*hadamard_dof0 +
                                 gate2*sigma_x_0 +
                                 gate1*hadamard_dof0);
  // Hadamard -> Z -> Hadamard
  Operator Ham_HZH_0 = OMEGA_eg*(gate3*hadamard_dof0 +
                                 gate2*sigma_z_0 +
                                 gate1*hadamard_dof0);
    // Hadamard -> Z -> Z (full rotation around equator)
  Operator Ham_HZZ_0 = OMEGA_eg*(gate3*sigma_z_0 +
                                 gate2*sigma_z_0 +
                                 gate1*hadamard_dof0);
  // Hadamard -> Y -> Hadamard
  Operator Ham_HYH_0 = OMEGA_eg*(gate3*hadamard_dof0 +
                                 gate2*sigma_y_0 +
                                 gate1*hadamard_dof0);
  // Hadamard -> S -> Hadamard
  Operator Ham_HSH_0 = OMEGA_eg*(gate3*hadamard_dof0 +
                                 gate2*s_gate_dof0 +
                                 gate1*hadamard_dof0);
  // Hadamard -> T -> Hadamard. Note T is S applied at half speed.
  Operator Ham_HTH_0 = OMEGA_eg*(gate3*hadamard_dof0 +
                                 gate2*(0.5*s_gate_dof0) +
                                 gate1*hadamard_dof0);
  // DOF 0,1
  // CNOT, if ground, keep same, if excited, reverse
  Operator Ham_CNOT = OMEGA_eg*((sigma_minus_0*sigma_plus_0)*id_1 +
                           (sigma_plus_0*sigma_minus_0)*sigma_x_1);
  // DOF 0,1,2
  // CCNOT
  Operator q1_q2_zero = (sigma_minus_0*sigma_plus_0)*(sigma_minus_1*sigma_plus_1)*id_2;
  Operator q1_zero_q2_one = (sigma_minus_0*sigma_plus_0)*(sigma_plus_1*sigma_minus_1)*id_2;
  Operator q1_one_q2_zero = (sigma_plus_0*sigma_minus_0)*(sigma_minus_1*sigma_plus_1)*id_2;
  Operator q1_one_q2_one = (sigma_plus_0*sigma_minus_0)*(sigma_plus_1*sigma_minus_1)*sigma_x_2;
  Operator Ham_CCNOT = OMEGA_eg*(q1_q2_zero + q1_zero_q2_one +
                                 q1_one_q2_zero + q1_one_q2_one);

  // The Lindblad operators
  // 1 Qubit, DOF 0
  double GAMMA=1.0;
  // COMMENT OUT THE LINDBLADS NOT BEING USED
  // const int num_lindblad_ops = 1;
  // phase damping
  // Operator lindblad_op = GAMMA*sigma_z_0;
  // amplitude damping
  // Operator lindblad_op = GAMMA*sigma_minus_0;
  // Operator L[num_lindblad_ops] = {lindblad_op};

  // depolarizing noise
  // from: https://ocw.mit.edu/courses/nuclear-engineering/22-51-quantum-theory-of-radiation-interactions-fall-2012/lecture-notes/MIT22_51F12_Ch8.pdf
  // The depolarizing channel is a model of a decohering qubit that has
  // particularly nice symmetry properties. We can describe it by saying that,
  // with probability 1 - p the qubit remains intact, while with probability
  // p an “error” occurs. The error can be of any one of three types,
  // where each type of error is equally likely. If {|0)|1)} is an orthonormal
  // basis for the qubit, the three types of errors can be characterized as:
  // 1. Bit-flip error
  // 2. Phase-flip error
  // 3. Both errors
  // We are chosing p = 0.27
  const int num_lindblad_ops = 3;
  Operator lindblad_op_1 = GAMMA*0.3*sigma_z_0;
  Operator lindblad_op_2 = GAMMA*0.3*sigma_x_0;
  Operator lindblad_op_3 = GAMMA*0.3*sigma_y_0;
  Operator L[num_lindblad_ops] = {lindblad_op_1, lindblad_op_2, lindblad_op_3};


  // TODO: verify w/ Bibek that state is initialized up
  // COMMENT OUT ALL STATES NOT BEING USED
  // 1 Qubit State
  State psi0(2,SPIN); // 1-freedom state ("down")
  psi0 *= sigma_plus_0; // spin up
  // psi0 *= Ham_H_0; // superposition of spin up and down

  // 2 Qubit States
  // State psi0(2,SPIN);
  // State psi1(2,SPIN);
  // State psilist[2] = {psi0,psi1};
  // State psi(2,psilist); // Product state (down, down)

  // 3 Qubit States
  // State psi0(2,SPIN);
  // State psi1(2,SPIN);
  // State psi2(2,SPIN);
  // State psilist[3] = {psi0,psi1, psi2};       
  // State psi(3,psilist); // Product state (down, down, down)
  // psi *= sigma_plus_2; // (down, down, up)
  // psi *= sigma_plus_0; // (up, down, down)
  // psi *= sigma_plus_1; // (down, up, down)
  // psi *= sigma_plus_0*sigma_plus_1; // (up, up, down)
  // psi *= sigma_plus_0*sigma_plus_2; // (up, down, up)
  // psi *= sigma_plus_1*sigma_plus_2; // (down, up, up)
  // psi *= sigma_plus_0*sigma_plus_1*sigma_plus_2; // (up, up, up)

  // The random number generator
  int seed = 74298;
  ACG gen(seed,55);
  ComplexNormal rand1(&gen);

  // Stepsize and integration time  
  double dt=0.001; // basic time step
  int numdts=10; // time interval between outputs = numdts*dt
  int num_gates = 3; // SPECIFY THE NUMBER OF GATES TO SIMULATE
  int numsteps=std::round(num_gates*(HALF_PI)/(dt*numdts)); // total integration time = numsteps*numdts*dt

  double accuracy = 0.000001;

  // Step function:
  // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
  // stochastic part: fixed stepsize Euler

  // Plug the relavent:
  // 1. State
  // 2. Hamiltonian
  // 3. Num of Lindbladian Ops
  // 4. Lindbladians
  AdaptiveStep theStepper(psi0, Ham_id_0, num_lindblad_ops, L, accuracy);

  // Output
  // Each Operator specified to be output would be written to
  // a seperate file.
  const int nOfOut = 4;
  // operators to ouput expectations and variances for
  Operator outlist[nOfOut] = {(id_0+sigma_z_0)*0.5, sigma_x_0, sigma_y_0, sigma_z_0};
  // output files list
  const char *flist[nOfOut] = {"sigma_z.out", "x", "y", "z"}; // Output files
  if (GAMMA > 0) {
    flist[0] = "sigma_z_open.out";
    flist[1] = "x_open";
    flist[2] = "y_open";
    flist[3] = "z_open";
  }
  // Each operator has a set of 4 values associated with it
  // expectation (real, imag) and variances (real, imag).
  // For our experiments we are concerned with real
  // expectations and variances.
  int pipe[8] = {1, 3, 5, 7, 9, 11, 13, 15};
  int numTraj = GAMMA > 0 ? 1000 : 1;
  // Run the trajectory algorithm
  Trajectory theTraject(psi0, dt, theStepper, &rand1);
  // Last parameter specifies the number of trajectories to compute
  theTraject.sumExp(nOfOut, outlist, flist, numdts, numsteps, numTraj);
}
