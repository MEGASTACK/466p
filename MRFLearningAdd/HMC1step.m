function [x,A,U_x,dUdx] = HMCgrid1step(x,Nleapfrog,leapfrogstep,U_x,dUdx,D)

% This function implements the hybrid monte carlo sampling routine specifically to use with
% HMCwrapper.m for grid posteriors. For a general version see the original HMC.m.org 
%
% INPUT : x    =   DXN matrix of postive data vectors where D is dimensionality of the
%                  feature space and N the number of samples. They are the starting
%                  points for the N markov chains.
% Nleapfrog    =   number of leapfrog steps.
% leapfrogstep =   size of the leapfrogstep.
% U_x          =   1XN energy at input x for each chain
% dUdx         =   DXN derivative of energy at x for each chain
% D            =   data to compute likelihood term of energy (NOTE: this code has been hardcoded for posterior estimation) 
%
% OUTPUT:  x   =   new matrix of negative data obtained by running the HMC procedure.
%          A   =   acceptance vector. The average is taken over resampling steps
%          U_x =   1XN energy at output x for each chain
%          dUdx=   DXN derivative of energy at output x for each chain
%                  
% SUBROUTINES:
%  Energy(x)   =   Converts DXN matrix of data 'x' into 1XN array of energies.
% DEnergy(x)   =   Converts DXN matrix of data 'x' into DXN array derivatives of energies.
%                  These functions define the model from which we are sampling.
%
% AUTHOR: MODIFIED - Sridevi Parise


[D,N] = size(x);


t = x;
p = randn(D,N);   % initialize momenta
K_x = sum(p.^2)/2;  % compute kinetic energy
Dx = dUdx;

for i=1:Nleapfrog   % loop over leapfrog steps
    p = p - 0.5*leapfrogstep * Dx;
    t = t + leapfrogstep * p;
    Dt =  DEnergy(t);  % use the derivative of your own energy function here
    p = p - 0.5*leapfrogstep * Dt;
    Dx = Dt;
end

U_t = Energy(t); % use your own energy function here
K_t = sum(p.^2)/2;

DH = U_x + K_x - U_t - K_t; % compute difference total energy

P = min(1,exp(DH));
A = rand(1,N) < P;  % accept according to metroplis-hastings criterium

x(:,A) = t(:,A);

U_x(:,A) = U_t(:,A);
dUdx(:,A) = Dt(:,A);
