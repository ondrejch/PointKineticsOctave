#{
Point Kinetics Equation (PKE) solver for education and research purposes.

Adapted from https://github.com/ondrejch/MSBR-ORNL-4528/tree/master/dynamic_model/Octave

Author: Vikram Singh, viikraam@gmail.com and Ondrej Chvala, ochvala@utk.edu

This program takes an input file './reactivity.dat' with time (in s), reactivity,
and external neutron source as three columns, formatted as
< time reactivity ext-source >
per line. The code plots for n(t) and C_i(t) for delayed neutron groups i=1,2,...6
#}
################################################################################

clear all;        # clear workspace
% pkg load odepkg;  # load odepkg, not needed in modern Octave

# Delayed neutron lifetimes 'beta' for the six precursor groups depending on
# fuel type. B = sum(bet). Retrieved from ORNL-MSR-67-102

#  global bet = [0.00023,0.00079,0.00067,0.00073,0.00013,0.00009]';  # U233 fuel
global bet = [0.000215,0.00142,0.00127,0.00257,0.00075,0.00027]';   # U235 fuel
global B   = sum(bet);

# Decay constants for the corresponding decay groups
global lam = [0.0126,0.0337,0.139,0.325,1.13,2.50]';

# Mean neutron generation time
global L = 0.0005;

# Initial values for n(t) and C_i(t)
nt    = 1.0;
Ct(1) = (bet(1)/(L*lam(1)))*nt;
Ct(2) = (bet(2)/(L*lam(2)))*nt;
Ct(3) = (bet(3)/(L*lam(3)))*nt;
Ct(4) = (bet(4)/(L*lam(4)))*nt;
Ct(5) = (bet(5)/(L*lam(5)))*nt;
Ct(6) = (bet(6)/(L*lam(6)))*nt;

# Read input file. Formatted as < time reactivity ext-source > per line.
global input_data = dlmread('./reactivity.dat');
global nrows      = rows(input_data);
global tmax       = input_data(nrows,1); # length of time for which to evaluate the equations

# Initial y[] and t values
y0 = [nt,Ct(1),Ct(2),Ct(3),Ct(4),Ct(5),Ct(6)]';
t0 = 0;

# Get reactivity value from input file for some t.
# Not general purpose. Only compatible with this project.
# Could be improved by using interpolation tables.
function rho=react(t)
  rho=0;
  global nrows;
  global input_data;
  if (t > input_data(nrows,1))
    rho = input_data(nrows,2);
  else
    for i = 1:nrows-1
      if (t>=input_data(i,1) & t<=input_data(i+1,1))
        rho = input_data(i,2);
        break
        else
          continue
      endif
    endfor
  endif
endfunction

function S=source(t)
  S=0;
  global nrows;
  global input_data;
  if (t > input_data(nrows,1))
    S = input_data(nrows,3);
  else
    for i = 1:nrows-1
      if (t>=input_data(i,1) & t<=input_data(i+1,1))
        S=input_data(i,3);
        break
        else
          continue
      endif
    endfor
  endif
endfunction

# Time-Derivative Column [dn/dt, dC_1/dt, dC_2/dt, dC_3/dt, ..]
function ndot=neudens(t,y,react,source,bet,B,lam,L)
  ndot(1) = source(t) + (((react(t)-B)/L)*y(1)) + (lam(1)*y(2)) + ...
  (lam(2)*y(3)) + (lam(3)*y(4)) + (lam(4)*y(5)) + (lam(5)*y(6)) + (lam(6)*y(7));
  ndot(2) = ((bet(1)/L)*y(1)) - (lam(1)*y(2));
  ndot(3) = ((bet(2)/L)*y(1)) - (lam(2)*y(3));
  ndot(4) = ((bet(3)/L)*y(1)) - (lam(3)*y(4));
  ndot(5) = ((bet(4)/L)*y(1)) - (lam(4)*y(5));
  ndot(6) = ((bet(5)/L)*y(1)) - (lam(5)*y(6));
  ndot(7) = ((bet(6)/L)*y(1)) - (lam(6)*y(7));
endfunction


# ODE solution stored in a matrix of 7 x (tmax*dt)
# Solver setup
vopt = odeset("RelTol", 1e-5, "AbsTol", 1e-5, "NormControl","on", "InitialStep",0.001, "MaxStep",0.02);
# Solver call
sol = ode45(@(t,y) neudens(t,y,@react,@source,bet,B,lam,L),[t0 tmax],y0,vopt);

# Extract x and y into separate vectors
tsol = sol.x; ysol = sol.y;

# Plot Figure 1, n(t) vs t.
clf('reset');
figure(1);
F1 = plot(tsol,ysol(1,:));
X1 = xlabel('time (in s)');
set(X1,'FontName','Times New Roman','fontsize',14);
axis([t0,tmax]);
Y1 = ylabel('Reactor Power (rel.)');
set(Y1,'FontName','Times New Roman','fontsize',14);

# Plot Figure 2, C_i(t) vs t.
figure(2)
F2 = plot(tsol,ysol(2:7,:));
X2 = xlabel('time (in s)');
set(X2,'FontName','Times New Roman','fontsize',14);
axis([t0,tmax]);
Y2 = ylabel('C_i(t)');
set(Y2,'FontName','Times New Roman','fontsize',14);
legend('C_1','C_2','C_3','C_4','C_5','C_6');
