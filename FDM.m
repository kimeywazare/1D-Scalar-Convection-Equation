% 1D Scalar Convection Equation
% du/dt + a*du/dx = 0 
clc ;close all; clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRE-PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Input Variables (Space and Time)
a = 1;                              % Velocity (1, 0.75, 0.5 & 0.25)
del_x = 0.005;                      % Grid spacing
del_t = 0.005;                      % Time step
CFL = (a*del_t)/del_x;              % Courant-Friedrichs-Lewy (CFL)
T0 = cputime;                       % Initial Computational time
t = 0;                              % Initial Time
t_end = (1/a);                      % Final Time
Nsteps = t_end/del_t;               % Number of time steps
%Nsteps = 200;
N = 1/del_x + 1;               % Number of nodes (length of domain/del_x)
%N = 1001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain Discretization
x = linspace(0,1,N);                 % Value of x (x = linspace(0,0.05,N))
u_n = zeros(N,1);

% Initial Conditions
for i = 1 : 0.05/del_x
    u_n(i) = 0.5*(1-cos(40*pi*x(i)));   % Input Cosine Equation
end   
%plot(x,u_n)                            % To check equation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN (SIMULATION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Select a scheme \nUpwind Scheme: 1\nLax-Friedrichs Scheme: 2\n')   
Choice = input('Scheme = ');   % Select a scheme to run
u = zeros(N,1);

for n = 1:Nsteps

    c = del_t/del_x;                  % Simplified Variable for equation
    u(1) = 0;                          % Boundary Conditions

%Implementing Explicit Schemes

   if Choice==1                          % Upwind Scheme
      for j = 2:N
          u(j) = u_n(j) - a*c*(u_n(j) - u_n(j-1));
      end
   elseif Choice==2                   % Lax Friedrich Scheme (Lax Method)    
       for j = 2:N-1
           u(j) = (0.5*(u_n(j+1) + u_n(j-1))) - (0.5*a*c)*(u_n(j+1)...
           - u_n(j-1)); 
       end
   end
   
%Update t & u for next step 
  t = t + del_t;                          % Addidtion of time
  u_n = u;                                % Omits boundary conditions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POST-PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots

%Plot for time domain
figure(1)
plot(x,u,'r-','linewidth',2);%axis('manual'); 
axis([0,1,-0.5,1]); grid on
xlabel('Chord Length','fontsize',18); 
ylabel('Scalar Variable','fontsize',18);
title(['Time = ',num2str(t),' seconds'],'fontsize',18); shg;
lgnd=legend('Wave');
set(lgnd,'Fontsize',18)
pause(0.1) 

%Plot for capturing wave at time instances
  if n==1
      figure(2)
      plot(x,u,'r-','linewidth',2);
      hold on
      axis([0,1,-0.5,1]); grid on
      xlabel('Chord Length','fontsize',14); 
      ylabel('Scalar Variable','fontsize',14);
  elseif n==round(Nsteps*0.25)
      figure(2)
      plot(x,u,'b-','linewidth',2);
      axis([0,1,-0.5,1]); grid on
      xlabel('Chord Length','fontsize',14); 
      ylabel('Scalar Variable','fontsize',14);
  elseif n==round(Nsteps*0.50)
      figure(2)
      plot(x,u,'g-','linewidth',2);
      axis([0,1,-0.5,1]); grid on
      xlabel('Chord Length','fontsize',14); 
      ylabel('Scalar Variable','fontsize',14);
  elseif n==round(Nsteps*0.75)
      figure(2)
      plot(x,u,'y-','linewidth',2);
      axis([0,1,-0.5,1]); grid on
      xlabel('Chord Length','fontsize',14); 
      ylabel('Scalar Variable','fontsize',14);
      lgnd=legend('t = 0','t = 0.25','t = 0.5','t = 0.75');
      set(lgnd,'Fontsize',14)
      title('Wave at various Time Instances (CFL = 1)','fontsize',14); 
      shg;
  end 
    
end
Tcomp = cputime - T0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%