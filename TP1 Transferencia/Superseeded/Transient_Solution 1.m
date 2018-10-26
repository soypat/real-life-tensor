%% 1-D FVM code. Bar with conduction, heat generation and convection
close all
Nx=50; % Number of volumes
L=1; %metro
x = linspace(0,L,Nx+1); % Surface positions for volumes
dx = L/Nx;
xmid = 0.5*(x(1:Nx) + x(2:Nx+1)); %Volume midpoints

%Let's solve a transient heat conduction 1-D problem.
%simple!
%dT/dt - alpha*div(T) = 0
%For our problem. 1D: dT/dt=alpha*(d^2T/dx^2)
%alpha=k/(rho*cp)

% Initial conditions: Temperature along bar (T=25 constant)
% Border conditions: T(1)=0
%                    T(Nx-1)= T(Nx) (Adiabatic condition)
%                    
%% Bar dimensions
A=0.1; %m^2 Area
% L=L; %Meters, declared above
D=sqrt(4*A/pi);
%% Datum
qgen=10; %Watt per meter cubed
cp=2; % Joule per kilogram kelvin (steel)
k=1;
hamb=10; %Watt per meter squared kelvin
Tamb=30;
rho=2000; % kilogram per meter cubed
alpha=k/rho/cp;
%% Solver
T=25*ones(Nx,1);
permissible_Error=1e-4; %How close to the exact solution do we want to be?
Tlast=T;
aux=solucion_analitica1(x);
Texact=aux(2:end)';
% T=Texact
T0=0;%Border condition
TN=30;
Tbc = [T0; T ; T(end)]; %We have two border conditions
dt=.5; %seconds
tfinal=dt*200; %How many iterations do we want?
t=0;%Initialize start time 
%{
T0     T(1)                T(Nx)     TN
|-------|------- ... -------|-------|
     <dT[1]>             <dT[Nx]>


BAR VOLUME
            NORTH(Convection)
            |-----------|
            |           |
EAST (Cond.)|   Qgen    | WEST (Conduction)
            |           |
            |           |
            |-----------|
%}
vol=dx*A;
dS=pi*D*dx; %Surface differential (of bar exterior for convection calculation)
Q=zeros(Nx,1);
while (t<tfinal) || Tdifference>permissible_Error
    Tbc = [T0; T ; T(end)]; %
    for i = 1:Nx
        dTdxwest = Tbc(i+1)-Tbc(i);
        dTdxeast = Tbc(i+1)-Tbc(i+2);
        Qwesteast= -k*A*(dTdxwest+dTdxeast)/dx;%Frame of reference taken on the volume
        
        % Generated
        Qgen=qgen*A*dx;
        % Convection
        Qnorth=hamb*pi*D*dx*(Tamb-Tbc(i+1));
        
        Q(i)=Qwesteast+Qgen+Qnorth;
    end
    t=t+dt;
    T=T+dt*Q/(rho*cp*vol);
    %Error control
    Tdifference=max(T-Tlast);
    Tlast=T;
    if max(T)>35%Runaway temperature
        rho=rho+10;
        t=0;
        T=25*ones(Nx,1);
        cla
        continue
    end
    %Plot current solution
     stairs(x-dx/2,[T0;T]);
     axis([0, 1, min(T0,TN), max(T0,TN)]);
     grid on;
     drawnow;
end
hold on
plot(x,solucion_analitica1(x))

