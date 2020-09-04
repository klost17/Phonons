%% Specific heat calculator

% In this code we compute the specific heat for different temperatures
% We will be using silicon atoms for our chain (we take the distance between
% the nearest neighbour atoms in a silicon crystal as our lattice parameter a)
m=28*1.660538921*10^(-27);%mass of an atom (kg)
K=59.7939;%Hook's constant (kg/s^2)
a=5.43*10^-10;%distance between atoms (m)
Kb=1.38064852*10^(-23);%Boltzmann constant
hbar=1.054571800*10^(-34);%Dirac constant
wo=sqrt(K/m);
N=20; %Number of atoms in the chain
% Computing the possible values of p, depending on the
% parity of N, in order to impose periodicity and to work in the 1st
% Brillouin Zone
if mod(N,2)==0 %N is even
    p=-N/2+[0:1:N-1];
else           %N is odd
    p=-N/2-1/2+[1:1:N];
end
Cvec=[]; Tvec=[];Cvec2=[];
%Computing the specific heat for different values of the temperature
for T=[1:1:1500]
C=0;
%Summing the contribution to the specific heat of each mode (each one fixed
%by p). We are assuming that all modes are contributing to the movement
for i=1:N
    m=p(i);
    if m==0
        C=C+Kb;%for p=0 we have an indetermination, when we do the limit w--->0 we obtain Kb
    else %the following formula is obtained derivating the mean energy respect temperature
        C=C+hbar^2*wo^2*(sin(pi*m/N))^2/(Kb*T^2*(sinh(wo*hbar/(Kb*T)*abs(sin(pi*m/N))))^2);
    end
end
Cvec=[Cvec C];
Tvec=[Tvec T];
end
figure(1)
%Plotting the specific heat in function of temperature.
%We also plot the limit NKb, to verify Dulong Petit law 
nvec=N*Kb*ones(length(Tvec));
plot(Tvec,Cvec,'b');hold on
plot(Tvec,nvec,'-r');hold on
h=legend('Specific Heat (T)','Classical Limit= N*Kb');
set(h,'Location','best');
xlabel('Temperature (K)');
ylabel('Specific Heat (J/K)');

    