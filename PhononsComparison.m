%% Phonons simulator: Comparison Version

%There are only some minor differences between this code and the one
%horizontal one. So, to keep it simple, we got rid of all the comments in
%this code except for the ones that explain those differences

a=5.43*10^-10*sqrt(3)/4;
Kb=1.38064852*10^(-23);
hbar=1.054571800*10^(-34);
N=50;
P=[1:N];

%Inputs for the first case (Blue)
m1=28*1.660538921*10^(-27);
K1=59.7939;
T1=295;
B1=1/(Kb*T1);

%Inputs for the second case (Red)
m2=m1;
K2=K1;
T2=595;
B2=1/(Kb*T2);

ii=mod(N,2);
if ii == 1
    n=[-N/2+1/2:N/2-1/2];
else
    n=[-N/2:N/2-1];
end

k=(2*pi*n)/(N*a);

w1=2*sqrt(K1/m1)*abs(sin(k*a/2));

E1=hbar.*w1.*(1./(exp(B1*hbar*w1)-1)+1/2);
[o,u]=find(n==0);
E1(u)=Kb*T1;

A1=sqrt(2*E1/K1);

%We now compute all the variables needed for the second case
w2=2*sqrt(K2/m2)*abs(sin(k*a/2));
E2=hbar.*w2.*(1./(exp(B2*hbar*w2)-1)+1/2);
E2(u)=Kb*T2;
A2=sqrt(2*E2/K2);

x=[1:N]*a;
xfine=[1:0.01:N]*a;

if ii == 1
    nfine=[-N/2+1/2:0.01:N/2-1/2];
else
    nfine=[-N/2:0.01:N/2-1];
end
Nfine=length(nfine);


figure;
timestart=1200;
timelapse=1280;

for t = timestart:timelapse
    
    dx1=[];
    dx2=[];
    tt=t*(10^(-14));
    dx1= sum((A1(P)'*ones(1,N)).*cos((w1(P)'*ones(1,N)).*tt-(k(P)'*n).*a),1); %First case's position variation 
    dx2= sum((A2(P)'*ones(1,N)).*cos((w2(P)'*ones(1,N)).*tt-(k(P)'*n).*a),1); %Second case's position variation 
    dxwave1= sum((A1(P)'*ones(1,Nfine)).*cos((w1(P)'*ones(1,Nfine)).*tt-(k(P)'*nfine).*a),1); %First case's position variation in fine grid
    dxwave2= sum((A2(P)'*ones(1,Nfine)).*cos((w2(P)'*ones(1,Nfine)).*tt-(k(P)'*nfine).*a),1); %Second case's position variation in fine grid
    xsol1=x+dx1; %First case's position
    xsol2=x+dx2; %Second case's position
    y=zeros(length(xsol1),1);
    
    if N > 6
        if ii == 1
            e=(N+1)/2;
        else
            e=N/2;
        end
        
        plot(x(e-2:e+2),y(e-2:e+2),'o','Color',[0 0 0],'MarkerSize',1,'LineWidth',1); %Equilibrium position
        hold on
        axis([x(e-3) x(e+3) -0.4*10^(-9) 0.4*10^(-9)]);
        plot(xsol1(e-2:e+2),y(e-2:e+2),'o','Color',[0 0 1],'MarkerSize',10,'LineWidth',1); %First case's atoms
        plot(xsol2(e-2:e+2),y(e-2:e+2),'o','Color',[1 0 0],'MarkerSize',10,'LineWidth',1); %Second case's Atoms
        hold on
        plot(xfine,dxwave1,'b') %First case's wave
        plot(xfine,dxwave2,'r') %Second case's wave
        
    else %case 2
        plot(x,y,'o','Color',[0 0 0],'MarkerSize',1,'LineWidth',2);
        hold on
        axis([0 x(N)+a -0.4*10^(-9) 0.4*10^(-9)]);
        plot(xsol1,y,'o','Color',[0 0 1],'MarkerSize',10,'LineWidth',1);
        plot(xsol2,y,'o','Color',[1 0 0],'MarkerSize',10,'LineWidth',1);
        hold on
        plot(xfine,dxwave1,'b')
        plot(xfine,dxwave2,'r')
    end
    frame= getframe(1);
    im{t} = frame2im(frame);
    hold off
    
end

filename = 'Monoatomic_chain_comparation.gif';

for i = timestart:timelapse
    [Ag,map] = rgb2ind(im{i},256);
    if i == timestart
        imwrite(Ag,map,filename,'gif','LoopCount',Inf,'DelayTime',.1);
    else
        imwrite(Ag,map,filename,'gif','WriteMode','append','DelayTime',.1);
    end
end


