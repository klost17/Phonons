%% Phonons simulator: Vertical Version

%There are only some minor differences between this code and the one
%horizontal one. So, to keep it simple, we got rid of all the comments in
%this code except for the ones that explain those differences

m=28*1.660538921*10^(-27);
K=59.7939;
a=5.43*10^-10*sqrt(3)/4;
N=50;
T=293;
P=[19:N-18];
Kb=1.38064852*10^(-23);
hbar=1.054571800*10^(-34);
B=1/(Kb*T);

ii=mod(N,2);
if ii == 1
    n=[-N/2+1/2:N/2-1/2];
else
    n=[-N/2:N/2-1];
end
k=(2*pi*n)/(N*a);

w=2*sqrt(K/m)*abs(sin(k*a/2));

E=hbar.*w.*(1./(exp(B*hbar*w)-1)+1/2);
[o,u]=find(n==0);
E(u)=Kb*T;

A=sqrt(2*E/K);

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
    dx=[];
    tt=t*(10^(-14));
    dx= sum((A(P)'*ones(1,N)).*cos((w(P)'*ones(1,N)).*tt-(k(P)'*n).*a),1);
    dxonda= sum((A(P)'*ones(1,Nfine)).*cos((w(P)'*ones(1,Nfine)).*tt-(k(P)'*nfine).*a),1);
    xsol=x+dx;
    y=zeros(length(xsol),1); %To get the equilibrium position mark of every atom in the x axis(y=0)
    if N > 6
        if ii == 1
             e=(N+1)/2;
        else
            e=N/2;
        end
        plot(x(e-2:e+2),y(e-2:e+2),'o','Color',[1 0 0],'MarkerSize',1,'LineWidth',1); % Equilibrium position
        %title(['movimiento cadena monoatómica 1D, modos normales: 1:' num2str( N ) ]) 
        hold on
        axis([x(e-3) x(e+3) -0.4*10^(-9) 0.4*10^(-9)]); 
        % Atoms: in this case we represent the movement of the atoms in the y axis
        plot(x(e-2:e+2),dx(e-2:e+2),'o','Color',[0 0 0],'MarkerSize',10,'LineWidth',1);
        hold on
        plot(xfine,dxonda) % Wave
        
    else 
        plot(x,y,'o','Color',[1 0 0],'MarkerSize',1,'LineWidth',2);
        %title(['movimiento cadena monoatómica 1D, modos normales: 1:' num2str( N ) ]) 
        hold on
        axis([0 x(N)+a -0.4*10^(-9) 0.4*10^(-9)]);
        plot(x,dx,'o','Color',[0 0 0],'MarkerSize',10,'LineWidth',1);
        hold on
        plot(xfine,dxonda)
    end
    frame= getframe(1);
    im{t} = frame2im(frame);
    hold off
    
end

% Creation of GIF

filename = 'Monoatomic_chain_vertical.gif';

for i = timestart:timelapse
    [Ag,map] = rgb2ind(im{i},256);
    if i == timestart
        imwrite(Ag,map,filename,'gif','LoopCount',Inf,'DelayTime',.1);
    else
        imwrite(Ag,map,filename,'gif','WriteMode','append','DelayTime',.1);
    end
end


