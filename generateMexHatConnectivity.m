function WS = generateMexHatConnectivity(N)

% Coupling Strength function
%Biological params for a network size 50-200 
CE = 21; 
CI = 21/2;
dE = 10;
dI = 21; 

%Making connectivity range smaller for computational speed with plasticity
%dM = round(1.5.*sqrt(N)); 
dM = 20;

% Distances
xd = -dM:dM;    yd = xd;
[X,Y] = meshgrid(xd,yd);
dXY = sqrt(X.^2 + Y.^2);

% Coupling Strengths
W = CE .* exp(-dXY.^2 ./ dE) - CI .* exp(-dXY.^2 ./ dI);

%Check ~balanced
%trapz(yd,trapz(xd,W,2))
% imagesc(W)

%Computing Sparse connectivity matrix
numConnect =(((2*dM+1)^2))*N;

I=zeros(numConnect,1);
J=zeros(numConnect,1);
V=zeros(numConnect,1);
h=1;

for n1=1:N^2 % for every neuron
    x=floor((n1-1) ./ N)+1;
    y=mod(n1-1, N)+1;
    
    for ii=-dM:dM
        for jj=-dM:dM
            x1=mod(x+ii-1,N)+1;
            y1=mod(y+jj-1,N)+1;
            I(h)=n1;
            J(h)=(x1-1)*N+y1;
            V(h)=W(dM+1+ii,dM+1+jj);
            h=h+1;
        end
    end
end

WS=sparse(I,J,V);

end
