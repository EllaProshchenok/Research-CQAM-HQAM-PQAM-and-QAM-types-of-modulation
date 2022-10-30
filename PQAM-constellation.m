clear variables
close all

M = 64;
Theta = pi/3;
MinDistance = 1.15;
BaseShift = MinDistance*cos(Theta)/2;
Shifts = [-BaseShift BaseShift];
LineSpacing = MinDistance*sin(Theta);

BaseMatrix = zeros(sqrt(M));
BaseStructure = [1 3; 0 2]; %Base Grey QPSK constellation

for j = 1:log2(M)/2-1
    
    if j == 1
        BaseMatrix(1:2^j, 1:2^j) = BaseStructure;
    end
    
    BaseMatrix(1:2^j, 2^j+1:2^(j+1)) = fliplr(BaseMatrix(1:2^j, 1:2^j)) + 4^(j);
    BaseMatrix(2^j+1:2^(j+1),1:2^(j+1)) = flipud(BaseMatrix(1:2^j, 1:2^(j+1))) + 2*4^(j);
    
end

for j=1:length(BaseMatrix)
    
    MatrixLinear(1,(j-1)*sqrt(M)+1:j*sqrt(M)) = BaseMatrix(j,:);
    
end

%Creating matrix of constellaion points
PointDistributionVer = ones(1,log2(M)/2) * LineSpacing;   %QAM

DistanceVer = zeros(1,2.^((log(M)/log(4))-1));
DistanceVer(1) = PointDistributionVer(1)/2;

for j = 1:(log(M)/log(4))-1
    
    init (j) = 2 + 2^(j-1) -1;
    
    step (j) = 2.^j;
    
    if j == 1
        ends (j) = 2.^((log(M)/log(4))-1);
    else
        ends (j) = ends(j-1) -  2.^(j-2);
    end
    
    DistanceVer((init(j):step(j):ends(j))) = PointDistributionVer(end - (j - 1));
    
end

DistanceVer = cumsum(DistanceVer);

PointDistributionHor = ones(1,log2(M)/2);   %QAM

DistanceHor = zeros(1,2.^((log(M)/log(4))-1));
DistanceHor(1) = PointDistributionHor(1)/2;

for j = 1:(log(M)/log(4))-1
    
    init (j) = 2 + 2^(j-1) -1;
    
    step (j) = 2.^j;
    
    if j == 1
        ends (j) = 2.^((log(M)/log(4))-1);
    else
        ends (j) = ends(j-1) -  2.^(j-2);
    end
    
    DistanceHor((init(j):step(j):ends(j))) = PointDistributionHor(end - (j - 1));
    
end

DistanceHor = cumsum(DistanceHor);

ConstImag = [-fliplr(DistanceVer) DistanceVer];
ConstImag = flipud(repmat(ConstImag, sqrt(M),1).');

ConstReal = [-fliplr(DistanceHor) DistanceHor];
ConstReal = repmat(ConstReal,sqrt(M),1);

for i = 1 : sqrt(M)
   ConstReal(i,:) =  ConstReal(i,:) + Shifts(rem(i,2) + 1);
end

ConstSquare = complex(ConstReal,ConstImag);

for j = 1 : length(ConstSquare)
    ConstLin(1,(j-1)*sqrt(M)+1:j*sqrt(M)) = ConstSquare(j,:);
end

ConstPrefinal = zeros(1,M);

for j = 0 : M-1
    ind = find(MatrixLinear == j);
    ConstPrefinal(1,j+1) = ConstLin(ind);
end

Constellation = ConstPrefinal;

coef = rms(Constellation);

save('16PQAM.mat','Constellation')

