clear variables
M = 64; %16 or 64
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

PointDistribution = ones(1,log2(M)/2);   %QAM


Distance = zeros(1,2.^((log(M)/log(4))-1));
Distance(1) = PointDistribution(1)/2;

for j = 1:(log(M)/log(4))-1
    
    init (j) = 2 + 2^(j-1) -1;
    
    step (j) = 2.^j;
    
    if j == 1
        ends (j) = 2.^((log(M)/log(4))-1);
    else
        ends (j) = ends(j-1) -  2.^(j-2);
    end
    
    Distance((init(j):step(j):ends(j))) = PointDistribution(end - (j - 1));
    
end

Distance = cumsum(Distance);

ConstReal = [-fliplr(Distance) Distance];
ConstReal = repmat(ConstReal, sqrt(M),1);

ConstImag = flipud(transpose(ConstReal(1,:)));
ConstImag = repmat(ConstImag,1,sqrt(M));

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

Filename = [sprintf( '%02d', M ) 'QAM.mat'];
save(Filename,'Constellation')
