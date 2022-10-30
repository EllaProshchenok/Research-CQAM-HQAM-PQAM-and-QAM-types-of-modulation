clear variables
close all

%Parameters
M = 64; %16 or 64
BasePSKConstM = 8;
CircleNum = M / BasePSKConstM;

BaseCircleRadius = 1;
InitPhase = 0;
BaseAngles = linspace(0,2*pi – 2*pi/BasePSKConstM,BasePSKConstM);
MinimalDistance = 2*BaseCircleRadius*sin((2*pi/BasePSKConstM)/2);

RadiusIncrease1 = (2*BaseCircleRadius*cos(pi/BasePSKConstM) + sqrt((2*BaseCircleRadius*cos(pi/BasePSKConstM)).^2 – 4*(BaseCircleRadius.^2 – MinimalDistance.^2)))/2 – BaseCircleRadius;
RadiusIncrease2 = MinimalDistance – RadiusIncrease1;
RadiusIncrease = [RadiusIncrease1 RadiusIncrease2];

RotationAngles = [pi/BasePSKConstM 0];

CircleRadius = BaseCircleRadius;


for I = 1 : CircleNum
    
    if I == 1
        CircleRadiusIncrease(i) = 0;
    else
        CircleRadiusIncrease(i) = RadiusIncrease(rem(I,2)+1);
        
    end
    
    RotationAngle = RotationAngles(rem(I,2)+1);
    CircleRadius = CircleRadius + CircleRadiusIncrease(i);
    InPhase = CircleRadius .* sin(BaseAngles + RotationAngle);
    Quadrature = CircleRadius .* cos(BaseAngles + RotationAngle);
    Points(i,:) = complex(InPhase, Quadrature);
    
end

Constellation = reshape(Points.’, [1 numel(Points)]);
