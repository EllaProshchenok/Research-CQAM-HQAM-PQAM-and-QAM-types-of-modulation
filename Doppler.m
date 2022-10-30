clearvars
close all
tic

%Generating binary sequence
Data = round(rand(1, 153600));

%Modulation (Manupulation)
ModOrder = [16 64]; %16 or 64
ModTypeOptions = {'CQAM' 'PQAM' 'HQAM' 'QAM'};
Fd = 384e3;
SymbolT = 0.1;

LowFreq = 10e3;
FreqBand = 8e3;
DeltaFreq = 1/SymbolT;
SubcarNum = FreqBand / DeltaFreq;
SubcarSpace = (1:SubcarNum);
MaxSpeed = 1;
MaxDopplerShift = floor((LowFreq + FreqBand/2)*(1+(MaxSpeed/1500)) - (LowFreq + FreqBand/2));

t = linspace(0,SymbolT,SymbolT*Fd);

%Noise Parameters
NoiseMagnitude = (0.2:0.05:1.5);

%Adding data
for a = 1 : length(ModTypeOptions)
    for b = 1 : length(ModOrder)
        
        clearvars  DataDec Subcar OFDMSymbol SNRSymbol RecConst DetectedSymbols
        
        M = ModOrder(b);
        ModType = cell2mat(ModTypeOptions(a));
        Filename = [sprintf( '%02d', M ) ModType '.mat'];
        load(Filename)
        %Decimal conversion
        parfor i = 1 : length(Data)/log2(M)
            
            DataDec(i) = bi2de(Data(log2(M)*(i-1)+1:log2(M)*i));
            
        end
        
        if length(ModType) == 3
            ModData = qammod(DataDec,M);
        else
            ModData = genqammod(DataDec,Constellation);
        end
        
        SymbolNum = length(ModData) / length(SubcarSpace);
        ModDataReshaped = reshape(ModData, [length(SubcarSpace) SymbolNum]).';
        
        %Quadrature Modulator
        
        for i = 1 : SymbolNum
            parfor j = 1 : length(SubcarSpace)
                
                ITX = real(ModDataReshaped(i,j)) .* real(exp(1i*2*pi*(LowFreq + SubcarSpace(j)/SymbolT).*t));
                QTX = imag(ModDataReshaped(i,j)) .* (-imag(exp(1i*2*pi*(LowFreq + SubcarSpace(j)/SymbolT).*t)));
                
                Subcar(j,:) = ITX + QTX;
                
            end
            PureOFDMSymbol = sum(Subcar);
            
            DopplerShift = (rand(1)-0.5)/2;
            DopplerFd = (Fd*(LowFreq + FreqBand/2+DopplerShift))/(LowFreq + FreqBand/2);
            Doppt = linspace(0,SymbolT,round(SymbolT*DopplerFd));
            DoppOFMDSymbol = interp1(t,PureOFDMSymbol,Doppt);
            
            if length(DoppOFMDSymbol) > Fd*SymbolT
                DoppOFMDSymbol(Fd*SymbolT+1:end) = [];
            else
                DoppOFMDSymbol = [DoppOFMDSymbol zeros(1,Fd*SymbolT-length(DoppOFMDSymbol))];
            end
            
            OFDMSymbol(i,:) = DoppOFMDSymbol ./ max(abs(DoppOFMDSymbol));
        end
        
        for  z = 1 : length(NoiseMagnitude)
            tic
            Noise = NoiseMagnitude(z) .* rand(SymbolNum, SymbolT*Fd);
            SymbolsRX = OFDMSymbol + NoiseMagnitude(z) .* Noise;
            
            for i = 1 : SymbolNum
                NoiseFiltered = bandpass(Noise(i,:),[LowFreq-500 LowFreq + FreqBand + 500],Fd);
                SNRSymbol(i) = snr(OFDMSymbol(i,:), NoiseFiltered);
                %Quadrature demodulator
                parfor j = 1 : SubcarNum
                    IRX = SymbolsRX(i,:) .* real(exp(1i*2*pi*(LowFreq + SubcarSpace(j)/SymbolT).*t));
                    QRX = SymbolsRX(i,:) .* (-imag(exp(1i*2*pi*(LowFreq + SubcarSpace(j)/SymbolT).*t)));
                    
                    IRXInt = cumtrapz(IRX);
                    QRXInt = cumtrapz(QRX);
                    
                    if z == 2 && a == 4 && b == 2
                    Check1(j,:) = IRXInt;
                    end
                    
                    Iconst = IRXInt(end);
                    Qconst = QRXInt(end);
                    
                    RecConst(i,j) = complex(Iconst,Qconst);
                end
                RealPart = real(RecConst(i,:));
                
                ImagPart = imag(RecConst(i,:));
                
                
                
                if length(ModType) == 3
                    RealPart = (mean(abs(real(ModDataReshaped(i,:))))/ mean(abs(real(RecConst(i,:))))) * RealPart;
                    ImagPart = (mean(abs(imag(ModDataReshaped(i,:))))/ mean(abs(imag(RecConst(i,:))))) * ImagPart;
                else
                    RealPart = (mean(abs(real(Constellation)))/ mean(abs(real(RecConst(i,:))))) * RealPart;
                    ImagPart = (mean(abs(imag(Constellation)))/ mean(abs(imag(RecConst(i,:))))) * ImagPart;
                end
                
                DetectedSymbols(i,:) =  complex(RealPart,ImagPart);
            end
            DetectedSymbolsVector = reshape(DetectedSymbols.',[1 numel(DetectedSymbols)]);
            
            if length(ModType) == 3
                DemodulatedDecData = qamdemod(DetectedSymbolsVector.', M);
            else
                DemodulatedDecData = genqamdemod(DetectedSymbolsVector.', Constellation);
            end
            
            
            DemodulatedData = de2bi(DemodulatedDecData,log2(M));
            ReceivedData = reshape(DemodulatedData.', [1 numel(DemodulatedData)]);
            
            [~,BER(z)] = biterr(Data, ReceivedData);
            SNR(z) = mean(SNRSymbol);
            toc
            
            if z == 2 && a == 4 && b == 2
                spectrogram(OFDMSymbol(1,:))
            end
            
        end
        Filename = ['BER_Doppler' sprintf( '%02d', M ) ModType '.mat'];
        save(Filename,'BER','SNR')
    end
end
toc
