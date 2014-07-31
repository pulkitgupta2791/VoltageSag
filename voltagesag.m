clear all
clc

% Signal Generation
   

    Ns=1*256;     % No of Samples per cycle, So sampling Freq = 256*50 = 12.8 kHZ
    k=[1:Ns*10];  %% 10 cycles total data to be analyzed
    freqs=50;   %% Frequency of the waveform
    spfreq=Ns*freqs;
    

%   Sag Generation  
    mag = 0.16 : .02 : 0.95    
    mm = length (mag);
    class=1;
    
    
    %%%%%% Parameters for Sag, swell and MI
es=2;
ed=[3 3.5 4 4.5 5];
k1=length(ed);

% %%%%%% Parameters for Harmonic
% harm2=0:0.005:0.1;
% harm3=0 :0.01:.3;
% harm5=0 :0.01:.3;
% harm7=0 :0.01:.3;
% harm9=0 :0.01:.3;
% mm=length(harm3);

for j=1:k1
    for i = 1: mm
      
      %%%%%%%%%%%%%%%   SAG GENERATION BLOCK     
      yp=1*sin(2*pi*freqs*k/spfreq);  % Normal Signal
      
      ef=(es+ed(j));
%    
%       %%%%% Sag / Swell / MIGeneration by Varying Sag Magnitude and Sag Duration
      yp(es*Ns+1:ef*Ns) = mag(i)* yp(es*Ns+1:ef*Ns);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
%       %%%%% Harmonic Generation
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm5(mm-i+1)*sin(2*pi*5*freqs*k/spfreq);  %Only 3rd and 5th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq);  %Only 3rd and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm9(mm-i+1)*sin(2*pi*9*freqs*k/spfreq);  %Only 3rd and 9th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm5(i)*sin(2*pi*5*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm5(i)*sin(2*pi*5*freqs*k/spfreq)+harm9(mm-i+1)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 9th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm7(i)*sin(2*pi*7*freqs*k/spfreq)+harm9(mm-i+1)*sin(2*pi*9*freqs*k/spfreq);  %Only 7th and 9th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm5(mm-i+1)*sin(2*pi*5*freqs*k/spfreq)+harm7(i)*sin(2*pi*7*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm5(mm-i+1)*sin(2*pi*5*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm5(i)*sin(2*pi*5*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq)+harm9(i)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm5(i)*sin(2*pi*5*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq)+harm9(mm-i+1)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm5(mm-i+1)*sin(2*pi*5*freqs*k/spfreq)+harm9(i)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm5(mm-i+1)*sin(2*pi*5*freqs*k/spfreq)+harm9(mm-i+1)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq)+harm9(i)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 7th
%         yp=1*sin(2*pi*freqs*k/spfreq)+harm3(i)*sin(2*pi*3*freqs*k/spfreq)+harm7(mm-i+1)*sin(2*pi*7*freqs*k/spfreq)+harm9(mm-i+1)*sin(2*pi*9*freqs*k/spfreq);  %Only 5th and 7th
  
       
       %%%%%%%%%%% SAG / SWELL With Harmonics
%     ef=(es+ed(j));
%     yp(es*Ns+1:ef*Ns) = mag(i)* yp(es*Ns+1:ef*Ns);
%     
    
subplot(2,2,1);
plot(yp);
title('original signal');
 y = yp;
        %%%%% Wavelet Decomposition %%%%%%%%%%%%%

        [C,L]=wavedec(y,11,'db4');

%         Finding cA 
       
        cA11=appcoef(C,L,'db4',11);A11=wrcoef('a',C,L,'db4',11);
        cA10=appcoef(C,L,'db4',10);A10=wrcoef('a',C,L,'db4',10);
        cA9=appcoef(C,L,'db4',9);A9=wrcoef('a',C,L,'db4',9);
        cA8=appcoef(C,L,'db4',8);A8=wrcoef('a',C,L,'db4',8);
        cA7=appcoef(C,L,'db4',7);A7=wrcoef('a',C,L,'db4',7);
        cA6=appcoef(C,L,'db4',6);A6=wrcoef('a',C,L,'db4',6);
        cA5=appcoef(C,L,'db4',5);A5=wrcoef('a',C,L,'db4',5);
        cA4=appcoef(C,L,'db4',4);A4=wrcoef('a',C,L,'db4',4);
        cA3=appcoef(C,L,'db4',3);A3=wrcoef('a',C,L,'db4',3);
        cA2=appcoef(C,L,'db4',2);A2=wrcoef('a',C,L,'db4',2);
        cA1=appcoef(C,L,'db4',1);A1=wrcoef('a',C,L,'db4',1);

%         fINDING cD
       
        cD11=detcoef(C,L,11);D11 = wrcoef('d',C,L,'db4',11);
        cD10=detcoef(C,L,10);sdv10=std(cD10);D10 = wrcoef('d',C,L,'db4',10);
        cD9=detcoef(C,L,9);sdv9=std(cD9);D9 = wrcoef('d',C,L,'db4',9);
        cD8=detcoef(C,L,8);sdv8=std(cD8);D8 = wrcoef('d',C,L,'db4',8);
        cD7=detcoef(C,L,7);sdv7=std(cD7);D7 = wrcoef('d',C,L,'db4',7);
        cD6=detcoef(C,L,6);sdv6=std(cD6);D6 = wrcoef('d',C,L,'db4',6);

        cD5=detcoef(C,L,5);sdv5=std(cD5);D5 = wrcoef('d',C,L,'db4',5);
        cD4=detcoef(C,L,4);sdv4=std(cD4);D4 = wrcoef('d',C,L,'db4',4);
        cD3=detcoef(C,L,3);sdv3=std(cD3);D3 = wrcoef('d',C,L,'db4',3);
        cD2=detcoef(C,L,2);sdv2=std(cD2);D2 = wrcoef('d',C,L,'db4',2);
        cD1=detcoef(C,L,1);sdv1=std(cD1);D1 = wrcoef('d',C,L,'db4',1);

        D=[D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11];
        ed1=sum(D1.*D1);ed2=sum(D2.*D2);ed3=sum(D3.*D3);ed4=sum(D4.*D4);ed5=sum(D5.*D5);ed6=sum(D6.*D6);
        ed7=sum(D7.*D7);ed8=sum(D8.*D8);ed9=sum(D9.*D9);ed10=sum(D10.*D10);ed11=sum(D11.*D11);
        
        edr1=sqrt(sum(D1.*D1)/length(D1));
        edr2=sqrt(sum(D2.*D2)/length(D2));
        edr3=sqrt(sum(D3.*D3)/length(D3));
        edr4=sqrt(sum(D4.*D4)/length(D4));
        edr5=sqrt(sum(D5.*D5)/length(D5));
        edr6=sqrt(sum(D6.*D6)/length(D6));
        edr7=sqrt(sum(D7.*D7)/length(D7));
        edr8=sqrt(sum(D8.*D8)/length(D8));
 edr9=sqrt(sum(D9.*D9)/length(D9));
 edr10=sqrt(sum(D10.*D10)/length(D10));
 edr11=sqrt(sum(D11.*D11)/length(D11));

 
%         Energy of Approximation
        A=[A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11];
        ea1=sqrt(sum(A1.*A1)/length(A1));
        ea2=sqrt(sum(A2.*A2)/length(A2));
        ea3=sqrt(sum(A3.*A3)/length(A3));
        ea4=sqrt(sum(A4.*A4)/length(A4));
        ea5=sqrt(sum(A5.*A5)/length(A5));
        ea6=sqrt(sum(A6.*A6)/length(A6));
        ea7=sqrt(sum(A7.*A7)/length(A7));
        ea8=sqrt(sum(A8.*A8)/length(A8));
        ea9=sqrt(sum(A9.*A9)/length(A9));
        ea10=sqrt(sum(A10.*A10)/length(A10));
        ea11=sqrt(sum(A11.*A11)/length(A11));
          EA=[ea1 ea2 ea3 ea4 ea5 ea6 ea7 ea8 ea9 ea10 ea11];
          subplot(2,2,2);
          
          plot(EA);
       title('energy of approxumation');
        
         ED=[ed1 ed2 ed3 ed4 ed5 ed6 ed7 ed8 ed9 ed10 ed11];
        EDR=[edr1 edr2 edr3 edr4 edr5 edr6 edr7 edr8 edr9 edr10 edr11];
        subplot(2,2,3);
        
        plot(EDR);
         title('energy of detail coefficients');
         STDV=[std(D1) std(D2) std(D3) std(D4) std(D5) std(D6) std(D7) std(D8) std(D9) std(D10) std(D11)];
        MN=[mean(D1) mean(D2) mean(D3) mean(D4) mean(D5) mean(D6) mean(D7) mean(D8) mean(D9) mean(D10) mean(D11)];
        KRT=[kurtosis(D1) kurtosis(D2) kurtosis(D3) kurtosis(D4) kurtosis(D5) kurtosis(D6) kurtosis(D7) kurtosis(D8) kurtosis(D9) kurtosis(D10) kurtosis(D11)];
        SKW=[skewness(D1) skewness(D2) skewness(D3) skewness(D4) skewness(D5) skewness(D6) skewness(D7) skewness(D8) skewness(D9) skewness(D10) skewness(D11)];
% %         ENTP =[entropy(D1) entropy(D2) entropy(D3) entropy(D4) entropy(D5) entropy(D6) entropy(D7) entropy(D8) entropy(D9) entropy(D10) entropy(D11) entropy(D12) entropy(D13)];
        ENTP =[wentropy(D1,'shannon') wentropy(D2,'shannon') wentropy(D3,'shannon') wentropy(D4,'shannon') wentropy(D5,'shannon') wentropy(D6,'shannon') wentropy(D7,'shannon') wentropy(D8,'shannon') wentropy(D9,'shannon') wentropy(D10,'shannon') wentropy(D11,'shannon')];
%         Feature of Wavelet Transform
        feature_pq_wave(j,i,:)=[EDR STDV MN KRT SKW ENTP];
        feature_pq_wave_new(mm*(j-1)+i,:)=[EDR STDV MN KRT SKW ENTP class];
        
    end
end