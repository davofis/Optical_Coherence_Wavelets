function [CapaReal, CapaVirtual]=EPMExacto(Nf,Lambda,Z,XiA,XA,TamanoApertura,GCoherencia,Intensidad)
% EPM calcula el espectro de potencia marginal correspondiente a la
% propagaci�n de campos en estados de coherencia arbitrarios, simulaci�n 
% basada en la teoria de onditas de coherencia espacial; es un algoritmo 
% enfocado a procesos de difracci�n unidimensionales bajo el enfoque del
% modelo num�rico sin aproximaciones, "c�lculo exacto" [1].
%
% ----------VARIABLES DE ENTRADA------------
% Nf:               N�mero de fuentes
% Lambda:           Longitud de onda
% Z:                Distancia del plano de entrada al plano de salida.
% XiA:              Coordenada XiA de entrada
% XA:               Coordenada XA de salida  
% TamanoApertura:   Tama�o de la apertura
% GCoherencia:      funci�n handle "Grado de coherencia" 
% Intensidad:       Distribuci�n de intensidad
%
% ----------VARIABLES DE SALIDA--------------
% CapaReal:         Capa real del espectro de potencia marginal
% CapaVirtual:      Capa virtual del espectro de potencia marginal
% EPMarginal:       Espectro de potencia marginal
%
%
% NOTA: TODAS LAS DISTANCIA EN MICROMETROS.
%
% AUTORES.  D. Vargas

%% MATRIZ BASE
n=length(XA);   % n�mero de pixeles en la ventana de observaci�n XA
m=2*Nf-1;       % n�mero de pixeles en la ventana de entrada XiA 
Malla=meshgrid(1:m,1:n);

%% MATRICES DE COORDENADAS DEL PLANO DE SALIDA XA Y ENTRADA XiA
[XAmatriz, XiAmatriz]=meshgrid(XA,XiA);
XiAmatriz=XiAmatriz';
XAmatriz=XAmatriz';

%% N�MERO DE ONDA, DISCRETIZACI�N DE RENDIJA
k=2*pi/Lambda;
TamanoPixel=TamanoApertura/(2*Nf-2);

%% CALCULAR CAPA REAL 
 LorentzReal=((Z + sqrt(Z^2 + (XAmatriz-XiAmatriz).^2))./(Z^2 + (XAmatriz-XiAmatriz).^2)).^2;
 CapaReal=(Intensidad/(4*Lambda^2))*LorentzReal.*mod(Malla,2);
 
%% CALCULAR CAPA VIRTUAL 
CapaVirtual=zeros(n,m);
Estructura=mod(Malla,2)<=0;
if Nf>=2
    for familia=2:Nf        
        Estructura(:,1:familia-1)=0;
        Estructura(:,2*Nf-1-(familia-2):end)=0;
        
        XiD=2*(familia-1)*TamanoPixel ;
        Fase1=(Z^2) + ((XAmatriz-XiAmatriz).^2) + ((XiD^2)/4) + (XiD*XiAmatriz) - (XiD*XAmatriz);
        Fase2=(Z^2) + ((XAmatriz-XiAmatriz).^2) + ((XiD^2)/4) - (XiD*XiAmatriz) + (XiD*XAmatriz);
       
        Factor1=(Z + sqrt(Fase1))./(Fase1);        
        Factor2=(Z + sqrt(Fase2))./(Fase2); 
        
        LorentzVirtual=(2*Intensidad/(4*Lambda^2))*Factor1.*Factor2;
        Coseno=cos(k*sqrt(Fase1) - k*sqrt(Fase2));
        
        FuenteVirtual=GCoherencia(XiD)*LorentzVirtual.*Coseno;
        
        Virtuales=Estructura.*FuenteVirtual;
        CapaVirtual=CapaVirtual+Virtuales;
        Estructura=Estructura==0; % revierte los 0�s por 1�s.
    end
end
    
 %% REFERENCIAS
    %[1] Casta�eda R, Sucerquia J. Non-approximated numerical modeling of 
    %    propagation of light in any state of spatial coherence

 