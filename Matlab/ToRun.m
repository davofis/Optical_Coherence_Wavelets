

%% PAR�METROS DE ENTRADA
Nf=20;              
TamanoApertura=10;   
Lambda=0.632; 	      
Z=100; 
Sigma=1000000;
Gamma=1000000;
Np=512; 

alph = 0.01;     
beta = (1 - sqrt(1 + 8*alph*Lambda*Z))/(4*alph*Lambda);
ExtensionXA =(sqrt(beta^2 + Z^2) + TamanoApertura);  
display(ExtensionXA);
%% GRADO COMPLEJO DE COHERENCIA
GCoherencia = @(d) exp(-d^2/(2*Sigma^2));
%% INTENSIDAD
Intensidad = 1;
%% EXTENCI�N DE LA VENTANA DE APERTURA COORDENA XiA
XiA=linspace(-0.5*TamanoApertura,0.5*TamanoApertura,2*Nf-1);
%% EXTENSI�N DE LA VENTANA DE OBSERVACI�N COORDENA XA
XA=linspace(-0.5*ExtensionXA,0.5*ExtensionXA,Np); 
 
%% LLAMAR EL ESPECTRO DE POTENCIA MARGINAL
[CapaReal, CapaVirtual]=EPMExacto(Nf,Lambda,Z,XiA,XA,TamanoApertura,GCoherencia,Intensidad);
EPMarginal=CapaReal + CapaVirtual;
 
%% CONTRIBUCIONES DE POTENCIA 
 % CAPA REAL 
 PotenciaEspectralReal=sum(CapaReal,2,'double');
 % CAPA VIRTUAL 
 PotenciaEspectralVirtual=sum(CapaVirtual,2,'double');
 % POTENCIAL SALIDA TOTAL 
 PotenciaEspectral = PotenciaEspectralReal + PotenciaEspectralVirtual;
 % POTENCIAL ENTRADA
  
%% SET GRAPHICS 
% ESTABLECER TAMA�O Y POSICI�N DE "FIGURE"
TamanoPantalla = get(0,'ScreenSize');
% figure('Position',[0, 0, TamanoPantalla(3), TamanoPantalla(4)])
figure('Position',[100, 50, 1200,642 ]) 

% Frente de onda discontinuo
subplot(1,3,1)
imagesc(XiA,XA,EPMarginal)
title({'(a)  Marginal power spectrum';''},'fontsize',9)
xlabel('\xi_A [um]','fontsize',10)
ylabel('X_A [um]','fontsize',10)
colormap gray
      
subplot(1,3,2)
area(XA,PotenciaEspectral,'FaceColor','y','EdgeColor','r')
title({'(b)  Power spectrum at the';'observation plane'},'fontsize',9)
xlabel('X_A  [um]','fontsize',10)
ylabel('S(X_A)','fontsize',10)
text(0.2*ExtensionXA,0.9,['Z=',num2str(Z),'um'],'fontsize',8)
    
subplot(1,3,3)
plot(XA,PotenciaEspectralReal,XA,PotenciaEspectralVirtual) 
title({'(c)  Power Contributions';''},'fontsize',9)
xlabel('X_A  [um]','fontsize',10)
ylabel('Arbitrary units','fontsize',10)
xlim([-ExtensionXA ExtensionXA])
hold
plot(XA,-PotenciaEspectralReal,'b')
hold
  
  
  
