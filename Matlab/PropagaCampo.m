% PropagaCampo  calcula el espectro de potencia marginal correspondiente
% a la propagaci�n de campos en estados de coherencia arbitrarios, es un
% algoritmo basado en el modelo de c�lculo exacto [1],el cual se ha dise�ado
% como una simulaci�n din�mica que permite estudiar el paso de interferencia
% a difraccci�n [2]. se debe resaltar la necesidad de definir el parametro 
% (b/Lambda) como la relaci�n entre la longitud de onda y el espaciamiento
% entre fuentes radiantes donde b=Tamano_Rendija/(Nf-1).Se muestra la
% propagaci�n de la potencia espectral a lo largo del eje z. 
%  
%
% VARIABLES PARA EL C�LCULO DEL ESPECTRO DE POTENCIA MARGINAL
% Nf:               N�mero de fuentes
% Np:               N�mero de pixeles de muetreo en XA
% Lambda:           Longitud de onda
% Zi:               Distancia de propagaci�n inicial.
% ZPaso:            Incremento en la distancia de propagaci�n.
% Zf:               Distancia de propagaci�n final.
% TamanoApertura:   Tama�o de la apertura
% ExtensionXiA:     Extension del plano de entrada
% Sigma:            Desviacion estandar distribuci�n de coherecia
%
% NOTA: TODAS LAS DISTANCIA EN MICROMETROS.
%
% AUTORES.  D. Vargas, E. Franco.

%% PAR�METROS DE ENTRADA
Nf=2; 
Np=5000; 
Lambda=0.632;
Zi=0.05; ZPaso=0.05; Zf=10;
TamanoApertura=5; 
ExtensionXiA=30;
ExtensionXA=ExtensionXiA;
Sigma=10000;

%% GRADO COMPLEJO DE COHERENCIA
GCoherencia = @(d) exp(-d^2/(2*Sigma^2));

%% INTENSIDAD
Intensidad = 1;

%% EXTENCI�N DE LA VENTANA DE APERTURA COORDENA XiA
XiA=linspace(-0.5*TamanoApertura,0.5*TamanoApertura,2*Nf-1);
%% EXTENSI�N DE LA VENTANA DE OBSERVACI�N COORDENA XA
XA=linspace(-0.5*ExtensionXA,0.5*ExtensionXA,Np); 
 
%% MATRIZ DE PERFILES DE POTENCIA ESPECTRAL
PerfilRadiante=[];
PerfilModulador=[];
PerfilPotenciaEspectral=[];
Zetas=[];

%% ESTABLECER TAMA�O Y POSICI�N DE "FIGURE"
TamanoPantalla = get(0,'ScreenSize');
%figure('Position',[100, 100, 0.9*TamanoPantalla(3), 0.6*TamanoPantalla(4)])

%% CORRER SIMULACI�N
tic
Nframe=1;
for Z=Zi:ZPaso:Zf
   
    %% LLAMAR EL ESPECTRO DE POTENCIA MARGINAL
    [CapaReal CapaVirtual]=EPMExacto(Nf,Lambda,Z,XiA,XA,TamanoApertura,GCoherencia,Intensidad);
    EPMarginal=CapaReal + CapaVirtual;
    
    %% POTENCIA ESPECTRAL NORMALIZADA
    PotenciaEspectralCapaRadianteNormailizada=sum(CapaReal,2,'double')/max(sum(CapaReal,2,'double'));
    PotenciaEspectralNormailizada=sum(EPMarginal,2,'double')/max(sum(EPMarginal,2,'double'));
    
    %% MATRIZ DE PERFILES 
    PerfilPotenciaEspectral=[PerfilPotenciaEspectral PotenciaEspectralNormailizada];
    Zetas=[Zetas Z];
    
    %% GRAFICAR DISTRIBUCI�N DE POTENCIA
    imagesc(Zetas,XA,PerfilPotenciaEspectral)
    title(sprintf('Potencia espectral a lo largo del eje z'))
    xlabel('Z  [um]')
    ylabel('\xi_A  [um]')
    colormap hot
%     subplot(1,2,1)
%     imagesc(Zetas,XA,PerfilPotenciaEspectral)
%     title(sprintf('Potencia espectral a lo largo del eje z'))
%     xlabel('Z  [um]')
%     ylabel('\xi_A  [um]')
%     colormap hot
%     
%     subplot(1,2,2)
%     plot(XA,PotenciaEspectralCapaRadianteNormailizada,'k','LineStyle',':')
%     hold
%     area(XA,PotenciaEspectralNormailizada,'FaceColor',[1,0.7,0.3],'EdgeColor',[1 0.1 0])
%     hold             
%     title('Power spectrum at the observation plane'),
%     xlabel('X_A        [m]');
%     ylabel('S(X_A)');  
%     
    %% guardar frame           
    M(Nframe)=getframe(gcf);
    Nframe=Nframe+1;      
    
end
movie2avi(M, 'VideoPropaga.avi','fps',10,'quality',100);
toc

%% REFERENCIAS
    %[1] Casta�eda R, Sucerquia J. Non-approximated numerical modeling of 
    %    propagation of light in any state of spatial coherence
    %[2] Casta�eda R, Vargas D,Franco E. New  criteria  for  non-paraxial
    %    phase-space  modeling  of  interference  and  diffraction.