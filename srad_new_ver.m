%
%
% Implementation of SRAD Filter
% Ref :Yongjian Yu and Scott T. Acton, "Speckle Reducing Anisotropic
% Diffusion",IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 11, NO. 11, NOVEMBER 2002
% Jeny Rajan, Chandrashekar P.S
% usage S=srad(I,T)
% I- noisy image
% T - threshold (greater than 0).
function S=srad_new_ver(I,T,niterations)
% tic
[x y]=size(I);
% I=double(I);
% Ic=double(I);
delta_t = 0.08;
% t=1;
eps=0.00000000001;
T_stop = T;
h = waitbar(0,'Filtering Image...');
%Initialize the picture Ic (new picture) with zeros
Ic=double(I);
%Apply niteration of the algorithm to the image
for i = 1:niterations            
  % fprintf('\rIteration %d',i);
  % if i >=2 
  %    I=Ic;
  % end
    for t = 1:T_stop     
        qt=exp(-t*.2);
        [Ix,Iy] = gradient(Ic);    
        di=sqrt(Ix.^2+Iy.^2);
        di2=del2(Ic);
        T1=0.5*((di./(Ic+eps)).^2);
        T2=0.0625*((di2./(Ic+eps)).^2);
        T3=(1+(0.25*(di2./(Ic+eps)))).^2;
        T=sqrt((T1-T2)./(T3+eps));
        dd=(T.^2-qt.^2)./((qt.^2*(1+qt.^2)+eps));
        cq=1./(1+dd);
        [D1,~]=gradient(cq.*Ix);
        [~,D4]=gradient(cq.*Iy);
        D=D1+D4;    
        Ic=real(Ic+delta_t .*D); 
        % update waitbar
        waitbar(((i-1)*T_stop+t)/(T_stop*niterations), h);
    end
end
% toc
S=uint8(Ic);
close(h)
return