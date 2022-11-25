function result=Trigrid_interpor_for_r_lambda(z,x,y,xout,yout)

%Nterms=N_elements(x)*N_elements(y)
Nterms=numel(x);
Xtemp=x' ;Xtemp= Xtemp(:);
Ytemp=repmat(y,1,Nterms/numel(y))'; 
SigTemp=z';SigTemp=SigTemp(:);
[xq,yq]=meshgrid(xout,yout);
%A=griddata(Xtemp,Ytemp,SigTemp,xq,yq,'cubic');
F=scatteredInterpolant(Xtemp,Ytemp,SigTemp,'natural','linear');
A = F(xq,yq);
% figure('Visible','on')
% plot(Xtemp,Ytemp,'.');
% hold on
% contourf(xq,yq,A)
% hold off


result=struct('Z',A,'X',xout,'Y',yout);
end

