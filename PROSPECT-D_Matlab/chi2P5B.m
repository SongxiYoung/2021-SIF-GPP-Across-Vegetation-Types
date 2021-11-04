% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2P5B(x,rmes,tmes)
N=x(1);
Cab=x(2);
Ccx=x(3);Cant=x(4);
Cbp=x(5);
Cw=x(6);
Cdm=x(7);
RT=prospect_DB(N,Cab,Ccx,Cant,Cbp,Cw,Cdm);
chi2=sqrt(sum((RT(:,2)-rmes).^2+(RT(:,3)-tmes).^2));
%chi2=sqrt(sum((RT(1:2051,2)-rmes).^2+(RT(1:2051,3)-tmes).^2));
%chi2=sqrt(sum((RT([174 290 725],2)-rmes).^2+(RT([174 290 725],3)-tmes).^2));


