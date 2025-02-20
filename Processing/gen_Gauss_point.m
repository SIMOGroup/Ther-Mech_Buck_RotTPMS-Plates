function [x_Gauss, w_Gauss] = gen_Gauss_point(x1,x2,nGauss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Gaussion integral points %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
eps = 1e-15;
m = (nGauss+1)/2;
xm = (1/2.)*(x2+x1);
xl = (1/2.)*(x2-x1);
pi = acos(-1.);
z1 = 0;

for i = 1:m
    z = cos(pi*(i-0.25)/(nGauss+0.5));
    p1 = 1.;
    p2 = 0.;
              
    while abs(z-z1) > eps
        p1 = 1.;
        p2 = 0.;
        for j = 1:nGauss
            p3 = p2;
            p2 = p1;
            p1 = ((2.*j-1.)*z*p2-(j-1.)*p3) / j;
        end
        pp = nGauss*(z*p1-p2)/(z^2-1.);
        z1 = z;
        z = z1-p1/pp;
    end
        x_Gauss(i) = xm-x1*z;
        x_Gauss(nGauss+1-i) = xm+x1*z;
        w_Gauss(i) = 2.*x1/((1.-z*z)*pp*pp);
        w_Gauss(nGauss+1-i) = w_Gauss(i);
end
x_Gauss = x_Gauss';
w_Gauss = w_Gauss';
end
