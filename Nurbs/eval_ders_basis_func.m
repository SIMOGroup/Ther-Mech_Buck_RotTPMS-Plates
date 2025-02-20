function dN = eval_ders_basis_func(i,p,u,Knot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluate all derivatives of BSpline basis finction at a natural knot point in 1D Knot vector %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Source: Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: Berlin 1995; pp. 72-73.
nders = 2;
left = zeros(1,p+1); right = zeros(1,p+1);
ndu = zeros(p+1,p+1); a = zeros(2,p+1);
dN = zeros(nders+1,p+1);

% --- Shape function ---
ndu(1,1) = 1;
for j = 1:p
    left(j+1) = u - Knot(i+1-j);
    right(j+1) = Knot(i+j) - u;
    saved = 0;
    for r = 0:j-1
        ndu(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = ndu(r+1,j) ./ ndu(j+1,r+1);   %%%
        ndu(r+1,j+1) = saved + right(r+2) .* temp;   
        saved = left(j-r+1) .* temp;
    end 
    ndu(j+1,j+1) = saved;
end 

for j = 0:p
    dN(1,j+1) = ndu(j+1,p+1);
end 

% --- Derivatives ---
for r = 0:p % loop over function index
s1 = 0;
% alternate rows in array a
s2 = 1;
a(1,1) = 1;

    for k = 1:nders % loop to compute kth derivative
        d = 0d+0;
        rk = fix(r-k);
        pk = fix(p-k);

        if r >= k
            a(s2+1,1) = a(s1+1,1) ./ ndu(pk+2,rk+1);
            d = a(s2+1,1) .* ndu(rk+1,pk+1);
        end

        if rk >= -1
            j1 = 1;
        else
            j1 = fix(-rk);
        end

        if (r-1) <= pk
            j2 = fix(k-1);
        else
            j2 = fix(p-r);
        end

        for j = j1:j2
            a(s2+1,j+1) =(a(s1+1,j+1) - a(s1+1,j))./ndu(pk+2,rk+j+1);
            d = d + a(s2+1,j+1).*ndu(rk+j+1,pk+1);
        end % j =fix(j2+1);

        if r <= pk
            a(s2+1,k+1) = -a(s1+1,k)./ndu(pk+2,r+1);
            d = d + a(s2+1,k+1).*ndu(r+1,pk+1);
        end

        dN(k+1,r+1) = d;
        j = fix(s1);
        s1 = fix(s2);
        s2 = fix(j);  % switch rows
    end 
end    

% Multiply through by the correct factors;
r = fix(p);
for k = 1:nders
   for j = 0:p
       dN(k+1,j+1) = dN(k+1,j+1).*r;
   end
r = fix(r.*(p-k));
end

end

