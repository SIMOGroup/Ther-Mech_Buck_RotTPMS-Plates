function [CP_ref, uKnot_ref, vKnot_ref] = refine_h_surf(p,q,uKnot,vKnot,CP,ref_u,ref_v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Refine NURBS surf by h-refinement %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Source: Adapted from algorithm in Piegl, Les. "The NURBS Book". Springer-Verlag: Berlin 1995; p. 167.
mcp = length(uKnot)-p-1; ncp = length(vKnot)-q-1;
ref_uKnot = refine_knot(uKnot,ref_u); ref_vKnot = refine_knot(vKnot,ref_v);  % Knot vector refinement

%% === Project 4D control points to 3D ===
for j = 1:ncp
    for i = 1:mcp
        Pw(i,j,1:3) = CP(i,j,1:3) * CP(i,j,4);
        Pw(i,j,4)   = CP(i,j,4);
    end
end

%% === Refine uKnot in u-direction ===
if isempty(ref_uKnot)
    uKnot_ref = uKnot;
    Qw = Pw;
else
    r = length(ref_uKnot); 
    mu = mcp+p+1;
    
    a = find_knot_span(ref_uKnot(1),uKnot,mcp);
    b = find_knot_span(ref_uKnot(r),uKnot,mcp) + 1;
    
    for col = 1:ncp
        for j = 1:a-p
            Qw(j,col,:) = Pw(j,col,:);
        end
        for j = b-1:mcp
            Qw(j+r,col,:) = Pw(j,col,:);
        end
    end

    for j = 1:a
        uKnot_ref(j) = uKnot(j);
    end
    for j = b+p:mu
        uKnot_ref(j+r) = uKnot(j);
    end

    i = b+p-1; k = i+r;
    for  j = r:-1:1
        while (ref_uKnot(j)<=uKnot(i) && i>a)
            for col = 1:ncp
                Qw(k-p-1,col,:) = Pw(i-p-1,col,:);
            end
            uKnot_ref(k) = uKnot(i);
            k = k-1; i = i-1;
        end

        for col = 1:ncp
            Qw(k-p-1,col,:) = Qw(k-p,col,:);
        end
        for l = 1:p
            ind = k-p+l;
            alfa = (ref_uKnot(j)-uKnot_ref(k+l)) / (uKnot(i-p+l)-uKnot_ref(k+l));
            for col = 1:ncp
                Qw(ind-1,col,:) = alfa*Qw(ind-1,col,:) + (1-alfa)*Qw(ind,col,:);
            end
        end
        uKnot_ref(k) = ref_uKnot(j);
        k = k-1;
    end
    Pw = Qw;
    mcp = mcp+r;
end

%% === Refine vKnot in v-direction ===
if isempty(ref_vKnot)
    vKnot_ref = vKnot;
else
    r = length(ref_vKnot); 
    mv = ncp+q+1;
    a = find_knot_span(ref_vKnot(1),vKnot,ncp);
    b = find_knot_span(ref_vKnot(r),vKnot,ncp) + 1;

    for row = 1:mcp
        for j = 1:a-q
            Qw(row,j,:) = Pw(row,j,:);
        end
        for j = b-1:ncp
            Qw(row,j+r,:) = Pw(row,j,:);
        end
    end 
    for j = 1:a
        vKnot_ref(j) = vKnot(j);
    end
    for j = b+q:mv
        vKnot_ref(j+r) = vKnot(j);
    end

    i = b+q-1; k = i+r;
    for j = r:-1:1
        while ((ref_vKnot(j)<=vKnot(i)) && (i>a))
            for row = 1:mcp
                Qw(row,k-q-1,:) = Pw(row,i-q-1,:);
            end
            vKnot_ref(k) = vKnot(i);
            k = k-1; i = i-1;
        end

        for row = 1:mcp
            Qw(row,k-q-1,:) = Qw(row,k-q,:);
        end
        for l = 1:q
            ind = k-q+l;
            alfa = (ref_vKnot(j)-vKnot_ref(k+l)) / (vKnot(i-q+l)-vKnot_ref(k+l));
            for row = 1:mcp
                Qw(row,ind-1,:) = alfa*Qw(row,ind-1,:) + (1-alfa)*Qw(row,ind,:);
            end
        end
        vKnot_ref(k) = ref_vKnot(j);
        k = k-1;
    end
end

%% === Project 3D control points to 4D ===
for j = 1:length(Qw(1,:,1))
    for i = 1:length(Qw(:,1,1))
        CP_ref(i,j,1:3) = Qw(i,j,1:3) / Qw(i,j,4);
        CP_ref(i,j,4) = Qw(i,j,4);
    end
end
end
