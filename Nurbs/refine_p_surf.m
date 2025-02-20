function [CP_ref, uKnot_ref, vKnot_ref, p_ref, q_ref] = refine_p_surf(p,q,uKnot,vKnot,CP,tp,tq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Refine NURBS surf by p-refinement %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Source: Adapted from algorithm in Piegl, Les. "The NURBS Book". Springer-Verlag: Berlin 1995; p. 206.
mcp = length(uKnot)-p-1; ncp = length(vKnot)-q-1;

%% === Project 4D control points to 3D ===
for i = 1:mcp
    for j = 1:ncp
        Pw(i,j,1:3) = CP(i,j,1:3)*CP(i,j,4);
        Pw(i,j,4)   = CP(i,j,4);
    end
end

%% === Refine uKnot in u-direction ===
if tp > 0
    mu = mcp+p+1;
    p_ref = p+tp;
    ku = 1;
    for i=1:mu-1
        if uKnot(i) ~= uKnot(i+1)
            ku = ku+1;
        end
    end
    uKnot_ref = zeros(1,mu+tp*ku);

    bezalfs(1,1) = 1.0;
    bezalfs(p_ref+1,p+1) = 1.0;

    for i = 1:fix(p_ref/2)   
        inv = 1/nchoosek(p_ref,i);
        mpi = min(p,i);
        for j = max(0,i-tp):mpi   
            bezalfs(i+1,j+1) = inv*nchoosek(p,j)*nchoosek(tp,i-j);
            bezalfs(p_ref-i+1,p-j+1) = bezalfs(i+1,j+1);    
        end
    end

    mh = p_ref;
    kind = p_ref+1;
    r = -1;
    a = p+1;
    b = p+2;
    cind = 1;
    ua = uKnot(a);
    Qw(1,:,:) = Pw(1,:,:);

    for i = 1:p_ref+1
        uKnot_ref(i) = ua;
    end

    for i = 1:p+1
        bpts(i,:,:) = Pw(i,:,:);
    end

    while b < mu
        i = b;
        while b<mu && uKnot(b) == uKnot(b+1)
            b = b+1;      
        end
        mult = b-i+1;   
        mh = mh+mult+tp;
        ub = uKnot(b);
        oldr = r;       
        r = p-mult;     

        if oldr > 0
            lbz = fix((oldr+2)/2);
        else
            lbz = 1;
        end

        if r > 0
            rbz = p_ref-fix((r+1)/2);
        else
            rbz = p_ref;
        end

        alfs = zeros(1,p-mult);
        if r > 0
            numer = ub - ua;
            for k = p:-1:mult+1
                alfs(k-mult) = numer/(uKnot(a+k)-ua);   
            end
            for j = 1:r   
                save = r-j+1;
                s = mult+j;
                for k = p+1:-1:s+1    
                    bpts(k,:,:) = alfs(k-s)*bpts(k,:,:)+(1-alfs(k-s))*bpts(k-1,:,:);
                end
                Nextbpts(save,:,:) = bpts(p+1,:,:);
            end
        end

        for i = lbz:p_ref
            ebpts(i+1,1:ncp,1:4) = 0;
            mpi = min(p,i);
            for j = max(0,i-tp):mpi   
                ebpts(i+1,:,:) = ebpts(i+1,:,:) + bezalfs(i+1,j+1)*bpts(j+1,:,:);
            end
        end

        if oldr > 1
            first = kind-2;
            last = kind;
            den = ub-ua;
            bet = (ub-uKnot_ref(kind))/den;
            for tr = 1:(oldr-1)
                i = first;
                j = last;
                kj = j-kind+1;
                while (j-i) > tr
                    if i < cind
                        alf = (ub-uKnot_ref(i+1))/(ua-uKnot_ref(i+1));
                        Qw(i+1,:,:) = (alf*Qw(i+1,:,:)+(1.0-alf)*Qw(i,:,:));
                    end
                    if j >= lbz
                        if (j-tr) <= (kind-p_ref+oldr)
                          gam = (ub-uKnot_ref(j-tr+1))/den;
                          ebpts(kj+1,:,:) = gam*ebpts(kj+1,:,:) + (1.0-gam)*ebpts(kj+2,:,:);
                        else
                          ebpts(kj+1,:,:) = bet*ebpts(kj+1,:,:) + (1.0-bet)*ebpts(kj+2,:,:);
                        end
                    end
                    i = i+1;
                    j = j-1;
                    kj = kj-1;
                end
                first = first-1;
                last = last+1;
            end
        end   

        if a ~= p+1
            for i = 0:(p_ref-oldr-1)
                uKnot_ref(kind+1) = ua;
                kind = kind+1;
            end
        end

        for j = lbz:rbz   
            Qw(cind+1,:,:) =  ebpts(j+1,:,:);
            cind = cind +1;
        end

        if b < mu 
            for j = 0:r-1
                bpts(j+1,:,:) = Nextbpts(j+1,:,:);
            end
            for j = r:p
                bpts(j+1,:,:) = Pw(b-p+j,:,:);
            end
            a = b;
            b = b+1;
            ua = ub;
        else   
            for i = 0:p_ref
                uKnot_ref(kind+i+1) = ub;
            end
        end
    end
elseif tp == 0
    uKot_ref = uKnot;
    p_ref = p;
    Qw = nPw;
end

%% === Refine vKnot in v-direction ===
if tq > 0
    clear ('Pw');
    Pw(:,:,:)=Qw(:,:,:);
    mcp = length(Pw(:,1,1));
    clear ('Qw');  
    clear ('bpts');
    clear ('Nextbpts');
    clear ('ebpts');
    clear ('bezalfs');

    mv = ncp+q+1;
    q_ref = q+tq;  
    kv = 1;
    for i=1:mv-1
        if vKnot(i) ~= vKnot(i+1) 
            kv = kv+1;   
        end
    end
    vKnot_ref = zeros(1,mv+tq*kv);

    bezalfs(1,1) = 1.0;
    bezalfs(q_ref+1,q+1) = 1.0;

    for i = 1:fix(q_ref/2)  
        inv = 1/nchoosek(q_ref,i);
        mpi = min(q,i);
        for j = max(0,i-tq):mpi   
            bezalfs(i+1,j+1) = inv*nchoosek(q,j)*nchoosek(tq,i-j);
            bezalfs(q_ref-i+1,q-j+1) = bezalfs(i+1,j+1);    
        end
    end

    mh = q_ref;
    kind = q_ref+1;
    r = -1;
    a = q+1;
    b = q+2;
    cind = 1;
    ua = vKnot(a);
    Qw(:,1,:) = Pw(:,1,:);

    for i = 1:q_ref+1
        vKnot_ref(i) = ua;
    end

    for i = 1:q+1
        bpts(:,i,:) = Pw(:,i,:);
    end

    while b < mv
        i = b;
        while b < mv && vKnot(b) == vKnot(b+1)
            b = b+1;      
        end
        mult = b-i+1;   
        mh = mh+mult+tq;
        ub = vKnot(b);
        oldr = r;       
        r = q-mult;     

        if oldr > 0
            lbz = fix((oldr+2)/2);
        else
            lbz = 1;
        end

        if r > 0
            rbz = q_ref-fix((r+1)/2);
        else
            rbz = q_ref;
        end

        if r > 0
            numer = ub - ua;
            for k = q:-1:mult+1
                alfs(k-mult) = numer/(vKnot(a+k)-ua);   
            end
            for j = 1:r   
                save = r-j+1;
                s = mult+j;
                for k = q+1:-1:s+1    
                    bpts(:,k,:) = alfs(k-s)*bpts(:,k,:)+(1-alfs(k-s))*bpts(:,k-1,:);
                end
                Nextbpts(:,save,:) = bpts(:,q+1,:);
            end
        end

        for i = lbz:q_ref
            ebpts(1:mcp,i+1,1:4) = 0;
            mpi = min(q,i);
            for j = max(0,i-tq):mpi   
                ebpts(:,i+1,:) = ebpts(:,i+1,:) + bezalfs(i+1,j+1)*bpts(:,j+1,:);
            end
        end

        if oldr > 1
            first = kind-2;
            last = kind;
            den = ub-ua;
            bet = (ub-vKnot_ref(kind))/den;
            for tr = 1:(oldr-1)
                i = first;
                j = last;
                kj = j-kind+1;
                while ((j-i)>tr)
                    if (i<cind)
                        alf = (ub-vKnot_ref(i+1))/(ua-vKnot_ref(i+1));
                        Qw(:,i+1,:) = (alf*Qw(:,i+1,:)+(1.0-alf)*Qw(:,i,:));
                    end
                    if j >= lbz
                        if (j-tr) <= (kind-q_ref+oldr)
                          gam = (ub-vKnot_ref(j-tr+1))/den;
                          ebpts(:,kj+1,:) = gam*ebpts(:,kj+1,:) + (1.0-gam)*ebpts(:,kj+2,:);
                        else
                          ebpts(:,kj+1,:) = bet*ebpts(:,kj+1,:) + (1.0-bet)*ebpts(:,kj+2,:);
                        end
                    end
                    i = i+1;
                    j = j-1;
                    kj = kj-1;
                end
                first = first-1;
                last = last+1;
            end
        end   

        if a ~= q+1 
            for i = 0:(q_ref-oldr-1)
                vKnot_ref(kind+1) = ua;
                kind = kind+1;
            end
        end

        for j = lbz:rbz   
            Qw(:,cind+1,:) =  ebpts(:,j+1,:);
            cind = cind +1;
        end

        if b < mv  
            for j = 0:r-1
                bpts(:,j+1,:) = Nextbpts(:,j+1,:);
            end
            for j = r:q
                bpts(:,j+1,:) = Pw(:,b-q+j,:);
            end
            a = b;
            b = b+1;
            ua = ub;
        else   
            for i = 0:q_ref
                vKnot_ref(kind+i+1) = ub;
            end
        end
    end
elseif tq == 0
    vKnot_ref = vKnot;
    q_ref = q;
end

%% === Project 3D control points to 4D ===
for i = 1:length(Qw(:,1,1))
    for j = 1:length(Qw(1,:,1))
        CP_ref(i,j,1:3) = Qw(i,j,1:3) / Qw(i,j,4);
        CP_ref(i,j,4) = Qw(i,j,4);
    end
end
