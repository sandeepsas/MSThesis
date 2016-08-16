function [T R TR, TT, ER, t] = GeometricICP(q,p,iter,wr)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual implementation

% Allocate vector for RMS of errors in every iteration.
t = zeros(iter+1,1);

% Start timer
tic;

Np = size(p,2);

% Transformed data point cloud
pt = p;

% Allocate vector for RMS of errors in every iteration.
% ER = zeros(iter+1,1);

% Initialize temporary transform vector and matrix.%
T = zeros(3,1);
% T=[0;-20;0];
R = eye(3,3);

% Initialize total transform vector(s) and rotation matric(es).
TT = zeros(3,1, iter+1);
TR = repmat(eye(3,3), [1,1, iter+1]);

Normals = lsqnormest(q,4);

kdOBJ = KDTreeSearcher(transpose(q));

t(1) = toc;
flag=0;
% Go into main iteration loop
for k=1:iter
%     if flag
%         break;
%     end
    
    disp(k);
    [match mindist] = match_kDtree(q,pt,kdOBJ);
    
    p_idx = true(1, Np);
    q_idx = match;
    
    if (wr>0 )
        edge = round((1-wr)*sum(p_idx));
        pairs = find(p_idx);
        [~, idx] = sort(mindist);
        p_idx(pairs(idx(edge:end))) = false;
        q_idx = match(p_idx);
        mindist = mindist(p_idx);
    end
    
    if k == 1
        ER(k) = sqrt(sum(mindist.^2)/length(mindist));
    else
        tau=abs(ER(k)-ER(k-1));
        disp(tau);
        if (tau<0.00001)
            flag=1;
        end
    end
    % Determine weight vector
    weights=ones(1,length(match));
    [R,T] = eq_plane(q(:,q_idx),pt(:,p_idx),Normals(:,q_idx),weights(p_idx));
    
    
    % Add to the total transformation
    TR(:,:,k+1) = R*TR(:,:,k);
    TT(:,:,k+1) = R*TT(:,:,k)+T;
    
    % Apply last transformation
    pt = TR(:,:,k+1) * p + repmat(TT(:,:,k+1), 1, Np);
    
    % Root mean of objective function
    ER(k+1)= rms_error(q(:,q_idx), pt(:,p_idx));
    disp(ER(k+1));
    % If Extrapolation, we might be able to move quicker
    t(k+1) = toc;
end


TR = TR(:,:,end);
TT = TT(:,:,end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match mindist] = match_kDtree(~, p, kdOBJ)
[match mindist] = knnsearch(kdOBJ,transpose(p));
match = transpose(match);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T] = eq_plane(q,p,n,weights)

n = n .* repmat(weights,3,1);

c = cross(p,n);

cn = vertcat(c,n);

C = cn*transpose(cn);

b = - [sum(sum((p-q).*repmat(cn(1,:),3,1).*n));
    sum(sum((p-q).*repmat(cn(2,:),3,1).*n));
    sum(sum((p-q).*repmat(cn(3,:),3,1).*n));
    sum(sum((p-q).*repmat(cn(4,:),3,1).*n));
    sum(sum((p-q).*repmat(cn(5,:),3,1).*n));
    sum(sum((p-q).*repmat(cn(6,:),3,1).*n))];

X = C\b;

cx = cos(X(1)); cy = cos(X(2)); cz = cos(X(3));
sx = sin(X(1)); sy = sin(X(2)); sz = sin(X(3));

R = [cy*cz cz*sx*sy-cx*sz cx*cz*sy+sx*sz;
    cy*sz cx*cz+sx*sy*sz cx*sy*sz-cz*sx;
    -sy cy*sx cx*cy];

T = X(4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the RMS error between two point equally sized point clouds with
% point correspondance.
% ER = rms_error(p1,p2) where p1 and p2 are 3xn matrices.

function [ER] = rms_error(p1,p2)
dsq = sum(power(p1 - p2, 2),1);
ER = sqrt(mean(dsq));
subd=p1-p2;
error_x=mean(subd(1,:));
error_y=mean(subd(2,:));
error_z=mean(subd(3,:));
error_sub=[error_x,error_y,error_z];
std_x=std(subd(1,:));
std_y=std(subd(2,:));
std_z=std(subd(3,:));
error_SD=[std_x,std_y,std_z];




function n = lsqnormest(p, k)
m = size(p,2);
n = zeros(3,m);

v = ver('stats');
if str2double(v.Version) >= 7.3
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));
else
    neighbors = kNearestNeighbors(p, p, k+1);
end

for i = 1:m
    x = p(:,neighbors(2:end, i));
    p_bar = 1/k * sum(x,2);
    
    P = (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %spd matrix P
    %P = 2*cov(x);
    
    [V,D] = eig(P);
    
    [~, idx] = min(diag(D)); % choses the smallest eigenvalue
    
    n(:,i) = V(:,idx);   % returns the corresponding eigenvector
end

