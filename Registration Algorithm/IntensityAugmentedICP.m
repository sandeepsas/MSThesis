function [TR, TT, ER] = IntensityAugmentedICP(q,p,iter,wr)

if nargin<3
    iter=10;
end

% Actual implementation
Np = size(p,2);

% Transformed data point cloud
pt = p;

% Allocate vector for RMS of errors in every iteration.
% ER = zeros(iter+1,1);

% Initialize temporary transform vector and matrix.
T = zeros(3,1);

% T=[-0.6;17;0.8];
R = eye(3,3);

% Initialize total transform vector(s) and rotation matric(es).
%tr1=[-0.6;17;0.8];
%TT=repmat(tr1,1,iter1);
TT = zeros(3,1,iter+1);

TR = repmat(eye(3,3), [1,1,iter+1]);

 kdOBJ = KDTreeSearcher(transpose(q));
flag=0;
% Go into main iteration loop
for k=1:iter
%     if flag==1
%         break;
%     end
    
     disp(['Iteration no. =' num2str(k)]);
    % [match MinDistance] =bruteForce(q,pt);
    [match MinDistance] = match_kDtree(q,pt,kdOBJ);
      
    p_idx = true(1, Np);
    q_idx = match;
    
     if (wr>0 )%&& k>10)
        edge = round((1-wr)*sum(p_idx));
        pairs = find(p_idx);
        [~, idx] = sort(MinDistance);
        p_idx(pairs(idx(edge:end))) = false;
        q_idx = match(p_idx);
        MinDistance = MinDistance(p_idx);
    end
    
    if k == 1
        ER(k) = sqrt(sum(MinDistance.^2)/length(MinDistance));
    else
        tau=(ER(k)-ER(k-1));
        disp(tau);
        disp(ER(k));
        if (tau<0.00001)
            flag=1;
        end
    end
    % Determine weight vector
    weights=ones(1,length(match));
    
    [R,T] = eq_point(q(:,q_idx),pt(:,p_idx), weights(p_idx));
    
    % Add to the total transformation
    TR(:,:,k+1) = R*TR(:,:,k);
    TT(:,:,k+1) = R*TT(:,:,k)+T;
    
    % Apply last transformation
    p1=p(1:3,:);
    pt1=pt(1:3,:);
    pt1 = TR(:,:,k+1) * p1 + repmat(TT(:,:,k+1), 1, Np);
    pt(1:3,:)=pt1;
    p(1:3,:)=p1;
    
    % Root mean of objective function
    ER(k+1)= rms_error(q(:,q_idx), pt(:,p_idx));
    
end


TR = TR(:,:,end);
TT = TT(:,:,end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [match mindist] = match_kDtree(~, p, kdOBJ)
	[match mindist] = knnsearch(kdOBJ,transpose(p));
    match = transpose(match);

function [R,T] = eq_point(q,p,weights)

m = size(p,2);
n = size(q,2);

p=p(1:3,:);
q=q(1:3,:);
% normalize weights
 weights = weights ./ sum(weights);

% find data centroid and deviations from centroid
q_bar = q * transpose(weights);
q_mark = q - repmat(q_bar, 1, n);
% Apply weights
q_mark = q_mark .* repmat(weights, 3, 1);

% find data centroid and deviations from centroid
p_bar = p * transpose(weights);
p_mark = p - repmat(p_bar, 1, m);

% Apply weights
%p_mark = p_mark .* repmat(weights, 3, 1);

N = p_mark*transpose(q_mark); % taking points of q in matched order

[U,~,V] = svd(N); % singular value decomposition

R = V*transpose(U);

T = q_bar - R*p_bar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Determine the RMS error between two point equally sized point clouds with
% point correspondance.
% ER = rms_error(p1,p2) where p1 and p2 are 3xn matrices.

function ER = rms_error(p1,p2)
dsq = sum(power(p1 - p2, 2),1);
ER = sqrt(mean(dsq));
% disp(['Error. =' num2str(ER)]);


