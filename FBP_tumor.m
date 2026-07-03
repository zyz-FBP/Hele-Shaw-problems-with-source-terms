clear
clc
format long

%% visualization setting
is_viz = true;

%% Function approx
R = 3;%radius of \bar{\Omega}
M = 200; % Number of basis functions

%% Gauss-Legendre nodes and weights on [0,R]
Nq = 2000;   % Gaussian quadrature points
[xg, wg] = lgwt(Nq, 0, R);   % xg = quadrature nodes, wg = weights
r = xg;

% exact source term f_0=f-\bar{f}
f = 10*(besseli(0,r)/besseli(0,R) - 2*besseli(1,3)/(3*besseli(0,3)));

%% zeros of J1 using Newton's method
beta = zeros(M,1);
maxIter = 50;

for m = 1:M
    Jx = (m + 1/4)*pi;
    for k = 1:maxIter
        J0 = besselj(0,Jx);
        J1 = besselj(1,Jx);
        J2 = besselj(2,Jx);
        dJ1 = 0.5*(J0 - J2);
        Jxnew = Jx - J1/dJ1;

        if abs(Jxnew - Jx) < 1e-12
           break;
        end
        Jx = Jxnew;
    end
    beta(m) = Jx;
end
lambda = (beta/R).^2;
    
%% Neumann Bessel basis
Phi = zeros(Nq,M);
for m = 1:M
    Phi(:,m) = besselj(0, beta(m)*r/R);
end

%% Gaussian quadrature projection coefficients
b = zeros(M,1);
for m = 1:M
    num = sum(wg .* f .* Phi(:,m) .* r);
    den = sum(wg .* Phi(:,m).^2 .* r);
    b(m) = num / den;
end
%% p1 coefficients
a = b ./ lambda;

%% initial settings
N = [400]; % number of boundary points
D1 = 0.1;
D2= 5;
for rr = 1:length(N)
[x,y,theta] = ellipse(N(rr),D1,D2);
theta0 = theta;
%inital boundary
curve0 = [x,y];
curve = curve0;
%analytical tangent vectors and normal vectors
if D1 == 0
    tvec = [-2*sin(theta) 2*cos(theta)];
    nvec = [-2*cos(theta) -2*sin(theta)];
else
    [tvec,nvec] = tvec_analytic(theta,curve0,D1,D2);
end

told = tvec;
nold = nvec;

%time step size
dt = 1.0e-5;
%start time
time = 0.;
%end time
time_end = 0.25;
%total time steps
T = time_end/dt;
%time step count
tdt = 0;

while time <= time_end+dt%(time_end+dt/2.)
    [A,b,t_gmls,n_gmls,center_x,center_y,p1_prime] = velocity(curve,N,told,nold,a,R, beta,M);
    told = t_gmls';
    nold = -n_gmls';
    % velocity in normal direction
    d_curve = -A\b - p1_prime';
    vn = n_gmls'.*d_curve;

    %radius
    radius = sqrt((curve(:,1)-center_x).^2 + (curve(:,2)-center_y).^2);
    max_radius_record(tdt+1,1) = max(radius);
    min_radius_record(tdt+1,1) = min(radius);
    
    curve_old = curve;
    %forward euler
    curve = curve + dt.*(n_gmls'.*d_curve);

    time = time + dt;
    tdt = tdt + 1;

    if is_viz
    figure (1)%plot boundary evolution and its initial boundary
    plot(curve0(:,1), curve0(:,2),'k','LineWidth',2,'MarkerSize',5);
    hold on
    plot(curve(:,1), curve(:,2),'b--','LineWidth',2,'MarkerSize',5);
    hold off
    axis equal
    title(['Boundary Moton at t=',num2str(tdt)])
    pause(0.1)
    end
end
end

function [A,b,t_gmls,n_gmls,center_x,center_y,p1_prime] = velocity(curve,N,told,nold,a,R, beta,M)
rightpt = zeros(N,1);
leftpt = zeros(N,1);
A = zeros(N,N);
b_sint = zeros(N,1);
x = curve(:,1);
y = curve(:,2);
% number of k-nearest points
if mod(ceil(sqrt(N)),2)==0
    k = ceil(sqrt(N))+1;
else
    k = ceil(sqrt(N));
end
epsilon = 1.0e-12;
trial=0;
%approximate tangent vector and normal vector by local svd
[t_rough,n_rough,inds] = localsvd(curve,k,1,told,nold);
n_rough(1,:) = -t_rough(2,:);
n_rough(2,:) = t_rough(1,:);
n_rough = -n_rough;
%approximate tangent vector and normal vector by gmls
[t_gmls6, alpha6] = localgmls(curve,inds,1,t_rough,n_rough,6);
t_gmls = t_gmls6;
n_gmls = zeros(size(t_gmls));
n_gmls(1,:) = -t_gmls(2,:);
n_gmls(2,:) = t_gmls(1,:);
for pp = 1:N
     if (n_rough(:,pp)'*n_gmls(:,pp)<0)
        n_gmls(:,pp) = -n_gmls(:,pp);
    end
end
% make \alpha_1 close to 0
while max(abs(alpha6(1,:)'))>epsilon
    [t_gmls6, alpha6] = localgmls(curve,inds,1,t_gmls,n_gmls,6);
t_gmls = t_gmls6;
n_gmls = zeros(size(t_gmls));
n_gmls(1,:) = -t_gmls(2,:);
n_gmls(2,:) = t_gmls(1,:);
for pp = 1:N
     if (n_rough(:,pp)'*n_gmls(:,pp)<0)
        n_gmls(:,pp) = -n_gmls(:,pp);
    end
end
trial = trial+1;
end
[~, alpha6] = localgmls(curve,inds,1,t_gmls,n_gmls,6);
alpha = alpha6;
alpha(1,:) = 0; 

%approximate curvature at each point
cur6_gmls =-2*alpha6(2,:)'./((1+alpha6(1,:)'.^2).^(3/2));
cur_gmls = cur6_gmls;

% pressure p1 with eigenfunction approximation
[p1_approx,grad_p1_approx] = func_approx(a,curve,R, beta,M);
p1 = p1_approx - 10/4*(2*besseli(1,3)/(3*besseli(0,3))-1/2)*(x.^2+y.^2);
grad_p1 = grad_p1_approx + [-10/4*(2*besseli(1,3)/(3*besseli(0,3))-0.1)*2*x -10/4*(2*besseli(1,3)/(3*besseli(0,3))-0.1)*2*y];
%update the boundary condition for the system of p2
p1_prime = dot(n_gmls,grad_p1');
cur_gmls = cur_gmls - p1;

for qq = 1:N-1
  rightpt(qq) = t_gmls(:,qq)'*[x(qq+1)-x(qq);y(qq+1)-y(qq)];
end
rightpt(N) = t_gmls(:,N)'*[x(1)-x(N);y(1)-y(N)];

for qq = 2:N
    leftpt(qq) = t_gmls(:,qq)'*[x(qq-1)-x(qq);y(qq-1)-y(qq)];
end
leftpt(1) = t_gmls(:,1)'*[x(N)-x(1);y(N)-y(1)];

%% polynomial
poly_matrix = [rightpt rightpt.^2 rightpt.^3 rightpt.^4 rightpt.^5 rightpt.^6];
poly = diag(poly_matrix*alpha,0);
%p^2(s)/s^2
fovers_lmatrix = [ones(N,1) leftpt leftpt.^2 leftpt.^3 leftpt.^4 leftpt.^5];
fovers_leftpt = diag(fovers_lmatrix*alpha,0);
fovers_rmatrix = [ones(N,1) rightpt rightpt.^2 rightpt.^3 rightpt.^4 rightpt.^5];
fovers_rightpt = diag(fovers_rmatrix*alpha,0);

%poly for midpoint
poly_half_matrix = [rightpt./2 (rightpt./2).^2 (rightpt./2).^3 (rightpt./2).^4 (rightpt./2).^5 (rightpt./2).^6];
poly_mid = diag(poly_half_matrix*alpha,0);

%% polynomial_prime
%sqrt(1+(p'(0))^2)=1 as alpha1=0
fprime0 = sqrt(1 + (alpha(1,:))'.^2);
%sqrt(1+(p'(delta s))^2)
poly_rprime_matrix = [ones(N,1) 2*rightpt 3*rightpt.^2 4*rightpt.^3 5*rightpt.^4 6*rightpt.^5];
% poly_rprime_matrix = [ones(N,1) 2*rightpt 3*rightpt.^2 4*rightpt.^3 5*rightpt.^4 6*rightpt.^5 7*rightpt.^6 8*rightpt.^7];
poly_rprime = diag(poly_rprime_matrix*alpha,0);
fprime_rightpt = sqrt(1 + poly_rprime.^2);
%sqrt(1+(p'(-delta s))^2)
poly_lprime_matrix = [ones(N,1) 2*leftpt 3*leftpt.^2 4*leftpt.^3 5*leftpt.^4 6*leftpt.^5];
% poly_lprime_matrix = [ones(N,1) 2*leftpt 3*leftpt.^2 4*leftpt.^3 5*leftpt.^4 6*leftpt.^5 7*leftpt.^6 8*leftpt.^7];
poly_lprime = diag(poly_lprime_matrix*alpha,0);
fprime_leftpt = sqrt(1 + poly_lprime.^2);

%centroid
int_x = (curve(:,1) + [curve(2:end,1);curve(1,1)].*fprime_rightpt).*rightpt./2;
int_s = (ones(N,1) + fprime_rightpt).*rightpt./2;
int_y = (curve(:,2) + [curve(2:end,2);curve(1,2)].*fprime_rightpt).*rightpt./2;
center_x = sum(int_x)/sum(int_s);
center_y = sum(int_y)/sum(int_s);

%% matrix A
for j=1:N
    %% Distance matrix and its norm
    D_matrix = [x(j),y(j)] - curve;
    norm_matrix = vecnorm(D_matrix,2,2);
   
%% Matrix A
norm_left = -1/(2*pi).*log(norm_matrix);
norm_right = [norm_left(2:end);norm_left(1)];
G_left = norm_left.*rightpt;
G_mid_right = norm_right.*fprime_rightpt.*rightpt;
G_right = [G_mid_right(end);G_mid_right(1:end-1)];
A(j,:) = (G_left + G_right)'./2;

if j==1
    A(j,N) = G_right(N)/2;
else
    A(j,j-1) = G_right(j-1)/2;
end
if j==N
    A(j,1) = G_left(1)/2;
else
    A(j,j+1) = G_left(j+1)/2;
end
fovers_lterm = - 1/(4*pi)*log(1+(fovers_leftpt(j))^2)*fprime_leftpt(j)*abs(leftpt(j))/2;
fovers_rterm = - 1/(4*pi)*log(1+(fovers_rightpt(j))^2)*fprime_rightpt(j)*rightpt(j)/2;
% singular terms
int_beta0 = -1/(2*pi)*(-leftpt(j)*log(-leftpt(j))+leftpt(j)+rightpt(j)*log(rightpt(j))-rightpt(j));
int_beta1 = -1/(2*pi)*(-1/2*leftpt(j)^2*log(-leftpt(j))+1/4*(leftpt(j))^2+1/2*rightpt(j)^2*log(rightpt(j))-1/4*rightpt(j)^2);
int_beta2 = -1/(2*pi)*(1/3*(-leftpt(j))^3*log(-leftpt(j))-1/9*(-leftpt(j))^3+1/3*rightpt(j)^3*log(rightpt(j))-1/9*rightpt(j)^3);

%psi_i(s) = beta0 + beta1*s + beta2*s^2
A(j,j) = int_beta0 - int_beta1*(rightpt(j)+leftpt(j))/(rightpt(j)*leftpt(j)) + int_beta2/(rightpt(j)*leftpt(j));
if j==1
    A(j,N) = A(j,N) + fovers_lterm - int_beta1*rightpt(j)^2*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    - int_beta2*rightpt(j)*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j)));
else
    A(j,j-1) = A(j,j-1) + fovers_lterm - int_beta1*rightpt(j)^2*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    - int_beta2*rightpt(j)*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j)));
end
if j==N
    A(j,1) = A(j,1) + fovers_rterm + int_beta1*leftpt(j)^2*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    + int_beta2*leftpt(j)*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j))) ;
else
    A(j,j+1) = A(j,j+1) + fovers_rterm + int_beta1*leftpt(j)^2*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    + int_beta2*leftpt(j)*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j))) ;
end

%% b
G_prime = 1/(2*pi)*D_matrix./(norm_matrix.^2);
G_prime(j,:) = 0;

norm_left_prime = G_prime;
norm_right_prime = [norm_left_prime(2:end,:);norm_left_prime(1,:)];

%simpson rule
norm_lprime = [norm_left_prime(end,:);norm_left_prime(1:end-1,:)];
norm_mprime = norm_left_prime;
norm_rprime = norm_right_prime;
G_lprime = norm_lprime.*abs(rightpt-leftpt).*fprime_leftpt.*[cur_gmls(end);cur_gmls(1:end-1)];
G_mprime = norm_mprime.*abs(rightpt-leftpt).*cur_gmls;
G_rprime = norm_rprime.*abs(rightpt-leftpt).*fprime_rightpt.*[cur_gmls(2:end);cur_gmls(1)];
if j==1
    b_sint(j) = 1/6*sum(dot(n_gmls,([G_lprime(2:end,:);G_lprime(1,:)]' + 4.*G_mprime' + [G_rprime(end,:);G_rprime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*(4/6*abs(rightpt(j)-leftpt(j))+1/6*fprime_rightpt(N)*abs(rightpt(N)-leftpt(N))+1/6*fprime_leftpt(j+1)*abs(rightpt(j+1)-leftpt(j+1)));
elseif j==N
    b_sint(j) = 1/6*sum(dot(n_gmls,([G_lprime(2:end,:);G_lprime(1,:)]' + 4.*G_mprime' + [G_rprime(end,:);G_rprime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*(4/6*abs(rightpt(j)-leftpt(j))+1/6*fprime_rightpt(j-1)*abs(rightpt(j-1)-leftpt(j-1))+1/6*fprime_leftpt(1)*abs(rightpt(1)-leftpt(1)));
else
    b_sint(j) = 1/6*sum(dot(n_gmls,([G_lprime(2:end,:);G_lprime(1,:)]' + 4.*G_mprime' + [G_rprime(end,:);G_rprime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*(4/6*abs(rightpt(j)-leftpt(j))+1/6*fprime_rightpt(j-1)*abs(rightpt(j-1)-leftpt(j-1))+1/6*fprime_leftpt(j+1)*abs(rightpt(j+1)-leftpt(j+1)));
end
end
%simpson rule
b = b_sint./2 + cur_gmls./2;
end

% initial point clouds
function [x,y,theta] = ellipse(n1,D1,D2)
%uniform random distributed angles
theta = [0:2*pi/n1:2*pi-2*pi/n1]';
if D1==0
    r = 2.0;
else
r =1.0+ D1.*cos(D2.*theta);
end
x = r .* cos(theta);
y = r .* sin(theta);
end

% approximate p1 and gradient of p1 at boundary
function [p1_approx,grad_p1_approx] = func_approx(a,curve,R, beta,M)
%% approximation
curve_x = curve(:,1);
curve_y = curve(:,2);
curve_r = sqrt(curve_x.^2 + curve_y.^2);

Phi_curve = zeros(length(curve_r),M);
grad_Phi_x = zeros(length(curve_r),M);
grad_Phi_y = zeros(length(curve_r),M);
for m = 1:M
    Phi_curve(:,m) = besselj(0, beta(m)*curve_r/R);
    grad_Phi_x(:,m) = -beta(m)/R*besselj(1,beta(m)*curve_r/R).*curve_x./curve_r;
    grad_Phi_y(:,m) = -beta(m)/R*besselj(1,beta(m)*curve_r/R).*curve_y./curve_r;
end
p1_approx = Phi_curve * a;
grad_p1_approx = [grad_Phi_x*a grad_Phi_y*a];
end

%% Compute the Legendre-Gauss nodes and weights on an interval [a,b] with
%% truncation order N (from matlab)
function [x,w]=lgwt(N,a,b)
N=N-1;
N1=N+1;
N2=N+2;
xu=linspace(-1,1,N1)';
%initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2)) +(0.27/N1)*sin(pi*xu*N/N2);
%Legendre-Gauss Vnadermonde Matrix
L=zeros(N1,N2);
%Derivative of LGVM
Lp=zeros(N1,1);
%Compute the zeros of the N+1 Legendre Poly using the recursion relation
%and the Newton-Raphson method
y0=2;
%Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0)) > eps
    L(:,1)=1;
    L(:,2)=y;
    for k=2:N1
        L(:,k+1)=((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end
    Lp=N2*(L(:,N1)-y.*L(:,N2))./(1-y.^2);
    y0=y;
    y=y0-L(:,N2)./Lp;
end
%Linear map from [-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;
%Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end

% analytical tangent vectors and normal vectors for initial boundary curve
function [tvec,nvec] =tvec_analytic(theta,x,D1,D2)
N = size(x,1);
n = size(x,2);
tvec = zeros(N,n);
nvec = zeros(N,n);
tvec(:,1) = - sin(theta).*(D1*cos(D2*theta) + 1) - D1*D2*sin(D2*theta).*cos(theta);
tvec(:,2) = cos(theta).*(D1*cos(D2*theta) + 1) - D1*D2*sin(D2*theta).*sin(theta);

normtvec = sqrt(tvec(:,1).^2+tvec(:,2).^2);

tvec = tvec./repmat(normtvec,1,n);
nvec(:,1) = -((cos(theta) + D1*cos(D2*theta).*cos(theta) - D1*D2*sin(D2*theta).*sin(theta)).*(2*D1*cos(D2*theta) + 2*D1^2*D2^2 + D1^2*cos(D2*theta).^2 - D1^2*D2^2*cos(D2*theta).^2 + D1*D2^2*cos(D2*theta) + 1))./(2*D1*cos(D2*theta) + D1^2*D2^2 + D1^2*cos(D2*theta).^2 - D1^2*D2^2*cos(D2*theta).^2 + 1).^(3/2);
nvec(:,2) = -((sin(theta) + D1*cos(D2*theta).*sin(theta) + D1*D2*sin(D2*theta).*cos(theta)).*(2*D1*cos(D2*theta) + 2*D1^2*D2^2 + D1^2*cos(D2*theta).^2 - D1^2*D2^2*cos(D2*theta).^2 + D1*D2^2*cos(D2*theta) + 1))./(2*D1*cos(D2*theta) + D1^2*D2^2 + D1^2*cos(D2*theta).^2 - D1^2*D2^2*cos(D2*theta).^2 + 1).^(3/2);
normnvec = sqrt(nvec(:,1).^2+nvec(:,2).^2);
nvec = nvec./repmat(normnvec,1,n);
end

% build multi-indices for multivariate polynomials if necessary
function index =  generatemultiindex(N,dim)
% input
% N : max degree of polynomials
% dim : dimension

P = nchoosek(N+dim,N);
index = zeros(dim,P);

Ntotal = (N+1)^dim;
allindex = zeros(dim,Ntotal);

for i=1:dim
    nskip = (N+1)^(dim-i);
    for k=1:Ntotal/nskip/(N+1)
        for j=1:N+1;
            allindex(i,(k-1)*nskip*(N+1)+(j-1)*nskip+1:(k-1)*nskip*(N+1)+j*nskip) = j-1;
        end
    end
end

index1 = find(sum(allindex)<=N);
index = allindex(:,index1);
end

% find k-nearest points (can use the built-in function (knn) directly)
function [ds,inds] = knnCPU(R,Q,k,maxMem)
% Find the k nearest neighbors for each element of the query set Q among
% the points in the reference set R
% Q is Nxd where N is the number of query points and d is the dimension
% R is Mxd where M is the number of reference points and d is the dimension
% k is the number of nearest neighbors desired
% maxMem is the number of GB of ram knnCPU can use
% ds is Nxk and contains the distances to the k nearest neighbors.
% inds is Nxk and contains the indices of the k nearest neighbors.

    if (nargin<4)
        maxMem = 2;
    end

    M = size(R,1);
    N = size(Q,1);    
    
    maxArray = (maxMem*2500)^2;
    blockSize = floor(maxArray/M);
    blocks = floor(N/blockSize);

    ds = zeros(N,k);
    inds = zeros(N,k);

    Nr = sum(R.^2,2);
    Nq = sum(Q.^2,2);

    for b = 1:blocks
        dtemp = -2*R*Q((b-1)*blockSize + 1:b*blockSize,:)';
        dtemp = bsxfun(@plus,dtemp,Nr);
        dtemp = bsxfun(@plus,dtemp,Nq((b-1)*blockSize + 1:b*blockSize)');
        [dst,indst] = sort(dtemp,1);
        ds((b-1)*blockSize + 1:b*blockSize,:) = dst(1:k,:)';
        inds((b-1)*blockSize + 1:b*blockSize,:) = indst(1:k,:)';
    end

    if (blocks*blockSize < N)
        dtemp = -2*R*Q(blocks*blockSize + 1:N,:)';
        dtemp = bsxfun(@plus,dtemp,Nr);
        dtemp = bsxfun(@plus,dtemp,Nq(blocks*blockSize + 1:N)');
        [dst,indst] = sort(dtemp,1);
        ds(blocks*blockSize + 1:N,:) = dst(1:k,:)';
        inds(blocks*blockSize + 1:N,:) = indst(1:k,:)';
    end

    ds = real(sqrt(ds));
end

% GMLS approximation of tangent vectors
function [t_tilde,alpha] = localgmls(x,inds,d,t_rough,n_rough,ell)
N = size(x,1);
n = size(x,2);
k0 = size(inds,2);

if d > 1
    index =  generatemultiindex(ell,d); % index is d*term    
    index = index(:,2:end);
else
    index = 1:ell;
end
term = size(index,2);

ind_alin = zeros(d,1);
for jj = 1:d
    ind_alin(jj) = find((sum(index,1)==1)&(index(jj,:)==1));
end

for pp=1:N
    dx = x(inds(pp,:),:)';
    dx = dx - repmat(x(pp,:)',1,k0); 

     lcx = zeros(k0,n);
     lcx = dx'*[t_rough(:,pp) n_rough(:,pp)];
                
     %%% Phi = [1 yy yy.^2 yy.^3 ...]
     Phi = ones(k0,term); % k*term
        
     %%% constant has 1, deg 1 has d, deg 2 has (d+1)*d/2
     for ss = 1:term
        for rr = 1:d
            Phi(:,ss) = Phi(:,ss).*lcx(:,rr).^index(rr,ss);
        end
     end
     % coefficients in polynomial by least squares
     a1 =(Phi'*Phi)\(Phi'*lcx(:,n));
     alpha(:,pp) = a1;
     % improved tangent vectors
     t_tilde(:,pp) = t_rough(:,pp) + a1(ind_alin)*n_rough(:,pp);
     t_tilde(:,pp) = t_tilde(:,pp)/vecnorm(t_tilde(:,pp));     
end
end

% Get approximate tangent vectors and normal vectors by local SVD method
function [t_rough,n_rough,inds] = localsvd(x,k0,d,tvec,nvec)
[~,inds] = knnCPU(x,x,k0);
N = size(x,1);
n = size(x,2);
   
for pp = 1:N 
    dx = x(inds(pp,:),:)';
    dx = dx - repmat(x(pp,:)',1,k0); 
    [U,~,~] = svd(dx); 
    
    t_rough(:,pp) = U(:,d);        
    n_rough(:,pp) = U(:,d+1:n);
   
    if (tvec(pp,:)*t_rough(:,pp)<0)
        t_rough(:,pp) = -t_rough(:,pp);
    end

    if (nvec(pp,:)*n_rough(:,pp)<0)
        n_rough(:,pp) = -n_rough(:,pp);
    end
end
end
