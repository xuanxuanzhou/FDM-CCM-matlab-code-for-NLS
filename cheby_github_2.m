%chebyshev collcation method for the schrodinger equation with delta potential and the directly boundary condition
%[-a,a]->[-a,0]+[0,a]->[-1,1]+[-1,1]
%循环内收敛过程三步到位
%需改正的地方：1，求逆矩阵那里，2，D矩阵的求法，3，非线性项的离散
clc;clear;
tic
Tt=1;a=20;
%% the controller of time and space step size
for fM=1:1
    fM
    % time
    N=100000;
    %N=10*2^(fM-1);
    tau=Tt/N;
    % space
    M=20*2^(fM-1);
    %M=400;
    %the value of numerical and exact function u,ue
    u=zeros(2*M-1,2);
    ue=zeros(2*M-1,1);u0=ue;eu=ue;
    % the coefficient of the nonlinear term and delta term
    mu=-1;lamda=-1;
    % spatial discretization
    h=zeros(M+1,1);hl=zeros(1,M);hr=zeros(1,M-1);hlx=hl;hrx=hr;
    x=zeros(2*M-1,1);
    for k=1:M+1
        %-1~1
        h(k)=-cos(pi*(k-1)/M);
    end
    hl=h(2:M+1)-1;
    hr=h(2:M)+1;
    hlx=a/2*hl;
    hrx=a/2*hr;
    x(1:M,1)=hlx;
    x(M+1:2*M-1,1)=hrx;
    % the first-derivation matrix
    c=zeros(M+1,1);omega=zeros(M+1,1);
    D=zeros(M+1,M+1);D_0=D;
    D_2=zeros(M-1,M-1);p1=zeros(M-1,1);p2=p1;
    for k=1:M-1
        c(k+1)=1;
    end
    c(1)=2;c(M+1)=2;
    
    for k=1:M+1
        omega(k)=pi/(c(k)*M);
        for kk=1:M+1
            if k~=kk
                %D(k,kk)=c(k)*(-1)^(k+kk)/(c(kk)*(h(k)-h(kk)));
                D(k,kk)=-c(k)*(-1)^(k+kk)/(2*c(kk)*sin((k+kk-2)*pi/(2*M))*sin((kk-k)*pi/(2*M)));
            elseif k~=1 && k~=M+1
                %D(k,kk)=-h(k)/(2*(1-h(k)^2));
                D(k,kk)=-h(k)/(2*sin((k-1)*pi/M)^2);
            elseif k==1
                D(k,kk)=-(2*M^2+1)/6;
            elseif k==M+1
                D(k,kk)=(2*M^2+1)/6;
            else
                disp('0');
            end
        end
    end
    D=2/a*D;
    D_0=D*D;D_2=D_0(2:M,2:M);
    p1(:,1)=1/2*D_0(2:M,M+1);p2(:,1)=1/2*D_0(2:M,1);
    %边界条件转化为未知向量的向量 inner boundary condition
    aa=zeros(2*M-1,1);
    aa(1:M-1,1)=-D(M+1,2:M);aa(M+1:2*M-1,1)=D(1,2:M);
    aa(M,1)=-D(M+1,M+1)+D(1,1)-lamda;
    % matrix
    E=eye(M-1);A=zeros(2*M-1,2*M-1);B=A;CC=A;
    A(1:M-1,1:M-1)=1i/tau*E+1/2*D_2;
    A(1:M-1,M)=p1;
    A(M,:)=aa;
    A(M+1:2*M-1,M)=p2;
    A(M+1:2*M-1,M+1:2*M-1)=1i/tau*E+1/2*D_2;
    B(1:M-1,1:M-1)=1i/tau*E-1/2*D_2;
    B(1:M-1,M)=-p1;
    B(M,:)=-aa;
    B(M+1:2*M-1,M)=-p2;
    B(M+1:2*M-1,M+1:2*M-1)=1i/tau*E-1/2*D_2;
    CC=inv(A);
    %%  the exponent of the nonlinear term
    for ip=1:4%
        ip
        p=1.5+(ip-1)*1.5;
        %% the exponent of the perturbation
        for iper=1:1
            omga=0.4*2^(1-iper)*0;
            muu=lamda^2/4+1;
            xi=-lamda/(2*sqrt(muu));
            % initial value
            u(:,1)=( (p+1)*muu/2*sech((p-1)/2*sqrt(muu)*abs(x)+ 1/2*logm( (1-lamda/(2*sqrt(muu)))/(1+lamda/(2*sqrt(muu)))  )  ).^2 ).^(1/(p-1));
            ue(:,1)=u(:,1);
            u(:,1)=u(:,1)+ omga*(1+0.5i)*exp(-(x+1).^2);
            %% computing
            for n=1:N
                %n
                ee=1;
                %iteration value
                non=zeros(2*M-1,1);%take the last level as iterative initial value
                u0=u(:,1);%*exp(1i);
                while ee > 10^(-12)
                    %nonlinear term
                    %slowly
                    %                     for k=1:2*M-1
                    %                         if abs(u(k,1))==abs(u0(k))
                    %                             %F(k,1)=0;
                    %                             non(k,1)=mu/2*(abs(u0(k)))^(p-1)*(u0(k)+u(k,1));
                    %                         else
                    %                             non(k,1)=mu/(p+1)*(-abs(u(k,1))^(p+1)+abs(u0(k))^(p+1))/(-abs(u(k,1))^2+abs(u0(k))^2)*(u0(k)+u(k,1));
                    %                         end
                    %                     end
                    F=(@(th)(th*(abs(u0)).^2+(1-th)*(abs(u(:,1))).^2).^(p/2-1/2));
                    non(:,1)=mu/2*quadv(F,0,1).*(u0+u(:,1));
                    non(M,1)=0;
                    non=non+ B*u(:,1);
                    u(:,2)=CC*non;
                    ee=norm(abs(u(:,2)-u0),'inf');
                    u0=u(:,2);
                end%while
                u(:,1)=u(:,2);
            end%time
            %% error analysis
            %chebyshev 为非均匀网格剖分，范数定义与以往不同
            eu=u(:,2)-exp(1i*muu*1)*ue;
            err_h1(ip,fM)=sqrt( sum( ((D(2:M+1,2:M+1)*eu(1:M,1)).*conj(D(2:M+1,2:M+1)*eu(1:M,1))).*omega(2:M+1,1))...
                +sum( ((D(2:M,2:M)*eu(M+1:2*M-1,1)).*conj(D(2:M,2:M)*eu(M+1:2*M-1,1))).*omega(2:M,1)) );
            err_l2(ip,fM)=sqrt(  sum( (eu(1:M,1).*conj(eu(1:M,1))).*omega(2:M+1,1))...
                +sum( (eu(M+1:2*M-1,1).*conj(eu(M+1:2*M-1,1))).*omega(2:M,1))  );
            err_infty(ip,fM)=norm(eu,'inf')
        end%iper
    end%ip
end%fM
toc


