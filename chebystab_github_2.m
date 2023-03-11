%chebyshev collcation method for the schrodinger equation with delta potential and the directly boundary condition
%[-a,a]->[-a,0]+[0,a]->[-1,1]+[-1,1]
%需改正的地方：1，求逆矩阵那里，2，D矩阵的求法，3，非线性项的离散
%orbital stability-2(nonlinear term-2)
clc;clear;
tic
Tt=10;a=20;
%% the controller of time and space step size
for fM=1:1
    fM
    % time
    N=10000;
    %N=10*2^(fN-1);
    tau=Tt/N;
    % space
    %M=40*2^(fM-1);
    M=200;
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
    for ip=1:1%
        ip
        %p=1.5+(ip-1)*1.5;
        p=ip+4;
        %% the exponent of the perturbation
        for iper=1:4
            omga=10^(-iper)*4;
           % omga=0.0001;
            muu=lamda^2/4+100;%2,5
            %muu=lamda^2/4*(1/(0.2922+ (iper-1)*0.0002) )^2;%6
            %muu=lamda^2/4*(1/(0.4407+ (iper-1)*0.0004) )^2;%7
            %muu=lamda^2/4*(1/(0.5666+ (iper-1)*0.0004) )^2;%8
            xi=-lamda/(2*sqrt(muu));
            % initial value
            u(:,1)=( (p+1)*muu/2*sech((p-1)/2*sqrt(muu)*abs(x)+ 1/2*logm( (1-lamda/(2*sqrt(muu)))/(1+lamda/(2*sqrt(muu)))  )  ).^2 ).^(1/(p-1));
            ue(:,1)=u(:,1);
            u(:,1)=u(:,1)+ omga*(1+0.5i)*exp(-(x+1).^2);
            %% orbital stability_1
            uer=u(:,1)-exp(1i*muu*0*tau)*ue(:,1);
            unum(iper,1)=sqrt( sum( ((D(2:M+1,2:M+1)*uer(1:M,1)).*conj(D(2:M+1,2:M+1)*uer(1:M,1))).*omega(2:M+1,1))...
                +sum( ((D(2:M,2:M)*uer(M+1:2*M-1,1)).*conj(D(2:M,2:M)*uer(M+1:2*M-1,1))).*omega(2:M,1)) );%一阶差商随时间的变化情况
            % ||u||_H^1
            ustab(iper,1)=sqrt( sum( ((D(2:M+1,2:M+1)*u(1:M,1)).*conj(D(2:M+1,2:M+1)*u(1:M,1))).*omega(2:M+1,1))...
                +sum( ((D(2:M,2:M)*u(M+1:2*M-1,1)).*conj(D(2:M,2:M)*u(M+1:2*M-1,1))).*omega(2:M,1)) );%一阶差商随时间的变化情况
            %  orbital stability_2
            ap=10000;
            for k=1:100
                uer2=u(:,1)-exp(1i*(k-1)/100*2*pi)*ue(:,1);
                ac=sqrt( sum( ((D(2:M+1,2:M+1)*uer2(1:M,1)).*conj(D(2:M+1,2:M+1)*uer2(1:M,1))).*omega(2:M+1,1))...
                    +sum( ((D(2:M,2:M)*uer2(M+1:2*M-1,1)).*conj(D(2:M,2:M)*uer2(M+1:2*M-1,1))).*omega(2:M,1)) );
                if  ac<ap
                    ap=ac;
                end
            end%
            ustab2(iper,1)=ap;
            
            %% computing
            for n=1:N
                %n
                ee=1;
                %iteration value
                non=zeros(2*M-1,1);%take the last level as iterative initial value
                u0=u(:,1);%*exp(1i);
                
                while ee > 10^(-12)
                    %nonlinear term
                    F=(@(th)(th*(abs(u0)).^2+(1-th)*abs(u(:,1)).^2).^(p/2-1/2));
                    non(:,1)=mu/2*quadv(F,0,1).*(u0+u(:,1));
                    non(M,1)=0;
                    non=non+ B*u(:,1);
                    u(:,2)=CC*non;
                    ee=norm(abs(u(:,2)-u0),'inf');
                    u0=u(:,2);
                end%while
                u(:,1)=u(:,2);
                %% orbital stability_1
                uer=u(:,1)-exp(1i*muu*n*tau)*ue(:,1);
                unum(iper,n+1)=sqrt( sum( ((D(2:M+1,2:M+1)*uer(1:M,1)).*conj(D(2:M+1,2:M+1)*uer(1:M,1))).*omega(2:M+1,1))...
                    +sum( ((D(2:M,2:M)*uer(M+1:2*M-1,1)).*conj(D(2:M,2:M)*uer(M+1:2*M-1,1))).*omega(2:M,1)) );%一阶差商随时间的变化情况
                % ||u||_H^1
                ustab(iper,n+1)=sqrt( sum( ((D(2:M+1,2:M+1)*u(1:M,1)).*conj(D(2:M+1,2:M+1)*u(1:M,1))).*omega(2:M+1,1))...
                    +sum( ((D(2:M,2:M)*u(M+1:2*M-1,1)).*conj(D(2:M,2:M)*u(M+1:2*M-1,1))).*omega(2:M,1)) );%一阶差商随时间的变化情况
                %%  orbital stability_2
                ap=10000;
                for k=1:100
                    uer2=u(:,1)-exp(1i*(k-1)/100*2*pi)*ue(:,1);
                    ac=sqrt( sum( ((D(2:M+1,2:M+1)*uer2(1:M,1)).*conj(D(2:M+1,2:M+1)*uer2(1:M,1))).*omega(2:M+1,1))...
                        +sum( ((D(2:M,2:M)*uer2(M+1:2*M-1,1)).*conj(D(2:M,2:M)*uer2(M+1:2*M-1,1))).*omega(2:M,1)) );
                    if  ac<ap
                        ap=ac;
                    end
                end%
                if ap<10000;
                    ustab2(iper,n+1)=ap;
                else
                    ustab2(iper,n+1)=NaN;
                end
                
            end%time
        end%iper
    end%ip
end%fM
toc


