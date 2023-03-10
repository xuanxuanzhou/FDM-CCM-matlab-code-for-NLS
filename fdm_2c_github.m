%time-dependent nonlinear schrodinger equation with dirac delta potential
%second-order in time and space
%ghost-point method
%参数与文章保持一致
%nonlinear term: numerical integration
%upload to github
clc;clear;
tic
Tt=1;a=20;
%% the controller of time and space step size
for fM=1:7
    fM
    % space
    M=800*20;
    %M=40*2^(fM-1)
    %time
    N=2*2^(fM-1);
    % N=1000;
    tau=Tt/N;h=a/M;
    % save the conservations( mass, energy)
    %mas=zeros(N+1,4);ene=mas;
    unum_0=zeros(2*M-1,1);
    %the value of numerical and exact function u,ue
    u=zeros(2*M-1,2);
    ue=zeros(2*M-1,1);non=ue;u0=ue;x=ue;eu=ue;
    x(:,1)=linspace(-a+h,a-h,2*M-1);
    % the coefficient of the nonlinear term and delta term
    mu=-1;lambda=-1;
    %  spatial matrix A,B Homogeneous Dirichlet Boundary
    E=ones(2*M-1,1);
    fm=zeros(2*M-1,1);fl=zeros(2*M-2,1);fr=fl;
    dm=fm;dl=fl;dr=fl;
    for k=1:2*M-1
        fm(k,1)=-2;
        if k<2*M-1
            fl(k,1)=1;
            fr(k,1)=1;
        end
    end
    fm(M,1)=fm(M,1)-lambda*h;
    fm=fm/h^2;fl=fl/h^2;fr=fr/h^2;
    dl=fl/2;dr=fr/2;dm=1i/tau*E+fm/2;
    fl=-fl/2;fr=-fr/2;fm=1i/tau*E-fm/2;
    %% the exponent of the nonlinear term
    for ip=1:4 %index_p
        p=1.5+(ip-1)*1.5;
        ip
        %% the exponent of the perturbation
        for iper=1:1  %index_perturbation
            %iper
            yta=0.4*2^(1-iper)*0;
            omga=lambda^2/4+1;
            % initial value
            u(:,1)=( (p+1)*omga/2*sech((p-1)/2*sqrt(omga)*abs(x)+ 1/2*logm( (1-lambda/(2*sqrt(omga)))/(1+lambda/(2*sqrt(omga)))  )  ).^2 ).^(1/(p-1));
            ue(:,1)=exp(1i*omga*Tt)*u(:,1);
            u(:,1)=u(:,1)+ yta*(1+0.5i)*exp(-(x+1).^2);
            u0=u(:,1);
            uxx=zeros(2*M-1,1);
            % mass and energy n=1
            %             for k=1:2*M-1
            %                 if k==1
            %                     uxx(k,1)=-2/h^2*u(k,1)+1/h^2*u(k+1,1);
            %                 elseif k==2*M-1
            %                     uxx(k,1)=1/h^2*u(k-1,1)-2/h^2*u(k,1);
            %                 else
            %                     uxx(k,1)=1/h^2*u(k-1,1)-2/h^2*u(k,1)+1/h^2*u(k+1,1);
            %                 end
            %             end
            %             for k=1:2*M-1
            %                 mas(1,ip)=mas(1,ip)+h*abs(u(k,1))^2;
            %                 ene(1,ip)=ene(1,ip)- h*uxx(k,1)*conj(u(k,1))+2*mu/(p+1)*h*abs(u(k,1))^(p+1) ;
            %             end
            %             ene(1,ip)=ene(1,ip)+lambda*abs(u(M,1))^2;
            %% computing
            for kk=1:N
                % kk
                non1=zeros(2*M-1,1);
                for k=1:2*M-1
                    if k==1
                        non1(k,1)=non1(k,1)+fm(k)*u(k,1)+fr(k)*u(k+1,1);
                    elseif k==2*M-1
                        non1(k,1)=non1(k,1)+fl(k-1)*u(k-1,1)+fm(k)*u(k,1);
                    else
                        non1(k,1)=non1(k,1)+fl(k-1)*u(k-1,1)+fm(k)*u(k,1)+fr(k)*u(k+1,1);
                    end
                end
                ee=1;
                %u0: Save the values from the previous time layer
                while ee>10^(-12)
                    %nonlinear term
                    F=(@(th)(th*(abs(u0)).^2+(1-th)*abs(u(:,1)).^2).^(p/2-1/2));
                    non(:,1)=mu/2*quadv(F,0,1).*(u0+u(:,1));
                    u(:,2)=trim1f(dl,dm,dr,non1+non);
                    ee=norm(u(:,2)-u0,'inf');
                    u0=u(:,2);
                end%while
                u(:,1)=u(:,2);
                % mass and energy n=kk
                %                 for k=1:2*M-1
                %                     if k==1
                %                         uxx(k,1)=-2/h^2*u(k,1)+1/h^2*u(k+1,1);
                %                     elseif k==2*M-1
                %                         uxx(k,1)=1/h^2*u(k-1,1)-2/h^2*u(k,1);
                %                     else
                %                         uxx(k,1)=1/h^2*u(k-1,1)-2/h^2*u(k,1)+1/h^2*u(k+1,1);
                %                     end
                %                 end
                %                 for k=1:2*M-1
                %                     mas(kk+1,ip)=mas(kk+1,ip)+h*abs(u(k,1))^2;
                %                     ene(kk+1,ip)=ene(kk+1,ip)- h*uxx(k,1)*conj(u(k,1)) +2*mu/(p+1)*h*abs(u(k,1))^(p+1) ;
                %                 end
                %                 ene(kk+1,ip)=ene(kk+1,ip)+lambda*abs(u(M,1))^2;
            end%kk
            %% error analysis
            eu=u(:,1)-ue;
            for k=1:2*M-1
                if k==1
                    uxx(k,1)=-2/h^2*eu(k,1)+1/h^2*eu(k+1,1);
                elseif k==2*M-1
                    uxx(k,1)=1/h^2*eu(k-1,1)-2/h^2*eu(k,1);
                else
                    uxx(k,1)=1/h^2*eu(k-1,1)-2/h^2*eu(k,1)+1/h^2*eu(k+1,1);
                end
            end
            ee1(fM,ip)=sqrt(h*sum(-uxx.*conj(eu)))%H_1
            ee2(fM,ip)=sqrt(h)*norm(eu,'fro')%L_2
            eeinfty(fM,ip)=norm(eu,'inf')%L_{\infty}
        end%iper
    end%ip
end%fM
toc


