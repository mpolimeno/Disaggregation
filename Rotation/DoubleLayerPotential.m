function DL = DoubleLayerPotential(ai,bj,myposint,myndir,myori,pos0,u_infty,M,flow)

    % k_hat
    k_hat = [0;0;1];

    DL = [0;0;0];
    
    xc = myposint; % center of each external face
    n_dir = myndir; % 1,2,3 for i,j,k
    n_ori = myori; % 1 or -1
    x0 = pos0;

    % get Q for ff face
    Q = rotation(n_dir,n_ori);
    Qt = Q';

    n_hat = Qt*k_hat; % normal vector for the mapped system
    
    % integrand
    Intg_ij = [0,0,0];
    for ii=1:length(ai)
        for jj=1:length(bj)

            % build gamma
            gammaij = [ai(ii,jj);bj(ii,jj);0];

            % shift to origin
            shift = xc;

            % build Tensor T
            Tl = (Qt*gammaij+shift')-x0';
            Ti = Tl';
            Ts = (shift'-x0');

            % background flow:
            if flow==1
                Ubg = u_infty;
            elseif flow==2
                Ubg = (cross(u_infty',(Qt*gammaij+shift')))';
            elseif flow==3
                x_ext = Qt*gammaij+shift';
                Ubg = M*x_ext;
            else
                x_shear = Qt*gammaij+shift';
                Ubg = (cross(u_infty',(Qt*gammaij+shift')))' + M*x_shear;
            end
            
            % numerator
            num = (Ubg'*Tl)*Ti*(Ts'*n_hat);

            % denominator
            den = norm((Qt*gammaij+shift')-x0')^5;

            % evaluate integrand
            Integrand = num/den;

            % Riemann sum
            Intg_ij = Intg_ij + Integrand;
        end
    end
    Integral = Intg_ij;
    DL = DL + Integral';
end
