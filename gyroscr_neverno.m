function [OUTF, OUTJ, p, Eff, Omega, ConLow, jout] = gyroscr(Nz, Nzi, Nt, Ne, ZAxis, ZAxisi, TAxis, Delta, Ic, dt, dz, dzi, tol, kpar2, INTT, INTZ, OUTNz, OUTNt, InitialField) %#codegen

WR = complex(zeros(Nt,1));
FNz = complex(zeros(Nt,1));
FNzm1 = complex(zeros(Nt,1));
JNz = complex(zeros(Nt,1));
JNzm1 = complex(zeros(Nt,1));
SigmaNz = complex(zeros(Nt,1));
SigmaNzm1 = complex(zeros(Nt,1));
WL = complex(zeros(Nt,1));
F0 = complex(zeros(Nt,1));
F1 = complex(zeros(Nt,1));
J0 = complex(zeros(Nt,1));
J1 = complex(zeros(Nt,1));
Sigma0 = complex(zeros(Nt,1));
Sigma1 = complex(zeros(Nt,1));
fmax = zeros(Nt, 1);
jmax = zeros(Nt, 1);
field = complex(zeros(Nz,1));
testfield = complex(zeros(Nz,1));
field_p = complex(zeros(Nz,1));
rfield_p = complex(zeros(Nz,1));
lfield_p = complex(zeros(Nz,1));
cu_p = complex(zeros(Nz,1));
cu = complex(zeros(Nz,1));
p_p = complex(zeros(Nzi,Ne));
J_p = complex(zeros(Nz,1));
J = complex(zeros(Nz,1));
OUTF = complex(zeros(OUTNz, OUTNt));
OUTJ = complex(zeros(OUTNz, OUTNt));
Eff = zeros(Nt,1);
Omega = zeros(Nt,1);
ConLow = zeros(Nt,1);
ZAxis_ipart = ZAxis(ZAxis <= ZAxisi(end));
ZAxis_ipart_end = length(ZAxis_ipart);
% theta = zeros(Nz, Ne);
% p = zeros(Nz, Ne);
% pv = zeros(Nz, 2*Ne);
% p0 = zeros(Ne,1);
% p0v = zeros(2*Ne,1);
% reidx = zeros(1,Ne);
% imidx = zeros(1,Ne);

if INTZ > 1
    IZ = 0:INTZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);
else
    IZ = 1:INTZ:length(ZAxis);
end

% kpar2 = zeros(length(ZAxis),1);
% N = length(ZAxis);
A = complex(zeros(Nz,1));
B = complex(zeros(Nz-1,1));
C = complex(zeros(Nz-1,1));
D = complex(zeros(Nz,1));

% SQR2 = sqrt(2.0D0);
SQRT2M2 = 2.828427124746190;
SQRT2D2 = 0.707106781186548;
SQRTDT = sqrt(dt);
SQRDZ = dz*dz;
SQRT2M2minus2p5 = SQRT2M2 - 2.5;

Nzm1 = Nz - 1;
C0 = 1.0D0;
CL = -1i*kpar2(1);
CR = -1i*kpar2(Nz);
% C1 = 1.0D0/sqrt(1i*pi);
C1 = 0;
C2 = 1.0D0/sqrt(1i*pi);
W0 = ((-1i*2.0D0/3.0D0*C0*dz/dt) - 1.0D0/dz);
W1 = ((-1i*C0/3.0D0*dz/dt) + 1.0D0/dz);
WNz = -((-1i*2.0D0/3.0D0*C0*dz/dt) - 1.0D0/dz);
WNzm1 = -((-1i*C0/3.0D0*dz/dt) + 1.0D0/dz);

A(1) =  1.0D0 - 4.0D0/3.0D0*C1*W0*SQRTDT;
A(2:Nzm1) = -2.0D0*(1.0D0 + 1i * SQRDZ/dt*C0);
A(Nz) = 1.0D0 + 4.0D0/3.0D0*C2*WNz*SQRTDT;
B(1) =   -4.0D0/3.0D0*C1*W1*SQRTDT;
B(2:Nzm1) = 1.0D0;
C(1:Nzm1 - 1) = 1.0D0;
C(Nzm1) = 4.0D0/3.0D0*C2*WNzm1*SQRTDT;

M = spdiags([[C; 0] A [0 ;B]], -1:1, Nz, Nz);

% Initial values
jout = 1;
field(:,1) = InitialField;
Sf = griddedInterpolant(ZAxis, field,'spline');
FforP = Sf(ZAxisi(1:Nzi-1));
OUTF(:, jout) = field(IZ,1);
th0 = 2.0D0*pi*(0:Ne-1)/Ne;
p0 = exp(1i*th0)';
p0v = [real(p0); imag(p0)];
reidx = 1:Ne;
imidx = Ne+1:2*Ne;

PforF   = complex(zeros(length(ZAxis_ipart), Ne)); % то что попадает на точки
PforF_p = complex(zeros(length(ZAxis_ipart), Ne)); % уравнения возбуждения

optsp = odeset('RelTol',1e-8,'AbsTol',1e-10);
Sp = cell(Ne,1);
p = oscill_ode(Sf, Nzi, ZAxisi, Delta, p0v, reidx, imidx, optsp);
% p = oscill_reim(field, Nzi, ZAxis, ZAxisi, Delta, p0v, reidx, imidx, opts);
% p = oscill_cmplx(field(1:Nzi), ZAxis(1:Nzi), Delta, p0);
for i=1:Ne
    Sp{i} = griddedInterpolant(ZAxisi, p(:,i),'spline');
    PforF(:,i) = Sp{i}(ZAxis_ipart);
    % plot(ZAxisi,abs(Sp{i}(ZAxisi)))
    % hold on
end
PforF_p = PforF;
J(1:ZAxis_ipart_end,1) = Ic * sum(PforF, 2)/Ne;
cu(:,1) = J(:) - 1i*kpar2(:).*field(:);
OUTJ(:,jout) = J(IZ,1);

% IDX = @(j) (j + 1);

fmax(IDX(0)) = max(abs(field(:,1)));
jmax(IDX(0)) = max(abs(cu(:,1)));
F0(IDX(0)) = field(1);
F1(IDX(0)) = field(2);
J0(IDX(0)) = cu(1);
J1(IDX(0)) = cu(2);
Sigma0(IDX(0)) = 0;
Sigma1(IDX(0)) = 0;
FNz(IDX(0)) = field(Nz);
FNzm1(IDX(0)) = field(Nzm1);
JNz(IDX(0)) = cu(Nz);
JNzm1(IDX(0)) = cu(Nzm1);
SigmaNz(IDX(0)) = 0;
SigmaNzm1(IDX(0)) = 0;
Eff(IDX(0)) = 1 - sum(abs(p(Nzi,:)).^2)/Ne;
Omega(IDX(0)) = 0;

WL(IDX(0)) = -dz * (-1i/6 * (2 * J0(IDX(0)) + J1(IDX(0))));
WR(IDX(0)) =  dz * (-1i/6 * (2 * JNz(IDX(0)) + JNzm1(IDX(0))));

SHOW = 1;
if SHOW == 1
    [lhfmax, lhfabs, lhjmax, lhjabs, hFig] = makeFig(ZAxis, ZAxis_ipart , TAxis);
end

%Coefficients
coeff_1i_m_C0_m_2_d_3_d_dt = 1i*C0*2.0D0/3.0D0/dt;
coeff_1i_m_C0_d_3_d_dt = 1i*C0/3.0D0/dt;
coeff_1i_d_6 = 1i/6.0D0;
coeff_4_d_3_m_SQRDT = 4.0D0/3.0D0*SQRTDT;
coeff_2_d_3_m_SQRDT = 2.0D0/3.0D0*SQRTDT;
coeff_1i_m_SQRDZ = 1i*SQRDZ;
coeff_1i_m_C0_m_SQRDZ_d_dt = 1i*C0*SQRDZ/dt;
coeff_C1_m_coeff_4_d_3_m_SQRDT = C1*coeff_4_d_3_m_SQRDT;
coeff_C2_m_coeff_4_d_3_m_SQRDT = C2*coeff_4_d_3_m_SQRDT;
coeff_exp_CL_m_dt = exp(CL*dt);
coeff_CL_m_dt = CL*dt;
coeff_exp_CR_m_dt = exp(CR*dt);
coeff_CR_m_dt = CR*dt;
coeff_dz_m_coeff_1i_d_6 = dz*coeff_1i_d_6;

num_st_test_iter = 0;
fmax_glob_old = max(abs(field(:,1)));

fprintf('\n');
timerVal = tic;
for step=1:Nt-1

    % step

    if SHOW == 1
        lhfmax.YData(1:step) = Omega(1:step);
        lhfmax.XData(1:step) = TAxis(1:step);
        lhfabs.YData = abs(field);

        lhjmax.YData(1:step) = Eff(1:step);
        lhjmax.XData(1:step) = TAxis(1:step);
        lhjabs.YData = abs(J);

        drawnow
    end    

    num_insteps = 0;
    % maxfield = max(abs(field(:,1)));
    testfield = field;
    % nst = 0;
    while 1        
        % nst = nst + 1
        Sf = griddedInterpolant(ZAxis,field_p,'spline');
        FforP_p = Sf(ZAxisi(2:Nzi));

        rhsP1 = -FforP' - 1i*p(1:Nzi-1,:).*(Delta - 1.0D0 + abs(p(1:Nzi-1,:)).^2);
        rhsP2 = -FforP_p' - 1i*p_p(2:Nzi,:).*(Delta - 1.0D0 + abs(p_p(2:Nzi,:)).^2);
        
        p_p(2:Nzi,:) = p(2:Nzi,:) + dzi/2*(rhsP1 + rhsP2); % предиктрор момента
        p_p(1,:) = p0'; % при z=0

        for i=1:Ne
            Sp{i} = griddedInterpolant(ZAxisi, p_p(:,i),'spline');
            PforF_p(:,i) = Sp{i}(ZAxis_ipart);
            % plot(ZAxisi(1:10),abs(Sp{i}(ZAxisi(1:10))))
            % hold on
        end

        J_p(1:ZAxis_ipart_end,1) = Ic * sum(PforF_p, 2)/Ne;
        cu_p(:,1) = J_p(:) - 1i*kpar2(:).*field_p(:);
   
        % plot(ZAxis, abs(J_p))
        % drawnow

        WL_PART = -dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt) * field(1)...
            + (-coeff_1i_m_C0_d_3_d_dt) * field(2)...
            - (2.0D0 * Sigma0(IDX(step-1)) + Sigma1(IDX(step-1))));
        
        WL(IDX(step)) = coeff_dz_m_coeff_1i_d_6*(2.0D0 * cu_p(1)  + 2.0D0 * cu(1)  + cu_p(2)  + cu(2)) + WL_PART;

        WR_PART =  dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt) * field(Nz)...
            + (-coeff_1i_m_C0_d_3_d_dt) * field(Nzm1)...
            - (2.0D0 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));

        WR(IDX(step)) = - coeff_dz_m_coeff_1i_d_6*(2.0D0 * cu_p(Nz) + 2.0D0 * cu(Nz) + cu_p(Nzm1) + cu(Nzm1)) + WR_PART;

        if step == 1
            IL = 0;
        elseif step == 2
            IL = coeff_4_d_3_m_SQRDT * (ul(0)*(1 - SQRT2D2) + ul(1)*(SQRT2M2minus2p5));
        else
            j = 1:step-2;
            IL = coeff_4_d_3_m_SQRDT * (ul(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
                + sum(ul(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
                + ul(step - 1)*(SQRT2M2minus2p5));
            % IL = coeff_4_d_3_m_SQRDT * (ul(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + ul(step - 1)*(SQRT2M2minus2p5));
            % for j = 1:step-2
            %     IL = IL + coeff_4_d_3_m_SQRDT * (ul(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
            % end
        end

        if step == 1
            IR = 0;
        elseif step == 2
            IR = coeff_4_d_3_m_SQRDT * (ur(0)*(1 - SQRT2D2) + ur(1)*(SQRT2M2minus2p5));
        else
            j = 1:step-2;
            IR = coeff_4_d_3_m_SQRDT * (ur(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
                + sum(ur(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
                + ur(step - 1)*(SQRT2M2minus2p5));
            % IR = coeff_4_d_3_m_SQRDT * (ur(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + ur(step - 1)*(SQRT2M2minus2p5));
            % for j = 1:step-2
            %     IR = IR + coeff_4_d_3_m_SQRDT * (ur(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
            % end
        end

        % D(1) = 0;
        %         D(1) = IN.TimeAmp * exp(1i * IN.TimeFreq * AxisTau(step));

        D_MIDDLE_PART = 2.0D0 * (1.0D0 - coeff_1i_m_C0_m_SQRDZ_d_dt) .* field(2:Nzm1)...
            - (field(1:Nz - 2) + field(3:Nz));

        D(2:Nzm1) = -coeff_1i_m_SQRDZ * (cu_p(2:Nzm1) + cu(2:Nzm1)) + D_MIDDLE_PART;

        D_0_PART =   C1 * (IL + coeff_2_d_3_m_SQRDT * (W0 * field(1)...
            + W1 * field(2) + WL(IDX(step-1))) * coeff_exp_CL_m_dt);

        D(1)  =   coeff_C1_m_coeff_4_d_3_m_SQRDT * WL(IDX(step)) + D_0_PART;

        D_END_PART = - C2 * (IR + coeff_2_d_3_m_SQRDT * (WNzm1 * field(Nzm1)...
            + WNz * field(Nz) + WR(IDX(step-1))) * coeff_exp_CR_m_dt);

        D(Nz) = - coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;

        % предиктро поля
        field_p = M \ D;
        % rfield_p = rtridag(C,A,B,D);
        % lfield_p = ltridag(C,A,B,D);
        % field_p = (rfield_p + lfield_p)/2.0D0;

        maxfield = max(abs(field_p(:,1)));
        if isnan(maxfield)
            print('Пиздец.')
            pause
        end
        maxdiff = max(abs(testfield - field_p));
        err = maxdiff/maxfield;
        if err < tol
            break
        end
        testfield = field_p;
        if num_insteps > 10
            error('\nToo many inner steps!\n');
            % fprintf('\nToo many inner steps!\n');
            % pause
        end
    end

    % num_insteps = 0;
    % maxfield = max(abs(field_p(:,1)));
    % testfield = field_p;
    % while 1
    %     num_insteps = num_insteps + 1;
    %     Sf = griddedInterpolant(ZAxis, field_p,'spline');
    %     p = oscill_ode(Sf, Nzi, ZAxisi, Delta, p0v, reidx, imidx, optsp);
    %     % p = oscill_reim(field, Nzi, ZAxis, ZAxisi, Delta, p0v, reidx, imidx, opts);
    %     % p = oscill_cmplx(field_p(1:Nzi), ZAxis(1:Nzi), Delta, p0);
    %     J_p(1:ZAxis_ipart_end,1) = Ic * sum(p, 2)/ Ne;
    %     cu_p(:,1) = J_p(:) - 1i*kpar2(:).*field_p(:);
    % 
    %     WL(IDX(step)) =  coeff_dz_m_coeff_1i_d_6 * (2.0D0 * cu_p(1)  + 2.0D0 * cu(1)  + cu_p(2)    + cu(2))    + WL_PART;
    %     WR(IDX(step)) = -coeff_dz_m_coeff_1i_d_6 * (2.0D0 * cu_p(Nz) + 2.0D0 * cu(Nz) + cu_p(Nzm1) + cu(Nzm1)) + WR_PART;
    % 
    %     D(2:Nzm1) = -coeff_1i_m_SQRDZ * (cu_p(2:Nzm1) + cu(2:Nzm1)) + D_MIDDLE_PART;
    % 
    %     D(1) =   coeff_C1_m_coeff_4_d_3_m_SQRDT * WL(IDX(step)) + D_0_PART;
    %     D(Nz) = -coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;
    % 
    % 
    %     % samosoglasovannoe pole
    %     field_p(:,1) = M \ D;
    %     % rfield_p(:,1) = rtridag(C,A,B,D);
    %     % lfield_p(:,1) = ltridag(C,A,B,D);
    %     % field_p = (rfield_p + lfield_p)/2.0D0;
    % 
    %     maxdiff = max(abs(testfield - field_p));
    %     err = maxdiff/maxfield;
    %     if err < tol
    %         break
    %     end
    %     testfield = field_p;
    %     if num_insteps > 1000
    %         fprintf('\nToo many inner steps!\n');
    %         pause;
    %     end
    % end    

    % rhsP1 = -FforP - 1i*p(1:Nzi-1,:).*(Delta - 1.0D0 + abs(p(1:Nzi-1,:)).^2);
    % rhsP2 = -Sf(ZAxisi(2:Nzi)) - 1i*p_p(2:Nzi,:).*(Delta - 1.0D0 + abs(p_p(2:Nzi,:)).^2);
    % 
    % p(2:Nzi,:) = p(2:Nzi,:) + dzi/2*(rhsP1 + rhsP2); % предиктрор момента
    % p(1,:) = p0'; % при z=0

    p = p_p;

    field(:,1) = field_p(:,1);
    Sf = griddedInterpolant(ZAxis, field,'spline');

    FforP = Sf(ZAxisi(1:Nzi-1)); % начальлное значение для следующего шага  

    for i=1:Ne
        Sp{i} = griddedInterpolant(ZAxisi, p(:,i),'spline');
        PforF(:,i) = Sp{i}(ZAxis_ipart);
        % plot(ZAxisi(1:10),abs(Sp{i}(ZAxisi(1:10))))
        % hold on
    end

    J(1:ZAxis_ipart_end,1) = Ic * sum(PforF, 2)/Ne;
    cu(:,1) = J(:) - 1i*kpar2(:).*field(:);
    fmax(IDX(step)) = max(abs(field(:,1)));
    jmax(IDX(step)) = max(abs(cu(:,1)));

    F0(IDX(step)) =  field(1);
    F1(IDX(step)) = field(2);
    FNz(IDX(step)) =  field(Nz);
    FNzm1(IDX(step)) = field(Nzm1);
    J0(IDX(step)) = cu(1);
    J1(IDX(step)) = cu(2);
    JNz(IDX(step)) = cu(Nz);
    JNzm1(IDX(step)) = cu(Nzm1);

    %     Omega(IDX(step)) = (angle(field(Nz)) - angle(FNz(IDX(step-1))))/dt;
    Omega(IDX(step)) = imag(log(field(Nz)/FNz(IDX(step-1))))/dt;
    Eff(IDX(step)) = 1 -  sum(abs(p(Nzi,:)).^2)/Ne;

    if (mod(num_st_test_iter,1000))
        fmax_glob_new = max(abs(field(:,1)));
        if abs(fmax_glob_new - fmax_glob_old)/fmax_glob_old < tol
            jout = jout + 1;
            OUTF(:, jout) = field(IZ,1);
            OUTJ(:, jout) = J(IZ,1);
            fprintf('Emergency exit!\n');
            return;
        end
        num_st_test_iter = num_st_test_iter + 1;
        fmax_glob_old = fmax_glob_new;
    end

    if mod(step,INTT) == 0
        jout = jout + 1;
        OUTF(:, jout) = field(IZ,1);
        OUTJ(:, jout) = J(IZ,1);
    end

    Sigma0(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(1) ...
        + (-coeff_1i_m_C0_d_3_d_dt) * F0(IDX(step - 1)) ...
        -coeff_1i_d_6*(cu(1) + J0(IDX(step - 1))) - Sigma0(IDX(step - 1));
    Sigma1(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(2) ...
        + (-coeff_1i_m_C0_d_3_d_dt) * F1(IDX(step - 1)) ...
        -coeff_1i_d_6*(cu(2) + J1(IDX(step - 1))) - Sigma1(IDX(step - 1));

    SigmaNz(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(Nz) ...
        + (-coeff_1i_m_C0_d_3_d_dt) * FNz(IDX(step - 1)) ...
        -coeff_1i_d_6*(cu(Nz) + JNz(IDX(step - 1))) - SigmaNz(IDX(step - 1));
    SigmaNzm1(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(Nzm1) ...
        + (-coeff_1i_m_C0_d_3_d_dt) * FNzm1(IDX(step - 1)) ...
        -coeff_1i_d_6*(cu(Nzm1) + JNzm1(IDX(step - 1))) - SigmaNzm1(IDX(step - 1));

    k = step + 1;

    ConLow(IDX(step)) = (2*imag(field(Nz)*conj(dfdzNz(step)) - field(1)*conj(dfdz0(step))) - Ic*Eff(IDX(step)))/(Ic*Eff(IDX(step)))*100;

    fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        'Step = %8i   Time = %10.4f   Bmax = %+15.10e   Jmax = %+15.10e   W = %+15.10e   E = %+15.10e   CL = %+10.5f %%'],...
        int64(step), TAxis(k), fmax(k), max(abs(cu(:,1))), Omega(IDX(step)), Eff(IDX(step)), ConLow(IDX(step)));
end

OUTJ(:,jout) = J(IZ,1);

fprintf("\n\n\n");

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

fprintf(" \n\n");

    function  f = ur(j)
        coder.inline("always");
        f = (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(coeff_CR_m_dt * (step - j));
    end

    function  f = ul(j)
        coder.inline("always");
        f = (W1 * F1(IDX(j)) + W0 * F0(IDX(j)) + WL(IDX(j))).' .* exp(coeff_CL_m_dt * (step - j));
    end

    function j = IDX(j)
        coder.inline("always");
        j = j + 1;
    end

    function  dfdz = dfdzNz(j)
        coder.inline("always");
        dfdz = (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j)));
    end

    function  dfdz = dfdz0(j)
        coder.inline("always");
        dfdz = (W1 * F1(IDX(j)) + W0 * F0(IDX(j)) + WL(IDX(j)));
    end
end




