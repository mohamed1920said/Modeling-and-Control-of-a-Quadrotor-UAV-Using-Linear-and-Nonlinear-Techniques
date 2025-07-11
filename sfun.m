function [sys,x0,str,ts,simStateCompliance] = sfun(t,x,u,flag,m,l,Ix,Iy,Iz,Jr,b,d,K1,K2,K3,K4,K5,K6,g)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u,m,l,Ix,Iy,Iz,Jr,b,d,K1,K2,K3,K4,K5,K6,g);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  =[0 0 0 0 0 0 0 0 0 0 0 0]';

%
% str is always an empty matrix
%
str = [];
%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%

function sys=mdlDerivatives(t,x,u,m,l,Ix,Iy,Iz,Jr,b,d,K1,K2,K3,K4,K5,K6,g)

ux=cos(x(1))*sin(x(3))*cos(x(5))+sin(x(1))*sin(x(5));
uy=cos(x(1))*sin(x(3))*sin(x(5))-sin(x(1))*sin(x(5));


a1=(Iy-Iz)/Ix;
a2=(Iz-Ix)/Iy;
a3=(Ix-Iy)/Iz;
%vecteur d'état 
phi=x(1);
dphi=x(2);
teta=x(3);
dteta=x(4);
psi=x(5);
dpsi=x(6);
xp=x(7);
dx=x(8);
y=x(9);
dy=x(10);
z=x(11);
dz=x(12);
%dynamique de drone
ddx=((cos(phi)*cos(psi)*sin(teta)+sin(phi)*sin(psi))*u(1)-K1*dx)/m;
ddy=((cos(phi)*sin(psi)*sin(teta)-sin(phi)*cos(psi))*u(1)-K2*dy)/m;
ddz=(cos(phi)*cos(teta)*u(1)+g-K3*dz)/m;
ddphi = a1 * dpsi * dteta - (Jr / Ix) * u(5) * dteta - (K4 / Ix) * dphi^2 + (1 / Ix) * u(2);
ddteta = a2 * dpsi * dphi - (Jr / Iy) * u(5) * dphi - (K5 / Iy) * dteta^2 + (1 / Iy) * u(3);
ddpsi = a3 * dphi * dteta - (K6 / Iz) * dpsi^2 + (l/ Iz) * u(4);
sys = [dx ddx dy ddy dz ddz dphi ddphi dteta ddteta dpsi ddpsi];


%%% 2eme METHODE 
%%% x=[phi dphi teta dteta psi dpsi x dx y dy z dz]';
% 
% a1=-K4/Ix;a2=(Iy-Iz)/Ix;a3=-Jr/Ix;a4=-K5/Ix;
% a5=(Iz-Ix)/Iy;a6=-Jr/Iy;a7=-K6/Iz;a8=(Ix-Iy)/Iz;
% a9=-K1/m;a10=-K2/m;a11=-K3/m;b1=l/Ix;b2=l/Iy;b3=l/Iz;

%%% representation d'état
% dx1=x(2);
% dx2=a1*x(2)^2+a2*x(4)*x(6)+a3*u(5)*x(4)+b1*u(2);
% dx3=x(4);
% dx4=a4*x(2)^2+a5*x(2)*x(6)+a6*u(5)*x(2)+b2*u(3);
% dx5=x(6);
% dx6=a7*x(6)^2+a8*x(2)*x(4)+a3*x(4)+b3*u(4);
% dx7=x(8);
% dx8=a9*x(8)+(ux*u(1))/m;
% dx9=x(10);
% dx10=a10*x(10)+(uy*u(1))/m;
% dx11=x(12);
% dx12=a11*x(12)+(cos(x(7))*cos(x(9))*u(1)-g)/m;
% sys =[dx1 dx2 dx3 dx4 dx5 dx6 dx7 dx8 dx9 dx10 dx11 dx12]';


% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,m,l,Ix,Iy,Iz,Jr,b,d,K1,K2,K3,K4,K5,K6,g)
phi=x(1);
dphi=x(2);
teta=x(3);
dteta=x(4);
psi=x(5);
dpsi=x(6);
xp=x(7);
dx=x(8);
y=x(9);
dy=x(10);
z=x(11);
dz=x(12);
sys = x;

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
