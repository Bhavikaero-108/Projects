 % Centrifugal Specific Speed:
                
% d31     = Density 
% T35     = Static Temperature 
% To35    = Stagnation Temperature 
% P35     = Static Pressure
% Po35    = Stagnation Pressure
% m3      = mass flow rate
% Ra      = Gas Constant
% gamaa   = Cp/Cv for air 
% cpa     = Specific Heat of air at constant pressure (Cp)
% cva     = Specific Heat of air at constant volume   (Cv)
% nhpc    = Efficiency of High Pressure Compressure
% r       = radius
% C       = velocity
% M       = Mach Number
% m       = mass flow rate
% A       = Area
% N       = Rotations Per Minute
% Nb      = Number of Blades
% alpha_r32 = metal angle (theoretical)
% alphar32  = flow angle (rel)
%etaitt  = Impeller Adiabatic Efficiency (total-to-total)
% ratio_r33_r32 = Radious ratio of r33 / r32
% p33rcy =  Static Pressure Recovery

% Subscript " h " = Hub
% Subscript " sh " = Shroud
% Subscript " o " = Stagnation
% Subscript " _ " = " * " Superscript

%                                   INPUT VARIABLES 

m31    = 4.082    ; % kg/s
m32    = 4.082    ; % kg/s
To31   = 416.4124 ; % K
Po31   = 201335   ; % Pa
To35   = 597.056  ; % K
Po35   = 608243.5 ; % Pa
T35    = 595.8643 ; % K
P35    = 604005  ; % Pa
Ra     = 287      ; % J/kgK
gamaa  = 1.4      ;
Nb     = 32       ; % Number of Blades
Nv     = 17       ; % Number of Blades
Whpc   = 741.074  ; % kW
nhpc   = 0.85     ; % 85 Percentage
N      = 29400    ; % Rotations Per Minute
alpha31    = 20   ; % Degrees
alpha_r32  =  -25       ; % Degrees
alpha_r32sh = alpha_r32 ; % Degrees
alphar31sh = alpha31 ; 
alpha_r31 = -25   ; % alpha STAR r 32    (range is -40 to 0)
B_2 = 0.19        ;   % B star 
B_3 =  0.19           ;
% V32 =   V31sh * 0.5; % This is assumed, Can be in range of  (V31sh*0.5) < V32 <   (V31sh*0.6);
cpa    = gamaa * Ra / (gamaa - 1) ;
etaitt  =    0.89 ;   % Adiabatic Efficiency (total-to-total) Eta i t-t
ratio_r33_r32 = 1.075 ; % Radious ratio of r33 / r32
                % Inlet Impeller Hub Radius r31h:

r51h = 0.09;
r31h = ( 0.722 * r51h);     % [1]
do31 = Po31 / (Ra * To31);  % [2]

% Z is a temporary variable
% d31new is a temporary variable


            % Inlet Impeller Shroud Radius r31sh
            
% Initialize variables
% minMr31sh = Inf;
% minMr31shIndex = 0;
                         % Initialize loop index
                         
Mr31sh_values = [];
Z_values = [];
r31sh_values = [];
A31_values = [];
d31_values = [];
Ca31_values = [];
C31_values = [];
M31_values = [];
a1_values = [];
d31_values = [];
U31sh_values = [];
Cw31_values = [];
Vw31sh_values = [];
V31sh_values = [];
alphar31sh_values = [];
Mr31sh_values = [];
r31m_values = [];
U31m_values = [];
Z_values = [];


            

for r31sh = 0.09301: 0.0001 : 0.190010 ;
    
    A31 = pi * ( ((r31sh)^2) - ((r31h)^2) );
    
    d31 = do31;
    Ca31 = m31 / ( d31 * A31 ) ; % [3] axial velocity
    C31 = Ca31 / cosd(alpha31) ; % [4] velocity
    M31 = sqrt( ((C31)^2) / ( (To31 * gamaa * Ra) - (((gamaa - 1)/2)*((C31)^2)) ) ); % [5] Mach Number
    a1 = C31 / M31;
    d31 = do31 / ( ( 1 + (( (gamaa - 1) / 2 )*((M31)^(2))) )^(1/(gamaa - 1))) ;
    
        % Define the tolerance value
        tolerance = 0.00001;

        % Initialize d31 with the given value
        d31 = do31;

        % Initialize the loop variable
        diff_d31 = tolerance + 1;

        % Execute the loop body at least once
        while diff_d31 > tolerance
            % Store the previous value of d31
            d31_prev = d31;

            Ca31 = m31 / ( d31 * A31 ) ; % [3] axial velocity

            C31 = Ca31 / cosd(alpha31) ; % [4] velocity

            M31 = sqrt( ((C31)^2) / ( (To31 * gamaa * Ra) - (((gamaa - 1)/2)*((C31)^2)) ) ); % [5] Mach Number

            a1 = C31 / M31; % [6]
            
            d31 = do31 / ( ( 1 + (( (gamaa - 1) / 2 )*((M31)^(2))) )^(1/(gamaa - 1))) ; % Calculate the new value of d31

            diff_d31 = abs(d31 - d31_prev); % Calculate the difference between the previous and current values of d31

        
        end
        
        U31sh = r31sh * (( N * pi) / 30) ; % [8]

        Cw31 = Ca31 * tand( alpha31 );

        Vw31sh = U31sh - Cw31 ;

        V31sh = sqrt((Ca31 ^ 2)+(Vw31sh ^ 2)); % [9]

        alphar31sh = atand(Vw31sh / Ca31);
        
        Mr31sh = V31sh / a1 ; % [10]

        r31m = (r31h + r31sh) / 2 ;

        U31m = r31m * (( N * pi) / 30) ;
        
        Z = Ca31 / U31m;
        
         % Stores the values of the variables in a array so that the
         % smallest of them can be called out later

            r31sh_values(end+1) = r31sh;
            A31_values(end+1)  = A31;
            Ca31_values(end+1)  = Ca31;
            C31_values(end+1)   = C31;
            M31_values(end+1)   = M31;
            a1_values(end+1)    = a1;
            d31_values(end+1)   = d31;
            U31sh_values(end+1) = U31sh;
            Cw31_values(end+1)  = Cw31;
            Vw31sh_values(end+1) = Vw31sh;
            V31sh_values(end+1) = V31sh;
            alphar31sh_values(end+1) = alphar31sh;
            Mr31sh_values(end+1) = Mr31sh;
            r31m_values(end+1)  = r31m;
            U31m_values(end+1)   = U31m;
            Z_values(end+1)     = Z;
            
            
end

% Find minimum value of Mr31sh
[minMr31sh, minMr31shIndex] = min(Mr31sh);

% Find minimum value of Mr31sh_values and its corresponding r31sh_value
[min_Mr31sh, min_Mr31sh_index] = min(Mr31sh_values);
min_r31sh = r31sh_values(min_Mr31sh_index);

% Plot graph of values
plot(Z_values, Mr31sh_values);
xlabel('Z');
ylabel('Mr_{3-1sh}');
title('Values of Mr_{3-1sh} and Z');
hold on;

% Plot second graph
plot(Z_values, r31sh_values);
xlabel('C_a_3_1 / U_3_1_m');
ylabel('Mr_3_1_s_h & r_3_1_s_h values');
title('Mr_3_1_s_h & r_3_1_s_h vs (C_a_3_1 / U_3_1_m)');

% Label minimum value on graph
text(Z_values(min_Mr31sh_index), min_Mr31sh, sprintf('(%.4f, %.4f)', Z_values(min_Mr31sh_index), min_Mr31sh) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Label minimum value on graph
text(Z_values(min_Mr31sh_index), min_r31sh, sprintf('(%.4f, %.4f)', Z_values(min_Mr31sh_index), min_r31sh), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Use index of minimum Mr3l value to find corresponding values of other variables

	r31sh = r31sh_values(min_Mr31sh_index);
	A31   =   A31_values(min_Mr31sh_index);
	d31   =   d31_values(min_Mr31sh_index);  
	Ca31  = Ca31_values(min_Mr31sh_index);
	C31   = C31_values(min_Mr31sh_index);
	M31   = M31_values(min_Mr31sh_index); 
	a1    = a1_values(min_Mr31sh_index);
	d31   = d31_values(min_Mr31sh_index);
	U31sh      = U31sh_values(min_Mr31sh_index);
	Cw31       = Cw31_values(min_Mr31sh_index);
	Vw31sh     = Vw31sh_values(min_Mr31sh_index);
	V31sh      = V31sh_values(min_Mr31sh_index);
	alphar31sh = alphar31sh_values(min_Mr31sh_index);
	Mr31sh     = Mr31sh_values(min_Mr31sh_index);
	r31m       =  r31m_values(min_Mr31sh_index);
	U31m       =  U31m_values(min_Mr31sh_index); 
	Z          = Z_values(min_Mr31sh_index); 

% Plot graph of values
plot(Z_values, Mr31sh_values);
xlabel('Z');
ylabel('Mr_{3-1sh}');
title('Values of Mr_{3-1sh} and Z');

% Plot graph of values
plot(Z_values, r31sh_values);
xlabel('Z');
ylabel('r31sh_values');
title('Values of r31sh_values and Z');


plot(Mr31sh, Z); % `plot(Z, Mr31sh);` plots a graph of `Mr31sh` vs `Z`.
        
        [minMr, minMrIndex] = min(Mr31sh); %`[minMr, minMrIndex] = min(Mr31sh);` finds the minimum value of `Mr31sh` and its index.
        minZ = Z(minMrIndex); %`minZ = Z(minMrIndex);` finds the value of `Z` corresponding to the minimum value of `Mr31sh`.
        fprintf('Minimum value of Mr31sh is %f at Z=%f\n', minMr, minZ); %`fprintf('Minimum value of Mr is %f at Z=%f\n', minMr, minZ);` displays the minimum value of `Mr31sh` along with its corresponding value of `Z`.
        
T31 = To31 / ( 1 + (( (gamaa - 1) / 2 )*((M31)^(2))) ) ;

P31 = Po31 * (( T31 / To31 )^((gamaa)/(gamaa - 1))) ;

d31 = P31 / (Ra * T31) ; 
        
Q31 = (m31 / d31) * 35.3147 ; % [ft^3/s]  
        
dhois = nhpc * cpa * ( To35 - To31 ) * 0.429923 ; % [Btu/lb]
        
Ns = (N * sqrt( Q31 )) / ((778.26 * dhois )^(3/4)); 

Vw31m = U31m - Cw31;

V32 =   V31sh * 0.5; % This is assumed, Can be in range of  (V31sh*0.5) < V32 <   (V31sh*0.6);

V31m = sqrt((Ca31 ^2)+(Vw31m ^2));

alphar31m = atand(Vw31m / Ca31);

Mr31m = V31m / sqrt( gamaa * Ra * T31) ;

U31h = (r31h * N * pi) / 30 ;

Vw31h = U31h - (Ca31 * tand(alpha31)) ;

V31h = sqrt((Ca31 ^ 2) + (Vw31h ^ 2));

alphar31h = atand(Vw31h / Ca31);

Mr31h = V31h / sqrt( gamaa * Ra * T31);

To32 = To35;



                                   % Inlet Impeller Incidence Angle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERIFY THIS PART %%%%%%%%%%%%%%%%%%
% 
% alpha_r31h = alphar31h - i31h;
% 
% alpha_r31m = alphar31m - i31m;
% 
% alpha_r31sh = alphar31sh - i31sh;


if ( alpha_r32 == (-30) ) || ( alpha_r32 > (-30) ) || ( alpha_r32 < 0) || ( alpha_r32 == 0) ;  
    
    U2corr1 = ( ( (alpha_r32sh) * ( 1348.8295 - 1226.8381 )) / (-30)) + 1226.8381 ;
    
elseif ( alpha_r32 == (-45) ) || ( alpha_r32 > (-45) ) || ( alpha_r32 < (-30)) ;
    
    U2corr1 = ( ( (alpha_r32sh - (-45) ) * ( 1348.8295 - 1422.189 ) ) / (-30 - (-45) )) + 1422.189 ; 

else
    
    disp('the entered alpha_r32 is out of range');
    
end
    
if ( alpha31 == (50) ) || ( alpha31 < (50) ) || ( alpha31 > 0) || ( alpha31 == 0);
    
    U2corr2 = ( ((alpha31) * ( 1241.1689 - 1190.26974 )) / (50)) + 1190.26975 ;
    
else
    
    disp('the entered alpha31 is out of range');
    
end

% r32 = (30 * U32 )/( N * pi );
% PRts = 3; 

U32corr = U2corr1 *(U2corr2 / 1190.269751);

                    % Impeller Blade Count Nb:
                    
% Nb = impeller blade count (Defined above)

                       % Impeller Exit Geometry: 

% V32 TO BE ASSUMES
                       
theta1 = (1.8 * To31)/ 518.7;
                       
U32 = (U32corr * sqrt(theta1) ) * 0.3048; % [1]

Cw32 = ( (cpa * (To32 - To31)) + (U31m * Cw31) )/(U32); 

sigma = ( 1 - ( (063 * pi) / Nb) );

Cw32in = Cw32 / sigma ; % [3]

Vw32in = - ( U32 - Cw32in );

V32in = Vw32in / sind( alpha_r32 ) ; % [4]

% V32 assumed and entered above

Vw32 = - ( U32 - Cw32 ) ;

alphar32 = asind( Vw32 / V32 );

Cr32 = Vw32 / ( tand(alphar32)); % [6]

if ( (Cr32 / Ca31) == (0.85) ) || (Cr32 / Ca31) == 0.9 || (Cr32 / Ca31) < 0.9 || (Cr32 / Ca31) > 0.85;
        disp('choose another value of alpha_r32');
elseif  alphar32 > 40 ;
        disp('choose another value of alpha_r32');
else
end

C32 = sqrt( (Cr32 ^ 2) + (Cw32 ^ 2) );

alphar32 = atand(Vw32 / Cr32) ;

M32 = sqrt( (C32 ^2) / ( (To32 * gamaa * Ra) - (((gamaa - 1)/2)*(C32 ^2)) ) ) ; 

Po32 = Po31 * (( 1 + ( (etaitt * (To32 - To31)) / To31))^(gamaa / (gamaa - 1)));

T32 = (C32 ^ 2) / ( gamaa * Ra * M32 * M32);

P32 = Po32 * (( T32 / To32)^(gamaa/(gamaa - 1))) ;

d32 = P32 / (Ra * T32);

Mr32 = V32 / (gamaa * Ra * T32) ;

r32 = (30 * U32)/(N * pi);

delta = alphar32 - alpha_r32 ;

                % Impeller Axial Length L:

% Input L12
L12 = ((r32 - r31h)/ 1.05) ;   %%%%%%%%%%%%%%%%% INPUT VARIABLE
                
% disp(['(r32 - r31h) / L12 = ' num2str(((r32 - r31h) / L12)),]);  
                
                % Channel Height at Impeller Exit b32:

 epsilonax = 0.000254; % m Tip Clearence
 
 % Cd2 = 1 - B_2 ;
 
 b32 = m32 / ((2* pi * r32 * d32 * Cr32)*(1 - B_2) ) ;
 
 h32 = b32 - epsilonax;
 
                 % Channel Height at Vaneless Diffuser b323:
                 
  b323 = b32;     %%%%%%%%%%%%%%% NEEDS TO BE VERIFIED and CHANGED
        
  b33  = b32;     %%%%%%%%%%%%%%% NEEDS TO BE VERIFIED and CHANGED
        
  disp(['b323 / r32 = ' num2str((b323 / r32)),]);
  
                  % Radius Ratio at Vaneless Diffuser 
 
% ratio_r33_r32 = (r33 / r32);
% 
% disp(['r33 / r32 = ' num2str((r33 / r32)),]);  

% Vaneless diffuser loss 1%
Po33 = Po32 * (1-0.01);
        
        
             % Flow Properties at Vaned Diffuser Leading Edge: 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 C33 = C32;
 
 r33 = ( ratio_r33_r32 * r32);
 
 Cw33 = Cw32 * (r32 / r33);
 
 Cr33t = sqrt((C33 ^2) - (Cw33 ^2)) ;
 
To33 = To32;

C33_previous = C33;
tolerance = 1e-6; % Set a small tolerance value

while true
    Cr33t = sqrt((C33_previous^2) - (Cw33^2));
    M33 = sqrt( (C33_previous^2) / ((To33 * gamaa * Ra) - (((gamaa - 1) / 2)*(C33_previous^2))) ) ; 
    T33 = (C33_previous^2) / (gamaa * Ra * M33 * M33);
    P33 = Po33 * ( ( T33 / To33) ^ (gamaa / (gamaa - 1)) );
    d33 = P33 /(Ra * T33);
    A33 = 2* pi * r33 * b33 ;
    m33 = d33 * Cr33t * A33 ;
    alpha33 = atand(Cr33t / Cw33 ); % Flow angle
    C33t  = sqrt((Cw33^2) + (Cr33t^2));
    
    % Check if the calculated value of C33t is equal to C33
    if abs(C33t - C33) < tolerance
        disp(['C33t = C33 = ' num2str(C33t),]);
        break;
    end
    
    % Update C33_previous with the calculated value of C33t
    C33_previous = C33t;
end

% C33t = sqrt( (Cw33 ^2) + (Cr33t ^2) );


                                 % Flow Properties at Vaned Diffuser Throat:
% Input B_3 i.e. B star 3


                
               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Needs to be Verified

To33cho = To33;
Po33cho = Po33;
M33cho = 1;

T33cho = To33 / ( 1 + (((gamaa - 1)/2) * M33cho * M33cho)) ;

P33cho = Po33 * ((T33cho / To33)^(gamaa / (gamaa - 1)));

p33rcy = (P33cho - P33) / (Po33 - P33) ; % Static Pressure Recovery 

disp(['p33rcy = ' num2str(p33rcy),]);

C33cho = sqrt(gamaa * Ra * T33cho);

m33cho = m33 * (1 + 0.01) ; 
        

    
    A33cho = (m33cho * Ra * T33cho )/ ((1 - B_3)* P33cho * C33cho) ;
    
    w_33 = A33cho / (Nv * b33) ;
    
%     if (( B_3 / w_33) == 0.8) || ((B_3 / w_33) > 0.8 ) || ((B_3 / w_33) == 1.2) || ((B_3 / w_33) < 1.2);
%         
%         break;     
%     end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             
%                                                 % Diffuser Vane 
%                                                 
% phi = 360 / Nv;
% 
% a2cl = - tand(90 - alpha_33 + phi );
% 
% b2ss = (r33 * cosd(phi)) + (r33 * sind(phi) * tand(90 - alpha_33 + phi - (omega / 2))) ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Something is missing here    %%%%%%%%%%%%%%%%%%%%%
% 
% omega = 2 * ( 90 - alpha_33 + phi - atand((( w_33 * sqrt( 1 + ( (tand( 90 - alpha_33 + phi)) ^ 2))) - (r33 * cosd( phi ))) / (r33 * sind(phi)))) ;
% 
%     
% 
%                     % Static Pressure Recovery Coefficient from Diffuser Throat to Exit:
% 
% M35 = 0.1;
% 
% T35 = To35 / (1 + (((gamaa - 1) / 2) * M35 * M35));
% 
% P35 = Po35 * ( (T35 / To35) ^ (gamaa / (gamaa - 1)) );
% 
% p35rcy = (P35 - P33cho) / (Po33 - P33cho) ;