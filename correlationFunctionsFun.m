function out = correlationFunctionsFun(p,t,options)
%return the correlation function corresponding to the paramters p and the
%input correlation function functional form in options.damping, which
%should correspond to the lineshape functions functions in
%analyticalResponseFunctionsFun

damping = options.damping;

switch damping,
  case {'overdamped','1exp'}
    %overdamped exp(-t/tau)

    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    tau1 = p(2); %first timescale (ps)

    %g = @(t) Delta^2/Lambda^2.*(exp(-Lambda.*t)-1+Lambda.*t);
    
    c2 = @(t) Delta1_cm^2.*exp(-t./tau1)
  case {'critical','1expcrit'}
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    tau1 = p(2); %first timescale (ps)
     
    %critically damped (1+2t/tau)exp(-2t/tau)
%    g = @(t) Delta^2/4/Lambda^2.*exp(-2.*Lambda.*t) ...
%      .*(3 + 2*Lambda*t + exp(2.*Lambda.*t).*(4*Lambda.*t - 3));
    c2 = @(t) Delta1_cm^2.*(1+2*t./tau1).*exp(-2*t./tau1);
    
  case '2exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(3); %first timescale (ps)
    tau2 = p(4); %second timescale (ps)

%     Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%     Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
%     Lambda1 = 1/tau1;
%     Lambda2 = 1/tau2;
%     anh = anh_cm*wavenumbersToInvPs*2*pi;
%     g = @(t) Delta1^2/4/Lambda1^2 ...
%       .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
%       + Delta2^2/4/Lambda2^2 ...
%       .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3));

    c2 = @(t) Delta1_cm^2.*exp(-t./tau1)+...
        Delta2_cm^2.*exp(-t./tau2);

    case '3exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(4); %first timescale (ps)
    tau2 = p(5); %second timescale (ps)
    tau3 = p(6); %second timescale (ps)

%     Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%     Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
%     Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
%     Lambda1 = 1/tau1;
%     Lambda2 = 1/tau2;
%     Lambda3 = 1/tau3;
%     anh = anh_cm*wavenumbersToInvPs*2*pi;
%     g = @(t) Delta1^2/4/Lambda1^2 ...
%       .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
%       + Delta2^2/4/Lambda2^2 ...
%       .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
%       + Delta3^2/4/Lambda3^2 ...
%       .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3));
    
    c2 = @(t) Delta1_cm^2.*exp(-t./tau1) + ...
        Delta2_cm^2.*exp(-t./tau2) + ...
        Delta3_cm^2.*exp(-t./tau3);

  case '3exp1off'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);%linewidth (sigma) in wavenumbers of static component
    tau1 = p(5); %first timescale (ps)
    tau2 = p(6); %second timescale (ps)
    tau3 = p(7); %second timescale (ps)

%     Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%     Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
%     Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
%     Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
%     Lambda1 = 1/tau1;
%     Lambda2 = 1/tau2;
%     Lambda3 = 1/tau3;
% 
%     g = @(t) Delta1^2/4/Lambda1^2 ...
%       .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
%       + Delta2^2/4/Lambda2^2 ...
%       .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
%       + Delta3^2/4/Lambda3^2 ...
%       .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3)) ...
%       + Delta4^2.*t.^2/2;


  c2 = @(t) Delta1_cm^2.*exp(-t./tau1) + ...
        Delta2_cm^2.*exp(-t./tau2) + ...
        Delta3_cm^2.*exp(-t./tau3) + ...
        Delta4_cm^2;
  case '1expcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    tau1 = p(3); %first timescale (ps)
      c2 = @(t) Delta1_cm^2.*(1+2*t./tau1).*exp(-2*t./tau1);

    case '2expcrit'
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
        tau1 = p(3); %first timescale (ps)
        tau2 = p(4); %second timescale (ps)
        c2 = @(t) Delta1_cm^2.*(1+2*t./tau1).*exp(-2*t./tau1) + ...
            Delta2_cm^2.*(1+2*t./tau2).*exp(-2*t./tau2);

    case '3expcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(4); %first timescale (ps)
    tau2 = p(5); %second timescale (ps)
    tau3 = p(6); %second timescale (ps)

%     Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%     Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
%     Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
%     Lambda1 = 1/tau1;
%     Lambda2 = 1/tau2;
%     Lambda3 = 1/tau3;
%     anh = anh_cm*wavenumbersToInvPs*2*pi;
%     g = @(t) Delta1^2/4/Lambda1^2 ...
%       .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
%       + Delta2^2/4/Lambda2^2 ...
%       .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
%       + Delta3^2/4/Lambda3^2 ...
%       .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3));
    
        c2 = @(t) Delta1_cm^2.*(1+2*t./tau1).*exp(-2*t./tau1) + ...
            Delta2_cm^2.*(1+2*t./tau2).*exp(-2*t./tau2) + ...
            Delta3_cm^2.*(1+2*t./tau3).*exp(-2*t./tau3);

  case '3exp1offcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);%linewidth (sigma) in wavenumbers of static component
    tau1 = p(5); %first timescale (ps)
    tau2 = p(6); %second timescale (ps)
    tau3 = p(7); %second timescale (ps)

%     Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%     Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
%     Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
%     Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
%     Lambda1 = 1/tau1;
%     Lambda2 = 1/tau2;
%     Lambda3 = 1/tau3;
% 
%     g = @(t) Delta1^2/4/Lambda1^2 ...
%       .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
%       + Delta2^2/4/Lambda2^2 ...
%       .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
%       + Delta3^2/4/Lambda3^2 ...
%       .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3)) ...
%       + Delta4^2.*t.^2/2;

     c2 = @(t) Delta1_cm^2.*(1+2*t./tau1).*exp(-2*t./tau1) + ...
        Delta2_cm^2.*(1+2*t./tau2).*exp(-2*t./tau2) + ...
        Delta3_cm^2.*(1+2*t./tau3).*exp(-2*t./tau3) + ...
        Delta4_cm^2;

  case 'hynesform'
      error('not implemented yet')
  %fit c2 to hynes type function, double integrate with Mathematica
  %to find g(t)
  a1 = 0.3232;
  k11 = 30.01; %ps-1
  k12 = 17.41; %ps-1
  a2 = 0.3378;
  k2 = 8.270; %ps-1
  a3 = 0.3455;
  k3 = 1.897; %ps-1
  
  %yuck
  g = @(t) Delta^2*exp(-t.*(k12+k2+k3))/(k2^2*k3^2*(k11^2+k12^2)^2) ...
      .*(a1*k2^2*k3^2.*exp((k2+k3).*t).*(cos(k11.*t).*(k12^2-k11^2)-2*k11*k12*sin(k11*t) ...
					 + exp(k12.*t).*(k11^2*(k12.*t+1)+k12^2*(k12.*t-1))) ... 
	 + (k11^2+k12^2)^2.*exp(k12.*t).*(a3*k2^2*exp(k2.*t).*(exp(k3.*t).*(k3.*t-1)+1) ...
					  +a2*k3^2*exp(k3.*t).*(exp(k2.*t).*(k2*t-1)+1)));
  otherwise
    error('damping value is unknown');
end

out = c2(t);