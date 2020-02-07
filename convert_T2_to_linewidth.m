%% T2_HOMO - Convert a T2 time to a homogeneous linewidth
%
%   This function will calculate the linewidth of a Lorenztian with a given
%   dephasing time (T2). The linewidth is 1/(pi*T2), followed by a
%   conversion from ps-1 to cm-1.
%
%   Input can be a single number, though it will also work for a vector or
%   a matrix. If you need to access the results, T2_cm will give you your
%   results in whatever format you put in.

invPsToWavenumbers = 33.3567;
T2 = input('Please enter the T2 time (in ps): ');
T2_cm = invPsToWavenumbers*1./(pi*T2);
if size(T2) == [1 1];
    Prompt1 = sprintf('The homogeneous linewidth is %.2f cm-1. \n',T2_cm);
    disp(Prompt1);
else
    Prompt2 = 'The homogeneous linewidths (in cm-1) are:';
    disp(Prompt2);
    disp(T2_cm);
end