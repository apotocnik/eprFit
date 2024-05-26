%--------------------------------------------------------------------------
% fun_lib - fit function library
%
% Author: Anton Potocnik, F5, IJS
% Date:   28.01.2009 - 22.03.2010
% Arguments:
%       funstr = equ_lib(name, number)
% Input:    name    ... function name ('dLorentz', 'dDyson', 'dGauss') no
% numbers or other characters, only letters
%           number  ... parameters number A# w# xc# alpha#
% Output:   funstr  ... string fitting function
%--------------------------------------------------------------------------
% NOTE!
% names must consist solely on letters, no numbers or other signs like *
%
function funstr = fun_lib(name, number)
    
% Sign needed!!!
switch name
    case 'dLorentz'     % w ... w(fwhm)
        funstr = ['-16*a#/pi*w#*(x-xc#)/(4*(x-xc#)^2+w#^2)^2' ...
                  '-16*a#/pi*w#*(x+xc#)/(4*(x+xc#)^2+w#^2)^2'];
              
    case 'dDyson'       % w ... 2*w(pp) alpha = 0abs...1dis 
        funstr = ['-4/pi*a#*(4*(x-xc#)*(w#*cos(alpha#)-2*(x-xc#)*sin(alpha#))/(4*(x-xc#)*(x-xc#)+w#*w#)/(4*(x-xc#)*(x-xc#)+w#*w#)+ sin(alpha#)/(4*(x-xc#)*(x-xc#)+w#*w#))' ... 
                  '-4/pi*a#*(4*(x+xc#)*(w#*cos(alpha#)-2*(x+xc#)*sin(alpha#))/(4*(x+xc#)*(x+xc#)+w#*w#)/(4*(x+xc#)*(x+xc#)+w#*w#)+ sin(alpha#)/(4*(x+xc#)*(x+xc#)+w#*w#))'];
              
              % old case 'dDyson'       % w ... w(fwhm) alpha = 0abs...pi/2 dis 
    %    funstr = ['+a#/w#/w#/pi*(-2*(x-xc#)/w#*(1-alpha#)+(1-(x-xc#)^2/w#^2)*alpha#)/(1+(x-xc#)^2/w#^2)^2' ...
    %              '+a#/w#/w#/pi*(-2*(x+xc#)/w#*(1-alpha#)+(1-(x+xc#)^2/w#^2)*alpha#)/(1+(x+xc#)^2/w#^2)^2'];
              
    case 'dGauss'       % w ... 2*w(pp)
       funstr =  ['-a#/sqrt(2*pi)/w#^3*(x-xc#)*exp(-(x-xc#)^2/2/w#^2)' ...
                  '-a#/sqrt(2*pi)/w#^3*(x+xc#)*exp(-(x+xc#)^2/2/w#^2)'];
                           
    case 'Lorentz'      % w ... w(fwhm)
        funstr = ['+2*a#/pi*w#/(4*(x-xc#)^2+w#^2)' ...
                  '+2*a#/pi*w#/(4*(x+xc#)^2+w#^2)'];
              
    case 'Gauss'        % w ... 2*w(pp)
       funstr =  ['-a#/sqrt(2*pi)/w#*exp(-(x-xc#)^2/2/w#^2)' ...
                  '-a#/sqrt(2*pi)/w#*exp(-(x+xc#)^2/2/w#^2)'];
                  
    case 'Dyson'       % w ... 2*w(pp)  alpha = 0abs...1dis 
        funstr = ['+a#*2/pi*(w#*cos(alpha#)+2*(x-xc#)*sin(alpha#))/(4*(x-xc#)^2+w#^2)' ...
                  '+a#*2/pi*(w#*cos(alpha#)+2*(x+xc#)*sin(alpha#))/(4*(x+xc#)^2+w#^2)'];
    
    case 'pVoigt'       % w ... 2*w(pp)  alpha = 0lor...1gauss 
        funstr = ['+a#*(alpha#*2/pi*w#/(4*(x-xc#)^2+w#^2)+(1-alpha#)*sqrt(4*log(2)/pi)/w#*exp(-(x-xc#)^2*4*log(2)/w#^2))' ];
              
    % old case 'Dyson'       % w ... 2*w(pp)  alpha = 0abs...1dis 
    %    funstr = ['+a#/pi*(w#*(1-alpha#)+(x-xc#)*alpha#)/(w#^2+(x-xc#)^2)' ...
    %              '+a#/pi*(w#*(1-alpha#)+(x+xc#)*alpha#)/(w#^2+(x+xc#)^2)'];
    
%     case 'Quadratic'
%         funstr = ['+q#*x*x'];
    
    case 'CoupledLorentzians'
%         funstr = {'a*((W2-2*(x-x0)*J)*(W1*W1+W2*W2)-4*((x-x0)*W2-(G0+2*J)*W1)*((x-x0)*W1+(G0+J)*W2))/(W1*W1+W2*W2)/(W1*W1+W2*W2)'
%                   'W1=((x-xca)*(x-xcb)-(wa+J)*(wb+J)+J*J)'
%                   'W2=((x-xca)*(wb-J)+(x-xcb)*(wa+J))'
%                   'x0 = 0.5*(xca+xcb)'
%                   'G0 = 0.5*(wa+wb)'
%                   };
%         funstr = {'a*((W2-2*(x-0.5*(xca+xcb))*J)*(W1*W1+W2*W2)-4*((x-0.5*(xca+xcb))*W2-(0.5*(wa+wb)+2*J)*W1)*((x-0.5*(xca+xcb))*W1+(0.5*(wa+wb)+J)*W2))/(W1*W1+W2*W2)/(W1*W1+W2*W2)'
%                   'W1=((x-xca)*(x-xcb)-(wa+J)*(wb+J)+J*J)'
%                   'W2=((x-xca)*(wb-J)+(x-xcb)*(wa+J))'
%                   'x0='
%                   'G0='
%                   }; 
% First Article corrected!!!
        funstr = 'a#*((((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#))-2*(x-0.5*(xca#+xcb#))*J#)*(((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)*((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)+((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#))*((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#)))-4*((x-0.5*(xca#+xcb#))*((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#))-(0.5*(wa#+wb#)+2*J#)*((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#))*((x-0.5*(xca#+xcb#))*((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)+(0.5*(wa#+wb#)+J#)*((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#))))/(((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)*((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)+((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#))*((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#)))/(((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)*((x-xca#)*(x-xcb#)-(wa#+J#)*(wb#+J#)+J#*J#)+((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#))*((x-xca#)*(wb#-J#)+(x-xcb#)*(wa#+J#)))';
        
    case 'twoDLorDiffW'     % w ... w(fwhm)
        funstr = ['-16*a#/pi*wa#*(x-xc#)/(4*(x-xc#)^2+wa#^2)^2' ...
                  '-16*a#/pi*wb#*(x-xc#)/(4*(x-xc#)^2+wb#^2)^2'];
              
    case 'rubyRA'      % w ... w(fwhm)
        funstr = ['+2/pi*( Ia#*Wa#/(4*(x-xc#+3/2*dx#)^2+Wa#^2) + Ib#*Wb#/(4*(x-xc#+1/2*dx#)^2+Wb#^2) ' ...
                        '+ Ic#*Wc#/(4*(x-xc#-1/2*dx#)^2+Wc#^2) + Id#*Wd#/(4*(x-xc#-3/2*dx#)^2+Wd#^2) )'];          
    case 'rubyR'      % w ... w(fwhm)
        funstr = ['+2/pi*( Ib#*335.28*Wb#/(4*(x-xxb#+3/2*dx#)^2+Wb#^2) + Ib#*230.85*Wb#/(4*(x-xxb#+1/2*1.1*dx#)^2+Wb#^2) ' ...
                        '+ Ib#*245.86*Wb#/(4*(x-xxb#-1/2*dx#)^2+Wb#^2) + Ib#*456.75*Wb#/(4*(x-xxb#-3/2*dx#)^2+Wb#^2) ' ...  
                        '+ Ia#*1512.8*Wa#/(4*(x-xxa#+3/2*dx#)^2+Wa#^2) + Ia#*1324.6*Wa#/(4*(x-xxa#+1/2*1.1*dx#)^2+Wa#^2) ' ...   
                        '+ Ia#*1197.3*Wa#/(4*(x-xxa#-1/2*dx#)^2+Wa#^2) + Ia#*908.05*Wa#/(4*(x-xxa#-3/2*dx#)^2+Wa#^2) )'];
                    
    case 'qxx'
        funstr = '+q#*x*x';
    
    case 's*x'
        funstr = '+s*x';
    
    case 'y0'
        funstr = '0*x+y0';
        
    otherwise
        funstr = {'dLorentz','dDyson','dGauss','Lorentz','Gauss','Dyson','pVoigt','CoupledLorentzians','twoDLorDiffW','rubyRA','rubyR','qxx','s*x','y0'}; % !!! s*x and y0 must be last two!!!
              
end

% Replace # with number
if ~iscell(funstr) 
    if nargin == 2
        funstr = strrep(funstr,'#',num2str(number));
    else
        funstr = strrep(funstr,'#','');
    end
end
