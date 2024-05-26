%--------------------------------------------------------------------------
% fun_lib_sim - simulation functions library
%
% Author: Anton Potocnik, F5, IJS
% Date:   28.01.2009 - 09.12.2011
% Arguments:
%       funstr = equ_lib(name, number)
% Input:    name    ... function name ('dLorentz', 'dDyson', 'dGauss') no
% numbers or other characters, only letters
%           number  ... parameters number A# w# xc# alpha#
% Output:   funstr  ... string fitting function
%           coefs   .... {'a','dH',...}
%--------------------------------------------------------------------------
% NOTE!
% function names must consist solely of letters, no numbers or other signs 
% like * are allowed
%

function [funstr coefs startVal]= fun_lib_sim(name, number)
    
% Sign needed!!!
switch name
%     case 'dLorentz'     % w ... w(fwhm)
%         funstr = ['-16*a#/pi*w#.*(x-xc#)./(4*(x-xc#).^2+w#^2).^2' ...
%                   '-16*a#/pi*w#.*(x+xc#)./(4*(x+xc#).^2+w#^2).^2'];
%         coefs = {'a','w','xc'};
%         startVal = [1e-3 1 334];
              
    case 'dLorentz'       % w ... 2*w(pp) phase = 0:abs...pi/2:dis 
        funstr = ['-4/pi*a#.*(4*(x-xc#).*(w#*cos(phase#)-2*(x-xc#)*sin(phase#))./(4*(x-xc#).*(x-xc#)+w#*w#)./(4*(x-xc#).*(x-xc#)+w#*w#)+ sin(phase#)./(4*(x-xc#).*(x-xc#)+w#*w#))' ... 
                  '-4/pi*a#.*(4*(x+xc#).*(w#*cos(phase#)-2*(x+xc#)*sin(phase#))./(4*(x+xc#).*(x+xc#)+w#*w#)./(4*(x+xc#).*(x+xc#)+w#*w#)+ sin(phase#)./(4*(x+xc#).*(x+xc#)+w#*w#))'];
        coefs = {'a','w','xc','phase'}; 
        startVal = [1e-3 1 334 0];
        
    case 'dGaussDisp'       % w ... 2*w(pp) no other component. Only works on 64bit processors !!!
       funstr =  ['-a#/sqrt(2*pi)./w#^3.*(x-xc#).*exp(-(x-xc#).^2/2/w#^2).*(cos(phase#)-sin(phase#).*erfi((x-xc#)/sqrt(2)/w#))-0.1*sin(phase#)/pi/w#^2'];% ...
%                   '-a#/sqrt(2*pi)./w#^3.*(x+xc#).*exp(-(x+xc#).^2/2/w#^2).*(cos(phase#)-sin(phase#).*erfi((x+xc#)/sqrt(2)/w#))-sin(phase#)/pi/w#^2'];
       coefs = {'a','w','xc','phase'}; 
       startVal = [1e-3 1 334 0];
    
    case 'dGauss'       % w ... 2*w(pp) no other component
       funstr =  ['-a#/sqrt(2*pi)./w#^3.*(x-xc#).*exp(-(x-xc#).^2/2/w#^2)' ...
                  '-a#/sqrt(2*pi)./w#^3.*(x+xc#).*exp(-(x+xc#).^2/2/w#^2)'];
       coefs = {'a','w','xc'}; 
       startVal = [1e-3 1 334];
       
    case 'dLorentzPowd'       % Nstep and freq(GHz) must be defined, dH ... ?2*w(pp) phase = 0abs...1dis 
        funstr = ['+1e-6*A#*powdersimu(x/1000,Nstep,freq,gx#,gy#,gz#,dHx#/1000,dHy#/1000,dHz#/1000,phase#,1)'];
        coefs = {'A','gx','gy','gz','dHx','dHy','dHz','phase'};
        startVal = [1 2.0023 2.0023 2.0023 1 1 1 0];
        
    case 'dGaussPowd'       % Nstep and freq(GHz) must be defined, dH ... ?2*w(pp) phase = 0abs...1dis 
        funstr = ['+1e-6*A#*powdersimu(x/1000,Nstep,freq,gx#,gy#,gz#,dHx#/1000,dHy#/1000,dHz#/1000,phase#,2)'];
        coefs = {'A','gx','gy','gz','dHx','dHy','dHz','phase'};
        startVal = [1 2.0023 2.0023 2.0023 1 1 1 0];
        
    case 'LorentzPowd'       % Nstep and freq(GHz) must be defined, dH ... ?2*w(pp) phase = 0abs...1dis 
        funstr = ['+1e-6*A#*powdersimu(x/1000,Nstep,freq,gx#,gy#,gz#,dHx#/1000,dHy#/1000,dHz#/1000,phase#,3)'];
        coefs = {'A','gx','gy','gz','dHx','dHy','dHz','phase'};
        startVal = [1 2.0023 2.0023 2.0023 1 1 1 0];
        
    case 'GaussPowd'       % Nstep and freq(GHz) must be defined, dH ... ?2*w(pp) phase = 0abs...1dis 
        funstr = ['+1e-6*A#*powdersimu(x/1000,Nstep,freq,gx#,gy#,gz#,dHx#/1000,dHy#/1000,dHz#/1000,phase#,4)'];
        coefs = {'A','gx','gy','gz','dHx','dHy','dHz','phase'};
        startVal = [1 2.0023 2.0023 2.0023 1 1 1 0];
                           
%     case 'Lorentz'      % w ... w(fwhm)
%         funstr = ['+2*a#/pi*w#./(4*(x-xc#).^2+w#^2)' ...
%                   '+2*a#/pi*w#./(4*(x+xc#).^2+w#^2)'];
%         coefs = {'a','w','xc'};
%         startVal = [1e-3 1 334];
              
%     case 'Gauss'        % w ... 2*w(pp)
%        funstr =  ['+a#/sqrt(2*pi)/w#.*exp(-(x-xc#).^2/2/w#^2)' ...
%                   '+a#/sqrt(2*pi)/w#.*exp(-(x+xc#).^2/2/w#^2)'];
%        coefs = {'a','w','xc'};
%        startVal = [1e-3 1 334];
    
    case 'GaussDisp'        % w ... 2*w(pp). Only works on 64bit processors !!!
       funstr =  ['+a#/sqrt(2*pi)/w#.*exp(-(x-xc#).^2/2/w#^2).*(cos(phase#)+sin(phase#).*erfi((x-xc#)/sqrt(2)/w#))'];% ...
       %           '+a#/sqrt(2*pi)/w#.*exp(-(x+xc#).^2/2/w#^2)'];
       coefs = {'a','w','xc','phase'};
       startVal = [1e-3 1 334 0];
       
    case 'Gauss'        % w ... 2*w(pp)
       funstr =  ['+a#/sqrt(2*pi)/w#.*exp(-(x-xc#).^2/2/w#^2)' ...
                  '+a#/sqrt(2*pi)/w#.*exp(-(x+xc#).^2/2/w#^2)'];
       coefs = {'a','w','xc'};
       startVal = [1e-3 1 334];
                  
    case 'Lorentz'       % w ... 2*w(pp)  alpha = 0abs...1dis 
        funstr = ['+a#*2/pi*(w#*cos(phase#)+2*(x-xc#)*sin(phase#))./(4*(x-xc#).^2+w#^2)' ...
                  '+a#*2/pi*(w#*cos(phase#)+2*(x+xc#)*sin(phase#))./(4*(x+xc#).^2+w#^2)'];
        coefs = {'a','w','xc','phase'};
        startVal = [1e-3 1 334 0];      
                    
    case 'qxx'
        funstr = '+q#*x.*x';
        coefs = {'q'};
        startVal = [0];
    
    case 'sx'
        funstr = '+s#*x';
        coefs = {'s'};
        startVal = [0];
    
    case 'offset'
        funstr = '+y#';
        coefs = {'y'};
        startVal = [0];
        
    otherwise
        funstr = {'dLorentz','dGauss','dLorentzPowd','dGaussPowd','LorentzPowd','GaussPowd','dGaussDisp','Lorentz','Gauss','GaussDisp','qxx','sx','offset'}; 
              
end

% Replace # with number
if ~iscell(funstr) 
    ind = strfind(funstr,'#');
    if nargin == 2
        funstr(ind) = num2str(number);
    else
        funstr(ind) = '';
    end
end
