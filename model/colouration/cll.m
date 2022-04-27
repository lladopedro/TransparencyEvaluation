function [CLL,freqs, ILD, ITD, Xcorr ] = cll(signall,signalr,srate)

% version 0.1

% Ville Pulkki 1999, Helsinki University of Technology
%
% Calculates the main spatialization cues using a simple 
% auditory model. Arguments signall and signalr are the
% sound signals that enter to each ear canal, srate is the
% sample rate. ILD is in phons, ITD in seconds, Xcorr is the cross
% correlation functions calculated, CLL is the composite loudness level  
% in phons, freqs are the frequency points. 
%
% You have to simulate the signal appearing in left and right ear 
% canals. 
% The coherent signals should enter ears within a 1ms 
% time window. This model is a simplification of the human auditory system.
% Precedence effect and other time-variant mechanisms are not implemented.
% 
% Roughly, this function can be used to simulate situations in which the
% signal is presented via headphones or via equivalently distant
% loudspeakers in an anechoic chamber. 
%
% Requires Auditory toolbox 
% path(path,'/work/matlab5/toolbox/AuditoryToolbox');
% or equivalent.
%


% number of bands
N=42;
lowFreq=200;
number_of_samples = max(size(signall));

%cochlea bpfiltering 
[f,b] = MakeERBFiltersB(srate,N,lowFreq);

%calculating the center frequencies of gammatone filters

EarQ = 9.26449;               %  Glasberg and Moore Parameters
minBW = 24.7;
order = 1;
cf = -(EarQ*minBW) + exp((1:N)'*(-log(srate/2 + EarQ*minBW) + ...
    log(lowFreq + EarQ*minBW))/N) ...
    *(srate/2 + EarQ*minBW);
freqs_rev = cf;

Xltmp=FilterBank(f,b,signall)'; 
Xrtmp=FilterBank(f,b,signalr)';

%turning the order, low frequencies first
for i=1:1:N
  Xl(N+1-i,:)=Xltmp(i,:);
  Xr(N+1-i,:)=Xrtmp(i,:);
  freqs(N+1-i)=freqs_rev(i);
end


[B,A]=butter(1,800/(srate/2));
delays=[-50:1:50];
for g=N:-1:1
	xl=(Xl(g,:));
	xr=(Xr(g,:));
	%rectification
	xl=xl.*(sign(xl)+1)/2;
	xr=xr.*(sign(xr)+1)/2;
	%filtering with lowpass filter
	xl=filter(B,A,xl);	
	xr=filter(B,A,xr);	
	% energies
	xautocorr=sqrt(xl * xl' * xr * xr');
	%cross-correlation
	Xcorr_unsc(g,:)=xcorr(xl, xr, 50);
	%energy-scaling
	Xcorr(g,:)=Xcorr_unsc(g,:)/xautocorr;
	%loudnesses
	louddenl(g)=sqrt(sqrt(sum(xl .* xl)/number_of_samples));
	louddenr(g)=sqrt(sqrt(sum(xr .* xr)/number_of_samples));
% 	g
end

% Cues
[qq,ITD]=max(Xcorr_unsc');
ITD=(ITD-51)/srate;
CLL=log2(louddenl+louddenr)*10+40;
ILD=(log2(louddenl)-log2(louddenr))*10;
