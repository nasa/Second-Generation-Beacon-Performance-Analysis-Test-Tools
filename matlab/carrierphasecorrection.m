% ///  Copyright(c) 2017 by Reese Bovard, all rights reserved
% ///  Proprietary Information of Concentric Real Time, LLC
% ///  Subject to Non-Disclosure.  Do Not Copy or Distribute
% ///  
% ///     Project: T21 Software Tools for Cospas-Sarsat DSSS Beacons
% ///     
% ///     Module: carrierphasecorrection 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 12 $ $Date: 2020-01-21 08:10:30 -0500 (Tue, 21 Jan 2020) $
% ///				$Id: carrierphasecorrection.m 12 2020-01-21 13:10:30Z reesebo $
% ///            

function [y,lamsav,t,freqsav,ensav,initialPhaseRad]=carrierphasecorrection(x,fs,bn,fchip,fig)
%CARRIERPHASECORRECTION Second-order QPSK Costas carrier-tracking loop.
%
%   [Y,LAMSAV,T,FREQSAV,ENSAV,INITIALPHASERAD] = ...
%       CARRIERPHASECORRECTION(...)
%
%   LAMSAV  Estimated input carrier phase, in radians.  The phase actually
%           applied to X is -LAMSAV.
%   FREQSAV Integral-path carrier-frequency estimate, in Hz.  This is the
%           appropriate diagnostic frequency trace; differentiating
%           LAMSAV sample by sample also differentiates the proportional
%           detector-error term and greatly amplifies its jitter.
%   ENSAV   Amplitude-normalized phase-detector error.
%   INITIALPHASERAD  Feed-forward initial carrier-phase estimate.  The
%           fourth-power estimate is inherently ambiguous by pi/2, which
%           is also the normal QPSK Costas-loop phase ambiguity and does
%           not affect carrier-frequency measurement.
%
%   The first three outputs retain their previous meanings, so existing
%   callers requesting [Y,LAMSAV,T] remain compatible.

if ~isvector(x) || isempty(x) || any(~isfinite(x))
    error('x must be a nonempty finite complex vector.');
end
if ~isscalar(fs) || ~isfinite(fs) || fs <= 0
    error('fs must be a positive finite scalar in Hz.');
end
if ~isscalar(fchip) || ~isfinite(fchip) || fchip <= 0 || fchip > fs
    error('fchip must be a positive finite scalar no greater than fs.');
end
if ~isscalar(bn) || ~isfinite(bn) || bn <= 0 || bn >= 0.5
    error('bn must be the normalized loop bandwidth relative to fchip (0 < bn < 0.5).');
end
if ~isscalar(fig) || ~isfinite(fig)
    error('fig must be a finite scalar.');
end

sps=round(fs/fchip);
if abs(fs/fchip-sps) > 1e-6
    warning('fs/fchip is not an integer; loop bandwidth uses sps=round(fs/fchip).');
end

% Second-order Costas loop with damping factor 1/sqrt(2).  bn is normalized to the
% chip rate, so convert it to normalized bandwidth at the sample rate.
damping=1/sqrt(2);
theta=(bn/sps)/(damping+1/(4*damping));
d=1+2*damping*theta+theta^2;
gi=4*theta^2/d;
gp=4*damping*theta/d;

y=zeros(size(x),'like',x);
ensav=zeros(size(x));
lamsav=zeros(size(x));               % estimated carrier phase in radians
freqsav=zeros(size(x));               % integral-path estimate in Hz
t=reshape((0:numel(x)-1)/fs,size(x));

% Estimate the phase before entering the feedback loop.  Normalizing each
% sample prevents the strongest samples from dominating.  For the diagonal
% QPSK constellation used by this detector, ideal symbols have x^4=-1;
% therefore -mean(x^4) contains four times the carrier phase.  A short
% observation (nominally 1 ms, at least eight chips) limits smearing from
% residual frequency while averaging over many modulation symbols.
initialPhaseSamples=min(numel(x),max(8*sps,round(1e-3*fs)));
initialBlock=x(1:initialPhaseSamples);
initialMagnitude=abs(initialBlock);
magnitudeFloor=max(initialMagnitude)*1e-3;
validInitial=initialMagnitude > magnitudeFloor;

initialPhaseRad=0;
if any(validInitial)
    normalizedInitial=initialBlock(validInitial)./initialMagnitude(validInitial);
    fourthMoment=mean(normalizedInitial.^4);
    if abs(fourthMoment) >= 0.05
        initialPhaseRad=angle(-fourthMoment)/4;
    else
        warning(['Initial fourth-power phase estimate is weak; ', ...
                 'starting the Costas loop at zero phase.']);
    end
else
    warning('Initial samples have insufficient magnitude; starting at zero phase.');
end

phaseState=initialPhaseRad;           % radians
freqState=0;                          % radians/sample

for n=1:numel(x)
    lamsav(n)=phaseState;
    freqsav(n)=freqState*fs/(2*pi);
    y(n)=x(n)*exp(-1j*phaseState);

    yi=real(y(n));
    yq=imag(y(n));
    den=abs(yi)+abs(yq);
    if den > eps(class(den))
        % Amplitude-normalized QPSK phase detector.  Its small-error gain
        % is approximately one radian/radian for any signal amplitude.
        en=(sign(yi)*yq-sign(yq)*yi)/den;
    else
        en=0;
    end
    ensav(n)=en;

    freqState=freqState+gi*en;
    phaseState=phaseState+freqState+gp*en;
end

if(fig>0)
    % The raw detector error contains modulation and chip-transition terms.
    % Overlay a short moving mean to expose the component that drives the
    % loop without changing the loop itself.
    errorMeanSamples=max(1,round(0.5e-3*fs));
    errorMean=movmean(ensav,errorMeanSamples);

    sfigure(fig);clf;
    subplot(411);
    showiq(y,fs);
    title('corrected signal');
    ylabel('I/Q');
    subplot(412);
    plot(t,ensav,'Color',[0.75 0.85 1]); hold on;
    plot(t,errorMean,'b','LineWidth',1.25); hold off;
    title('phase-detector error');
    ylabel('\phi error');
    legend('raw','0.5 ms mean','Location','best');
    subplot(413);
    plot(t,-lamsav);
    title(sprintf('applied phase correction (initial estimate %.3f rad)', ...
                  initialPhaseRad));
    ylabel('phase (rad)');
    subplot(414);
    plot(t,-freqsav);
    title('integral-path frequency correction');
    xlabel('time (s)');
    ylabel('frequency (Hz)');


%     h=scatterplot(x,1,0,'b');hold on;
%     h=scatterplot(y,1,0,'yx',h);


    %hi=eyediagram(y(N:N+75*sps).',2*sps,sps,0,'b');
end

