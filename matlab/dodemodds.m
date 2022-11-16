% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: dodemodsimple 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: dodemodds.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
%
function [bits,symbols,phii,py]=dodemodds(x,xs,spc,fig)

clf;
SAMPLES_PER_CHIP = spc;
NUM_SYMBOLS = 150;
NUM_PREAMBLE_SYMBOLS = 25;
CHIPS_PER_SYMBOL = 256;

pad=10;
x = [ zeros(1,pad) x zeros(1,pad) ];

b = [ 1 2 3 4 3 2 1 ]; b = b/sum(b);

start_idx = round(pad*3/4+length(b)/2);

x = filter( b, 1, x );
LX = length(x);
    
% loop over 150 correlation results
idx1 = [ start_idx+1 : start_idx+SAMPLES_PER_CHIP*CHIPS_PER_SYMBOL]; % offset data, 256 chips per bit, oversampled by four
idx2 = [ 1 : SAMPLES_PER_CHIP*CHIPS_PER_SYMBOL ];
result_sum = zeros(NUM_SYMBOLS,1);
result_idx = zeros(NUM_SYMBOLS,1);
for k=1:NUM_SYMBOLS

    % correlate over conjugatged and non-conjugated spreading sequences

    corr1 = dot(x(idx1),xs(idx2) ); 
    corr2 = dot(x(idx1),conj(xs(idx2))); 

    % use the sum over the correlations
    if ( abs(corr1)  > abs(corr2) )
        result_idx(k) = 1;
        result_sum(k) = corr1;
    else
        result_idx(k) = 2;
        result_sum(k) = corr2;
    end;

    % advance indices to next symbol
    idx1 = idx1 + SAMPLES_PER_CHIP*CHIPS_PER_SYMBOL;
    idx2 = idx2 + SAMPLES_PER_CHIP*CHIPS_PER_SYMBOL;

end;

symbols = result_sum;

% make initial phase correction using first five symbols
phi = angle( sum( symbols(1:5) ) );
symbols = symbols * exp( -j*phi );
phii=phi;

% % make final phase trend removal
py = zeros(1,NUM_SYMBOLS/5);
r = 1:5;
for k=1:NUM_SYMBOLS/5
    
    % calculate phase of next block of 5 symbols
    if ( r(end) <= NUM_PREAMBLE_SYMBOLS )
        sym_sum = sum( symbols(r) );
    else
        sym_sum = 0;
        for n=r
            if ( real( symbols(n) ) > 0 )
                sym_sum = sym_sum + symbols(n);
            else
                sym_sum = sym_sum - symbols(n);
            end;
        end;
    end;
	phi = angle( sym_sum );    
   
        % update cumulative phase correction
    for n=k:NUM_SYMBOLS/5
        py(n) = py(n) + phi;
    end
    
    % apply correction cumulatively to symbols
	symbols(r(1):end) = symbols(r(1):end) * exp( -j*phi );
    
    % increment range vector
    r = r + 5;
    
end;


% make bit decisions
bits = zeros(1,2*NUM_SYMBOLS);
for k=1:NUM_SYMBOLS
    if ( (real(symbols(k)) > 0) && (result_idx(k) == 1) )
        bits( 2*(k-1) + 1 : 2*k ) = [ 0 0 ];
    elseif ( (real(symbols(k)) > 0) && (result_idx(k) == 2) )
        bits( 2*(k-1) + 1 : 2*k ) = [ 0 1 ];
    elseif ( (real(symbols(k)) < 0) && (result_idx(k) == 1) )
        bits( 2*(k-1) + 1 : 2*k ) = [ 1 1 ];
    else
        bits( 2*(k-1) + 1 : 2*k ) = [ 1 0 ];
    end;
end;

%normalize and plot symbols

symbols = 0.7071*(symbols/mean(abs(symbols))+(-2j*result_idx+3j*ones(length(result_idx),1)));

sfigure(fig);
clf; hold on;
plot(real(symbols),imag(symbols),'.');
plot([0,0], [-1.5,1.5], 'k');
plot([-1.5,1.5], [0,0], 'k');
Constellation=.7071*[1+1i,1-1i,-1+1i,-1-1i];
SymbolMapping=[0,1,2,3];
for jj=1:4
        text(real(Constellation(jj))-0.15,...,
        imag(Constellation(jj))+0.15,...
        dec2base(SymbolMapping(jj),2,2));
    plot(real(Constellation(jj)),...,
        imag(Constellation(jj)),'rd');
end
title('symbol mapping');
drawnow;

