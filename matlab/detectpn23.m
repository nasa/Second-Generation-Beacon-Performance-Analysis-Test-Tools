% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: detectpn23 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: detectpn23.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
function crossings=detectpn23(x,maxshift,mode,sps,fig)

x=x(:).';
% get sequence
[~, preamble]=getwfPN23_upsample(sps, '', 'none', mode); %getwavPN23('',mode);
fchip = 38400;
fs = sps*fchip;

NALPHA = 500;
NFFT = 16384;

% faxis vector for displays
faxis = fs*[-NFFT/2:NFFT/2-1]/NFFT;

% create a vector for the average power
expavg = zeros(1,NFFT);
tnr    = 10;       % todo, calculate realistic value
buffer = zeros(1,NFFT);
threshold = zeros(1,NFFT);
bw = 1000;
nbins = fix(bw/(fs/NFFT));
nbins = 2*fix(nbins/2)+1;

faxis=[-(NFFT/2) : (NFFT/2-1)] *fs/NFFT;

doppler = (faxis>-40000 & faxis<40000);

% truncate and conjugate the preamble
ref = conj( preamble(1:NFFT) );

% suppress threshold crossing detection until exp average has time to
% converge
init_count = 1;

% set maximum values
maxval = 0;
maxX = zeros(1,NFFT);

disp('init exp avg');
crossings=[];
crossctr=1;

first=1;
quiettime=0;
try
for m=1:2
   buffer=x(1:NFFT);
   samples = [x(NFFT+1:end)];
    for n=1:max(NALPHA,maxshift)
        
        % shift a sample into the buffer
        buffer(1:end-1) = buffer(2:end);
        buffer(end) = samples(n);
        
        % do fft
        X = fftshift( abs( fft( buffer .* ref ) ) );
        
        % update exponential average
        alpha  = 1/init_count;
        expavg = (1-alpha)*expavg + alpha*X;
        
        % look for threshold crossings
        if ( init_count ~= NALPHA )
            init_count = init_count + 1;
        else
            if(first)
                first=0;
                break;
            end
            % calculate detection threshold
            tx = zeros(1,NFFT+(nbins-1)/2);
            tx(1:NFFT) = expavg;
            tx = filter( ones(1,nbins)/nbins, 1, tx );
            threshold = tnr * tx( (nbins-1)/2 + 1 : (nbins-1)/2 + NFFT );
            
            % look for crossings withing valid doppler region
            crossing = (X-threshold).*doppler;
            crossindex=find(crossing>0);
            
            if( ~isempty(crossindex))
                quiettime=0;
                if(fig)
                gcf;clf;
                plot( faxis, 20*log10(maxX), 'c', faxis, 20*log10(X), 'b', faxis, 20*log10(expavg), 'm', faxis, 20*log10(threshold), 'g'); grid; drawnow;
                hold on;
                plot( [faxis(1) faxis(end)], [20*log10(mean(maxX)) 20*log10(mean(maxX))],'k');
                plot( faxis(crossindex), 20*log10(X(crossindex)), 'ro');
                %title(sprintf('%s', datestr(dn+m/(24*3600))));
                hold off;
                drawnow;
                end
                for(izz=1:length(crossindex))
                                        % sample, magnitude, frequency bin
                %crossings(crossctr,:) = [(m)*fs+n-NFFT X(crossindex(izz)) crossindex(izz)];
                crossings(crossctr,:) = [ n X(crossindex(izz)) faxis(crossindex(izz)) 20*log10(X(crossindex(izz))/mean(expavg(crossindex(izz)-10:crossindex(izz)+10)))];
                crossctr=crossctr+1;
 
                end
                %crossings
                %save(matfname);

            else
                quiettime = quiettime+1;
                if(crossctr>1 && quiettime>100)
                    break;
                end
            end
            % update maximum
            [mv,mi] = max(X);
            if ( mv > maxval )
                if(fig)
                    plot( faxis, 20*log10(maxX), 'c', faxis, 20*log10(X), 'b', faxis, 20*log10(expavg), 'm', faxis, 20*log10(threshold), 'g'); grid; drawnow;
                    hold on; %plot( [faxis(1) faxis(end)], [20*log10(mean(maxX)) 20*log10(mean(maxX))],'k');
%                     Xsm = smoothnoise(20*log10(X),20);%
%                     plot(faxis,Xsm,'k');drawnow;
                    hold off;
                end
                maxval = mv;
                maxX   = X;
            end
            

        end


        
    end

end
catch e
    %getReport(e)
end




