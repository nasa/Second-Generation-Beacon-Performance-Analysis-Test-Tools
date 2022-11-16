% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
%      
%       Module: mesageStruct.m 
%      
%       Author:     Reese Bovard
%                   Concentric Real Time, LLC
%
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: messageStruct.m 15 2022-09-29 15:45:13Z reesebo $
% T.x01 Message Structure
% Basic structure: 202 information bits
% 48 BCH FEC bits BCH(250,202)

classdef messageStruct < handle
    properties (SetAccess = private)
                      
        msg=[];
        sfact=[];
        hex=[];                
        msgvec=[];
        id=[];

    end
    properties (Hidden = true)
                msgbit=[];
                encodedData=[];
                N;
                K;
                enc;
    end % properties
    methods (Access=public)
        % constructor
        function MS=messageStruct

             MS.msg=struct(  ...
                    'TACNumber',        {str2bin(MS,dec2bin(1010,16))}, ... %%ones(1,16)}, ...
                    'SerialNumber',     {str2bin(MS,dec2bin(101,14))}, ...
                    'CountryCode',      {str2bin(MS,dec2bin(303,10))}, ...
                    'StatusOfHomingDev',{1}, ...
                    'RLSFunction',      {0}, ...
                    'SelfTest',         {0}, ...        % 1=normal, 0=self test
                    'EncodedLocation',  {struct( 'NS',           {1}, ...            % North=0, South=1
                        'DegreesN',     {ones(1,7)}, ...    % Degrees (0 to 90) in 1 degree increments
                        'FractionN',     {str2bin(MS,'000001111100000')}, ...   % Decimal parts of a degree (0.5 to 0.00003 – 3.4metres resolution)
                        'EW',           {1}, ...            % East=0, West=1
                        'DegreesE',     {ones(1,8)}, ...    % Degrees (0 to 180) in 1 degree increments
                        'FractionE',     {str2bin(MS,'111110000011111')} ...   % Decimal parts of a degree (0.5 to 0.00003 – 3.4metres resolution)                        
                        )}, ...
                    'VesselIDType',         {zeros(1,3)}, ...
                    'VesselID',         {zeros(1,44)}, ...
                    'BeaconType',       {zeros(1,3)}, ...
                    'CancellationMessageStatus',       {zeros(1,14)},...
                    'rfRotatingFieldIdentfier', {zeros(1,4)}, ...              % 0000 – G.008 Objective Requirements (this table) 0001 - Inflight Emergency (TBD) 0010 – Return Link Service (TBD) 0011 - National Use (TBD) 0100 to 1111 - Spares ( for future use)
                    'rfElapsedTimeSinceActivation', {zeros(1,6)}, ...                %0 to 63 hours in one hour steps (actual time since activation shall be truncated, not rounded e.g. between 1 hour and 2 hours after activation shall be encoded as 1 hour).   
                    'rfTimeFromLastEncodedLocation', {zeros(1,11)}, ...  % 0 to 2047 minutes (34 hours and 7 minutes) in one minute steps
                    'rfAltitudeofEncodedLocation', {zeros(1,10)}, ...     % -400 metres to 15968 metres in steps of 16 metres             
                    'rfDilutionofPrecision', {zeros(1,8)}, ...           % The value of HDOP of the encoded location shall be reported followed by the value of VDOP
                    'rfActivationMethod', {zeros(1,2)}, ...                   
                    'rfRemainingBatteryCapacity', {zeros(1,3)}, ...
                    'rfGNSSStatus', {zeros(1,2)}, ...                                   
                    'rfFill',           {zeros(1,2)} ...                     %'pad',              {zeros(1,22)} ...
                    );

       MS.sfact=struct('FractionN',                  {struct('premult', 2^15, 'posmult', 1, 'preadd', 0, 'postadd', 0)}, ...
                    'FractionE',                  {struct('premult', 2^15, 'posmult', 1, 'preadd', 0, 'postadd', 0)}, ...
                    'rfAltitudeofEncodedLocation',{struct('premult', 1, 'posmult', 1/16, 'preadd', 400, 'postadd', 0)} ...
                    );
                
                %                     'deltaTfromLastEncodedLocation', {zeros(1,11)}, ...
%                     'AltitudeofEncodedLocation', {zeros(1,9)}, ...
%                     'DOP',              {ones(1,7)}, ...
%                     'EmergencyType',    {ones(1,4)}, ...
%                     'ActivationType',   {zeros(1,2)}, ...
%                     'BatteryLevel',     {ones(1,3)}, ...
%                     'CancellationMessageStatus',        {zeros(1,4)}, ...
%                     'Reserved',         {zeros(1,5)}, ...
%                     'Messaging',        {zeros(1,10)},...
%                     'GNSSStatus',      {zeros(1,4)}, ...

            getall(MS);
            MS.hex=bin2hexi(MS,MS.msgvec,1);

            MS.N=250;
            MS.K=202;
            MS.id=getBeaconID(MS);
            MS.bchencode;
            %MS.enc = comm.BCHEncoder(MS.N,MS.K);
        end
       function o=getMessagelength(MS)

          o=length(MS.msgvec);
          fprintf('length of msgvec = %d bits, length of hex=%d chars\n', o, length(MS.hex));
       end
       function o=getBeaconID(ms)
            o=getHex23(ms);
       end
       function o=getbin(ms)
           o=str2bin(ms,hex2bin(ms,ms.hex,4*length(ms.hex)-2,1));
       end
       function o=getFullMessageBits(ms)
           ms.bchencode;
           o=str2bin(ms,hex2bin(ms,ms.hex,252,1));
           o=o(3:end);
       end
        % show: tabular output of all fields with corresponding bit 
        % positions, lengths and values.
        function o=show(MS)
                fn = fieldnames(MS.msg);
                cnt=1;
                MS.msgvec=[];
                
                colnames={'startbit', 'stopbit', 'length', 'hex', 'dec', 'scaled', 'bits'};
                content={};
                rownames={};
                for(ix=1:numel(fn))

                    if(~isstruct(MS.msg.(fn{ix})))
                        ln=length(MS.msg.(fn{ix}));
                        content{cnt,1} = MS.msgbit(cnt);
                        content{cnt,2} = MS.msgbit(cnt)+ ln-1;
                        content{cnt,3} = ln;
                        content{cnt,4} = ['x.' bin2hexi(MS, MS.msg.(fn{ix}))];
                        content{cnt,5} = hex2dec(content{cnt,4}(3:end));
                        content{cnt,6} = unscale(MS,fn{ix},content{cnt,5});
                        content{cnt,7} = [ 'b.' hex2bin(MS,content{cnt,4}(3:end),ln)];
                        rownames{cnt}=fn{ix};
                        cnt=cnt+1;
                        
                    else
                        sfn=fieldnames(MS.msg.(fn{ix}));
                        for(iy=1:numel(sfn))
                            ln=length(MS.msg.(fn{ix}).(sfn{iy}));
                            content{cnt,1} = MS.msgbit(cnt);
                            content{cnt,2} = MS.msgbit(cnt)+ ln-1;
                            content{cnt,3} = ln;
                            content{cnt,4} = ['x.' bin2hexi(MS, MS.msg.(fn{ix}).(sfn{iy}))];
                            content{cnt,5} = hex2dec(content{cnt,4}(3:end));
                            content{cnt,6} = unscale(MS,sfn{iy},content{cnt,5});                            
                            content{cnt,7} = [ 'b.' hex2bin(MS,content{cnt,4}(3:end),ln)];


                            rownames{cnt}= [fn{ix} '.' sfn{iy}];
                            cnt=cnt+1;
                        end
                    end

                end
                o=cell2table(content,'rownames', rownames, 'variablenames', colnames)
        end

       %%% Methods specific to fields
       function setTAC(MS, in)
           MS.setbits('TACNumber', str2bin(dec2bin(in,16)) );
       end
       function setSerialNumber(MS, in)
           MS.setbits('SerialNumber', str2bin(dec2bin(in,14)) );
       end
      function setCountryCode(MS, in)
           MS.setbits('CountryCode', str2bin(dec2bin(in,10)) );
      end
      
      % ****** added by JC July 24, 2015 ******
      function setHoming(MS, in)
           MS.setbits('StatusOfHomingDev', str2bin(dec2bin(in,1)) );
      end
      function setSelfTest(MS, in)
           MS.setbits('SelfTest', str2bin(dec2bin(in,1)) );
      end
      
      function setBeaconType(MS, in)
           MS.setbits('BeaconType', str2bin(dec2bin(in,3)) );
      end
      function setRLS(MS, in)
           MS.setbits('RLSFunction', str2bin(dec2bin(in,1)) );
      end  
      function setCancellationMessageStatus(MS, in)
           MS.setbits('CancellationMessageStatus', str2bin(dec2bin(in,14)) );
      end
      
      % ****** modified by JC
%       function setRotate(MS, in)
%            MS.setbits('Rotate', str2bin(dec2bin(in,48)) );
%       end
      function setrfRotatingFieldIdentfier(MS, in)
           MS.setbits('rfRotatingFieldIdentfier', str2bin(dec2bin(in,4)) );
      end
      
      function setrfElapsedTimeSinceActivation(MS, in)
           MS.setbits('rfElapsedTimeSinceActivation', str2bin(dec2bin(in,6)) );
      end
      
      function setrfTimeFromLastEncodedLocation(MS, in)
           MS.setbits('rfTimeFromLastEncodedLocation', str2bin(dec2bin(in,11)) );
      end
 
      
      function setrfAltitudeofEncodedLocation(MS, in)
          alt_scale = scale(MS,'rfAltitudeofEncodedLocation',in); %floor((in+400)/16);
                    
          %ret = unscale(MS,'rfAltitudeofEncodedLocation',alt_scale); %floor((in+400)/16);

          if(alt_scale<0)
              alt_scale = 1;
          end
          if(alt_scale>2^10)
              alt_scale = 2^10;
          end
           MS.setbits('rfAltitudeofEncodedLocation', str2bin(dec2bin(alt_scale,10)) );
      end
      function out=getrfAltitudeofEncodedLocation(MS)
          alt=MS.msg.('rfAltitudeofEncodedLocation');
          out = unscale(MS,'rfAltitudeofEncodedLocation',hex2dec(bin2hexi(MS,alt,0))); %floor((in+400)/16);
      end
      function setrfDilutionofPrecision(MS, in)
           MS.setbits('rfDilutionofPrecision', str2bin(dec2bin(in,8)) );
      end     

      function setrfActivationMethod(MS, in)
           MS.setbits('rfActivationMethod', str2bin(dec2bin(in,2)) );
      end  
 
      function setrfRemainingBatteryCapacity(MS, in)
           MS.setbits('rfRemainingBatteryCapacity', str2bin(dec2bin(in,3)) );
      end  
      
      function setrfGNSSStatus(MS, in)
           MS.setbits('rfGNSSStatus', str2bin(dec2bin(in,2)) );
      end   
      
      function setrfFill(MS, in)
           MS.setbits('rfFill', str2bin(dec2bin(in,2)) );
      end          
      %***************************************
      
      
      function setVesselID(MS, in)
           MS.setbits('VesselID', str2bin(dec2bin(in,44)) );
      end
      function setVesselIDType(MS, in)
           MS.setbits('VesselIDType', str2bin(dec2bin(in,3)) );
      end
      function o=getVesselID(ms)
          msg=ms.msg;
          o=bin2hexi(ms,[msg.VesselID]);
      end
      function setEncodedLongitude(ms, in)
          
          if( abs(in) > 180)
              disp('error: please maintain -180< Input < 180');
              return;
          end
          
          if(in<0)
              ew=1;
          else
              ew=0;
          end
          ms.setbits('EncodedLocation.EW', ew);
          in=abs(in);
          base = floor(in);
          frac = min(2^15-1,round((in-base) * (2^15)) );
          
          ms.setbits('EncodedLocation.DegreesE', str2bin(dec2bin(base,8)));
          ms.setbits('EncodedLocation.FractionE', str2bin(dec2bin(frac,15)));
          %disp(['long output: ' dec2bin(ew) ' ' dec2bin(base,8) ' ' dec2bin(frac,15)]);

      end
      function setEncodedLatitude(ms, in)
          if( abs(in) > 90)
              disp('error: please maintain -90< Input < 90');
              return;
          end
          
          if(in<0)
              ns=1;
          else
              ns=0;
          end
          ms.setbits('EncodedLocation.NS', ns);
          
          in=abs(in);
          base = floor(in);
          frac = min(2^15-1,round((in-base) * (2^15)) );
          
          ms.setbits('EncodedLocation.DegreesN', str2bin(dec2bin(base,7)))
          ms.setbits('EncodedLocation.FractionN', str2bin(dec2bin(frac,15)))
          
          %disp(['long output: ' dec2bin(ns) ' ' dec2bin(base,7) ' ' dec2bin(frac,15)]);


      end
      function setvec(MS,vec)
        setvecb(MS,MS.str2bin(MS.hex2bin(vec,202,1,-2)));
      end
      function setvecb(MS,vec)
        
        MS.msgvec=vec;
        fn = fieldnames(MS.msg);

        cnt=1;
        for(ix=1:numel(fn))

            if(~isstruct(MS.msg.(fn{ix})))
                ln=length(MS.msg.(fn{ix}));
                MS.msg.(fn{ix}) = MS.msgvec(MS.msgbit(cnt):MS.msgbit(cnt)+ln-1);
                cnt=cnt+1;
            else
                sfn=fieldnames(MS.msg.(fn{ix}));
                for(iy=1:numel(sfn))
                    ln= length(MS.msg.(fn{ix}).(sfn{iy}));
                    MS.msg.(fn{ix}).(sfn{iy}) = MS.msgvec(MS.msgbit(cnt):MS.msgbit(cnt)+ln-1);
                    cnt=cnt+1;
                end
            end

        end
        getall(MS);
        bchencode(MS);
        MS.id=getBeaconID(MS);

      end
      function bchencode(ms)
          %ms.encodedData = step(ms.enc,ms.msgvec(:))';
          ms.encodedData =ms.sgbbchencode( ms.msgvec(:));
          ms.hex=bin2hexi(ms,[0 0 ms.encodedData],1);
          %disp(['updated hex:' ms.hex])
      end
      function showall(ms)
          disp(ms.msg);
          disp(ms.msg.EncodedLocation);
          disp(['Hex63:' ms.hex])
          disp(['BeaconID: ' getBeaconID(ms)]);
          show(ms);
      end
      function showhex(ms)
          disp(['Hex63:' ms.hex])
          disp(['BeaconID: ' getBeaconID(ms)]);
      end
      function o=getHex23(ms)
          msg=ms.msg;
          o=bin2hexi(ms,[ 1 msg.CountryCode 1 0 1 msg.TACNumber msg.SerialNumber msg.SelfTest msg.VesselIDType msg.VesselID]);
      end

     function o=hex2bin(ms,x,numbits,varargin)

        if(nargin>3)
            front=varargin{1};
        else
            front=0;
        end

       
        o=[];

        for i=1:length(x)
            temp = sprintf('f%s', x(i));
            temp = dec2bin(hex2dec(temp));
            o = [ o temp(5:end)];
        end

        if(nargin>4)
            rightshift=varargin{2};
            if(rightshift>0)
                o=o(1:end-rightshift);
            elseif(rightshift<0)                  
                o=[o(abs(rightshift)+1:end) 0 0];
            end
        end
        
        if(~front)
            outdelta = length(o) - numbits+1;
            if(outdelta > 0)
                o = o(outdelta:end);
            end
        else
            if(length(o) >numbits)
                o=o(1:numbits);
            end
        end
     end
      function o=str2bin(ms,x)

        for i=1:length(x)
            if(x(i) == '1')
                o(i) = 1;
            else
                o(i) = 0;
            end
        end
      end
      function o = sgbbchencode( ms, m )
            % m is the binary message vector
            % note n,k are not used here...the generator polynomial is created
            % elsewhere using g=bchgenpoly(n,k);
            % g is the generator polynomial
            g=[1 1 1 0 0 0 1 1 1 1 1 1 0 1 0 1 1 1 0 0 0 0 1 0 1 1 1 0 1 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 1 0 1 1 1];

            % get vector lengths
            LM = length(m);
            LG = length(g);
            LY = LM + LG - 1;

            % initialize vectors
            y = zeros(1,LY);
            y(1:LM) = m;
            x = y;

            % find the remainder
            idx = 1:LG;
            for k=1:LM
                if ( x(k) == 1 )
                  x(idx) = xor( x(idx), g );
                end;
                idx = idx + 1;
            end;

            % attach remainder
            y(LM+1:end) = x(LM+1:end);  
            o=y;
      end
    
          
      function help(ms)
fprintf('Method								Description                                                                                 \n');
fprintf('_______                            _____________                                                                               \n');
fprintf('bchencode                           Generate the BCH code of beacon                                                            \n');
fprintf('getrfAltitudeofEncodedLocation      Returns decimal altitude                                                                   \n');
fprintf('messageStruct                       Constructor                                                                                \n');
fprintf('getMessageLength                    Returns the length in bits of message                                                      \n');
fprintf('setCountryCode                      Set country code, input decimal value                                                      \n');
fprintf('setEncodedLatitude                  Floating point lat in                                                                      \n');
fprintf('setEncodedLongitude                 Floating point lng in                                                                      \n');
fprintf('setHoming                           Homing bit 1=on                                                                            \n');
fprintf('setSelfTest                         Self Test bit 1=on                                                                         \n');
fprintf('setSerialNumber                     Set serial number, input decimal number                                                    \n');
fprintf('setCancellationMessageStatus                        Input decimal number, converts to 17 bits                                                  \n');
fprintf('setTAC                              Input decimal TAC number                                                                   \n');
fprintf('setUserCancel                       Input user cancel bit 1=cancel                                                             \n');
fprintf('setVesselIDType                     Input decimal vessel ID Type                                                               \n');
fprintf('setVesselID                         Input decimal vessel ID                                                                    \n');
fprintf('setrfActivationMethod               Input decimal activation method                                                            \n');
fprintf('setrfAltitudeofEncodedLocation      Input altitude integer, 16m resolution                                                     \n');
fprintf('setrfDilutionofPrecision            Input DOP from GPS, input decimal number                                                   \n');
fprintf('setrfElapsedTimeSinceActivation     Input elapsed time in seconds                                                              \n');
fprintf('setrfFill                           Input fill, decimal                                                                        \n');
fprintf('setrfGNSSStatus                     Input GNSS Status, decimal                                                                 \n');
fprintf('setrfRemainingBatteryCapacity       Input battery capacity, decimal                                                            \n');
fprintf('setrfRotatingFieldIdentfier         Input rotating field id, decimal                                                           \n');
fprintf('setrfTimeFromLastEncodedLocation    Input elapsed time from last encoded location, decimal seconds                             \n');
fprintf('setvec                              Insert a full set of message bits into the structure                                       \n');
fprintf('						for example: ms.setvec(''FFFFF195528FE0F81FFF07C0000000000000000000000000000'') \n');
fprintf('show                           		Show a table of the structure                                                           \n');
fprintf('showall                      		Show even more than show                                                                    \n');
      end
    end % public methods
    methods (Access=private)

     function out=scale(MS, fld, in)
          try
            out = round(((MS.sfact.(fld).premult*in)+MS.sfact.(fld).preadd)*MS.sfact.(fld).posmult+MS.sfact.(fld).postadd);
          catch
            out=in;
          end
          
      end
      function out=unscale(MS,fld,in)
        try
            out = ((((in-MS.sfact.(fld).postadd)/MS.sfact.(fld).posmult)-MS.sfact.(fld).preadd)/MS.sfact.(fld).premult);
        catch
            out=in;
        end
      end
            % getall: rebuild the msgvec from the msg structure.
        function out = get(ms,fld)
          out=pullbits(ms,fld);
        end
        function getall(MS)
                fn = fieldnames(MS.msg);
                cnt=1;
                MS.msgvec=[];
                for(ix=1:numel(fn))

                    if(~isstruct(MS.msg.(fn{ix})))
                        MS.msgbit(cnt) = length(MS.msgvec)+1;
                        cnt=cnt+1;
                        MS.msgvec=[MS.msgvec MS.msg.(fn{ix})];
                    else
                        sfn=fieldnames(MS.msg.(fn{ix}));
                        for(iy=1:numel(sfn))
                            MS.msgbit(cnt) = length(MS.msgvec)+1;
                            cnt=cnt+1;
                            MS.msgvec=[MS.msgvec MS.msg.(fn{ix}).(sfn{iy})];
                        end
                    end

                end
        end
        function setbits(MS,ifn,bits)

            MS.encodedData =[]; % clear out previously encoded
            
            ptr=strfind(ifn,'.');
            if(isempty(ptr))
                MS.msg.(ifn) = bits;
            else
                ifn2=ifn(ptr+1:end);
                ifn=ifn(1:ptr-1);
                MS.msg.(ifn).(ifn2) = bits;
            end
            vc = struct2cell(MS.msg);
            ovc={};
            for(ix=1:length(vc))
                if(isstruct(vc{ix}))
                    ovc = [ovc struct2cell(vc{ix})'];
                else
                    ovc = [ovc vc(ix)];
                end                               
            end
            vc=ovc;
            vm=[];
            for(ix=1:length(vc))
                fd = cell2mat(vc(ix));
                vm=[vm ;fd(:)];
            end
            MS.msgvec=vm';
            %MS.hex=bin2hexi(MS,MS.msgvec,1);
            MS.bchencode;
            MS.id=getBeaconID(MS);
            %disp(MS.msg)
            %disp(['updated hex:' MS.hex])
        end
        function o=pullbits(MS,ifn)
               ptr=strfind(ifn,'.');
            if(isempty(ptr))
                o=MS.msg.(ifn);
            else
                ifn2=ifn(ptr+1:end);
                ifn=ifn(1:ptr-1);
                o=MS.msg.(ifn).(ifn2);
            end

        end
        function o=pulldec(MS,ifn)
            bits=pullbits(MS,ifn);
            o=bin2dec(hex2bin(bin2hexi(MS,bits),length(bits)));
        end
       function o=pullhex(MS,ifn)
            bits=pullbits(MS,ifn);
            o=bin2hexi(MS,bits);
       end
       function fn=getFields(MS)
           fn = fieldnames(MS.msg);
       end

       function o=bin2hexi(MS, msg, varargin)

        if(nargin>2)
            front=varargin{1};
        else
            front=0;
        end
        
        sz=size(msg);

        ln=length(msg);
        
        if(ln<4)
            o=dec2hex(bin2dec(num2str(msg)));
            return;
        end
        
        for(ix=1:sz)
            message=msg(ix,:);
            begin=0;
            m = mod(length(message),4);
            endcap=[];
            if(m>0)
                if(~front)
                    x=message(1:m);
                    message = message(m+1:end);
                    s=sprintf('%d',x);
                    oo(1)=dec2hex(bin2dec(s));
                    begin=1;
                else
                    message=[message zeros(1,m)];
                    begin=0;
                end
            end
            y=reshape(message,4,length(message)/4);

            for(i=1:length(message)/4)
                x=y(:,i);
                s=sprintf('%d%d%d%d',x(1),x(2),x(3),x(4));
                oo(begin+i)=dec2hex(bin2dec(s));
            end
            o(ix,:) = oo;
        end
        
       end
               
    end % private methods
end
        