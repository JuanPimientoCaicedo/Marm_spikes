


function DAQ = mReadNIDAQ(nidaqBin)
%------------------------------------------
% DAQ = mReadNIDAQ(nidaqBin)
%------------------------------------------
% reads NIDAQ file from SpikeGLX as I have it setup in lab 0267
% *** not a universal function for every setup ***
% *** it depends on how the digital and analog ***
% *** are wired to the PXIe Chassis IO card ***
%------------------------------------------
% INPUTS:
%   nidaqBin - full filepath (including file extension) 
%              to the nidaq binary file created by SpikeGLX
%              file extension should be something like [...]_g0_t0.nidq.bin
% 
% OUTPUTS:
%    DAQ - structure with fields
%
% Samples and sampling rate 
%       .fs       - sampling rate 
%       .nSamples - numbner fo samples in NIDAQ bin file 
%
% Eventmarkers 
%       .eventMarker  - [Mx1] event markers (16-bit)
%       .eventSample  - [Mx1] sample of events marker (in NIDAQ bin file) 
%       .eventTimeSec - [Mx1] time (sec) in recording 
%    -> where M = number of eventmarkers recorded 
%
% Sync bit (Sync bit is channel 385 in SpikeGLX neural recordings)
%       .syncStream  - [Nx1] stream of sync bit (length of recording)
%       .syncEdge    - [Kx1] sync edges (signed)
%       .syncSample  - [Kx1] sample of sync edge 
%       .syncTimeSec - [Kx1] time (sec) of sync edge 
%     -> where N = number of samples in NIDAQ bin file  
%     ->   and K = number of sync edges (lo->hi and hi->lo)  
%
% Analog eye data from Eyelink  
%       .eyeData    - [Nx3] x-pos, y-pos, pupil 
%     -> where N = number of samples in NIDAQ bin file  
%
% Photodiode     
%       .photodiode - [Nx1] photodiode stream (length of recording)
%     -> where N = number of samples in NIDAQ bin file  
%
% *** Note, to match the NIDAQ data to the Neural Data you need to 
%      1. align sync edges to account for drift in the clocks 
%      2. convert time in sec to samples at the fs of the Neural Data 
%        (fs ~ 30e3)
%
%------------------------------------------
% Written by Jarrod R. Dowdall (2023)
%------------------------------------------
% function version 1.0.1 - November 17, 2023
%------------------------------------------

DAQ = [];

[projDir,niFile,niExt] = fileparts(nidaqBin); %#ok<ASGLU>
meta  = readmeta(fullfile(projDir,[niFile,'.meta']));

fs = str2double(meta.niSampRate);
DAQ.fs = fs;

nSamp = str2double(meta.fileTimeSecs) * fs;
DAQ.nSamples = nSamp;

% 4 * 16-bit, 16-bit, 16-bit 
nChan = str2double(meta.nSavedChans);

%--------------------------------------
fid = fopen(nidaqBin, 'r', 'ieee-le');
binArray = fread(fid, Inf, '*int16');
fclose(fid);
%--------------------------------------
binArray = reshape(binArray,nChan,[]);

%--------------------------------------
% analog, eyes (x,y), pupil, photo-diode 
analogScaleFactor = str2double(meta.niAiRangeMax) ./ str2double(meta.niMaxInt);
XA = single(binArray(1:4,:)) * analogScaleFactor;
%-----------------------------------------
% digitial pulses
%DW = reshape(logical(typecast(binArray(5,:),'uint8')),8,[]);
bits = logical(dec2bin(typecast(binArray(5,:),'uint16')) - '0');
%-----------------------------------------
% eventmarkers 
EVTS = typecast(binArray(6,:),'uint16');
%-----------------------------------------
% 1,4
% 1 - Sync bit, 
% 4 - event trigger 
%-----------------------------------------
% DW -------------------------------------
% Event markers 
eventSample = find(diff([0; bits(:,4)])==1);
eventMarker = double(EVTS(eventSample));
eventTime   = eventSample ./ fs;

DAQ.eventMarker  = eventMarker(:);
DAQ.eventSample  = eventSample(:);
DAQ.eventTimeSec = eventTime(:);

%-----------------------------------------
% sync bit, find edges 
synstream  = bits(:,1);
[synidx,~,synedge] = find(diff([0; synstream]));

DAQ.syncStream  = synstream;
DAQ.syncEdge    = synedge;
DAQ.syncSample  = synidx;
DAQ.syncTimeSec = synidx/fs;

% XA -------------------------------------
% Analog 

% Eye Data
DAQ.eyeData    = XA(1:3,:)';
% photodiode     
DAQ.photodiode = XA(4,:)';

%-----------------------------------------

end

% =========================================================
% Parse ini file returning a structure whose field names
% are the metadata left-hand-side tags, and whose right-
% hand-side values are MATLAB strings. We remove any
% leading '~' characters from tags because MATLAB uses
% '~' as an operator.
%
% If you're unfamiliar with structures, the benefit
% is that after calling this function you can refer
% to metafile items by name. For example:
%
%   meta.fileCreateTime  // file create date and time
%   meta.nSavedChans     // channels per timepoint
%
% All of the values are MATLAB strings, but you can
% obtain a numeric value using str2double(meta.nSavedChans).
% More complicated parsing of values is demonstrated in the
% utility functions below.
%
%--------------------------------
% using 'path' overloads matlab's 'path'
% resulting in the message "Unrecognized field name "internal"."

function meta = readmeta(metapath)

% Create the matching metafile name
%[dumPath,name,dumExt] = fileparts(binName);
%metaName = strcat(name, '.meta');
%fullfile(dataPath, metaName)

% Parse ini file into cell entries C{1}{i} = C{2}{i}
fid = fopen(metapath, 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
fclose(fid);

% create structure
% fieldnames
fields = cellfun(@(str) matlab.lang.makeValidName(str), C{1}, 'unif',0);
% values 
values = C{2};

%isnum  = cellfun(@(str) isfinite(str2double(str)), values);
%values(isnum) = cellfun(@(str) eval(sprintf('[%s]',str)), values(isnum), 'unif',0);
meta = cell2struct(values,fields,1);

end




