function [out_header,out_data]=EEG_import_MEGA(input_folder,session_number);
% - Imports MEGA data = EEG signal recorded using NeurOne system (Bittium)
% - Creates files compatible with the letswave EEG processing toolbox
% 
% Input:
%   input_folder - path to the MEGA folder with subfolders from one session
%   session_number - numbers of subfolders that are to be imported
%
% Original author Andre Mouraux
% Modified by Domi for GABA-AD project (2020)

out_header=[];
out_data=[];

%recording
recording=module_read_neurone(input_folder,session_number);

%prepare header
out_header.filetype='time_amplitude';
out_header.name=[input_folder,'_',num2str(session_number)];
out_header.tags={};
out_header.history(1).configuration=[];
out_header.datasize=double([1 length(recording.signalTypes) 1 1 1 recording.properties.length*recording.properties.samplingRate]);
out_header.xstart=1/recording.properties.samplingRate;
out_header.ystart=0;
out_header.zstart=0;
out_header.xstep=1/recording.properties.samplingRate;
out_header.ystep=1;
out_header.zstep=1;

%dummy chanloc
chanloc.labels='';
chanloc.topo_enabled=0;
chanloc.SEEG_enabled=0;
%set chanlocs
for chanpos=1:length(recording.signalTypes);
    chanloc.labels=recording.signalTypes{chanpos};
    out_header.chanlocs(chanpos)=chanloc;
end;
            
%set events
out_header.events=[];
if isempty(recording.markers.index);
else
    numevents=length(recording.markers.index);
    for eventpos=1:numevents;
        event.code='unknown';
        if ~isempty(recording.markers.type(eventpos));
            event.code=recording.markers.type{eventpos};
        end;
        if isnumeric(event.code);
            event.code=num2str(event.code);
        end;
        event.latency=recording.markers.time(eventpos);
        event.epoch=1;
        events(eventpos)=event;
    end;
    out_header.events=events;
end;

%data
out_data=zeros(round(out_header.datasize));
for k=1:length(recording.signalTypes)
    eval(['out_data(1,k,1,1,:)=squeeze(recording.signal.',recording.signalTypes{k},'.data);']);
end

out_header.datasize=size(out_data);
