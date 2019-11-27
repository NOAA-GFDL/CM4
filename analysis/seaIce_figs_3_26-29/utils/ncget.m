function out = ncget(infile,query,r1,r2,r3,r4) 
% This function reads variables from a netcdf file.  
% usage: variable_out = getnc('input_file','variable_name',r1,r2,r3,r4)
%
% See below for ways to get information about the netcdf file!
%
% Inputs 
%   required:
%     input_file = the name of the file to open; can include path 
%        (e.g. '/home/gerber/sample.nc' or 'sample.nc')
%     variable_name = the variable to read (e.g. 'lat')
%
%   optional:
%     r1,r2,r3,r4 the range of the variable to get. If no range is
%     given, you'll get the entire variable.  Make sure to
%     ask only for dimensions available; otherwise you'll get an error.
%     Specify the range like any matlab array.  Some examples.  
%        Suppose 'ps' is 128 x 64 x 100 (lon, lat, time) and 
%     1) You'd like just half of it (one hemisphere)
%        >  ps_nh=ncget(input_file,'ps',:,33:64,:);
%     2) You'd like just date 27
%        >  ps27=ncget(input_file,'ps',:,:,27);
%     3) You'd like every other longitude and latitude
%        >  ps_course=ncget(infile,'ps',1:2:127,1:2:63,:);
%     4) There are many ways to ask for the entire variable ps:
%        >  ps=ncget(input_file,'ps');
%        >  ps=ncget(input_file,'ps',:,:,:);
%        >  ps=ncget(input_file,'ps',1:128,1:64,1:100); etc.
%
%   Unfortunately this script can only handle up to 4 dimensional
%   arrays.
%
% Output:
%    The variable (and range) you requested.
%
% To get information about a file, just give the input_file name
%       > ncget(input_file)
%   or  
%       > ncget(input_file,'info')
% It will write out all the variables, giving the size of
% dimensions, and the dimensions of the input variables.
% 
% For information on a specific variable, try
%       > out=ncget(input_file,variable,'info')
% It will right out the size of dimensions and the dimensions of
% variables, and out will contain the size of the dimension/variable.


% usage and information:
if nargin == 0
  disp('You must at least specify an input file.  Type help ncget')
  disp('for usage information.')
  out=0; return
elseif  nargin == 1 | (nargin==2 & strcmp(query,'info'))
  % If only the file name is specified, or info requested, describe
  % what's in the file!  
  if exist(infile,'file')
    disp('The dimensions (lengths) and variables (dimensions) available in:')
    disp(infile)
    % open the file
    fin=netcdf.open(infile,'NC_NOWRITE');
    % obtain some general information on the file
    [ndims,nvars,natts,unlimdimID]= netcdf.inq(fin);
    % step through the variables
    for vstep=0:nvars-1
      [varname,xtype,dimids,natts] = netcdf.inqVar(fin,vstep);
      if size(dimids)==1
	vec=length(netcdf.getVar(fin,vstep));
	disp(strcat(varname, ' (',num2str(vec),')'))
      else
	dinfo='';
	for step=1:length(dimids),
	  dinfo=strcat(dinfo,netcdf.inqDim(fin,dimids(step)));
	  if step~=length(dimids)
	    dinfo=strcat(dinfo,', ');
	  end
	end
       disp(strcat(varname,' (',dinfo,')'))
      end
    end 
    % close up the file
    netcdf.close(fin);
    out=0; return
  else
    disp('You must specify a valid input_file_name.')
    out=0; return
  end
  
elseif ( ~exist(infile,'file') | ~ischar(query) )
  % chec for usage
  disp('First two inputs must be an input_file_name and a variable_name.')
  out=0; return 

elseif nargin==3 & strcmp(r1,'info')
  % Lastly, one is just looking for information on a single variable
  % open the file
  fin=netcdf.open(infile,'NC_NOWRITE');
  % obtain some general information on the file
  [ndims,nvars,natts,unlimdimID]= netcdf.inq(fin);
  
  % this "toggle" is used in case it can't find the variable
  toggle=0;
  % step through the variables, pick out query, as defined above
  for vstep=0:nvars-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(fin,vstep);
    if(strcmp(varname,query)),
      if size(dimids)==1
	vec=length(netcdf.getVar(fin,vstep));
	disp(strcat(varname, ' (',num2str(vec),')'))
      else
	dinfo='';
	for step=1:length(dimids),
	  dinfo=strcat(dinfo,netcdf.inqDim(fin,dimids(step)));
	  if step~=length(dimids)
	    dinfo=strcat(dinfo,', ');
	  end
	  [dummy dimlen]=netcdf.inqDim(fin,dimids(step));
	  vec(step)=dimlen;
	end
	disp(strcat(varname,' (',dinfo,')'))
      end
      toggle=1;
      out=vec;
    end
  end 
  
  % close up the file
  netcdf.close(fin);
  
  if toggle==0
    disp(strcat('Requested variable (',query,') not found!'))
    out=0;
  end 
 
  return
end


% okay -- we want more than just information.  Let's get the variable!

ndim=nargin-2;

% if no range arguments are specified, take everything
if ndim==0

  % open the file
  fin=netcdf.open(infile,'NC_NOWRITE');
  % obtain some general information on the file
  [ndims,nvars,natts,unlimdimID]= netcdf.inq(fin);
  
  % this "toggle" is used in case it can't find the variable
  toggle=0;
  % step through the variables, pick out query, as defined above
  for vstep=0:nvars-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(fin,vstep);
    if(strcmp(varname,query)),
      out=squeeze(netcdf.getVar(fin,vstep));
      toggle=1;
    end
  end 
  
  % close up the file
  netcdf.close(fin);
  
  if toggle==0
    disp(strcat('Requested variable (',query,') not found!'))
    out=0;
  end 
  
else
 
  % here it's a bit more complicated, as we're looking for only
  % part of the variable.  First, let's collect the information
  % we'll need for using netcdf.getVar.  If ':' is specified, that
  % means the user wants everything from that dimension.
  
  take_all=zeros(1,ndim);
  start=zeros(1,ndim);
  count=zeros(1,ndim);
  stride=ones(1,ndim);
  
  for step=1:ndim
    check=eval(strcat('r',num2str(step)));
    if ischar(check) 
      if strcmp(check,':')
	take_all(step)=1;
      else
	disp('Range inputs not valid. Must be vectors or :, to get all of a given dimension.');
	out=0; return
      end
    elseif isvector(r1)
      % find out where to start, recalling that nc arrays start with 0!
      start(step)=check(1)-1;
      % next, find out how much to take
      count(step)=length(check);
      % lastly, what is the stride
      if length(check)>1
	stride(step)=check(2)-check(1);
      end
    else
      disp('Range inputs not valid. Must be vectors or :, to get all of a given dimension.');
      out=0; return
    end
  end
  

  % open the file
  fin=netcdf.open(infile,'NC_NOWRITE');
  % obtain some general information on the file
  [ndims,nvars,natts,unlimdimID]= netcdf.inq(fin);
  
  % this "toggle" is used in case it can't find the variable
  toggle=0;
  % step through the variables, pick out query, as defined above
  for vstep=0:nvars-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(fin,vstep);
    if(strcmp(varname,query)),
      
      % first check to make sure we have the right size
      if length(dimids)~=ndim
	disp(sprintf('Problem: the variable you requested has %d dimensions, but you specified only %d.', length(dimids),ndim))
	netcdf.close(fin);
	out=0; return
      end
      
      % next, any colons in the range info mean the user wants all
      % of a dimension.  
      for step=1:ndim
	
	[dummy dimlen]=netcdf.inqDim(fin,dimids(step));
	
	if take_all(step)
	  count(step)=dimlen;
	  % note that start(step) will already be 0 and stride(step) 1
	else 
	  % If a range value was specified, make sure it's plausible.
	  check=eval(strcat('r',num2str(step)));
	  if check(1)<1
	    disp('Range values cannot be 0 or negative!')
	    netcdf.close(fin);
	    out=0; return
	  elseif max(check) > dimlen;
	    disp(sprintf('You requested a larger range for dimension %d than possible.',step)) 
	    disp(sprintf('Its size is %d.',dimlen))
	    netcdf.close(fin);
	    out=0; return
	  end
	end
      end
	
      % okay, with hope this should work - get the variable!
      out=squeeze(netcdf.getVar(fin,vstep,start,count,stride));
      toggle=1;
    end 
  end  
  % close up the file
  netcdf.close(fin);
  
  if toggle==0
    disp(strcat('Requested variable (',query,') not found!'))
    out=0;
  end 
  
end

% we've had trouble when the netcdf packs the variable as an integer
% This automatically converts the output to double precision real.
out=double(out);

% written by Edwin Gerber.  Please send any comments to: epg2@nyu.edu



