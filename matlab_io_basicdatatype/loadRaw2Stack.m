function img = loadRaw2Stack(filename, dfactor, b_readThumbNail, b_readMiddleFrameOnly)
%function img = loadRaw2Stack(filename, dfactor, b_readThumbNail, b_readMiddleFrameOnly)
%
% Load raw format generated by Hanchuan Peng as a stack.
% This function is a wrapper to call the actual loadRaw2Stack_c function.
%
% Note: the old .m function with the same name is no longer used. It was
% stored as the .nouse file.
%
% by Hanchuan Peng
% 2006-05-10
% 2007-08-19: add two new parameters about thumb nails and
% middle-frame-only
%
% Note: as of 07/08/19, the b_readThumbNail and b_readMiddleFrameOnly are
% only valid for LSM files.
%
% 2008-08-22: add dfactor to downscale the intensity for 16 bits *tif* images

if nargin<1,
    help loadRaw2Stack;
    error('Usage: img = loadRaw2Stack(filename, dfactor, b_readThumbNail, b_readMiddleFrameOnly)');
end;

if exist(filename, 'file')==0,
    error('You have specified a file which does not exist.');
end;

if ~exist('dfactor', 'var'),
    dfactor = 0;
end;

if ~exist('b_readThumbNail', 'var'),
    b_readThumbNail = 0;
end;

if ~exist('b_readMiddleFrameOnly', 'var'),
    b_readMiddleFrameOnly = 0;
end;


img = squeeze(loadRaw2Stack_c(filename, dfactor, b_readThumbNail, b_readMiddleFrameOnly));
