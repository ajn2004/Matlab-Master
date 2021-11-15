function view_parsing(varargin)

[args,pvpairs] = parseparams(varargin);
assignin('base','args',args);
assignin('base','pvpairs',pvpairs);