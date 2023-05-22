classdef usysfitOptions
	% USYSFITOPTIONS Create options set for the usysfit command.
	%
	% opt = USYSFITOPTIONS returns the default options for the usysfit
	% command.
properties
	Display = true;
	Balreal = false;
	Prescale = false;
	Ssbal = false;
	ReduceUnc = true;
	Unames = {};
end

methods
	function obj = usysfitOptions(varargin)
		cnt = 1;
		isarg = @(str)(~isempty(find(strcmp(varargin, str))));
		while cnt < nargin
			opt_name = varargin{cnt};
			try                            
				obj.(opt_name) = varargin{cnt + 1};
			catch er
				if strcmp(er.identifier, 'MATLAB:noPublicFieldForClass')
						error('Unknown oprion "%s" for the command "usysfit".', opt_name);
				else
						throw(er);
				end
			end
			cnt = cnt + 2;
		end        
	end
end
end