function [usys, info] = usysfit(sys, base, opt)
	%USYSFIT Construct uncertain model from state-space array
	%Usage
	%  usys = USYSFIT(sys, base)
	%  usys = USYSFIT(sys, base, opt)
	%Inputs:
	%  sys: State-space array.
	%  base: Cell array of function pointers of the form @(del)(...). The
	%  argument 'del' is assumed to be a vector of parameter values. The
	%  state-space matrices of usys are of the form A = A_1 * base{1} + A_2 *
	%  base{2} + ...
	%  opt: Option set created by usysfitOptions.
	%Outputs:
	%  usys: uss with parametric uncertainty. The name of the uncertain
	%  parameters are taken from the parameter names in sys.SamplingGrid if
	%  defined.
	%Example:
	% Load the example system.
	% 	load('usysfit_example.mat');
	% 	sys.SamplingGrid
	% 		ans = 
	% 		struct with fields:
	% 			del_1: [7×5×3×9 double]
	% 			del_2: [7×5×3×9 double]
	% 			del_3: [7×5×3×9 double]
	% 			del_4: [7×5×3×9 double]
	%
	% Fit an uncertain system using 1 as the only basis function.
	% 	base = {@(del)(1)};
	% 	usys = usysfit(sys_grid, base)
	% 	usys =
	% 		Uncertain continuous-time state-space model with 2 outputs, 3 inputs,
	%			10 states, and no uncertain blocks.
	%
	% To improve accuracy, add del_1, del_2, del_3, and del_4 to the basis functions.
	% 	base = {@(del)(1), @(del)(del(1)), @(del)(del(2)), @(del)(del(3)), @(del)(del(4))};
	% 	usys = usysfit(sys_grid, base)
	% 	usys = 
	% 		Uncertain continuous-time state-space model with 2 outputs, 3 inputs, 10 states.
	% 		The model uncertainty consists of the following blocks:
	% 		del_1: Uncertain real, nominal = 0, range = [-1,1], 11 occurrences
	% 		del_2: Uncertain real, nominal = 0, range = [-1,1], 10 occurrences
	% 		del_3: Uncertain real, nominal = 0, range = [-1,1], 11 occurrences
	% 		del_4: Uncertain real, nominal = 0, range = [-1,1], 11 occurrences
	%
	% See also usysfitOptions
	if nargin < 3
		opt = usysfitOptions;
	end
	sys = squeeze(sys);
	dims = size(sys);
	dims = dims(3 : end);
	if isempty(dims)
		usys = uss(sys);
		return;
	end
	unames = {};
	uvars = {};
	dels = {};
	delgrids = cell(1, numel(dims));
	for kk = 1 : numel(dims)
		if numel(opt.Unames) >= kk
			uname = opt.Unames{kk};
		else
			uname = sprintf('del_%d', kk);
		end
		unames{kk} = uname;
		uvars{kk} = ureal(uname, 0);
		dels{kk} = linspace(-1, 1, dims(kk));
	end
	uvars = [uvars{:}]';
	[delgrids{:}] = ndgrid(dels{:});
	if numel(dims) == 1 || numel(dims) == 2
		dims2permute = [2, 1];
	else
		dims2permute = flip(1 : numel(dims));
	end
	tran = @(ar)(permute(ar, dims2permute));
	delgrids = cellfun(tran, delgrids, 'UniformOutput', false);
	del = cellfun(@(c)(c(:)), delgrids, 'UniformOutput', false);
	del = [del{:}]';
	
	sys = numcond(sys, opt);
	
	fit = @(m_grid)(fit_matrix(m_grid, del, uvars, unames, base, opt.Display));
	[uared, ua] = fit(sys.a);
	[ubred, ub] = fit(sys.b);
	[ucred, uc] = fit(sys.c);
	[udred, ud] = fit(sys.d);
	info = struct;
	info.ss.a = ua;
	info.ss.b = ub;
	info.ss.c = uc;
	info.ss.d = ud;
	info.ssred.a = uared;
	info.ssred.b = ubred;
	info.ssred.c = ucred;
	info.ssred.d = udred;
	if opt.ReduceUnc
		usys = ss(uared, ubred, ucred, udred);
	else
		usys = ss(ua, ub, uc, ud);
	end
	usys.y = sys.y;
	usys.u = sys.u;
	
	info = check_result(sys, usys, unames, info, opt);
end

function sys_cond = numcond(sys, opt)
	if ~opt.Balreal && ~opt.Prescale
		sys_cond = sys;
		return;
	end
	s = size(sys);
	ind_ref = cell(1, numel(s));
	ind_ref(1 : 2) = {':'};
	ind_ref(3 : end) = arrayfun(@(m)(num2cell(m)), ceil(s(3 : end) / 2));
	sys_cond = sys;
	% iteration variable
	ndim = ndims(sys) - 2;
	indmax = size(sys);
	indmax = indmax(3 : end);
	
	% ssbal
	if opt.Ssbal
		sys_cond = ssbal(sys_cond);
	end
	
	% balreal
	if opt.Balreal
		[~, ~, Tbal] = balreal(sys_cond(ind_ref{:}));
		ind = num2cell(ones(1, ndim));
		done = false;
		while ~done
			sys_cond(:, :, ind{:}) = ss2ss(sys(:, :, ind{:}), Tbal);
			done = true;
			for kk = 1 : ndim
				ind{kk} = ind{kk} + 1;
				if ind{kk} <= indmax(kk)
					done = false;
					break;
				end
				ind{kk} = 1;
			end
		end
	end
	
	% prescale
	if opt.Prescale
		[~, info] = prescale(sys_cond(ind_ref{:}));
		Tleft = diag(info.SL);
		Tright = diag(info.SR);
		ind = num2cell(ones(1, ndim));
		done = false;
		while ~done
			[a, b, c, d] = ssdata(sys(:, :, ind{:}));
			sys_cond(:, :, ind{:}) = ss(Tleft * a * Tright, Tleft * b, c * Tright, d);
			done = true;
			for kk = 1 : ndim
				ind{kk} = ind{kk} + 1;
				if ind{kk} <= indmax(kk)
					done = false;
					break;
				end
				ind{kk} = 1;
			end
		end
	end
end

function info = check_result(sys, usys, unames, info, opt)
	if ~isempty(unames)
		for kk = 1 : numel(unames)
			uname = unames{kk};
			if kk == 1
				sys_fit = usubs(usys, uname, linspace(-1, 1, size(sys, 2 + kk)));
			else
				sys_fit = squeeze(usubs(sys_fit, uname, linspace(-1, 1, size(sys, 2 + kk)), '-batch'));
			end
		end
	else
		sys_fit = usys;
	end
	mnames = {'a', 'b', 'c', 'd'};
	err_abs = inf(4, 1);
	err_rel = inf(4, 1);
	sqsum = inf(4, 1);
	for kk = 1 : numel(mnames)
		m = sys.(mnames{kk});
		m_fit = sys_fit.(mnames{kk});
		if ~isempty(m)
			m_err = m - m_fit;
			err_abs(kk) = max(abs(m_err(:)));
			err_rel(kk) = err_abs(kk) / mean(abs(m(:)));
			sqsum(kk) = sqrt(sum(abs(m_err(:).^2)));
		else
			err_abs(kk) = 0;
			err_rel(kk) = 0;
			sqsum(kk) = 0;
		end
	end
	% check dc gains
	dc = dcgain(sys);
	dc_fit = dcgain(sys_fit);
	dc_err = abs(dc - dc_fit);
	err_abs(end + 1) = max(abs(dc_err(:)));
	err_rel(end + 1) = err_abs(end) / mean(abs(dc_err(:)));
	sqsum(end + 1) = sqrt(sum(abs(dc_err(:).^2)));
	mnames = [mnames, 'dc'];
	for kk = 1 : numel(mnames)
		info.(mnames{kk}).abs_err = err_abs(kk);
		info.(mnames{kk}).rel_err = err_rel(kk);
		info.(mnames{kk}).sqsum = sqsum(kk);
		if opt.Display
			if kk <= 4
				n = 'matrix';
			else
				n = 'gain';
			end
			fprintf('%s %s\n\tabs\t\trel\t\tsqrsum\n', mnames{kk}, n);
			fprintf('\t%.4f\t', err_abs(kk));
			if err_abs(kk) < 100
				fprintf('\t');
			end
			fprintf('%.4f\t', err_rel(kk));
			if err_rel(kk) < 100
				fprintf('\t');
			end
			fprintf('%.4f\n', sqsum(kk));
		end
	end
end

function [m_unc_red, m_unc] = fit_matrix(m, del, uvars, unames, base, disp)
	if isempty(m)
		m_unc_red = [];
		m_unc = [];
		return;
	end
	warning('off', 'MATLAB:nargchk:deprecated');
	straighten = @(m)(reshape(permute(m, [1, 2, flip(3 : numel(size(m)))]),...
		size(m, 1), size(m, 2), []));
	if nargin < 4
		disp = false;
	end
	% fit
	m_grid = straighten(m);
	[coeffs, eval_fit] = lsfitmx(m_grid, del, base);
	m_unc = eval_fit(coeffs, uvars, base);
	% reduce
	unc_lfr = cellfun(@(un)(lfr(un)), unames, 'UniformOutput', false);
	unc_lfr = [unc_lfr{:}];
	base_subs_lfr = cellfun(@(b)(minlfr(b(unc_lfr))), base, 'UniformOutput', false);
	base_subs_lfr = [base_subs_lfr{:}];
	mats = num2cell(coeffs, [1, 2]);
	m_lfr = gmorton(mats, base_subs_lfr);
	m_unc_red = lfr2rob(m_lfr);
	m_unc_red = umat(m_unc_red.d);
end