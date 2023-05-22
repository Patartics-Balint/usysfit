function [coeffs, eval_fit, base] = lsfitmx(A, rho, base)
	[nr, nc, nrho, nb, base] = parse_input(A, rho, base);		
	X = [];
	for krho = 1 : nrho
		X_row = [];
		for kb = 1 : nb
			X_row = [X_row, base{kb}(rho(:, krho)) * eye(nr)];
		end
		X = [X; X_row];
	end

	XX = (X' * X) \ X';
	AA = reshape(permute(A, [2, 1, 3]), [nc, nr * nrho])';
	Y = XX * AA;
	coeffs  = permute(reshape(Y', [nc, nr, nb]), [2, 1, 3]);
	eval_fit = @eval_fit_;
end

function [nr, nc, nrho, nb, base] = parse_input(A, rho, base)
	nr = size(A, 1);
	nc = size(A, 2);
	nrho = size(rho, 2);
	if nrho ~= size(A, 3)
		error('The number of matrix and parameter samples must be the same.');
	end
	if nargin < 2
		base = 0;
	end
	if isa(base, 'double') && base >= 0 && floor(base) == base
		d = base;
		base = {};
		for kk = 0 : d
			base{end + 1} = @(rho)(rho .^ kk);
		end
	end
	if isa(base, 'function_handle')
		base = {base};
	end
	nb = numel(base);
end

function A_fit = eval_fit_(coeffs, rho, base)
	if isa(base, 'function_handle')
		base = {base};
	end
	nb_ = numel(base);
	nrho_ = size(rho, 2);
	for krho_ = 1 : nrho_
		A_fit_new = zeros(size(coeffs(:, :, 1)));
		for kb_ = 1 : nb_
			A_fit_new = A_fit_new + coeffs(:, :, kb_) * base{kb_}(rho(:, krho_));
		end
		A_fit(:, :, krho_) = A_fit_new;
	end
end