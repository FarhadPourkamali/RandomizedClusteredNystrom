function [ err ] = errFrobNorm( K , K_F , L )
% [ err ] = errFrobNorm( K , K_F , L ) computes normalized approximation error in terms of the
% Frobenius norm:
%                     | K - L*L'|_F / |K|_F
% where K is nxn and L is nxr
% K_F is |K|_F

sqnorm = K_F^2 + norm(L'*L,'fro').^2 - (2 * trace((L'*K)*L));
err = sqrt(sqnorm)./K_F;

end

