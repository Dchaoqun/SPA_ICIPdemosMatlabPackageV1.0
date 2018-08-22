% [VEC] = vector(MTX)
% 
% Pack elements of MTX into a column vector.  Same as VEC = MTX(:)
% Previously named "vectorize" (changed to avoid overlap with Matlab's
% "vectorize" function).

function vec = vector(mtx)

vec = mtx(:);

if iscell(mtx),
    L = length(mtx);
    vec = [];
    for n = 1:L,
        vec = [vec; vector(mtx{n})];
    end
end

