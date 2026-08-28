function mask = to_mask(S)
%TO_MASK Occupancy lattice as a logical lesion mask.
%
%   mask = TO_MASK(S) returns S.occ directly. It exists so callers treat the
%   CA's output through a named interface rather than reaching into S, and so
%   a future off-lattice representation can change the internals without
%   touching the measurement bridge.

mask = S.occ;
end
