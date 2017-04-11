function [] = MultiMatchDemo()
% Demo of MultiMatch:
% Multiple object matching using BBS for initialization and NCC for searching comparison
% Algorithm outline:
% Input: I, input image; T, template for matching
% Output: AA, array of Affine transformations
% 1. Downsample _I_ and _T_ w.r.t. some rules (same scale or to a proper scale and then adjust the scale range).
% 2. Use BBS to approximately locate possible matches. Generate a mask _M_, and only search where _M_==1 in the following steps.
% 2*. Run multiple BBS rounds with different scaling of _T_ may help constraint scale range.
% 3. Generate configs for searching.
% 4. Search for the 1st match _A1_.
% 4*. If _M_ consists of disjoint regions, then search each region for a best match.
% 5. Set _M(A1(T))=0_ and search for the 2nd match, and so on.
% 6. Stop when score is below a criterion.
% 7. Output the result.


end