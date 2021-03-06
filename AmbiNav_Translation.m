function T = AmbiNav_Translation(Li, Lo, d, kVec)
%AMBINAV_TRANSLATION Ambisonics translation coefficients matrix.
%   T = AMBINAV_TRANSLATION(LI,LO,D,K) computes the ambisonic translation
%   coefficients matrix T, for input ambisonics order LI, output order LO,
%   translation position vector D (given in Cartesian coordinates), and for
%   angular wavenumber K.
%
%   K may be a vector, in which case T is (LI+1)^2-by-(LO+1)^2-by-LENGTH(K).
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also GUMEROV2005, AMBINAV_INTERPOLATIONFILTERS.

%   ==============================================================================
%   This file is part of the 3D3A AmbiNav Toolkit.
%   
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2018 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

%   References:
%     [1] Gumerov and Duraiswami (2005) Fast Multipole Methods for the
%         Helmholtz Equation in Three Dimensions.
%     [2] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

narginchk(4,4);

if numel(d) == 3
    [AZIM,ELEV,R] = cart2sph(d(1),d(2),d(3));
else
    error('Translation vector D should have three elements.');
end

kLen = length(kVec);

Ni = (Li + 1)^2;
No = (Lo + 1)^2;

if all(kVec*R < AmbiNav_KDThreshold())
    Qi = eye(Ni);
    Qo = eye(No);
else
    Qi = AmbiNav_ZRotation(AZIM, ELEV, Li);
    Qo = AmbiNav_ZRotation(AZIM, ELEV, Lo).';
end

T = zeros(Ni,No,kLen);
Tz = AmbiNav_ZTranslation(kVec*R, max([Li, Lo]));
for kk = 1:kLen
    if kVec(kk)*R < AmbiNav_KDThreshold()
        T(:,:,kk) = eye(Ni,No);
    else
        T(:,:,kk) = Qi * Tz(1:Ni,1:No,kk) * Qo;
    end
end

end