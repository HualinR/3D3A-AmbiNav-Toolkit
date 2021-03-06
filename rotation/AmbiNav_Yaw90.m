function Qy = AmbiNav_Yaw90(L)
%AMBINAV_YAW90 Ambisonics rotation of 90 degrees yaw.
%   Q = AMBINAV_YAW90(L) computes the ambisonic rotation coefficients
%   matrix Q, up to ambisonics order L, for a rotation of 90 degrees yaw.
%
%   See also AMBINAV_YAW.

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

N = (L + 1)^2;
[Lvec, Mvec] = getAmbOrder(0:N-1);

Qy = zeros(N);
for no = 1:N
    lo = Lvec(no);
    mo = Mvec(no);
    for ni = 1:N
        li = Lvec(ni);
        mi = Mvec(ni);
        if (lo == li) && (abs(mo) == abs(mi))
            if mo*mi >= 0
                if ~mod(mi,2)
                    Qy(no,ni) = (-1)^(mi/2);
                end
            else
                if mod(mi,2)
                    Qy(no,ni) = (-1)^((mi-1)/2);
                end
            end
        end
    end
end

end