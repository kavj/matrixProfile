% build file

% This won't cover all of them, but Clang supports most of the basic gcc
% flags

cc = mex.getCompilerConfigurations('C++');
devel = cc.Manufacturer;

if strcmp(devel, 'GNU')
    mex('-v', '-R2018a', 'mpx_cimpl.cpp', 'pearson_split.cpp', 'CXXFLAGS="-O3 -march=native -funroll-loops -fwrapv -ffp-contract=fast"');
    mex('-v', '-R2018a', 'cross_cov_cimpl.cpp', 'cross_cov.cpp ', 'CXXFLAGS="-O3 -march=native -funroll-loops -fwrapv -ffp-contract=fast"');

elseif strcmp(devel, 'Microsoft')
   mex('-v', '-R2018a', 'mpx_cimpl.cpp', 'pearson_split.cpp', 'CXXFLAGS="/O2 /GL /arch:AVX2"');
   mex('-v', '-R2018a', 'cross_cov_cimpl.cpp', 'cross_cov.cpp ', 'CXXFLAGS="/O2 /GL /arch:AVX2"');

end

