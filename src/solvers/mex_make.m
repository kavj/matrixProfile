% build file

% This won't cover all of them, but Clang supports most of the basic gcc
% flags

cc = mex.getCompilerConfigurations('C++');
devel = cc.Manufacturer;
mpxfiles = 'mpx_cimpl.cpp pearson_split.cpp ';
crosscovfiles = 'cross_cov_cimpl.cpp cross_cov.cpp ';
edfiles = 'euc_dist_cimpl.cpp euc_dist.cpp ';
if strcmp(devel, 'GNU')
   flags = 'CXXFLAGS="-O3 -march=native -funroll-loops -fwrapv -ffp-contract=fast"';
    
elseif strcmp(devel, 'Microsoft')
   % Matlab makes it inconvenient to check enable architectural features
   % Can't easily access setbv or mayIuse if supported. I may revisit this
   % if I add explicit vectorization Auto vectorization and pragma simd
   % probably won't do much for dtw, particularly without some manual
   % unrolling in key areas.
   
   % side side note, MSVC does not have a specific flag for floating point
   % contraction
   
   flags = '/O2 /GL /arch:AVX2';
    
end

eval(['mex ' mpxfiles '-R2018a']);
eval(['mex ' crosscovfiles '-R2018a']);
eval(['mex ' edfiles '-R2018a']);
%mex('-v', '-R2018a ', files, flags);

