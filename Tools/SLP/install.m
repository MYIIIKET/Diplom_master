% This script compiles *.cpp with mex tool.
% You should run `mex -setup C++` firstly to configure the 
% C++ language compilation.

% mex -setup C++

mex -largeArrayDims iteration.cpp

mex -largeArrayDims find2cells.cpp