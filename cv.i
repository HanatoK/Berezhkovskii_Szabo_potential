%module cv
%{
#include "cv.h"
%}
%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
};
%include "cv.h"
