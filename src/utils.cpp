#include "utils.h"

mata1 JuppRmat(const veca1 & y1, const veca1 & y2){
  veca1 sum = y1 + y2;
  a1type denom = 1.0 + y1.dot(y2);//(y1.transpose() * y2).coeff(0,0)
  mata1 ident = mata1::Identity(y1.size(), y1.size());
  mata1 out = ((sum * sum.transpose()).array() / denom).matrix() - ident;
  return out;
}