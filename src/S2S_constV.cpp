#include "S2S_constV.h"
#include "mnlink_cpp.h"
#include "uldSvMF.h"

mata1 JuppRmat(const veca1 & y1, const veca1 & y2){
  veca1 sum = y1 + y2;
  a1type denom = 1.0 + y1.dot(y2);//(y1.transpose() * y2).coeff(0,0)
  mata1 ident = mata1::Identity(y1.size(), y1.size());
  mata1 out = ((sum * sum.transpose()).array() / denom).matrix() - ident;
  return out;
}


veca1 ull_S2S_constV(mata1 y, mata1 xs, mata1 xe, mnlink_Omega_cpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 Kstar){
  int p = om.p1.size();
  //check that ncol(y) == p
  if (y.cols() != p){Rcpp::stop("width of y does not equal length of p1");}

  // project Omega matrix to satisfy orthogonality to p1 and q1
  mnlink_Omega_cpp<a1type> om_projected = Omega_proj_cpp(om);
  veca1 omvec_projected = mnlink_Omega_cpp_vec(om_projected);

  //get mean
  mata1 ypred;
  ypred = mnlink_cpp(xs, xe, omvec_projected, p); 

  //evaluate SvMF density of each observation
  veca1 ld(y.rows());
  mata1 Gstar = getHstar(om_projected.p1) * Kstar;
  mata1 G(p, p);
  veca1 a(p);
  a(0) = a1;
  if (aremaining.size() != p - 1){ Rcpp::stop("aremaining must have length p - 1."); }
  a.segment(1, p-1) = aremaining;
  for (int i = 0; i < y.rows(); ++i){
     G.col(0) = ypred.row(i);
    G.block(0, 1, p, p-1) = JuppRmat(om_projected.p1, ypred.row(i)) * Gstar;
    ld(i) = uldSvMF_cann(y.row(i), k, a, G)(0);
  }
  return ld;
}



veca1 ull_S2S_constV_forR(mata1 y, mata1 xs, mata1 xe, veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 Kstar){
   mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, y.cols(), xe.cols());
   veca1 ld = ull_S2S_constV(y, xs, xe, om, k, a1, aremaining, Kstar);
   return ld;
}


// Cayley Transform: C(A) = (I - A)^{-1} * (I + A)
// A is skew symmetric
// [[Rcpp::export]]
mata1 cayleyTransform(const mata1 &A) {
  int n = A.rows();  // Size of the matrix
  
  // Identity matrix
  mata1 I = mata1::Identity(n, n);

  // Calculate (I - A) and (I + A)
  mata1 I_minus_A = I - A;
  mata1 I_plus_A = I + A;

  // Cayley transform: (I - A)^{-1} * (I + A)
  mata1 result = I_minus_A.inverse() * I_plus_A;

  return result;
}

// Inverse Cayley Transform: C^{-1}(M) = (M - I) * (M + I)^{-1}
// M must be orthogonal and determinant of not -1
// [[Rcpp::export]]
mata1 inverseCayleyTransform(const mata1 &M) {
  int n = M.rows();  // Size of the matrix
  
  // Identity matrix
  mata1 I = mata1::Identity(n, n);

  // Calculate (M - I) and (M + I)
  mata1 M_minus_I = M - I;
  mata1 M_plus_I = M + I;

  // Inverse Cayley transform: (M - I) * (M + I)^{-1}
  mata1 result = M_minus_I * M_plus_I.inverse();

  return result;
}


// [[Rcpp::export]]
veca1 vectorizeLowerTriangle(const mata1 &A) {
  int n = A.rows();  // Size of the matrix
  int num_elements = (n * (n - 1)) / 2;  // Number of elements in the lower triangle

  veca1 result(num_elements);  // Vector to store the lower triangular elements

  int idx = 0;  // Index to keep track of the vector position

  // Loop through the lower triangular part of the matrix
  // like stacking the columns on top each other omitting the diagonal and upper triangle
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      result(idx) = A(i, j);  // Extract element from lower triangle
      ++idx;
    }
  }

  return result;
}

// [[Rcpp::export]]
mata1 inverseVectorizeLowerTriangle(const veca1 &vec) {
  // n is the size of the skew-symmetric matrix to reconstruct
  int n = (1 + std::sqrt(8*vec.size() + 1))/2;
  mata1 A = mata1::Zero(n, n);  // Initialize an n x n matrix with zeros

  int idx = 0;  // Index to keep track of the vector position

  // Fill the lower triangular part of the matrix
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      A(i, j) = vec(idx);   // Fill element in lower triangle
      A(j, i) = -vec(idx);  // Use skew-symmetric property: A(j, i) = -A(i, j)
      ++idx;
    }
  }

  return A;  // Return the reconstructed skew-symmetric matrix
}

// aremaining becomes log(1/a) with the first element dropped
// [[Rcpp::export]]
veca1 S2S_constV_nota1_tovecparams(veca1 & omvec, a1type k, veca1 aremaining, mata1 Kstar){
  int p = aremaining.size() + 1;
  if (Kstar.cols() != Kstar.rows()){Rcpp::stop("Kstar should be square");}
  if (Kstar.cols() != p-1){Rcpp::stop("Kstar should have dimension p-1");}
  a1type detKstar = Kstar.determinant();
  if (std::abs(CppAD::Value(detKstar) + 1.0) < 1e-8) {Rcpp::stop("Kstar has a determinant very close to -1, please change the sign of one of Kstar's columns");}
  
  veca1 vecCayaxes = vectorizeLowerTriangle(inverseCayleyTransform(Kstar)); 
  veca1 result(omvec.size() + 1 + aremaining.size() - 1 + vecCayaxes.size());
  result.segment(0, omvec.size()) = omvec;
  result(omvec.size()) = k;
  result.segment(omvec.size() + 1, aremaining.size() - 1) = aremaining.tail(aremaining.size() - 1).array().log();
  result.segment(omvec.size() + 1 + aremaining.size() - 1, vecCayaxes.size()) = vecCayaxes;
  return result;
}

// reverse function
// without ce, would only need q = qs + qe to be passed in
std::tuple<veca1, a1type, veca1, mata1> S2S_constV_nota1_fromvecparams(const veca1 & mainvec, int p, int qs, int qe) {
  if (mainvec.size() != p + (qs + qe) + p*(qs + qe) + (qe>0) * (1 + p) + 1 + (p-2) + ((p-2) * (p-1) / 2)  ) {
    Rcpp::stop("Input vector size does not match expected dimensions.");
  }
  
  veca1 omvec = mainvec.segment(0,  p + (qs + qe) + p*(qs + qe) + (qe>0) * (1 + p));
  a1type k = mainvec(omvec.size());
  veca1 laremaining_m1(p - 2); //convert log remaining to full aremaining
  laremaining_m1 = mainvec.segment(omvec.size() + 1, p - 2);
  veca1 aremaining(p-1);
  aremaining[0] = CppAD::exp(-laremaining_m1.sum());
  aremaining.tail(aremaining.size() - 1) = laremaining_m1.array().exp();
  veca1 vecCayaxes = mainvec.segment(omvec.size() + 1 + aremaining.size() - 1, (p-1) * (p-2)/2);
  mata1 Kstar = cayleyTransform(inverseVectorizeLowerTriangle(vecCayaxes));
  
  return std::make_tuple(omvec, k, aremaining, Kstar);
}

//export the reverse function
// [[Rcpp::export]]
Rcpp::List S2S_constV_nota1_fromvecparamsR(const veca1 & mainvec, int p, int qs, int qe) {
  auto result = S2S_constV_nota1_fromvecparams(mainvec, p, qs, qe);
  veca1 omvec = std::get<0>(result);
  a1type k = std::get<1>(result);
  veca1 aremaining = std::get<2>(result);
  mata1 Kstar = std::get<3>(result);
  
  // Return as a list
  return Rcpp::List::create(
    Rcpp::Named("omvec") = omvec,
    Rcpp::Named("k") = k,
    Rcpp::Named("aremaining") = aremaining,
    Rcpp::Named("Kstar") = Kstar
  );
}


pADFun tape_ull_S2S_constV_nota1(veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 Kstar, vecd & p_in, vecd & qe_in, matd & yx){
  int p = int(p_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  int qe = int(qe_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  int qs = yx.cols() - qe - p; 
  if (p!=3){Rcpp::warning("This function approximates the vMF normalising constant when p!=3.");}
  if (qs < 0){
    Rcpp::stop("yx has insufficient columns to match p and qe.");
  }
  // separate the response the covariates
  mata1 y = yx.leftCols(p);
  mata1 xs = yx.rightCols(qs + qe).leftCols(qs);
  mata1 xe = yx.rightCols(qe);

  // Get all parameters except a1 into a vector
  veca1 mainvec = S2S_constV_nota1_tovecparams(omvec, k, aremaining, Kstar);
  veca1 a1vec(1);
  a1vec(0) = a1;

// tape with main vector and a1 as a dynamic
  CppAD::Independent(mainvec, a1vec);
  // split mainvec into parts, overwriting passed arguments
  auto result = S2S_constV_nota1_fromvecparams(mainvec, p, qs, qe);
  omvec = std::get<0>(result);
  k = std::get<1>(result);
  aremaining = std::get<2>(result);
  Kstar = std::get<3>(result);
  
  mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, p, qe);

  veca1 ld = ull_S2S_constV(y, xs, xe, om, k, a1, aremaining, Kstar);

  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(mainvec, ld);
  tape.check_for_nan(false);

  pADFun out(tape, mainvec, a1vec, "ull_S2S_constV_nota1");
  return(out);
}


