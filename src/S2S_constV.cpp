#include "S2S_constV.h"
#include "mnlink_cpp.h"
#include "uldSvMF.h"
#include "utils.h"

//rG0 specifies the axes of the SvMF according to the reference coords.
//referencecoords * rG0 are back in the same coords as everything else
//(referencecoords * rG0) are pararallel tranported so that the first column equals p1
veca1 ull_S2S_constV(mata1 y, mata1 xs, mata1 xe, mnlink_Omega_cpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 rG0, matd referencecoords){
  int p = om.p1.size();
  //check that ncol(y) == p
  if (y.cols() != p){Rcpp::stop("width of y does not equal length of p1");}
  //check that referencecoords are orthonormal (wont be taped as uses purely matd objects)
  if ((referencecoords.transpose() * referencecoords - matd::Identity(referencecoords.rows(), referencecoords.rows())).norm() > 1E-8){
    Rcpp::stop("referencecoords columns are not an orthonormal basis");
  }

  // project Omega matrix to satisfy orthogonality to p1 and q1
  mnlink_Omega_cpp<a1type> om_projected = Omega_proj_cpp(om);
  veca1 omvec_projected = mnlink_Omega_cpp_vec(om_projected);

  //get mean
  mata1 ypred;
  ypred = mnlink_cpp(xs, xe, omvec_projected, p); 

  //evaluate SvMF density of each observation
  veca1 ld(y.rows());
  //rG0 is provided in coordinate system given by referencecoords
  //to get G0 in the canonical reference system of p1 etc, use referencecoords * G0
  //then to parallel transport G0 to p1, use JuppRmat
  mata1 G0 = referencecoords.cast<a1type>() * rG0;
  mata1 G0star = JuppRmat(G0.col(0), om_projected.p1) * G0.rightCols(G0.cols() - 1);
  mata1 G(p, p);
  veca1 a(p);
  a(0) = a1;
  if (aremaining.size() != p - 1){ Rcpp::stop("aremaining must have length p - 1."); }
  a.segment(1, p-1) = aremaining;
  for (int i = 0; i < y.rows(); ++i){
    G.col(0) = ypred.row(i);
    G.block(0, 1, p, p-1) = JuppRmat(om_projected.p1, ypred.row(i)) * G0star;
    ld(i) = uldSvMF_cann(y.row(i), k, a, G)(0);
  }
  return ld;
}



veca1 ull_S2S_constV_forR(mata1 y, mata1 xs, mata1 xe, veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 rG0, matd referencecoords){
   mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, y.cols(), xe.cols());
   veca1 ld = ull_S2S_constV(y, xs, xe, om, k, a1, aremaining, rG0, referencecoords);
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
// G0star are the orientation axes of SvMF in cannonical coordinate (p x p-1 matrix). These axes must be orthogonal to p1. It is converted to a p-1 x p-1 skew-symmetrix matrix
// ideally G0star is close to the referencecoords axes
// [[Rcpp::export]]
veca1 S2S_constV_nota1_tovecparams(veca1 & omvec, a1type k, veca1 aremaining, mata1 G0star, matd referencecoords){
  int p = aremaining.size() + 1;

  //check G0star
  if (G0star.cols() != p-1){Rcpp::stop("G0star should have p-1 columns");}
  if (G0star.rows() != p){Rcpp::stop("G0star should have p columns");}
  if ((omvec.segment(0, p).transpose() * G0star).norm() > 1e-8){Rcpp::stop("G0star columns should be orthogonal to p1.");}

  //convert G0star to referencecoords
  G0star = referencecoords.cast<a1type>().transpose() * G0star;
  Rcpp::Rcout << "rG0star:" << std::endl << G0star << std::endl;
  //parallel transport along p1 to referencecoords[,1] so that first row of G0star is zeros
  G0star = JuppRmat(referencecoords.transpose() * omvec.segment(0,p), veca1::Unit(p,0)) * G0star;
  Rcpp::Rcout << "rrG0star:" << std::endl << G0star << std::endl;
  //drop first row of zeros
  mata1 Kstar(p-1,p-1);
  Kstar = G0star.bottomRows(p-1);
  Rcpp::Rcout << "Kstar:" << std::endl << Kstar << std::endl;
  a1type detKstar = Kstar.determinant();
  if (std::abs(CppAD::Value(detKstar) + 1.0) < 1e-8) {Rcpp::stop("Kstar has a determinant very close to -1, please change the sign of one of G0star's columns");}
  veca1 vecCayaxes = vectorizeLowerTriangle(inverseCayleyTransform(Kstar)); 

  // put everything into a vector
  veca1 result(omvec.size() + 1 + aremaining.size() - 1 + vecCayaxes.size());
  result.segment(0, omvec.size()) = omvec;
  result(omvec.size()) = k;
  result.segment(omvec.size() + 1, aremaining.size() - 1) = aremaining.tail(aremaining.size() - 1).array().log(); //log the a
  result.segment(omvec.size() + 1 + aremaining.size() - 1, vecCayaxes.size()) = vecCayaxes;
  return result;
}

// reverse function
// without ce, would only need q = qs + qe to be passed in
std::tuple<veca1, a1type, veca1, mata1> S2S_constV_nota1_fromvecparams(const veca1 & mainvec, int p, int qs, int qe) {
  if (mainvec.size() != Omega_veclength(p, qs, qe) + 1 + (p-2) + ((p-2) * (p-1) / 2)  ) {
    Rcpp::stop("Input vector size does not match expected dimensions.");
  }
  
  veca1 omvec = mainvec.segment(0,  Omega_veclength(p, qs, qe));
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


//G0star are the orientation axes of SvMF in cannonical coordinate (p x p-1) matrix
pADFun tape_ull_S2S_constV_nota1(veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0star, vecd & p_in, vecd & qe_in, matd & yx, matd referencecoords){
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
  veca1 mainvec = S2S_constV_nota1_tovecparams(omvec, k, aremaining, G0star, referencecoords);
  veca1 a1vec(1);
  a1vec(0) = a1;

  // tape with main vector and a1 as a dynamic
  CppAD::Independent(mainvec, a1vec);
  // split mainvec into parts, overwriting passed arguments
  auto result = S2S_constV_nota1_fromvecparams(mainvec, p, qs, qe);
  omvec = std::get<0>(result);
  k = std::get<1>(result);
  aremaining = std::get<2>(result);
  mata1 Kstar = std::get<3>(result);
  
  mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, p, qe);

  //identify G0 with p1
  mata1 rG0(p,p);
  rG0.col(0) = referencecoords.transpose() * om.p1; //rG0 because it is G0 expressed in the referencecoords coordinate system
  //Kstar should be rG0[,-1] parallel transported so that rG0[,1] = nthpole, then with the first row dropped
  //rG0[,-1] back add a row of zeros, then use reverse parallel transport to get full axes
  mata1 rrG0star(p,p-1);
  rrG0star << mata1::Zero(1,p-1),Kstar;
  rG0.block(0, 1, p, p-1) = JuppRmat(rG0.col(0),veca1::Unit(p,0)).transpose() * rrG0star;

  veca1 ld = ull_S2S_constV(y, xs, xe, om, k, a1, aremaining, rG0, referencecoords);

  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(mainvec, ld);
  tape.check_for_nan(false);

  pADFun out(tape, mainvec, a1vec, "ull_S2S_constV_nota1");
  return(out);
}


