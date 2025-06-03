#include "S2S_constV.h"
#include "mnlink_cpp.h"
#include "uldSvMF.h"
#include "utils.h"

//G0 specifies the axes of the SvMF.
//G0star is pararallel transported along G01 -> predicted mean
veca1 ull_S2S_constV(mata1 y, mata1 xs, mata1 xe, mnlink_Omega_cpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 G0){
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
  mata1 G0star = G0.rightCols(G0.cols() - 1);
  mata1 G(p, p);
  veca1 a(p);
  a(0) = a1;
  if (aremaining.size() != p - 1){ Rcpp::stop("aremaining must have length p - 1."); }
  a.segment(1, p-1) = aremaining;
  for (int i = 0; i < y.rows(); ++i){
    G.col(0) = ypred.row(i);
    G.block(0, 1, p, p-1) = JuppRmat(G0.col(0), ypred.row(i)) * G0star;
    ld(i) = uldSvMF_cann(y.row(i), k, a, G)(0);
  }
  return ld;
}



veca1 ull_S2S_constV_forR(mata1 y, mata1 xs, mata1 xe, veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0){
   mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, y.cols(), xe.cols());
   veca1 ld = ull_S2S_constV(y, xs, xe, om, k, a1, aremaining, G0);
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
//' @param G0 are the orientation axes of SvMF in cannonical coordinate (p x p matrix). Ideally G0 is close to the referencecoords axes. G0 must be a rotation matrix (det > 0) so that the Cayley transform representation works.
//' @param referencecoords is a p x p orthonormal matrix specifying the reference coordinates for the Cayley transforms. It is best if referencecoords is close to the best G0 (so rG0 is close the identity) and it will fail if `G01` is the antepode of `referencoords[,1]`.
// [[Rcpp::export]]
veca1 S2S_constV_nota1_tovecparams(veca1 & omvec, a1type k, veca1 aremaining, mata1 G0, matd referencecoords, std::string G01behaviour){
  int p = aremaining.size() + 1;
  // check G01behaviour
  if ((G01behaviour != "p1") && (G01behaviour != "fixed") && (G01behaviour != "free")){Rcpp::stop("G01behaviour not understood");}

  //check G0
  if (G0.cols() != p){Rcpp::stop("G0 should have p columns");}
  if (G0.rows() != p){Rcpp::stop("G0 should have p columns");}
  if ((G0.transpose() * G0 - mata1::Identity(p,p)).norm() > 1e-8){Rcpp::stop("G0 columns should be orthonormal.");}
  if (CppAD::Value(G0.determinant()) < 0.) {Rcpp::stop("G0 has a negative determinant, please change the sign of one of G0's columns so that it is a rotation");}
  if (G01behaviour == "p1"){ if ((G0.col(0) - omvec.segment(0,p)).norm() > 1e-8){Rcpp::stop("G01 should be equal to p1");} }
  
  //check that referencecoords are orthonormal (wont be taped as uses purely matd objects)
  if ((referencecoords.transpose() * referencecoords - matd::Identity(referencecoords.rows(), referencecoords.rows())).norm() > 1E-8){
    Rcpp::stop("referencecoords columns are not an orthonormal basis");
  }

  //convert G0 to referencecoords
  G0 = referencecoords.cast<a1type>().transpose() * G0;
  mata1 rotmat; //rot matrix to represent using a Cayley transform
  //if G01 is fixed or p1 then get a p-1 x p-1 matrix reprensenting the remaining free columns
  if ((G01behaviour == "p1") || (G01behaviour == "fixed")){
    //parallel transport along p1 to referencecoords[,1] so that first row of G0star is zeros
    mata1 G0star = JuppRmat(G0.col(0), veca1::Unit(p,0)) * G0.rightCols(p-1);
    //drop first row of zeros
    rotmat = G0star.bottomRows(p-1);
  } else if (G01behaviour == "free") {
    rotmat = G0;
  }
  if (CppAD::Value(rotmat.determinant()) < 0.) {Rcpp::stop("Determinant is negative. Please change the sign of one of G0's columns.");}
  veca1 vecCayaxes = vectorizeLowerTriangle(inverseCayleyTransform(rotmat)); 

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
//' @param G01 only for the G01behaviour == "fixed" situation is G01 needed to recover full parameter set
std::tuple<veca1, a1type, veca1, mata1> S2S_constV_nota1_fromvecparams(const veca1 & mainvec, int p, int qs, int qe, matd referencecoords, std::string G01behaviour, vecd G01 = vecd(0)) {
  // check G01behaviour
  if ((G01behaviour != "p1") && (G01behaviour != "fixed") && (G01behaviour != "free")){Rcpp::stop("G01behaviour not understood");}

  //check length
  int vecCay_length = 0;
  if (G01behaviour == "free"){vecCay_length = ((p-1) * p / 2);}
  else {vecCay_length = ((p-2) * (p-1) / 2);}
  if (mainvec.size() != Omega_veclength(p, qs, qe) + 1 + (p-2) + vecCay_length) {
    Rcpp::stop("Input vector size does not match expected dimensions.");
  }
  
  veca1 omvec = mainvec.segment(0,  Omega_veclength(p, qs, qe));
  a1type k = mainvec(omvec.size());
  veca1 laremaining_m1(p - 2); //convert log remaining to full aremaining
  laremaining_m1 = mainvec.segment(omvec.size() + 1, p - 2);
  veca1 aremaining(p-1);
  aremaining[0] = CppAD::exp(-laremaining_m1.sum());
  aremaining.tail(aremaining.size() - 1) = laremaining_m1.array().exp();

  //get back rotmat from vecCay
  veca1 vecCayaxes = mainvec.tail(vecCay_length);
  mata1 rotmat = cayleyTransform(inverseVectorizeLowerTriangle(vecCayaxes));

  // get back G0
  mata1 G0;
  //if G01 is fixed or p1 then get a p-1 x p-1 matrix reprensenting the remaining free columns
  if (G01behaviour == "p1"){G0.col(0) = referencecoords.transpose() * omvec.segment(0,p);}
  if (G01behaviour == "fixed"){G0.col(0) = referencecoords.transpose() * G01;}
  if ((G01behaviour == "p1") || (G01behaviour == "fixed")){
    mata1 G0star = mata1::Zero(p, p-1);
    G0star.bottomRows(p-1) = rotmat;
    //undo: parallel transport along p1 to referencecoords[,1] so that first row of G0star is zeros
    G0star = JuppRmat(G0.col(0), veca1::Unit(p,0)).transpose() * G0star;
    G0.rightCols(p-1) = G0star;
  } else if (G01behaviour == "free") {
    G0 = rotmat;
  }
  //undo: convert G0 to referencecoords
  G0 = referencecoords.cast<a1type>() * G0;

  return std::make_tuple(omvec, k, aremaining, G0);
}

//export the reverse function
// [[Rcpp::export]]
Rcpp::List S2S_constV_nota1_fromvecparamsR(const veca1 & mainvec, int p, int qs, int qe, matd referencecoords, std::string G01behaviour, Rcpp::Nullable<vecd> G01 = R_NilValue) {
  vecd G01_;
  if (G01.isNotNull()){G01_ = Rcpp::as<vecd>(G01);}
  else {G01_ = vecd();}

  auto result = S2S_constV_nota1_fromvecparams(mainvec, p, qs, qe, referencecoords, G01behaviour, G01_);
  veca1 omvec = std::get<0>(result);
  a1type k = std::get<1>(result);
  veca1 aremaining = std::get<2>(result);
  mata1 G0 = std::get<3>(result);
  
  // Return as a list
  return Rcpp::List::create(
    Rcpp::Named("omvec") = omvec,
    Rcpp::Named("k") = k,
    Rcpp::Named("aremaining") = aremaining,
    Rcpp::Named("G0") = G0
  );
}


//G0 are the orientation axes of SvMF in cannonical coordinate (p x p) matrix
pADFun tape_ull_S2S_constV_nota1(veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0, vecd & p_in, vecd & qe_in, matd & yx, matd referencecoords, std::string G01behaviour){
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

  //for G01 fixed prepare G01 from input G0
  vecd tapeG01(G0.rows());
  for (int i = 0; i < G0.rows(); ++i) {
    tapeG01(i) = CppAD::Value(G0(i, 0));
  }


  // Get all parameters except a1 into a vector
  veca1 mainvec = S2S_constV_nota1_tovecparams(omvec, k, aremaining, G0, referencecoords, G01behaviour);
  veca1 a1vec(1);
  a1vec(0) = a1;

  // tape with main vector and a1 as a dynamic
  CppAD::Independent(mainvec, a1vec);
  // split mainvec into parts, overwriting passed arguments
  auto result = S2S_constV_nota1_fromvecparams(mainvec, p, qs, qe, referencecoords, G01behaviour, tapeG01);//final argument only used if G01behaviour == "fixed"
  omvec = std::get<0>(result);
  k = std::get<1>(result);
  aremaining = std::get<2>(result);
  G0 = std::get<3>(result);
  
  mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, p, qe);

  veca1 ld = ull_S2S_constV(y, xs, xe, om, k, a1, aremaining, G0);

  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(mainvec, ld);
  tape.check_for_nan(false);

  pADFun out(tape, mainvec, a1vec, "ull_S2S_constV_nota1");
  return(out);
}


