#include <Rcpp.h>
using namespace Rcpp;
#include <Rcpp.h>
using namespace Rcpp;

//' C++ Kalman Smoother Engine
 //'
 //' Internal function for fast smoothing.
 //' @param y Numeric Vector of data
 //' @param degree Integer, 1 or 2
 //' @param H Double
 //' @param Q Numeric Vector
 //' @param init_state Double
 //' @param init_var Double
 //' @export
 // [[Rcpp::export]]
 NumericVector kalman_smoother_cpp(NumericVector y, int degree, double H, NumericVector Q, double init_state, double init_var) {
   int n = y.size();
   int m = degree;

   // -- Init Matrices --
   NumericMatrix a(n + 1, m);          // State Estimate
   NumericMatrix P_mat(n + 1, m * m);  // State Variance (flattened)

   std::fill(a.begin(), a.end(), 0.0);
   std::fill(P_mat.begin(), P_mat.end(), 0.0);

   // Set initial state
   a(0, 0) = init_state;
   P_mat(0, 0) = init_var;
   if (m == 2) P_mat(0, 3) = init_var / 10.0;

   NumericVector v(n);     // Prediction Error
   NumericVector F(n);     // Prediction Error Variance
   NumericMatrix K(n, m);  // Kalman Gain

   // --- 1. Forward Filter ---
   for (int t = 0; t < n; t++) {
     // Current State
     double a0 = a(t, 0);
     double a1 = (m == 2) ? a(t, 1) : 0.0;

     // Prediction (Time Update)
     double a_pred_0 = a0 + a1;
     double a_pred_1 = a1;

     double p00 = P_mat(t, 0);
     double p01 = (m == 2) ? P_mat(t, 1) : 0.0;
     double p10 = (m == 2) ? P_mat(t, 2) : 0.0;
     double p11 = (m == 2) ? P_mat(t, 3) : 0.0;

     double P_pred_00, P_pred_01, P_pred_10, P_pred_11;

     if (m == 1) {
       P_pred_00 = p00 + Q[0];
     } else {
       // T = [1 1; 0 1]
       P_pred_00 = p00 + p01 + p10 + p11 + Q[0];
       P_pred_01 = p01 + p11;
       P_pred_10 = p10 + p11;
       P_pred_11 = p11 + Q[1];
     }

     // Measurement Update
     if (NumericVector::is_na(y[t])) {
       // Missing Data: Propagate prediction
       a(t + 1, 0) = a_pred_0;
       P_mat(t + 1, 0) = P_pred_00;
       if (m == 2) {
         a(t + 1, 1) = a_pred_1;
         P_mat(t + 1, 1) = P_pred_01; P_mat(t + 1, 2) = P_pred_10; P_mat(t + 1, 3) = P_pred_11;
       }
       v[t] = 0.0; F[t] = 1.0; // Placeholders
     } else {
       // Standard Update
       v[t] = y[t] - a_pred_0;
       F[t] = P_pred_00 + H;
       if (F[t] < 1e-9) F[t] = 1e-9; // Stability Clip

       double invF = 1.0 / F[t];
       double K0 = P_pred_00 * invF;
       double K1 = (m == 2) ? P_pred_10 * invF : 0.0;

       a(t + 1, 0) = a_pred_0 + K0 * v[t];
       // P = (I - KH)P
       P_mat(t + 1, 0) = P_pred_00 * (1.0 - K0);

       if (m == 2) {
         a(t + 1, 1) = a_pred_1 + K1 * v[t];
         P_mat(t + 1, 1) = P_pred_01 * (1.0 - K0);
         P_mat(t + 1, 2) = -K1 * P_pred_00 + P_pred_10;
         P_mat(t + 1, 3) = -K1 * P_pred_01 + P_pred_11;
       }
       K(t, 0) = K0; if (m == 2) K(t, 1) = K1;
     }
   }

   // --- 2. Backward Smoother (RTS) ---
   NumericVector smoothed(n);
   NumericMatrix r(n, m);
   std::fill(r.begin(), r.end(), 0.0);

   for (int t = n - 1; t >= 0; t--) {
     double r0 = (t == n - 1) ? 0.0 : r(t, 0);
     double r1 = (t == n - 1) ? 0.0 : ((m==2) ? r(t, 1) : 0.0);
     double r0_new, r1_new;

     if (NumericVector::is_na(y[t])) {
       // L = T transpose
       r0_new = r0;
       r1_new = (m == 1) ? 0.0 : (r0 + r1);
     } else {
       // L = T - K Z
       double invF = 1.0 / F[t];
       double ZTv = v[t] * invF;

       if (m == 1) {
         r0_new = ZTv + (1.0 - K(t,0)) * r0;
       } else {
         r0_new = ZTv + (1.0 - K(t,0)) * r0 + (-K(t,1)) * r1;
         r1_new =       1.0 * r0 + 1.0 * r1;
       }
     }

     if (t > 0) {
       r(t-1, 0) = r0_new;
       if (m == 2) r(t-1, 1) = r1_new;
     }

     double sm = a(t, 0) + P_mat(t, 0) * r0_new;
     if (m == 2) sm += P_mat(t, 1) * r1_new;
     smoothed[t] = sm;
   }

   return smoothed;
 }
