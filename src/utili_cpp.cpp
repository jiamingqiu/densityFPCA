#include <Rcpp.h>
using namespace Rcpp;

// C++ funtion computing L2 inner product using trapezoidal rule via Kahan sum
// grid must be sorted
// [[Rcpp::export]]
double cpp_dotL2(NumericVector f1,
                  NumericVector f2,
                  NumericVector grid) {
  double res = 0;
  int n = grid.size();

  NumericVector tm(n);
  double tmSum = 0;
  double tmAdd = 0;
  double compens = 0;

  tm[0] = f1[0] * f2[0];
  // using Kahan Sum
  for(int i = 0; i < n - 1; i++){
    tm[i+1] = f1[i+1] * f2[i+1];
    tmAdd = (tm[i] + tm[i+1]) * (grid[i+1] - grid[i]) - compens;
    tmSum = res + tmAdd;
    compens = (tmSum - res) - tmAdd;
    res = tmSum;
  }
  return res/2;
}

// This is not very helpful in terms of speeding things up, rewrite semiDscWass
// C++ funtion computing cdf with pdf using trapezoidal rule via "Kahan" sum
// grid must be sorted
// [[Rcpp::export]]
NumericVector cpp_calCdf(NumericVector pdf, NumericVector grid){
  int n = grid.size();
  NumericVector res(n);

  double tmAdd = 0;
  double compens = 0;

  res[0] = 0;
  for(int i = 1; i < n; i++){
    tmAdd = 0.5 * (pdf[i-1] + pdf[i]) * (grid[i] - grid[i-1]) - compens;
    res[i] = res[i-1] + tmAdd;
    compens = (res[i] - res[i-1]) - tmAdd;
  }
  // for(int i = 0; i < n; i++){
  //   res[i] = res[i] / res[n-1];
  // }
  return res;
}


// C++ funtion computing wasserstein distance
// cdf1 and cdf2 are cdf, non-negative, s.t.
// on [gird[i], grid[i+1]), cdf1(x) = cdf1[i]
// cdfs must starts from 0
// grid must be sorted
// [[Rcpp::export]]
double cpp_distWass(NumericVector cdf1,
                    NumericVector cdf2,
                    NumericVector grid,
                    double p) {


  // getting ready for wass dist
  //re-unify
  // int n = grid.size();
  // NumericVector wcdf1(n+1);
  // NumericVector wcdf2(n+1);
  // NumericVector wgrid(n+1);
  // // padding so that cdf always starts with 0
  // wcdf1[0] = 0;
  // wcdf2[0] = 0;
  // wgrid[0] = grid[0] - 1;
  // for(int i = 1; i < n; i++){
  //   wcdf1[i] = cdf1[i-1] / cdf1[n-1];
  //   wcdf2[i] = cdf2[i-1] / cdf2[n-1];
  //   wgrid[i] = grid[i-1];
  // }
  // wcdf1[n] = 1;
  // wcdf2[n] = 1;

  double from = 0;
  double to = 0;
  int idx1 = 0;
  int idx2 = 0;
  double res = 0;

  double tmAdd = 0;
  double compens = 0;
  double tmSum = 0;

  // double diff = 0;
  // std::cout << "idx1\tidx2\tcdf1[idx1+1]\tcdf2[idx2+1]\tdiff\tfrom\tto\tres\n";
  // waste my time: use std::abs for primitive type, and Rcpp::abs, which is default
  // only works for *Vector, *Matrix e.t.c., c.f.
  // https://stackoverflow.com/questions/45713089/rcpp-no-matching-function-for-call-to-abs

  // until from = 1
  while(from < 1.0){
    // declare: cdf1 will jump next
    if(cdf1[idx1 + 1] < cdf2[idx2 + 1]){
      to = cdf1[idx1 + 1];

      // diff = std::abs(grid[idx1+1] - grid[idx2+1]);
      // res = res + (to - from) * std::pow(diff, p);
      // std::cout << idx1<<"\t" << idx2 << "\t"<<
      //   cdf1[idx1+1]<<"\t"<<cdf2[idx2+1]<< "\t"<< diff << "\t " << from << "\t" << to << "\t" << res << "\n";

      // balance during (from, to]
      tmAdd = (to - from) * std::pow(std::abs(grid[idx1+1] - grid[idx2+1]), p) - compens;
      tmSum = res + tmAdd;
      compens = (tmSum - res) - tmAdd;
      res = tmSum;

      // reset duration
      from = to;
      // cdf1 go up 1 step
      idx1 = idx1 + 1;
    }
    // declare: cdf2 will jump next
    if(cdf1[idx1 + 1] > cdf2[idx2 + 1]){
      to = cdf2[idx2 + 1];
      // diff = std::abs(grid[idx1+1] - grid[idx2+1]);
      // res = res + (to - from) * std::pow(diff, p);
      // std::cout << idx1<<"\t" << idx2 << "\t"<<
      //   cdf1[idx1+1]<<"\t"<<cdf2[idx2+1]<< "\t"<< diff << "\t " << from << "\t" << to << "\t" << res << "\n";

      tmAdd = (to - from) * std::pow(std::abs(grid[idx1+1] - grid[idx2+1]), p) - compens;
      tmSum = res + tmAdd;
      compens = (tmSum - res) - tmAdd;
      res = tmSum;

      from = to;
      idx2 = idx2 + 1;
    }
    // declare: both will jump next
    if(cdf1[idx1 + 1] == cdf2[idx2 + 1]){
      to = cdf2[idx2 + 1];
      // diff = std::abs(grid[idx1+1] - grid[idx2+1]);
      // res= res + (to - from) * std::pow(diff, p);
      // std::cout << idx1<<"\t" << idx2 << "\t"<<
      //   cdf1[idx1+1]<<"\t"<<cdf2[idx2+1]<< "\t"<< diff << "\t " << from << "\t" << to << "\t" << res << "\n";

      tmAdd = (to - from) * std::pow(std::abs(grid[idx1+1] - grid[idx2+1]), p) - compens;
      tmSum = res + tmAdd;
      compens = (tmSum - res) - tmAdd;
      res = tmSum;

      from = to;
      idx1 = idx1 + 1;
      idx2 = idx2 + 1;
    }
  }
  // remark: no need to worry ties, i.e., flat cdf1 or cdf2
  // if flat, to - from = 0

  return std::pow(res, 1/p);
}

// // Just a test for some basic math functions
// double cpp_pow(double x, double pwr){
//   return pow(x, pwr);
// }

// C++ funtion computing ISE of 2 array of function values,
// values may be evaluated at different grids
// using trapezoidal rule via Kahan sum
// grid must be sorted
// Do NOT FORGET [[Rcpp::export]]
// double cpp_ISE(NumericVector f1,
//                  NumericVector f2,
//                  NumericVector g1,
//                  NumericVector g2) {
//   double res = 0;
//   int n1 = g1.size();
//   int n2 = g2.size();
//
//   NumericVector cmGrid(n1 + n2);
//   NumericVector tm(n1 + n2);
//   int p1 = 0, p2 = 0, pcm = 0;
//
//   double tmSum = 0;
//   double tmAdd = 0;
//   double compens = 0;
//
//   for(int i = 0; i< n1 + n2; i++){
//     if(g1[p1] < g2[p2]){
//       cmGrid[pcm] = g1[p1];
//       pcm++;
//       p1++;
//     }else{
//       cmGrid[pcm] = g2[p2];
//       pcm++;
//       p2++;
//     }
//   }
//
//   tm[0] = f1[0] - f2[0];
//   // using Kahan Sum
//   for(int i = 0; i < n - 1; i++){
//     tm[i+1] = f1[i+1] * f2[i+1];
//     tmAdd = (tm[i] + tm[i+1]) * (grid[i+1] - grid[i]) - compens;
//     tmSum = res + tmAdd;
//     compens = (tmSum - res) - tmAdd;
//     res = tmSum;
//   }
//   return res/2;
// }
