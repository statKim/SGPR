#ifndef sgp_header_H
#define sgp_header_H

double mcp(double theta, double l, double gamma);

double scad(double theta, double l, double gamma);

double ep(double theta, double l, double tau);

double lasso(double theta, double l);

double approx(double l_v, double l_g, double g_tau, double v_tau, double g_gamma, double v_gamma,
              double old_b, double z_v, double g_norm, double g_norm_active,
              int method_v, int method_g);

double get_loss(Rcpp::DoubleVector r, int n);

double xty(Rcpp::DoubleVector x, Rcpp::DoubleVector y, int n, int v);

double get_norm(Rcpp::DoubleVector X);

#endif
