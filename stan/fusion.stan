// Spatiotemporal fusion of monthly NO2 and PM2.5 in Fortaleza. One program for
// every model of the paper; the models differ only by the data flags below.
// The marginal likelihood is computed with a 4-state Kalman filter whose static
// effects (PM2.5 deviations at the reference locations, low-cost device
// effects) are regression effects of an augmented filter, integrated in closed
// form at the end (Durbin and Koopman 2012, section 6.2.2; de Jong 1991).
//
// Latent field (log scale, month t = 1..T from 2020-01, pollutant p in {N, P}):
//   eta_p(s, t) = b_p + u_p(s) + v_p(t)
//   u_p(s)      = X_p(s)' theta_p + log R_p(s) + n_p(s)     n_P only, at the reference locations
//   v_p(t)      = Z_p[t] beta_p + Hm[t] h_p + e_p[t]
// Spatial profile:
//   a_p = U zc_p sd_p,  r_kp = exp(a_kp) / sum_j wbar_j exp(a_jp),  R_p(s) = sum_k w_k(s) r_kp,
//   so sum_k wbar_k r_kp = 1 and the SatPM2.5 municipal mean of
//   C_P(., t) = exp(b_P + v_P(t)) R_P(.) is exactly exp(b_P + v_P(t)).
//   zc_N = z_N; zc_P = rho_S z_N + sqrt(1 - rho_S^2) z_P (profile correlation, couple_space);
//   z_N, z_P ~ N(0, I); a pollutant without profile (field_p = 0) has log R_p = 0.
// Temporal components: e_t = diag(rho_N, rho_P) e_(t-1) + w_t, w_t ~ N(0, Qw),
//   Qw = [[iN^2, r iN iP], [r iN iP, iP^2]], iN = s_N sqrt(1 - rho_N^2), e_1 stationary
//   with c0 = r iN iP / (1 - rho_N rho_P) (coupled = 0 sets the innovation correlation r = 0).
// Observation rows (ly = log y), one per month and source:
//   NO2 background (cls 1):  b_N + u_N(s) + v_N(t) + eps,                eps ~ N(0, sig_N^2)
//   NO2 traffic (cls 2, 3):  b_N + u_N(s) + g[c] + kappa v_N(t) + eps,   eps ~ N(0, sig_NT^2)
//       amplitude kappa ~ U(0, 1) on the whole v_N(t) (use_kappa = 0: kappa = 0, static traffic);
//       n_g = 2 gives avenue and terminal their own offsets, n_g = 1 one shared offset.
//   Reference PM2.5 (type 2): b_P + u_P(s) + v_P(t) + eps_P,             eps_P ~ N(0, sig_P^2)
//   SatPM2.5 (type 3):       d_A + H1[t] gamma_A + b_P + b_A v_P(t) + eps_A,  b_A = 1 unless free_bA
//   Low-cost (type 4, device d): a_M + b_M (b_P + u_P(s_d) + v_P(t)) + v_d + m_t + eps_M
//       v_d ~ N(0, tau_M^2) static per device, m_t ~ N(0, tau_m^2) per month;
//       the deviation of a reference location enters a low-cost row there scaled by b_M.
//   Calibration pair k:      log moqa_k = a_M + v_C + b_M log bam_k + eps_C, eps_C ~ N(0, sig_C^2),
//       v_C ~ N(0, tau_M^2) shared by the n_C pairs, in closed form (Sherman-Morrison).
//   TROPOMI (use_column; cell i, month t; natural signed scale, unit 1e-5 mol/m2):
//     Y_it = alpha_T + H1[t] gamma_T + lambda_S H_i[h_S] + beta_cell[i]
//            + (lambda_T / s_N) e_N[t] + q_t + eps_it,       eps_it ~ N(0, d_it^2 / kappa_T)
//     H_i[h_S] = B_col z_N (cell mean of the standardized NO2 contrast, without sd_N;
//       absent without the NO2 profile), beta_cell = U_cell b_cell (sum zero),
//       b_cell ~ N(0, sigma_beta^2), q_t a stationary AR(1) (sigma_q, rho_q) in state 4.
//     col_temporal = 0 removes the temporal loading lambda_T (the final model);
//     sign_constraint = 1 restricts the spatial loading to lambda_S > 0, which
//     identifies the sign of the NO2 profile.
//     The month's n_t rows enter the filter collapsed into one pseudo-row (Durbin and
//     Koopman 2012, section 6.5; Jungbacker and Koopman 2015): with r_it = Y_it - mu_it,
//     w_it = 1 / d_it^2, S_t = sum_i w_it and rbar_t = sum_i w_it r_it / S_t,
//       sum_i log N(r_it | c_t, d_it^2 / kappa_T) = log N(rbar_t | c_t, 1 / (kappa_T S_t)) + R_t,
//       R_t = -((n_t - 1) / 2) log(2 pi) + ((n_t - 1) / 2) log kappa_T - sum_i log d_it
//             - (1/2) log S_t - (kappa_T / 2) sum_i w_it (r_it - rbar_t)^2,
//     where c_t = (lambda_T / s_N) e_N[t] + q_t is the one state combination of the month.
//
// Inference. The core state x_t = (e_N, e_P, m, q) has 4 dimensions. The static
// effects u = (n_P at the Jn training reference locations, v_d at the D training
// devices) ~ N(0, Dg), Dg = diag(tau_P^2 ..., tau_M^2 ...), enter row i as g_i' u
// (g_i: 1 for a reference PM2.5 row at its location, 1 for a low-cost row at its
// device, b_M for a low-cost row at a reference location). Given u the innovations
// are affine in u: the filter is run on the augmented observation rows
// [y_i - c_i, g_i'] with state means A (4 x (M + 1)), so that v_i(u) = v_i^0 - G_i' u
// with F_i free of u. With W_i = (v_i^0, G_i') / sqrt(F_i), q = sum_i (v_i^0)^2 / F_i,
// s = sum_i G_i v_i^0 / F_i and S = sum_i G_i G_i' / F_i (the blocks of W'W),
//   log p(y) = -(1/2) sum_i log(2 pi F_i) - q/2 + s' Aq^-1 s / 2 - log|Aq| / 2 - log|Dg| / 2,
//   Aq = S + Dg^-1
// (Durbin and Koopman 2012, section 6.2.2; de Jong 1991; Rasmussen and Williams
// 2006, A.2). Rows are processed one at a time in month order (univariate
// treatment, Durbin and Koopman 2012, sections 4.3, 6.4, 7.2). NUTS samples only
// the fixed effects, the hyperparameters, the profile coordinates and the
// n_cells - 1 cell contrasts.
//
// Priors (paper, Table 2; Supplementary Section S2 for the full model):
//   b_N, b_P ~ N(log 10, 1)        g ~ N(log 1.8, 0.35) each
//   beta, h, theta ~ N(0, 0.5)     d_A ~ N(0, 0.3)   gamma_A ~ N(0, 0.3)
//   b_A ~ N(1, 0.5) on (0.5, 2.5)  a_M ~ N(0, 1)     b_M ~ N(1, 0.3) on [0.5, 2]
//   alpha_T ~ N(0, 3)  gamma_T ~ N(0, 0.3)  lambda_S, lambda_T ~ N(0, 0.5)  log kappa_T ~ N(0, 1)
//   rho_N, rho_P, kappa ~ U(0, 1)
//   (r + 1) / 2, (rho_S + 1) / 2, (rho_q + 1) / 2 ~ Beta(2, 2)
//   every scale (s, tau, sig, sd, sigma_q, sigma_beta) ~ half exponential-power with k = 4,
//   p(sigma) = 4 / (a Gamma(1/4)) exp(-(sigma / a)^4), a = U / qgamma(0.05, 1/4, upper)^(1/4),
//   with U per parameter in prior_scale (quartic = 1), or half-N(0, a) (quartic = 0).
//   The `~` statements drop their Normal constants (b_cell keeps -log sigma_beta);
//   target += terms keep theirs.
functions {
  // Half exponential-power prior of a scale:
  //   log f(s | a) = log 4 - log a - lgamma(1/4) - (s / a)^4 on s >= 0,
  // or the half-Normal with SD a when quartic = 0, both with their constants.
  real scale_lpdf(real s, real a, int quartic) {
    if (quartic) return log(4) - log(a) - lgamma(0.25) - pow(s / a, 4);
    return log(2) + normal_lpdf(s | 0, a);
  }
  // log r_k of the normalized positive profile with coordinates z and SD sd:
  // a = U z sd, log r_k = a_k - log sum_j wbar_j exp(a_j), so sum_k wbar_k r_k = 1.
  vector log_profile(matrix U, vector wbar, vector z, real sd) {
    int K = rows(U);
    vector[K] a = U * z * sd;
    a = a - max(a);
    return a - log_sum_exp(log(wbar) + a);
  }
}
data {
  int<lower=2> T;
  int<lower=0> P_N;
  int<lower=0> P_P;
  matrix[T, P_N] Z_N;
  matrix[T, P_P] Z_P;
  int<lower=0> H;
  matrix[T, H] Hm;
  matrix[T, 2] H1;
  // Spatial basis shared by both pollutants: K knots, Helmert contrasts U (K x K-1),
  // municipal reference weights wbar (sum 1) and the support of each training site.
  int<lower=2> K;
  matrix[K, K - 1] U;
  vector<lower=0>[K] wbar;
  int<lower=0, upper=1> field_N;
  int<lower=0, upper=1> field_P;
  int<lower=0, upper=1> couple_space;
  int<lower=1> J;
  matrix<lower=0>[J, K] W_site;
  int<lower=0> Q_N;
  matrix[J, Q_N] Xs_N;
  int<lower=0> Q_P;
  matrix[J, Q_P] Xs_P;
  // PM2.5 deviation: static-effect index of each training reference location (0: none).
  int<lower=0, upper=1> use_nugget;
  int<lower=0> Jn;
  array[J] int<lower=0, upper=Jn> nug;
  int<lower=0, upper=1> coupled;
  int<lower=0, upper=1> use_kappa;
  int<lower=1, upper=2> n_g;
  int<lower=0, upper=1> acag_seasonal;
  int<lower=0, upper=1> free_bA;
  // Rows in month order: type 1 NO2, 2 reference PM2.5, 3 SatPM2.5, 4 low-cost.
  int<lower=1> n;
  vector[n] ly;
  array[n] int<lower=1, upper=4> type;
  array[n] int<lower=0, upper=3> cls;
  array[n] int<lower=0, upper=J> site;
  array[T] int<lower=1> start;
  array[T] int<lower=0> count;
  int<lower=0> D;
  array[n] int<lower=0, upper=D> dev;
  int<lower=0> n_C;
  vector[n_C] lmoqa_C;
  vector[n_C] lbam_C;
  // TROPOMI (R/column.R): rows in month
  // order with per-month start/count (an empty month has count 0); the cell
  // means B_col of the standardized contrast and the cell contrasts U_cell.
  // With use_column = 0 the rows are zero-size and n_cells = 1 is a dummy.
  int<lower=0, upper=1> use_column;
  // 1: TROPOMI loads e_N by lambda_T / s_N; 0: no temporal loading.
  int<lower=0, upper=1> col_temporal;
  // 1: lambda_S > 0; 0: lambda_S unconstrained.
  int<lower=0, upper=1> sign_constraint;
  int<lower=1> n_cells;
  int<lower=2> K_col;
  matrix[n_cells, K_col - 1] B_col;
  matrix[n_cells, n_cells - 1] U_cell;
  int<lower=0> n_col;
  vector[n_col] y_col;
  vector<lower=0>[n_col] d_col;
  array[n_col] int<lower=1, upper=n_cells> cell_col;
  array[n_col] int<lower=1, upper=T> t_col;
  array[T] int<lower=1> start_col;
  array[T] int<lower=0> count_col;
  // Scale priors: quartic tails (1) or half-Normal (0), with the scale a of each:
  // 1 s_N, 2 s_P, 3 tau_P, 4 tau_M, 5 tau_m, 6 sig_N, 7 sig_NT, 8 sig_P, 9 sig_A,
  // 10 sig_M, 11 sig_C, 12 sd_N, 13 sd_P, 14 sigma_q, 15 sigma_beta.
  int<lower=0, upper=1> quartic;
  vector<lower=0>[15] prior_scale;
}
transformed data {
  int has_a = 0;
  int has_m = 0;
  int use_cal;
  int n_rhoS = couple_space * field_N * field_P;
  int has_col = use_column;
  array[3] int g_of = {0, 1, n_g};
  for (i in 1:n) if (type[i] == 3) has_a = 1;
  for (i in 1:n) if (type[i] == 4) has_m = 1;
  use_cal = has_m && n_C > 0;
  if (has_m && n_C == 0) reject("Low-cost rows need calibration pairs (n_C > 0).");
  if (use_nugget == 0 && Jn > 0) reject("Jn must be 0 without the PM2.5 deviation.");
  if (use_column == 0 && n_col > 0) reject("Column rows need use_column = 1.");
  if (use_column == 1 && n_col == 0) reject("use_column = 1 needs column rows.");
  if (has_col && K_col != K) reject("K_col must equal the profile rank K.");
  for (t in 1:T) if (count_col[t] > 0) {
    for (i in start_col[t]:(start_col[t] + count_col[t] - 1)) {
      if (t_col[i] != t) reject("Column rows must be in month order with start_col/count_col.");
    }
  }
  // Per-month sums of the column rows: S_t = sum_i 1 / d_it^2 and sum_i log d_it,
  // and the state-free part of R_t.
  vector[n_col] inv_d2_col = inv_square(d_col);
  vector[n_col] log_d_col = log(d_col);
  vector[T] S_col = rep_vector(0, T);
  vector[T] sum_log_d_col = rep_vector(0, T);
  int n_months_col = 0;
  real col_constant = 0;
  for (t in 1:T) if (count_col[t] > 0) {
    S_col[t] = sum(segment(inv_d2_col, start_col[t], count_col[t]));
    sum_log_d_col[t] = sum(segment(log_d_col, start_col[t], count_col[t]));
    n_months_col += 1;
    col_constant += -0.5 * (count_col[t] - 1) * log(2 * pi()) - sum_log_d_col[t] - 0.5 * log(S_col[t]);
  }
  // Column rows by month and cell: precision weights w_tc = 1 / d_tc^2 (0 when
  // absent), w_tc y_tc and the month sums sum_c w y, sum_c w y^2, so that the
  // month's residual sum of squares is a function of the cell means alone;
  // 1 / S_t where the month has rows, 0 otherwise (never read there).
  matrix[T, n_cells] Wm = rep_matrix(0, T, n_cells);
  matrix[T, n_cells] WY = rep_matrix(0, T, n_cells);
  vector[T] WY1 = rep_vector(0, T);
  vector[T] WY2 = rep_vector(0, T);
  vector[T] inv_S_col = rep_vector(0, T);
  for (i in 1:n_col) {
    Wm[t_col[i], cell_col[i]] = inv_d2_col[i];
    WY[t_col[i], cell_col[i]] = inv_d2_col[i] * y_col[i];
    WY1[t_col[i]] += inv_d2_col[i] * y_col[i];
    WY2[t_col[i]] += inv_d2_col[i] * square(y_col[i]);
  }
  for (t in 1:T) if (count_col[t] > 0) inv_S_col[t] = 1 / S_col[t];
  // Static effects as regression columns: 1..Jn the deviations, Jn+1..Jn+D the devices.
  int M = Jn + has_m * D;
  // Filter rows: the ground rows and one collapsed TROPOMI row per month with rows.
  int n_rows = n + has_col * n_months_col;
  real log2pi = log(2 * pi());
}
parameters {
  real b_N;
  real b_P;
  vector[n_g] g;
  vector[P_N] beta_N;
  vector[P_P] beta_P;
  vector[H] h_N;
  vector[H] h_P;
  vector[Q_N] theta_N;
  vector[Q_P] theta_P;
  array[has_a] real d_A;
  vector[has_a * acag_seasonal * 2] gamma_A;
  array[has_a * free_bA] real<lower=0.5, upper=2.5> b_A;
  real<lower=0> s_N;
  real<lower=0> s_P;
  real<lower=0, upper=1> rho_N;
  real<lower=0, upper=1> rho_P;
  array[coupled] real<lower=-1, upper=1> r;
  array[use_kappa] real<lower=0, upper=1> kappa;
  real<lower=0> sig_N;
  real<lower=0> sig_NT;
  real<lower=0> sig_P;
  array[has_a] real<lower=0> sig_A;
  vector[field_N * (K - 1)] z_N;
  vector[field_P * (K - 1)] z_P;
  array[field_N] real<lower=0> sd_N;
  array[field_P] real<lower=0> sd_P;
  array[n_rhoS] real<lower=-1, upper=1> rho_S;
  array[use_nugget] real<lower=0> tau_P;
  array[has_m] real a_M;
  array[has_m] real<lower=0.5, upper=2> b_M;
  array[has_m] real<lower=0> tau_M;
  array[has_m] real<lower=0> sig_M;
  array[use_cal] real<lower=0> sig_C;
  array[has_m] real<lower=0> tau_m;
  // TROPOMI column: level, seasonal pair, spatial and temporal loadings, the
  // log precision factor, the common AR(1) error and the cell contrasts.
  array[has_col] real alpha_T;
  vector[has_col * 2] gamma_T;
  array[has_col * field_N] real<lower=(sign_constraint ? 0 : negative_infinity())> lambda_S;
  array[has_col * col_temporal] real lambda_T;
  array[has_col] real log_kappa_T;
  array[has_col] real<lower=0> sigma_q;
  array[has_col] real<lower=-1, upper=1> rho_q;
  array[has_col] real<lower=0> sigma_beta;
  vector[has_col * (n_cells - 1)] b_cell;
}
model {
  real k = 0;
  real bA = 1;
  real rr = 0;
  real tm2 = 0;
  real tq2 = 0;
  real rq = 0;
  real kappa_T = 1;
  real lamT = 0;
  real bM = 1;
  // Precision-weighted mean residual of the column rows of each month.
  vector[T] rbar_col = rep_vector(0, T);
  real iN;
  real iP;
  real c0;
  vector[T] f_N = rep_vector(0, T);
  vector[T] f_P = rep_vector(0, T);
  vector[K] lr_N = rep_vector(0, K);
  vector[K] lr_P = rep_vector(0, K);
  vector[J] u_N = rep_vector(0, J);
  vector[J] u_P = rep_vector(0, J);
  vector[4] phi4;
  matrix[4, 4] Q4;
  matrix[4, 4] P;
  // Augmented state means: column 1 the observations, columns 2..M+1 the
  // static-effect regressors; the scaled innovations W and log F of every row.
  matrix[4, M + 1] A = rep_matrix(0, 4, M + 1);
  matrix[n_rows, M + 1] W;
  vector[n_rows] logF;
  int irow = 0;

  b_N ~ normal(log(10), 1);
  b_P ~ normal(log(10), 1);
  // Traffic/background NO2 ratio: UK NO2 network, roadside vs urban background
  // national means 2000-2004, ratio 1.77-1.90.
  g ~ normal(log(1.8), 0.35);
  beta_N ~ normal(0, 0.5);
  beta_P ~ normal(0, 0.5);
  h_N ~ normal(0, 0.5);
  h_P ~ normal(0, 0.5);
  theta_N ~ normal(0, 0.5);
  theta_P ~ normal(0, 0.5);
  target += scale_lpdf(s_N | prior_scale[1], quartic);
  target += scale_lpdf(s_P | prior_scale[2], quartic);
  target += scale_lpdf(sig_N | prior_scale[6], quartic);
  target += scale_lpdf(sig_NT | prior_scale[7], quartic);
  target += scale_lpdf(sig_P | prior_scale[8], quartic);
  if (coupled) target += beta_lpdf((r[1] + 1) / 2 | 2, 2);
  if (has_a) {
    d_A ~ normal(0, 0.3);
    gamma_A ~ normal(0, 0.3);
    target += scale_lpdf(sig_A[1] | prior_scale[9], quartic);
    if (free_bA) b_A ~ normal(1, 0.5);
  }
  if (field_N) {
    z_N ~ std_normal();
    target += scale_lpdf(sd_N[1] | prior_scale[12], quartic);
  }
  if (field_P) {
    z_P ~ std_normal();
    target += scale_lpdf(sd_P[1] | prior_scale[13], quartic);
  }
  if (n_rhoS) target += beta_lpdf((rho_S[1] + 1) / 2 | 2, 2);
  if (use_nugget) target += scale_lpdf(tau_P[1] | prior_scale[3], quartic);
  if (has_m) {
    a_M ~ normal(0, 1);
    b_M ~ normal(1, 0.3);
    target += scale_lpdf(tau_M[1] | prior_scale[4], quartic);
    target += scale_lpdf(sig_M[1] | prior_scale[10], quartic);
    target += scale_lpdf(tau_m[1] | prior_scale[5], quartic);
    tm2 = square(tau_m[1]);
    bM = b_M[1];
    if (use_cal) {
      // Pairs jointly N(a_M + b_M log bam, s2 I + t2 11') with v_C integrated:
      // |S| = s2^(n_C - 1) (s2 + n_C t2), S^-1 = (I - t2 / (s2 + n_C t2) 11') / s2.
      vector[n_C] rC = lmoqa_C - a_M[1] - b_M[1] * lbam_C;
      real s2 = square(sig_C[1]);
      real t2 = square(tau_M[1]);
      real den = s2 + n_C * t2;
      target += scale_lpdf(sig_C[1] | prior_scale[11], quartic);
      target += -0.5 * (n_C * log2pi + (n_C - 1) * log(s2) + log(den)
        + (dot_self(rC) - t2 * square(sum(rC)) / den) / s2);
    }
  }
  if (has_col) {
    // TROPOMI priors and the state-free residuals of every row.
    vector[n_cells] mu_cell = alpha_T[1] + U_cell * b_cell;
    kappa_T = exp(log_kappa_T[1]);
    alpha_T ~ normal(0, 3);
    gamma_T ~ normal(0, 0.3);
    if (field_N) {
      lambda_S ~ normal(0, 0.5);
      mu_cell += lambda_S[1] * (B_col * z_N);
    }
    if (col_temporal) lambda_T ~ normal(0, 0.5);
    log_kappa_T ~ normal(0, 1);
    target += scale_lpdf(sigma_q[1] | prior_scale[14], quartic);
    target += scale_lpdf(sigma_beta[1] | prior_scale[15], quartic);
    target += beta_lpdf((rho_q[1] + 1) / 2 | 2, 2);
    b_cell ~ normal(0, sigma_beta[1]);
    // With a_c the cell mean (mu_cell) and h_t = H1[t] gamma_T, the month's
    // sums m1_t = sum_c w_tc (y_tc - a_c) and m2_t = sum_c w_tc (y_tc - a_c)^2
    // give rbar_t = m1_t / S_t - h_t and sum_c w_tc (r_tc - rbar_t)^2 =
    // m2_t - m1_t^2 / S_t (the seasonal term cancels within the month).
    {
      vector[T] m1 = WY1 - Wm * mu_cell;
      vector[T] m2 = WY2 - 2 * (WY * mu_cell) + Wm * square(mu_cell);
      rbar_col = m1 .* inv_S_col - H1 * gamma_T;
      target += -0.5 * kappa_T * sum(m2 - square(m1) .* inv_S_col);
    }
    tq2 = square(sigma_q[1]);
    rq = rho_q[1];
    if (col_temporal) lamT = lambda_T[1] / s_N;
    // The parts of R_t that do not depend on the residuals.
    target += col_constant + 0.5 * (n_col - n_months_col) * log_kappa_T[1];
  }

  if (use_kappa) k = kappa[1];
  if (has_a && free_bA) bA = b_A[1];
  if (coupled) rr = r[1];
  if (H > 0) {
    f_N += Hm * h_N;
    f_P += Hm * h_P;
  }
  if (P_N > 0) f_N += Z_N * beta_N;
  if (P_P > 0) f_P += Z_P * beta_P;
  if (field_N) lr_N = log_profile(U, wbar, z_N, sd_N[1]);
  if (field_P) {
    vector[K - 1] zc = z_P;
    if (n_rhoS) zc = rho_S[1] * z_N + sqrt(1 - square(rho_S[1])) * z_P;
    lr_P = log_profile(U, wbar, zc, sd_P[1]);
  }
  for (j in 1:J) {
    u_N[j] = log_sum_exp(log(W_site[j])' + lr_N);
    u_P[j] = log_sum_exp(log(W_site[j])' + lr_P);
  }
  if (Q_N > 0) u_N += Xs_N * theta_N;
  if (Q_P > 0) u_P += Xs_P * theta_P;
  iN = s_N * sqrt(1 - square(rho_N));
  iP = s_P * sqrt(1 - square(rho_P));
  c0 = rr * iN * iP / (1 - rho_N * rho_P);
  // Core transition diag(phi4): AR(1) for e and q, restart (0) for m.
  phi4 = [rho_N, rho_P, 0, rq]';
  Q4 = [[square(iN), rr * iN * iP, 0, 0], [rr * iN * iP, square(iP), 0, 0],
        [0, 0, tm2, 0], [0, 0, 0, tq2 * (1 - square(rq))]];
  P = [[square(s_N), c0, 0, 0], [c0, square(s_P), 0, 0], [0, 0, tm2, 0], [0, 0, 0, tq2]];
  for (t in 1:T) {
    if (t > 1) {
      A = diag_pre_multiply(phi4, A);
      P = quad_form_diag(P, phi4) + Q4;
    }
    for (q in 1:count[t]) {
      int i = start[t] + q - 1;
      real m;
      real noise2;
      real F;
      vector[4] Pz;
      row_vector[M + 1] vrow;
      if (type[i] == 1) {
        real amp = 1;
        real shift = 0;
        if (cls[i] > 1) {
          shift = g[g_of[cls[i]]];
          amp = k;
        }
        m = b_N + u_N[site[i]] + shift + amp * f_N[t];
        noise2 = cls[i] > 1 ? square(sig_NT) : square(sig_N);
        Pz = amp * col(P, 1);
        F = amp * Pz[1] + noise2;
        vrow = -amp * A[1];
      } else if (type[i] == 2) {
        m = b_P + u_P[site[i]] + f_P[t];
        noise2 = square(sig_P);
        Pz = col(P, 2);
        F = Pz[2] + noise2;
        vrow = -A[2];
        if (nug[site[i]] > 0) vrow[1 + nug[site[i]]] += 1;
      } else if (type[i] == 3) {
        m = d_A[1] + b_P + bA * f_P[t];
        if (acag_seasonal) m += H1[t] * gamma_A;
        noise2 = square(sig_A[1]);
        Pz = bA * col(P, 2);
        F = bA * Pz[2] + noise2;
        vrow = -bA * A[2];
      } else {
        m = a_M[1] + bM * (b_P + u_P[site[i]] + f_P[t]);
        noise2 = square(sig_M[1]);
        Pz = bM * col(P, 2) + col(P, 3);
        F = bM * Pz[2] + Pz[3] + noise2;
        vrow = -(bM * A[2] + A[3]);
        vrow[1 + Jn + dev[i]] += 1;
        if (nug[site[i]] > 0) vrow[1 + nug[site[i]]] += bM;
      }
      vrow[1] += ly[i] - m;
      irow += 1;
      W[irow] = vrow / sqrt(F);
      logF[irow] = log(F);
      A += (Pz / F) * vrow;
      P -= (Pz * Pz') / F;
    }
    if (has_col && count_col[t] > 0) {
      // The month's column rows collapsed into one pseudo-row: the
      // precision-weighted mean residual rbar_col[t] loads e_N by lambda_T / s_N
      // and q by 1 with variance 1 / (kappa_T S_t).
      real F;
      vector[4] Pz;
      row_vector[M + 1] vrow;
      Pz = lamT * col(P, 1) + col(P, 4);
      F = lamT * Pz[1] + Pz[4] + inv_S_col[t] / kappa_T;
      vrow = -(lamT * A[1] + A[4]);
      vrow[1] += rbar_col[t];
      irow += 1;
      W[irow] = vrow / sqrt(F);
      logF[irow] = log(F);
      A += (Pz / F) * vrow;
      P -= (Pz * Pz') / F;
    }
  }
  // Prediction-error decomposition over the rows, then the Gaussian integral
  // over the M static effects: W'W = [[q, s'], [s, S]].
  target += -0.5 * n_rows * log2pi - 0.5 * sum(logF);
  {
    matrix[M + 1, M + 1] G = crossprod(W);
    target += -0.5 * G[1, 1];
    if (M > 0) {
      vector[M] dinv;
      for (j in 1:Jn) dinv[j] = 1 / square(tau_P[1]);
      for (d in 1:(has_m * D)) dinv[Jn + d] = 1 / square(tau_M[1]);
      {
        matrix[M, M] L = cholesky_decompose(add_diag(G[2:(M + 1), 2:(M + 1)], dinv));
        vector[M] w = mdivide_left_tri_low(L, G[2:(M + 1), 1]);
        target += 0.5 * dot_self(w) - sum(log(diagonal(L))) + 0.5 * sum(log(dinv));
      }
    }
  }
}
