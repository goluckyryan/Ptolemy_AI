// Fast SixJ using Racah's sum formula (replaces 6-nested m-loop)
// Reference: Varshalovich, Moskalev, Khersonskii - "Quantum Theory of Angular Momentum"
// Sum over integer k, typically 5-15 terms regardless of J magnitude.

static double log_factorial_cache[10000] = {0.0};
static bool lf_init = false;

static void init_log_factorial() {
    if (lf_init) return;
    log_factorial_cache[0] = 0.0;
    for (int i = 1; i < 10000; i++)
        log_factorial_cache[i] = log_factorial_cache[i-1] + std::log((double)i);
    lf_init = true;
}

static double log_fact(int n) {
    if (n < 0) return -1e300; // signals invalid
    if (n < 10000) return log_factorial_cache[n];
    // Stirling for large n
    return n * std::log(n) - n + 0.5 * std::log(2.0 * M_PI * n);
}

// Triangle coefficient delta(a,b,c) = sqrt[(a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!]
// Returns log of delta^2 = log[(a+b-c)!(a-b+c)!(-a+b+c)!] - log[(a+b+c+1)!]
static double log_delta2(double a, double b, double c) {
    int ia = (int)std::round(2*a), ib = (int)std::round(2*b), ic = (int)std::round(2*c);
    int n1 = (ia+ib-ic)/2, n2 = (ia-ib+ic)/2, n3 = (-ia+ib+ic)/2, n4 = (ia+ib+ic)/2+1;
    if (n1 < 0 || n2 < 0 || n3 < 0) return -1e300;
    return log_fact(n1) + log_fact(n2) + log_fact(n3) - log_fact(n4);
}

static bool triangle(double a, double b, double c) {
    double sum = a + b + c;
    return c >= std::abs(a-b) - 1e-9 && c <= a+b + 1e-9 &&
           std::abs(std::round(sum) - sum) < 1e-9; // a+b+c must be integer
}

double SixJSymbol_Racah(double J1, double J2, double J3, double J4, double J5, double J6) {
    init_log_factorial();
    // Check triangle conditions for all 4 triads
    if (!triangle(J1,J2,J3)) return 0.0;
    if (!triangle(J1,J5,J6)) return 0.0;
    if (!triangle(J4,J2,J6)) return 0.0;
    if (!triangle(J4,J5,J3)) return 0.0;

    // Half-integer to integer (doubled)
    int tJ1=(int)round(2*J1), tJ2=(int)round(2*J2), tJ3=(int)round(2*J3);
    int tJ4=(int)round(2*J4), tJ5=(int)round(2*J5), tJ6=(int)round(2*J6);

    // Prefactor from 4 delta functions:
    // W = delta(J1,J2,J3)*delta(J1,J5,J6)*delta(J4,J2,J6)*delta(J4,J5,J3)
    //   * sum_k (-1)^k (k+1)! / [...]
    double log_pre = 0.5 * (log_delta2(J1,J2,J3) + log_delta2(J1,J5,J6) +
                             log_delta2(J4,J2,J6) + log_delta2(J4,J5,J3));
    if (log_pre < -200) return 0.0;

    // Summation limits (in half-integer units × 2):
    int k_min_2 = std::max({tJ1+tJ2+tJ3, tJ1+tJ5+tJ6, tJ4+tJ2+tJ6, tJ4+tJ5+tJ3});
    int k_max_2 = std::min({tJ1+tJ2+tJ4+tJ5, tJ1+tJ3+tJ4+tJ6, tJ2+tJ3+tJ5+tJ6});

    if (k_min_2 > k_max_2) return 0.0;

    double result = 0.0;
    for (int tk = k_min_2; tk <= k_max_2; tk += 2) {
        int k = tk / 2; // k is always integer (since all J's are half-integer or integer, sum is always integer)
        // Denominator factorials:
        int d1=(tJ1+tJ2+tJ3)/2, d2=(tJ1+tJ5+tJ6)/2;
        int d3=(tJ4+tJ2+tJ6)/2, d4=(tJ4+tJ5+tJ3)/2;
        int d5=(tJ1+tJ2+tJ4+tJ5)/2-k, d6=(tJ1+tJ3+tJ4+tJ6)/2-k;
        int d7=(tJ2+tJ3+tJ5+tJ6)/2-k;
        int n1=k-d1, n2=k-d2, n3=k-d3, n4=k-d4;

        if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0) continue;
        if (d5 < 0 || d6 < 0 || d7 < 0) continue;

        double log_term = log_fact(k+1)
            - log_fact(n1) - log_fact(n2) - log_fact(n3) - log_fact(n4)
            - log_fact(d5) - log_fact(d6) - log_fact(d7);

        double term = std::exp(log_term + log_pre);
        if (k % 2 == 0) result += term;
        else            result -= term;
    }

    return result;
}
