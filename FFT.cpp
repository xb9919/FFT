#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <complex>
#define N 1000100
#define PI acos(-1)
using namespace std;
typedef complex<double> cd;
typedef vector<cd> VC;

inline void FFT(VC& a, bool inv)
{
    int n = a.size();
    for (int i = 0, j = 0; i < n; ++i) {
        if (j > i) swap(a[i], a[j]);
        int k = n;
        for (; j & (k >>= 1); j &= ~k);
        j |= k;
    }
    double pi = inv ? -PI : PI;
    for (int step = 1; step < n; step <<= 1) {
        double alp = pi / step;
        for (int k = 0; k < step; ++k) {
            cd wk = exp(cd(0, alp * k));
            for (int e = k; e < n; e += step << 1) {
                int o = e + step;
                cd w = wk * a[o];
                a[o] = a[e] - w;
                a[e] += w;
            }
        }
    }
    if (inv) for (int i = 0; i < n; ++i) a[i] /= n;
}

inline VC operator * (const VC& x, const VC& y)
{
    int n = x.size(), m = y.size(), s = 2;
    for (; s < n + m; s <<= 1);
    VC a(s, 0), b(s, 0), c(s, 0);
    for (int i = 0; i < n; ++i) a[i] = x[i]; FFT(a, 0);
    for (int i = 0; i < m; ++i) b[i] = y[i]; FFT(b, 0);
    for (int i = 0; i < s; ++i) c[i] = a[i] * b[i]; FFT(c, 1);
    VC ret(n + m - 1);
    for (int i = 0; i < n + m - 1; ++i) ret[i] = fabs(c[i].real());
    return ret;
}

int n;
char s[N];

int main()
{
    //freopen("test.in", "r", stdin);
    int n, m; scanf("%d%d", &n, &m); ++n, ++m;
    VC a(n << 1), b(m << 1);
    for (int i = 0, x; i < n; ++i) scanf("%d", &x), a[i] = x;
    for (int i = 0, x; i < m; ++i) scanf("%d", &x), b[i] = x;
    for (int i = n; i < (n << 1); ++i) a[i] = 0;
    for (int i = m; i < (m << 1); ++i) b[i] = 0;
    VC ret = a * b;
    for (int i = 0; i < n + m - 1; ++i) printf("%.0f ", ret[i].real());
    return 0;
}