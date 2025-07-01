#include "kband.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
constexpr int INF = 1e9;
int mat = 1, mismat = 2, ogap = 3, egap = 1, k = 25;
int match(char si, char tj) {
  return si == tj ? mat : -mismat;
}
inline int convert(int j, int i) {
  return j - i + k;
}
int PSA_Kband(const std::string& _s, const std::string& _t, std::string* alignedS = nullptr, std::string* alignedT = nullptr) {
  std::string s = _s, t = _t;
  bool isSwap = 0;
  if (s.size() > t.size()) {
    std::swap(s, t);
    isSwap = 1;
  }
  int n = s.size();
  int m = t.size();
  int diff = m - n;
  // allocate
  std::vector<std::vector<int>> a(n + 1, std::vector<int>(diff + 2 * k + 1, -INF));
  std::vector<std::vector<int>> b(a), c(a);
  a[0][convert(0, 0)] = 0;
  auto isInsiderStrip = [&](int i, int j) {
    if (i < 0 || i > n + 1 || j < 0 || j > m + 1) return false;
    return (-k <= j - i && j - i <= k + diff);
  };
  // dp
  for (int i = 0; i <= n; i++) {
    for (int tj = 0; tj <= diff + 2 * k; tj++) {
      int j = tj + i - k;
      if (i - 1 >= 0 && j - 1 >= 0 && isInsiderStrip(i - 1, j - 1)) {
        a[i][convert(j, i)] = match(s[i - 1], t[j - 1]) + std::max({ a[i - 1][convert(j - 1, i - 1)], b[i - 1][convert(j - 1, i - 1)], c[i - 1][convert(j - 1, i - 1)] });
      }
      if (j - 1 >= 0 && isInsiderStrip(i, j - 1))
        b[i][convert(j, i)] = std::max(b[i][convert(j, i)], a[i][convert(j - 1, i)] - ogap);
      if (j - 1 >= 0 && isInsiderStrip(i, j - 1))
        b[i][convert(j, i)] = std::max(b[i][convert(j, i)], b[i][convert(j - 1, i)] - egap);
      if (i - 1 >= 0 && isInsiderStrip(i - 1, j))
        c[i][convert(j, i)] = std::max(c[i][convert(j, i)], a[i - 1][convert(j, i - 1)] - ogap);
      if (i - 1 >= 0 && isInsiderStrip(i - 1, j))
        c[i][convert(j, i)] = std::max(c[i][convert(j, i)], c[i - 1][convert(j, i - 1)] - egap);
    }
  }
  int ans = std::max({ a[n][convert(m, n)], b[n][convert(m, n)], c[n][convert(m, n)] });
  std::string S, T;

  //getWay
  int i = n, j = m;
  char op = a[n][convert(m, n)] > std::max(b[n][convert(m, n)], c[n][convert(m, n)]) ? 'a' : (b[n][convert(m, n)] > c[n][convert(m, n)] ? 'b' : 'c');

  while (i > 0 || j > 0) {
    // std::cout << op << " " << i << " " << j << "\n";
    if (op == 'a' && i - 1 >= 0 && j - 1 >= 0 && isInsiderStrip(i - 1, j - 1)) {
      S += s[i - 1];
      T += t[j - 1];
      op = a[i - 1][convert(j - 1, i - 1)] > std::max(b[i - 1][convert(j - 1, i - 1)], c[i - 1][convert(j - 1, i - 1)]) ? 'a' : (b[i - 1][convert(j - 1, i - 1)] > c[i - 1][convert(j - 1, i - 1)] ? 'b' : 'c');
      i--, j--;
    }
    else if (op == 'b' && j - 1 >= 0 && isInsiderStrip(i, j - 1)) {
      S += '-';
      T += t[j - 1];
      op = a[i][convert(j - 1, i)] - ogap > b[i][convert(j - 1, i)] - egap ? 'a' : 'b';
      j--;
    }
    else if (op == 'c' && i - 1 >= 0 && isInsiderStrip(i - 1, j)) {
      S += s[i - 1];
      T += '-';
      op = a[i - 1][convert(j, i - 1)] - ogap > c[i - 1][convert(j, i - 1)] - egap ? 'a' : 'c';
      i--;
    }
  }
  assert(op == 'a');
  reverse(S.begin(), S.end());
  reverse(T.begin(), T.end());
  if (isSwap) {
    // std::swap(s, t);
    std::swap(S, T);
  }
  if (alignedS != nullptr) {
    *alignedS = S;
  }
  if (alignedT != nullptr) {
    *alignedT = T;
  }
  // std::cout << "W: " << ans << "\n";
  // std::cout << "S: " << S << "\n";
  // std::cout << "T: " << T << "\n";
  return ans;
}


