#include "kband.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
constexpr int INF = 1e9;
int mat = 1, mismat = 2, ogap = 3, egap = 1;
int match(char si, char tj) {
  return si == tj ? mat : -mismat;
}
inline int convert(int j, int i, int k) {
  return j - i + k;
}
int PSA_Kband(const std::string& _s, const std::string& _t, std::string* alignedS = nullptr, std::string* alignedT = nullptr) {
  // 使用指针避免字符串复制
  const std::string* s_ptr = &_s;
  const std::string* t_ptr = &_t;
  bool isSwap = false;

  if (_s.size() > _t.size()) {
    std::swap(s_ptr, t_ptr);
    isSwap = true;
  }
  const std::string& s = *s_ptr;
  const std::string& t = *t_ptr;

  int n = s.size();
  int m = t.size();
  int diff = m - n;
  int k = 25;
  int w = diff + 2 * k + 1;
  // allocate
  std::vector<int> a_flat((n + 1) * w, -INF);
  std::vector<int> b_flat((n + 1) * w, -INF);
  std::vector<int> c_flat((n + 1) * w, -INF);

  auto get_row = [w](std::vector<int>& mat, int i) {
    return mat.data() + i * w;
  };
  int* a0 = get_row(a_flat, 0);
  a0[convert(0, 0, k)] = 0;
  auto isInsiderStrip = [&](int i, int j) {
    if (i < 0 || i > n + 1 || j < 0 || j > m + 1) return false;
    return (-k <= j - i && j - i <= k + diff);
  };
  // dp
  for (int i = 0; i <= n; i++) {
    int* a_pre = get_row(a_flat, i - 1);
    int* b_pre = get_row(b_flat, i - 1);
    int* c_pre = get_row(c_flat, i - 1);

    int* a_cur = get_row(a_flat, i);
    int* b_cur = get_row(b_flat, i);
    int* c_cur = get_row(c_flat, i);
    const int diag_min = std::max(-k, -i);
    const int diag_max = std::min(k + diff, m - i);
    for (int offset = diag_min; offset <= diag_max; offset++) {
      int j = i + offset;
      int curtj = convert(j, i, k);
      if (isInsiderStrip(i - 1, j - 1)) {
        int pretj = convert(j - 1, i - 1, k);
        a_cur[curtj] = match(s[i - 1], t[j - 1]) + std::max({ a_pre[pretj], b_pre[pretj], c_pre[pretj] });
      }
      if (isInsiderStrip(i, j - 1)) {
        int pretj = convert(j - 1, i, k);
        b_cur[curtj] = std::max(b_cur[curtj], a_cur[pretj] - ogap);
        b_cur[curtj] = std::max(b_cur[curtj], b_cur[pretj] - egap);
      }
      if (isInsiderStrip(i - 1, j)) {
        int pretj = convert(j, i - 1, k);
        c_cur[curtj] = std::max(c_cur[curtj], a_pre[pretj] - ogap);
        c_cur[curtj] = std::max(c_cur[curtj], c_pre[pretj] - egap);
      }
    }
  }
  int* an = get_row(a_flat, n);
  int* bn = get_row(b_flat, n);
  int* cn = get_row(c_flat, n);
  int curtj = convert(m, n, k);
  int ans = std::max({ an[curtj], bn[curtj], cn[curtj] });
  std::string S, T;

  //getWay
  int i = n, j = m;

  char op = an[curtj] > std::max(bn[curtj], cn[curtj]) ? 'a' : (bn[curtj] > cn[curtj] ? 'b' : 'c');

  while (i > 0 || j > 0) {
    // std::cout << op << " " << i << " " << j << "\n";
    int* a_nxt = get_row(a_flat, i - 1);
    int* b_nxt = get_row(b_flat, i - 1);
    int* c_nxt = get_row(c_flat, i - 1);

    int* a_cur = get_row(a_flat, i);
    int* b_cur = get_row(b_flat, i);
    int* c_cur = get_row(c_flat, i);
    if (op == 'a' && isInsiderStrip(i - 1, j - 1)) {
      S += s[i - 1];
      T += t[j - 1];
      int nxttj = convert(j - 1, i - 1, k);
      op = a_nxt[nxttj] > std::max(b_nxt[nxttj], c_nxt[nxttj]) ? 'a' : (b_nxt[nxttj] > c_nxt[nxttj] ? 'b' : 'c');
      i--, j--;
    }
    else if (op == 'b' && isInsiderStrip(i, j - 1)) {
      S += '-';
      T += t[j - 1];
      int nxttj = convert(j - 1, i, k);
      op = a_cur[nxttj] - ogap > b_cur[nxttj] - egap ? 'a' : 'b';
      j--;
    }
    else if (op == 'c' && isInsiderStrip(i - 1, j)) {
      S += s[i - 1];
      T += '-';
      int nxttj = convert(j, i - 1, k);
      op = a_nxt[nxttj] - ogap > c_nxt[nxttj] - egap ? 'a' : 'c';
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


