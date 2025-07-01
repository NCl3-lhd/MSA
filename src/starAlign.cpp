#include "starAlign.h"
#include "kband.h"
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
int findCenterSeq(const std::vector<std::string>& seqs) {
  std::vector<int> scores(seqs.size());
  int maxPos = -1;
  for (int i = 0; i < seqs.size(); i++) {
    int sumScore = 0;
    for (int j = 0; j < seqs.size(); j++) {
      if (j == i) continue;
      sumScore += PSA_Kband(seqs[i], seqs[j], nullptr, nullptr);
    }
    if (maxPos == -1 || sumScore > scores[maxPos]) {
      scores[i] = sumScore;
      maxPos = i;
    }
  }
  return maxPos;
}
std::vector<std::array<std::string, 2>> getAlignedSeqs(int idxC, const std::vector<std::string>& seqs) {
  std::vector<std::array<std::string, 2>> res;
  std::string s, t;
  for (int i = 0; i < seqs.size(); i++) {
    if (i == idxC) continue;
    PSA_Kband(seqs[idxC], seqs[i], &s, &t);
    res.push_back({ s, t });
  }
  return res;
}
void getMarkInsertion(const std::vector<std::array<std::string, 2>>& starsAligned, std::vector<int>& markInsertion) {
  for (const std::array<std::string, 2> &str : starsAligned) {
    int gapCount = 0;
    int pi = 0;
    for (int i = 0; i < str[0].size(); i++) {
      char ch = str[0][i];
      if (ch == '-') {
        gapCount++;
      }
      else {
        markInsertion[pi] = std::max(markInsertion[pi], gapCount);
        pi++;
        gapCount = 0;
      }
    }
    markInsertion[pi] = std::max(markInsertion[pi], gapCount);
  }
}
std::string insertGap(const std::vector<int>& mark, const std::string& str) {
  assert(mark.size() == str.size() + 1);
  std::string res;
  for (int i = 0; i < str.size(); i++) {
    for (int j = 0; j < mark[i]; j++) res += '-';
    res += str[i];
  }
  for (int j = 0; j < mark[str.size()]; j++) res += '-';
  return res;
}
std::vector<std::string> insertGapToAlignedSeqs(const std::vector<std::string>& seqs, int idxC, const std::vector<std::array<std::string, 2>>& starsAligned, const std::vector<int>& markInsertion) {
  std::vector<std::string> res(seqs.size());
  res[idxC] = insertGap(markInsertion, seqs[idxC]);
  for (int i = 0; i < starsAligned.size(); i++) {
    const std::array<std::string, 2>& str = starsAligned[i];
    std::vector<int> mark(str[0].size() + 1);
    int gapCount = 0;
    int pi = 0;
    for (int j = 0; j < str[0].size(); j++) {
      char ch = str[0][j];
      if (ch == '-') {
        gapCount++;
      }
      else {
        mark[j - gapCount] = markInsertion[pi++] - gapCount;
        gapCount = 0;
      }
    }
    mark[int(str[0].size()) - gapCount] = markInsertion[pi] - gapCount;
    res[i + (i >= idxC)] = insertGap(mark, str[1]);
  }
  return res;
}

void starAlign(const std::vector<std::string>& seqs, std::vector<std::string>& alignedSeqs) {
  int idxC = findCenterSeq(seqs);
  std::vector<std::array<std::string, 2>> starsAligned = getAlignedSeqs(idxC, seqs);
  std::vector<int> markInsertion(seqs[idxC].size() + 1);
  getMarkInsertion(starsAligned, markInsertion);
  alignedSeqs = insertGapToAlignedSeqs(seqs, idxC, starsAligned, markInsertion);
}