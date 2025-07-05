#include "starAlign.h"
#include "kband.h"
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include "ThreadPool.h"
int findCenterSeqByLenth(const std::vector<std::string>& seqs) {
  int maxPos = -1;
  for (int i = 0; i < seqs.size(); i++) {
    if (maxPos == -1 || seqs[i].size() > seqs[maxPos].size()) {
      maxPos = i;
    }
  }
  return maxPos;
}
int findCenterSeqByScore(const std::vector<std::string>& seqs) {
  std::vector<long long> scores(seqs.size());
  int maxPos = -1;
  for (int i = 0; i < seqs.size(); i++) {
    long long sumScore = 0;
    std::cout << i << "\n";
    for (int j = 0; j < seqs.size(); j++) {
      if (j == i) continue;
      sumScore += PSA_Kband(seqs[i], seqs[j], nullptr, nullptr);
    }
    if (maxPos == -1 || sumScore > scores[maxPos]) {
      scores[i] = sumScore;
      maxPos = i;
    }
    // std::cout << scores[i] << "\n";
  }
  // std::cout << scores[maxPos] << "\n";
  return maxPos;
}
std::vector<std::array<std::string, 2>> getAlignedSeqs(int idxC, const std::vector<std::string>& seqs, int thread_count) {
  std::vector<std::array<std::string, 2>> res;
  res.reserve(int(seqs.size()) - 1);
  if (thread_count > 1) {
    ThreadPool pool(thread_count);
    std::vector<std::future<std::array<std::string, 2>>> results;
    const std::string seqC = seqs[idxC];
    for (int i = 0; i < seqs.size(); i++) {
      if (i == idxC) continue;
      results.emplace_back(pool.enqueue([&seqC, seqI = seqs[i]] {
        std::array<std::string, 2> aligned;
        PSA_Kband(seqC, seqI, &aligned[0], &aligned[1]);
        return aligned;
      }));
      // std::cout << i << "\n";
    }
    std::cout << seqs.size() << "\n";
    int cur = 0;
    for (auto&& result : results) {
      res.emplace_back(result.get());
      if (cur++ % 1000 == 0) std::cout << cur << "\n";
    }
    std::cout << "finish" << "\n";
  }
  else {
    std::array<std::string, 2> aligned;
    // int cur = 0;
    for (int i = 0; i < seqs.size(); i++) {
      if (i == idxC) continue;
      PSA_Kband(seqs[idxC], seqs[i], &aligned[0], &aligned[1]);
      res.emplace_back(std::move(aligned));
      // std::cout << cur++ << "\n";
    }
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

void starAlign(const std::vector<std::string>& seqs, std::vector<std::string>& alignedSeqs, int thread_count) {
  std::cerr << "findCenterSeqByLenth: " << "\n";
  int idxC = findCenterSeqByLenth(seqs);
  std::cerr << "getAlignedSeqs: " << "\n";
  std::vector<std::array<std::string, 2>> starsAligned = getAlignedSeqs(idxC, seqs, thread_count);
  std::cerr << "markInsertion: " << "\n";
  std::vector<int> markInsertion(seqs[idxC].size() + 1);
  getMarkInsertion(starsAligned, markInsertion);
  std::cerr << "insertGapToAlignedSeqs: " << "\n";
  alignedSeqs = insertGapToAlignedSeqs(seqs, idxC, starsAligned, markInsertion);
}