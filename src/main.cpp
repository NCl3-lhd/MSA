#include <iostream>
#include <cstdlib> // 包含 EXIT_FAILURE 宏
#include "fasta.h"
#include "starAlign.h"
int main(int argc, char** argv) {

  if (argc < 4) {
    std::cerr << "Too few input parameters\n";
    std::cerr << "Usage: " << argv[0] << " <input_path> <output_path> <thread_count>\n";
    exit(EXIT_FAILURE);
  }

  std::string readPath = argv[1];
  std::string writePath = argv[2];
  int thread_count = std::stoi(std::string(argv[3]));
  // std::cout << readPath << "\n";
  std::vector<std::string> seqs, labels;
  readFasta(readPath, seqs, labels);
  // while (seqs.size() > 20) seqs.pop_back();
  // while (labels.size() > 20) labels.pop_back();
  std::vector<std::string> alignedSeqs;
  starAlign(seqs, alignedSeqs, thread_count);
  writeFasta(writePath, alignedSeqs, labels);
  // std::string s = "TGGGTACCACCCAAGTATTGA", t = "TGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATGT";
  // PSA(s, t);
  return 0;
}