#include<iostream>
#include <cstdlib> // 包含 EXIT_FAILURE 宏
#include "kband.h"
#include "fasta.h"
int main(int argc, char** argv) {

  if (argc < 3) {
    std::cerr << "Too few input parameters\n";
    std::cerr << "Usage: " << argv[0] << " <input_path> <output_path>\n";
    exit(EXIT_FAILURE);
  }

  std::string readPath = argv[1];
  std::string writePath = argv[2];
  // std::cout << readPath << "\n";
  std::vector<std::string> seqs, labels;
  readFasta(readPath, seqs, labels);
  writeFasta(writePath, seqs, labels);
  // std::string s = "TGGGTACCACCCAAGTATTGA", t = "TGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATGT";
  // PSA(s, t);
  return 0;
}