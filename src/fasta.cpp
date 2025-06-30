#include "fasta.h"
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream> 
#include <cassert>

void readFasta(const std::string& path, std::vector<std::string>& seqs, std::vector<std::string>& labels) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + path);
  }

  std::string line;
  std::string currentSeq;

  while (std::getline(file, line)) {
    // 跳过空行
    if (line.empty()) continue;

    // 处理标签行
    if (line[0] == '>') {
      // 保存上一个序列（如果有）
      if (!currentSeq.empty()) {
        seqs.push_back(currentSeq);
        currentSeq.clear();
      }
      // 添加新标签（去掉开头的'>'）
      labels.push_back(line.substr(1));
    }
    else {// 处理序列行
      currentSeq += line;
    }
  }
  // 保存最后一个序列
  if (!currentSeq.empty()) {
    seqs.push_back(currentSeq);
  }
  file.close();
  assert(seqs.size() == labels.size());
}

void writeFasta(const std::string& path, std::vector<std::string>& seqs, std::vector<std::string>& labels) {
  std::ofstream file(path);
  // 检查序列和标签数量是否匹配
  if (seqs.size() != labels.size()) {
    throw std::invalid_argument("序列和标签数量不匹配");
  }

  if (!file.is_open()) {
    throw std::runtime_error("无法打开文件: " + path);
  }

  const int line_width = 70; // 每行60~80个字符的标准FASTA格式

  for (int i = 0; i < seqs.size(); ++i) {
    // 写入标签行
    file << '>' << labels[i] << '\n';
    
    // 写入序列（按指定宽度换行）
    const std::string& sequence = seqs[i];
    int pos = 0;
    const int len = sequence.length();
    
    while (pos < len) {
      // 计算本行应写入的字符数
      int writeCount = std::min(line_width, len - pos);
      
      // 写入一行序列
      file << sequence.substr(pos, writeCount) << '\n';
      
      // 移动到下一行起始位置
      pos += writeCount;
    }
    file << '\n';
  }
  file.close();
} 