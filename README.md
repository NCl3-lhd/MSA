# MSA

> multiple sequence alignment for RNA/DNA sequences

## How to use?

1. git clone this repository
   ```bash
   git clone https://github.com/NCl3-lhd/MSA.git
   ```

2. make
   ```bash
   cd MSA
   cmake -S . -B build
   cmake --build build
   ```

3. run
   ```bash
   ./build/MSA ./test/mtDNA.fasta ./test/result.fasta
   ```

4. clean
   ```bash
   cmake --build build --target clean
   ```

## License

Staralign is released under the MIT license. See the file LICENSE for more details.

## Contact

If you have any problems, please email me.(3300236038@qq.com)