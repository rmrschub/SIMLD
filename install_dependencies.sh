 #!/usr/bin/env bash

git clone https://github.com/s-yata/marisa-trie.git extern/marisa-trie
cd extern/marisa-trie
 autoreconf -i
 ./configure --enable-native-code
 make
 make install

cd ..
cd ..

git clone https://github.com/powturbo/TurboPFor-Integer-Compression.git extern/TurboPFor-Integer-Compression
cd extern/TurboPFor-Integer-Compression
make