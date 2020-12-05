# CE Experimentos

Experimentos para a disciplina de Computação Experimental.
O objetivo desses experimentos é testar a performance de implementações de bitvectors.

As bibliotecas testadas são [sdsl-lite](https://github.com/simongog/sdsl-lite) e [Sux](https://github.com/vigna/sux).

# Criando arquivos de teste

Os arquivos de teste representam bitvectors de tamanhos e densidades diferentes. Cada caractere (`0` ou `1`) representa um bit do bitvector (sim, não é tão eficiente...)

O tamanho e densidade de cada bitvector é definido no arquivo `createWorkloads.cpp`:
```
const int64_t bitVectorSize[11] = {1, 16, 64, 128, 256, 512, 1024}; // Tamanho em MB do bitVector
const double densityPercentage[6] = {0.05, 0.1, 0.25, 0.5, 0.75, 0.9}; // ( Num of 1s / length)
```

1. Crie uma pasta `data/`
```bash
$ mkdir data
```

2. Compile e execute o arquivo `createWorkloads.cpp`
```bash
$ g++ --std=c++17 createWorkloads.cpp
$ ./a.out
```
Isso deve gerar os arquivos de teste.