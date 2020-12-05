# CE Experimentos

Experimentos para a disciplina de Computação Experimental.
O objetivo desses experimentos é testar a performance de implementações de bitvectors.

As bibliotecas testadas são [sdsl-lite](https://github.com/simongog/sdsl-lite) e [Sux](https://github.com/vigna/sux).

As estruturas que foram testadas:
* SDSL-LITE
    * rank_support_v (rank)
    * rank_support_v5 (rank)
    * rank_support_rrr (rank)
    * rank_support_sd (rank)
    * select_support_mcl (select)
    * select_support_rrr (select)
    * select_support_sd (select)

* SUX
    * Rank9Sel (rank + select)

Note que somente estruturas baseadas em bitvectors foram testadas.

Os resultados dos nossos testes podem ser vistos na pasta `results/`.
Os testes foram executados num ambiente WLS2 distro Ubuntu 20.4, numa máquina com Intel Core i7-7500U, com 16Gb de memória DDR4.

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

Para executar o arquivo de teste (`sdsl-lite.cpp` ou `sux.cpp`), basta compilá-lo (a primeira linha do arquivo contém a instrução de compilação, com todas as flags) e executar o arquivo gerado (default `a.out`);