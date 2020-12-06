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

Os resultados dos nossos testes podem ser vistos na pasta `statistics/`.
Os testes foram executados num ambiente WLS2 distro Ubuntu 20.4, numa máquina com Intel Core i7-7500U, com 16Gb de memória DDR4.

# Criando arquivos de teste

## Dados
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

## Índices
3. Gere o arquivo de índices

Para as comparações ficarem mais justas, os índices consultados por todas as estruturas são as mesmas. Eles são criados a partir do arquivo `createIndex.cpp`. As duas macros importantes nesse arquivo são:

```
#define RUN_SIZE 1000  // Quantidade de índices que serão consultados por rodada
#define SIZE 64        // Tamanho do vetor para o qual esse índice deve ser gerado
```

Os arquivos de índice estão na pasta `utils/`.    
Os índices seguem uma ideia de _sliding window_, para tentar diminuir a quantidade de cache misses durante os testes. Cada índice está a uma distância de, no máximo, 1024 posições do índice anterior. Esses offsets são gerados usando uma distribuição uniforme.     

Para compilar e executar o arquivo gerador de índices pode-se usar os comandos abaixo.
```bash
$ g++ --std=c++17 createIndex.cpp
$ ./a.out
```
Isso deve gerar os arquivos de teste.

Para executar o arquivo de teste (`experiments.cpp`), basta compilá-lo (a primeira linha do arquivo contém a instrução de compilação, com todas as flags) e executar o arquivo gerado (default `a.out`). Os resultados serão colocados na pasta `statistics/`.

# Customizando os Testes

Macros relevantes para customização dos testes. Lembre de gerar os bitvector e os arquivos de índices **antes** de executar os testes.
```c++
#define SAMPLE_SIZE 500               // Numero de rodadas por teste
#define RUN_SIZE 1000                 // Numero de operações por rodada

#define SIZE 64                       // Tamanho do bitvector em MB
#define DENSITIES {0.05, 0.5, 0.9}    // Densidades
```