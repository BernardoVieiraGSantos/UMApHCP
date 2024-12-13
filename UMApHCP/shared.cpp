#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shared.h"

int readInstance()
{
    // Inicialização padrão
    instanceEntries.nodeQuantity = 0;
    instanceEntries.coordinates = NULL;

    // Abrir o arquivo
    FILE* file = fopen(sourceFile, "r");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        printf("Caminho: %s\n", sourceFile);
        return 1;
    }

    if (showLogs) {
        printf("Arquivo aberto: %s\n", sourceFile);
    }

    // Ler a quantidade de hubs
    if (fscanf(file, "%d", &instanceEntries.nodeQuantity) != 1)
    {
        printf("Erro ao ler a quantidade de hubs no arquivo.\n");
        fclose(file);
        return 1;
    }

    if (showLogs) {
        printf("A quantidade de hub: %u\n", instanceEntries.nodeQuantity);
    }

    // Alocar memória para as coordenadas
    instanceEntries.coordinates = (Coordinate*)malloc(instanceEntries.nodeQuantity * sizeof(Coordinate));
    if (instanceEntries.coordinates == NULL)
    {
        printf("Erro ao alocar memória para coordenadas.\n");
        fclose(file);
        instanceEntries.nodeQuantity = 0;
        return 1;
    }

    if (showLogs) {
        printf("Memória alocada: %zu bytes\n", instanceEntries.nodeQuantity * sizeof(Coordinate));
    }
    // Ler as coordenadas
    for (int i = 0; i < instanceEntries.nodeQuantity; i++)
    {
        if (fscanf(file, "%lf %lf", &instanceEntries.coordinates[i].x, &instanceEntries.coordinates[i].y) != 2)
        {
            printf("Erro ao ler as coordenadas do hub %d.\n", i + 1);
            free(instanceEntries.coordinates);
            fclose(file);
            return 1;
        }
    }

    if (showLogs) {
        printf("Coordenadas lidas!!!\n");
    }

    fclose(file);
    return 0;
}

int freeCoordinatesInstance() {
    free(instanceEntries.coordinates);

    return 0;
}

int generateCostMatriz() {
    int size = instanceEntries.nodeQuantity; // Número de coordenadas
    Coordinate* coordinates = instanceEntries.coordinates;

    // Alocar a matriz de custos (matriz dinâmica)
    costMatrix = (double**)malloc(size * sizeof(double*));
    if (costMatrix == NULL) {
        printf("Erro ao alocar a matriz de custos.\n");
        return 1;
    }

    if (showLogs) {
        printf("Matriz de custo alocada\n");
    }

    for (int i = 0; i < size; i++) {
        costMatrix[i] = (double*)malloc(size * sizeof(double));
        if (costMatrix[i] == NULL) {
            printf("Erro ao alocar a linha %d da matriz de custos.\n", i);
            // Liberar linhas previamente alocadas antes de retornar
            for (int j = 0; j < i; j++) {
                free(costMatrix[j]);
            }
            free(costMatrix);
            return 1;
        }
    }

    if (showLogs) {
        printf("Matriz de custo lida\n");
    }

    // Preencher a matriz com os custos (distâncias)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                costMatrix[i][j] = 0.0; // Distância de um ponto para ele mesmo é 0
            }
            else {
                double dx = coordinates[i].x - coordinates[j].x;
                double dy = coordinates[i].y - coordinates[j].y;
                costMatrix[i][j] = sqrt(dx * dx + dy * dy); // Distância Euclidiana
            }
        }
    }

    if (showLogs) {
        printf("Matriz de custo preenchida\n");
    }

    return 0;
}

int freeCostMatriz() {
    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        free(costMatrix[i]);
    }
    free(costMatrix);

    return 0;
}

int generateCostAggMatriz() {
    costAggMatrix = (double**) calloc(sizeof(double*), instanceEntries.nodeQuantity);

    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        costAggMatrix[i] = (double*) calloc(sizeof(double), instanceEntries.nodeQuantity);
    }

    //preenchendo a matriz de custo agregada
    for (int i = 0; i < instanceEntries.nodeQuantity; i++)
    {
        for (int j = 0; j < instanceEntries.nodeQuantity; j++)
        {
            for (int m = 0; m < instanceEntries.nodeQuantity; m++)
            {
                for (int k = 0; k < instanceEntries.nodeQuantity; k++)
                {
                    costAggMatrix[i][k] += costMatrix[i][k] + alpha * costMatrix[k][m];
                    costAggMatrix[j][m] += costMatrix[m][j];
                }
            }
        }
    }

    return 0;
}

int freeCostAggMatriz() {
    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        free(costAggMatrix[i]);
    }

    free(costAggMatrix);

    return 0;
}