// Implementation of the Map ADT

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Map.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static char *myStrdup(char *s);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct map {
   int numCities; // number of cities
   int numRoads; // number of roads
   char **cities; // stores the names for all the cities in order
   double **edges; // adjacency matrix storing positive weights, 0 if node is not adjacent
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Creates a new map with the given number of cities
 * Assumes that `numCities` is positive
 */
Map MapNew(int numCities) {
    struct map *new = malloc(sizeof(*new));
    new->numRoads = 0;
    new->numCities = numCities;

    // creating array of strings with char limit 20
    new->cities = malloc(numCities * sizeof(char *));
    for (int i = 0; i < numCities; i++) {
        new->cities[i] = myStrdup("0");
    }

    // creating adjacencty matrix
    new->edges = malloc(numCities * sizeof(double *));
    for (int i = 0; i < numCities; i++) {
        new->edges[i] = calloc(numCities, sizeof(double));
    }
    return new;
}

/**
 * Frees all memory allocated to the given map
 */
void MapFree(Map m) {
    // freeing the adjacency matrix
    for (int i = 0; i < m->numCities; i++) {
        free(m->edges[i]);
    }
    free(m->edges);
    // freeing the array of strings
    for (int i = 0; i < m->numCities; i++) {
        free(m->cities[i]);
    }
    free(m->cities);
    free(m);
}

/**
 * Returns the number of cities on the given map
 */
int MapNumCities(Map m) {
    return m->numCities;
}

/**
 * Returns the number of roads on the given map
 */
int MapNumRoads(Map m) {
    return m->numRoads;
}

/**
 * Sets the name of the given city
 * If the city's name has already been set, renames it
 */
void MapSetName(Map m, int city, char *name) {
    free(m->cities[city]);
    m->cities[city] = myStrdup(name);
}

/**
 * Returns the name of the given city, or "unnamed" if the city's name
 * has not been set
 */
char *MapGetName(Map m, int city) {
    char* name = m->cities[city];

    int hasName = strcmp(name, "0");
    // strcmp returns 0 if the strings are matching
    if (hasName == 0) {
        return "unnamed";
    }
    return name;
}

/**
 * Inserts a road between two cities with the given length
 * Does nothing if there is already a road between the two cities
 * Assumes that the cities are valid and are not the same
 * Assumes that the length of the road is positive
 */
void MapInsertRoad(Map m, int city1, int city2, int length) {
    if (m->edges[city1][city2] == 0) {
        m->edges[city1][city2] = length;
        m->edges[city2][city1] = length;
        m->numRoads++;
    }
}

/**
 * Returns the length of the road between two cities, or 0 if no such
 * road exists
 */
int MapContainsRoad(Map m, int city1, int city2) {
    return m->edges[city1][city2];
}

/**
 * Stores the roads connected to the given city in the given `roads`
 * array in any order and returns the number of roads stored. The `from`
 * field should be equal to `city` for all roads in the array.
 * Assumes that the roads array is at least as large as the number of
 * cities on the map.
 */
int MapGetRoadsFrom(Map m, int city, struct road roads[]) {
    int counter = 0;
    for (int i = 0; i < m->numCities; i++) {
        if (m->edges[city][i] != 0) {
            roads[counter].from = city;
            roads[counter].to = i;
            roads[counter].length = m->edges[city][i];
            counter++;
        }
    }
    return counter;
}

/**
 * Displays the map
 * !!! DO NOT EDIT THIS FUNCTION !!!
 * This function will work once the other functions are working
 */
void MapShow(Map m) {
    printf("Number of cities: %d\n", MapNumCities(m));
    printf("Number of roads: %d\n", MapNumRoads(m));

    struct road *roads = malloc(MapNumRoads(m) * sizeof(struct road));
    if (roads == NULL) {
        fprintf(stderr, "error: out of memory\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < MapNumCities(m); i++) {
        printf("[%d] %s has roads to:", i, MapGetName(m, i));
        int numRoads = MapGetRoadsFrom(m, i, roads);
        for (int j = 0; j < numRoads; j++) {
            if (j > 0) {
                printf(",");
            }
            printf(" [%d] %s (%d)", roads[j].to, MapGetName(m, roads[j].to),
                   roads[j].length);
        }
        printf("\n");
    }

    free(roads);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static char *myStrdup(char *s) {
    char *copy = malloc((strlen(s) + 1) * sizeof(char));
    return strcpy(copy, s);
}
