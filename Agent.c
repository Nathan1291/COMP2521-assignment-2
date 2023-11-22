// Implementation of the Agent ADT

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Agent.h"
#include "Map.h"

// Stack Adt implementation
// Linked list without duplicates to suit algorithm
// if it is to insert a duplicate, insert into the start of list
// then remove the value that was already in the list
struct stack {
    struct node *head;
};

// Linked list structure to store ther path for the dfs and djikstras
struct node {
    int id;
    struct node *next;
};

struct linkedList {
    struct node *head;
    struct node *tail;
};


// This struct stores information about an individual agent and can be
// used to store information that the agent needs to remember.
struct agent {
    char *name;
    int startLocation;
    int location;
    int maxStamina; // max stamina
    int stamina;    // current stamina
    int strategy;
    Map map;

    // TODO: Add more fields here as needed
    int *visited; // array to keep track of how many times a city has been visited
    struct linkedList *dfs; // using a linked list to keep track of the dfs algorithm that must be traversed
    struct linkedList *djikstras; // using a linked list to keep track of the shortest path that must be traversed for stage 3
};

static struct move chooseRandomMove(Agent agent, Map m);
static struct move chooseCheapestLeastVisited(Agent agent, Map map);
static struct move CheapestLeastVisitedAlgorithm(Agent agent, Map map);
static struct move chooseDFS(Agent agent, Map map);
static struct move DFSAlgorithm(Agent agent, Map map);


static int filterRoads(Agent agent, struct road roads[], int numRoads,
                       struct road legalRoads[]);
static int findLeastVisited(Agent agent, int numLegalRoads, struct road legalRoads[],
                       struct road leastVisitedRoads[]);
static struct road leastStamina(int numLeastVisited, struct road leastVisitedRoads[]);
void createDFSPath(Agent agent, Map map);
int findValidPredecessor(Map map, int visited[], int predecessors[], struct linkedList *list, int id);
int getStaminaCost(Map map, int location, int next);
int leastTurnMostStaminaCity(int visited[], int stamina[], int citiesTraversedSize[], int numCities);

////////////////////////////////////////// Node Functions /////////////////////////////////////////
struct node *newNode(int id);
////////////////////////////////////////// Node Functions /////////////////////////////////////////

/////////////////////////////////////////// LL Functions //////////////////////////////////////////
struct linkedList *newLinkedList(void);
void linkedListInsert(struct linkedList *dfs, int id);
void linkedListPop(struct linkedList *dfs);
void freeLinkedList(struct linkedList *dfs);
/////////////////////////////////////////// LL Functions //////////////////////////////////////////

///////////////////////////////////////// Stack Functions /////////////////////////////////////////
struct stack *newStack(void);
void stackInsert(struct stack *stack, int next);
int stackPop(struct stack *stack);
int stackIsEmpty(struct stack *stack);
void freeStack(struct stack *stack);

///////////////////////////////////////// Stack Functions /////////////////////////////////////////


/**
 * Creates a new agent
 */
Agent AgentNew(int start, int stamina, int strategy, Map m, char *name) {
    if (start >= MapNumCities(m)) {
        fprintf(stderr, "error: starting city (%d) is invalid\n", start);
        exit(EXIT_FAILURE);
    }

    Agent agent = malloc(sizeof(struct agent));
    if (agent == NULL) {
        fprintf(stderr, "error: out of memory\n");
        exit(EXIT_FAILURE);
    }

    agent->startLocation = start;
    agent->location = start;
    agent->maxStamina = stamina;
    agent->stamina = stamina;
    agent->strategy = strategy;
    agent->map = m;
    agent->name = strdup(name);


    int numCities = MapNumCities(m);
    // Storing all values as 0 for the visited array and initialising the start location as visited
    agent->visited = calloc(numCities, sizeof(int));
    agent->visited[start] = 1;
    // Initialising the dfs variable as NULL
    agent->dfs = newLinkedList();
    // Initialising the djikstras variable as NULL
    agent->djikstras = newLinkedList();
    return agent;
}

/**
 * Frees all memory associated with the agent
 * NOTE: You should not free the map because the map is owned by the
 *       main program, and the main program will free it
 */
void AgentFree(Agent agent) {
    // TODO: You may need to update this to free any extra memory you use
    free(agent->name);
    // freeing the visited array
    free(agent->visited);
    // freeing the linked list
    freeLinkedList(agent->dfs);
    free(agent->dfs);
    freeLinkedList(agent->djikstras);
    free(agent->djikstras);
    free(agent);
}

////////////////////////////////////////////////////////////////////////
// Gets information about the agent
// NOTE: It is expected that these functions do not need to be modified

/**
 * Gets the name of the agent
 */
char *AgentName(Agent agent) {
    return agent->name;
}

/**
 * Gets the current location of the agent
 */
int AgentLocation(Agent agent) {
    return agent->location;
}

/**
 * Gets the amount of stamina the agent currently has
 */
int AgentStamina(Agent agent) {
    return agent->stamina;
}

////////////////////////////////////////////////////////////////////////
// Making moves

/**
 * Calculates the agent's next move
 * NOTE: Does NOT actually carry out the move
 */
struct move AgentGetNextMove(Agent agent, Map m) {
    switch (agent->strategy) {
        case STATIONARY:              return (struct move){agent->location, 0};
        case RANDOM:                  return chooseRandomMove(agent, m);
        // TODO: Implement more strategies here
        case CHEAPEST_LEAST_VISITED : return chooseCheapestLeastVisited(agent, m);
        case DFS :                    return chooseDFS(agent, m);

        default:
            printf("error: strategy not implemented yet\n");
            exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////// STAGE 0 ////////////////////////////////////////////////
static struct move chooseRandomMove(Agent agent, Map m) {
    struct road *roads = malloc(MapNumCities(m) * sizeof(struct road));
    struct road *legalRoads = malloc(MapNumCities(m) * sizeof(struct road));

    // Get all roads to adjacent cities
    int numRoads = MapGetRoadsFrom(m, agent->location, roads);

    // Filter out roads that the agent does not have enough stamina for
    int numLegalRoads = filterRoads(agent, roads, numRoads, legalRoads);

    struct move move;
    if (numLegalRoads > 0) {
        // Sort the roads using insertion sort
        for (int i = 1; i < numLegalRoads; i++) {
            struct road r = legalRoads[i];
            int j = i;
            while (j > 0 && r.to < legalRoads[j - 1].to) {
                legalRoads[j] = legalRoads[j - 1];
                j--;
            }
            legalRoads[j] = r;
        }

        // nextMove is randomly chosen from the legal roads
        int k = rand() % numLegalRoads;
        move = (struct move){legalRoads[k].to, legalRoads[k].length};
    } else {
        // The agent must stay in the same location
        move = (struct move){agent->location, 0};
    }

    free(legalRoads);
    free(roads);
    return move;
}
//////////////////////////////////////// STAGE 0 ////////////////////////////////////////////////


//////////////////////////////////////// STAGE 1 ////////////////////////////////////////////////=
// Algorithm Explanation
// Simple algorithm that simply narrows the choices which the detective can move through
// based on the given criteria

// main function that stores basic algorithm
// checks if there is an existing djikstras path otherwise just run the regular CheapestLeastVisited Algorithm
// agentNextMove automatically updates visited so the function does not need any extra implementation
static struct move chooseCheapestLeastVisited(Agent agent, Map map) {
    // If an informant path does not exist
    struct move move;
    if (agent->djikstras->head == NULL) {
        move = CheapestLeastVisitedAlgorithm(agent, map);
        return move;
    }
    // if an informant path exists
    // extracting what city is next in the djikstras path
    int next = agent->djikstras->head->id;
    // Checking if there is enough stamina for the move
    int cost = getStaminaCost(map, agent->location, next);
    // if there isnt enough stamina, stay at current location and end the function
    if (agent->stamina < cost) {
        move = (struct move){agent->location, 0};
        return move;
    }
    // if there is enough stamina, pop the first node of the dfs linked list and return the move
    linkedListPop(agent->djikstras);
    return (struct move){next, cost};
}

// finds the cheapest, least visited and legal road the agent can move to
// takes in the road with smallest id if there are overlaps
static struct move CheapestLeastVisitedAlgorithm(Agent agent, Map map) {
    struct road *roads = malloc(MapNumCities(map) * sizeof(struct road));
    struct road *legalRoads = malloc(MapNumCities(map) * sizeof(struct road));
    struct road *leastVisitedRoads = malloc(MapNumCities(map) * sizeof(struct road));
    struct move move;

    // Get all roads to adjacent cities
    int numRoads = MapGetRoadsFrom(map, agent->location, roads);

    // Filter out roads that the agent does not have enough stamina for
    int numLegalRoads = filterRoads(agent, roads, numRoads, legalRoads);

    // If there are no legal roads, make the agent stay and end the function
    if (numLegalRoads == 0) {
        move = (struct move){agent->location, 0};
        free(roads);
        free(legalRoads);
        free(leastVisitedRoads);
        return move;
    }

    // Filter out roads that are least visited so that we can find the right road
    int numLeastVisited = findLeastVisited(agent, numLegalRoads, legalRoads, leastVisitedRoads);

    // if there is only 1 least visited, simply use that road for move
    if (numLeastVisited == 1) {
        move = (struct move){leastVisitedRoads[0].to, leastVisitedRoads[0].length};
    }
    // if there is more than 1 least visited, find the one with the lowest stamina cost
    else {
        struct road leastStaminaRoad = leastStamina(numLeastVisited, leastVisitedRoads);
        move = (struct move){leastStaminaRoad.to, leastStaminaRoad.length};
    }

    // free all malloced items and return results
    free(roads);
    free(legalRoads);
    free(leastVisitedRoads);
    return move;
}

//////////////////////////////////////// STAGE 1 ////////////////////////////////////////////////


//////////////////////////////////////// STAGE 2 ////////////////////////////////////////////////

// main function that stores basic algorithm
// checks if there is an existing djikstras path otherwise just run the regular Dfs Algorithm
// The current dfs path is cleared if there is an existing djikstras path
static struct move chooseDFS(Agent agent, Map map) {
    // If an informant path does not exist
    struct move move;
    if (agent->djikstras->head == NULL) {
        move = DFSAlgorithm(agent, map);
        return move;
    }
    // if an informant path exists
    // Clearing the current dfs path to ensure a new one is developed when the djikstras algorithm completes
    freeLinkedList(agent->dfs);
    agent->dfs->head = NULL;
    agent->dfs->tail = NULL;

    // extracting what city is next in the djikstras path
    int next = agent->djikstras->head->id;
    // Checking if there is enough stamina for the move
    int cost = getStaminaCost(map, agent->location, next);
    // if there isnt enough stamina, stay at current location and end the function
    if (agent->stamina < cost) {
        move = (struct move){agent->location, 0};
        return move;
    }
    // if there is enough stamina, pop the first node of the dfs linked list and return the move
    linkedListPop(agent->djikstras);
    return (struct move){next, cost};
}

// Algorithm Explanation
// The Algorithm will check if the dfs property is NULL, if it is, this means that we need to
// calculate a dfs path. Once a dfs path is found, it will be stored within the dfs property in the
// form of a linked list. The movement logic will run regardless of if the dfs property is NULL or not,
// and will simply use the information provided in the linked list to check if it is a valid movement
static struct move DFSAlgorithm(Agent agent, Map map) {
    struct move move;
    // If the dfs path doesn't exist, we need to calculate the path before we move the agent

    if (agent->dfs->head == NULL) {
        // creates a path but removes the first variable since it is the current position
        createDFSPath(agent, map);
        linkedListPop(agent->dfs);
    }
    // Moving the agent
    // extracting what city is next in the dfs path
    int next = agent->dfs->head->id;
    // Checking if there is enough stamina for the move
    int cost = getStaminaCost(map, agent->location, next);
    // if there isnt enough stamina, stay at current location and end the function
    if (agent->stamina < cost) {
        move = (struct move){agent->location, 0};
        return move;
    }
    // if there is enough stamina, pop the first node of the dfs linked list and return the move
    linkedListPop(agent->dfs);
    return (struct move){next, cost};
}
//////////////////////////////////////// STAGE 2 ////////////////////////////////////////////////

//////////////////////////////////////// STAGE 3 ////////////////////////////////////////////////


//////////////////////////////////////// STAGE 3 ////////////////////////////////////////////////

//////////////////////////////////////// Helper Functions ////////////////////////////////////////////////

// Takes an array with all the possible roads and puts the ones the agent
// has enough stamina for into the legalRoads array
// Returns the number of legal roads
static int filterRoads(Agent agent, struct road roads[], int numRoads,
                       struct road legalRoads[]) {
    int numLegalRoads = 0;
    for (int i = 0; i < numRoads; i++) {
        if (roads[i].length <= agent->stamina) {
            legalRoads[numLegalRoads++] = roads[i];
        }
    }
    return numLegalRoads;
}

/**
 * Executes a given move by updating the agent's internal state
 */
void AgentMakeNextMove(Agent agent, struct move move) {
    if (move.to == agent->location) {
        agent->stamina = agent->maxStamina;
    } else {
        agent->stamina -= move.staminaCost;
    }
    agent->location = move.to;

    // TODO: You may need to add to this to handle different strategies
    // adds to the number of times the city has been visited
    agent->visited[move.to]++;
}

// Takes an array of possible roads and finds which cities has been least visited
// and adds them into the leastVisitedRoads array
// return the number of leastvisited cities there are
static int findLeastVisited(Agent agent, int numLegalRoads, struct road legalRoads[], struct road leastVisitedRoads[]) {
    // method is less efficient than upkeeping leastvisitedroads as you go through the
    // list in a single for loop but is easier to implement and is not much slower
    int min = agent->visited[legalRoads[0].to];
    // find the least number of visits in the given roads
    for (int i = 0; i < numLegalRoads; i++) {
        if (min > agent->visited[legalRoads[i].to]) {
            min = agent->visited[legalRoads[i].to];
        }
    }

    // keeps track of the number of leastvisited roads
    int numLeastVisited = 0;
    // adds all the roads with the given number of visits in the array
    for (int i = 0; i < numLegalRoads; i++) {
        if (min == agent->visited[legalRoads[i].to]) {
            leastVisitedRoads[numLeastVisited] = legalRoads[i];
            numLeastVisited++;
        }
    }
    return numLeastVisited;
}

// takes an array of possible roads with equal numbers of visits and returns
// which road requires the least stamina
// if there are multiple roads with equal stamina, return the road with the lowest id
static struct road leastStamina(int numLeastVisited, struct road leastVisitedRoads[]) {
    // goes through each road and finds the lowest stamina cost
    int min = leastVisitedRoads[0].length;
    for (int i = 0; i < numLeastVisited; i++) {
        if (min > leastVisitedRoads[i].length) {
            min = leastVisitedRoads[i].length;
        }
    }
    // find the road with lowest stamina cost, since roads are already in ascending order
    // automatically returns road with the lowest id
    for (int i = 0; i < numLeastVisited; i++) {
        if (min == leastVisitedRoads[i].length) {
            return leastVisitedRoads[i];
        }
    }
    return leastVisitedRoads[0];
}

// calculates a dfs path as a linked list and stores it within the agent->dfs variable
// stupid ass stack that isnt a stack (linked list that doesnt allow duplicates)
// keeps track of a predecessor for each node. if an endpoint is reached, access predecessor
// array until an available node is open and push it to the stack, adds all predecessors to dfs path
void createDFSPath(Agent agent, Map map) {
    int numCities = MapNumCities(map);
    int location = agent->location;
    // array to track if we have covered all of the nodes in the map
    int visitedCount = 0;
    int *visited = calloc(numCities, sizeof(int));

    // array to track the predecessor of each array
    // initialise all values as numCities, as it is not possible for an id to be equal to numCites
    int *predecessors = calloc(numCities, sizeof(int));
    for (int i = 0; i < numCities; i++) {
        predecessors[i] = numCities;
    }
    int currentPredecessor = numCities;

    struct stack *stack = newStack();
    stackInsert(stack, location);
    // loops through until we have visited all the cities
    while (!stackIsEmpty(stack) && visitedCount < numCities) {
        int id = stackPop(stack);
        // Visit the city, if it hasn't already
        if (visited[id] == 0) {
            visited[id] = 1;
            visitedCount++;
            if (currentPredecessor != numCities) {
                predecessors[id] = currentPredecessor;
            }
        }
        // Insert the node into our path following the visit
        linkedListInsert(agent->dfs, id);
        // if we are on our last node, dont check for further cities to travel to
        if (visitedCount != numCities) {
            // instructions
            // add all possible roads from the id that is popped in ascending order
            struct road roads[numCities];
            int numRoads = MapGetRoadsFrom(map, id, roads);

            // if we are in a dead end, backtrack using find valid predecessors
            // dead ends are cities with adjacent citiees that we have already all visited
            int availableRoads = 0;
            for (int i = 0; i < numRoads; i++) {
                if (visited[roads[i].to] == 0) {
                    availableRoads++;
                }
            }
            if (availableRoads == 0) {
                id = findValidPredecessor(map, visited, predecessors, agent->dfs, id);
                numRoads = MapGetRoadsFrom(map, id, roads);
            }

            // adding all the possible adjacent cities to the current stack
            for (int i = numRoads-1; i >= 0; i--) {
                if (visited[roads[i].to] == 0) {
                    stackInsert(stack, roads[i].to);
                }
            }
            // stores the predecessor as the city that we inserted roads into the stack with
            currentPredecessor = id;
        }
    }
    // freeing everything
    free(visited);
    free(predecessors);
    freeStack(stack);
}

// Helper function for createDFSPath
// Finds the predecessor city which has an available road we can travel to
// visits the cities as we go down the predecessor list
int findValidPredecessor(Map map, int visited[], int predecessors[], struct linkedList *dfs, int id) {
    // bool to track if a city has been found
    int hasAvailableCity = 0;
    // tracking the cities we are backtracking to, starting from the previous city
    int curr = predecessors[id];
    int validCity = 0;
    int numCities = MapNumCities(map);

    while (!hasAvailableCity) {
        // adding to the path as we go
        linkedListInsert(dfs, curr);
        // checking if the current road we are on has a possible branch off
        struct road roads[numCities];
        int numRoads = MapGetRoadsFrom(map, curr, roads);
        for (int i = 0; i < numRoads; i++) {
            if (visited[roads[i].to] == 0) {
                hasAvailableCity = 1;
                validCity = curr;
            }
        }
        if (curr != numCities && predecessors[curr] != numCities) {
            curr = predecessors[curr];
        }
    }
    return validCity;
}

// returns the stamina cost of the given road
int getStaminaCost(Map map, int location, int next) {
    int numCities = MapNumCities(map);
    struct road roads[numCities];
    int numRoads = MapGetRoadsFrom(map, location, roads);
    int staminaCost = 0;
    for (int i = 0; i < numRoads; i++) {
        if (roads[i].to == next) {
            staminaCost = roads[i].length;
        }
    }
    return staminaCost;
}


////////////////////////////////////////////////////////////////////////
// Learning information

/**
 * Tells the agent where the thief is
 */
// Implementation of djikstras algorithm
// Stores 4 given arrays and updates the path as they go along
// keeps track of the shortest path to each respective node in the array
void AgentTipOff(Agent agent, int thiefLocation) {
    int numCities = MapNumCities(agent->map);
    // Stores which cities we have calculated distances from i.e. visited
    // 0 if it hasnt been visited and 1 if it has
    int numVisited = 0;
    int visited[numCities];
    for (int i = 0; i < numCities; i++) visited[i] = 0;
    // Stores the current distance to the given node, setting initial values to -1 to indicate its uncalculated
    int stamina[numCities];
    for (int i = 0; i < numCities; i++) stamina[i] = 0;
    // Stores the size of the following array of the cities stored
    int citiesTraversedSize[numCities];
    for (int i = 0; i < numCities; i++) citiesTraversedSize[i] = 0;
    // Stores the cities travelled through for the shortest path to the given node
    // for each given row, it lists the city in order
    // maximum path length is 2*numCities, in the case where the cities are all linearly connected by 1 road
    // and must constantly stop
    int citiesTraversed[numCities][2*numCities];

    int initialStamina = agent->stamina;
    // begin the process for the current location, running djikstras for the position we are a
    visited[agent->location] = 1;
    numVisited++;
    struct road roads[numCities];
    int numRoads = MapGetRoadsFrom(agent->map, agent->location, roads);
    for (int i = numRoads-1; i >= 0; i--) {
        if (visited[roads[i].to] == 0) {
            if (initialStamina < roads[i].length) {
                printf("%d", roads[i].to);
                stamina[roads[i].to] = agent->maxStamina - roads[i].length;
                citiesTraversed[roads[i].to][citiesTraversedSize[0]] = agent->location;
                citiesTraversed[roads[i].to][citiesTraversedSize[1]] = agent->location;
                citiesTraversedSize[roads[i].to] = 2;

            }
            else {
                stamina[roads[i].to] = initialStamina - roads[i].length;
                citiesTraversed[roads[i].to][citiesTraversedSize[0]] = agent->location;
                citiesTraversedSize[roads[i].to] = 1;
            }
        }
    }

    // loop through until all nodes have been visited
    while (numVisited < numCities) {
        // find the least turn, most stamina road
        int next = leastTurnMostStaminaCity(visited, stamina, citiesTraversedSize, numCities);

        initialStamina = stamina[next];
        visited[next] = 1;
        numVisited++;
        numRoads = MapGetRoadsFrom(agent->map, next, roads);
        for (int i = numRoads-1; i >= 0; i--) {
            // if we havent visited yet
            if (visited[roads[i].to] == 0) {
                // if there city hasn't had any paths leading to it, update it to the current one
                if (citiesTraversedSize[roads[i].to] == 0) {
                    // if we need to wait at the current city before proceeding
                    if (initialStamina < roads[i].length) {
                        stamina[roads[i].to] = agent->maxStamina - roads[i].length;
                        // matching the next city with the shortest path to the current city
                        for (int j = 0; j < citiesTraversedSize[next]; j++) {
                            citiesTraversed[roads[i].to][j] = citiesTraversed[next][j];
                        }
                        citiesTraversed[roads[i].to][citiesTraversedSize[next]] = next;
                        citiesTraversed[roads[i].to][citiesTraversedSize[next]+1] = next;
                        citiesTraversedSize[roads[i].to] = citiesTraversedSize[next] + 2;
                    }
                    // if we don't need to wait at the current city before moving
                    else {
                        stamina[roads[i].to] = initialStamina - roads[i].length;
                        // matching the next city with the shortest path to the current city
                        for (int j = 0; j < citiesTraversedSize[next]; j++) {
                            citiesTraversed[roads[i].to][j] = citiesTraversed[next][j];
                        }
                        citiesTraversed[roads[i].to][citiesTraversedSize[next]] = next;
                        citiesTraversedSize[roads[i].to] = citiesTraversedSize[next] + 1;
                    }
                }
                // if the city already has a path, we need to figure out which one is more efficient
                // if it is not as efficient we do not need to do anything
                else {
                    if (initialStamina < roads[i].length) {
                        // case where it is more efficient
                        if (citiesTraversedSize[next] + 2 < citiesTraversedSize[roads[i].to]) {
                            stamina[roads[i].to] = agent->maxStamina - roads[i].length;
                            // matching the next city with the shortest path to the current city
                            for (int j = 0; j < citiesTraversedSize[next]; j++) {
                                citiesTraversed[roads[i].to][j] = citiesTraversed[next][j];
                            }
                            citiesTraversed[roads[i].to][citiesTraversedSize[next]] = next;
                            citiesTraversed[roads[i].to][citiesTraversedSize[next]+1] = next;
                            citiesTraversedSize[roads[i].to] = citiesTraversedSize[next] + 2;
                        }
                        // case where it is equally efficient, we need to check stamina
                        if (citiesTraversedSize[next] + 2 == citiesTraversedSize[roads[i].to]) {
                            if (agent->maxStamina - roads[i].length < stamina[roads[i].to]) {
                                stamina[roads[i].to] = agent->maxStamina - roads[i].length;
                                // matching the next city with the shortest path to the current city
                                for (int j = 0; j < citiesTraversedSize[next]; j++) {
                                    citiesTraversed[roads[i].to][j] = citiesTraversed[next][j];
                                }
                                citiesTraversed[roads[i].to][citiesTraversedSize[next]] = next;
                                citiesTraversed[roads[i].to][citiesTraversedSize[next]+1] = next;
                                citiesTraversedSize[roads[i].to] = citiesTraversedSize[next] + 2;
                            }
                        }
                    }
                    else {
                        // case where it is more efficient
                        if (citiesTraversedSize[next] + 1 < citiesTraversedSize[roads[i].to]) {
                            stamina[roads[i].to] = initialStamina - roads[i].length;
                            // matching the next city with the shortest path to the current city
                            for (int j = 0; j < citiesTraversedSize[next]; j++) {
                                citiesTraversed[roads[i].to][j] = citiesTraversed[next][j];
                            }
                            citiesTraversed[roads[i].to][citiesTraversedSize[next]] = next;
                            citiesTraversedSize[roads[i].to] = citiesTraversedSize[next] + 1;
                        }
                        // case where it is equally efficient, we need to check stamina
                        if (citiesTraversedSize[next] + 2 == citiesTraversedSize[roads[i].to]) {
                            if (agent->maxStamina - roads[i].length < stamina[roads[i].to]) {
                                stamina[roads[i].to] = initialStamina - roads[i].length;
                                // matching the next city with the shortest path to the current city
                                for (int j = 0; j < citiesTraversedSize[next]; j++) {
                                    citiesTraversed[roads[i].to][j] = citiesTraversed[next][j];
                                }
                                citiesTraversed[roads[i].to][citiesTraversedSize[next]] = next;
                                citiesTraversedSize[roads[i].to] = citiesTraversedSize[next] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // once all the shortest paths have been calculated, we simply need to add it through onto the linked list implementation
    for (int i = 0; i <citiesTraversedSize[thiefLocation]; i++) {
        linkedListInsert(agent->djikstras, citiesTraversed[thiefLocation][i]);
    }
    linkedListInsert(agent->djikstras, thiefLocation);

    linkedListPop(agent->djikstras);
}


// Helper function to find the next city to visit
// of the available cities to visit, we search for the one with the lowest
// turns to get to and if there are cities with equal turns we pick the one with the highest
// stamina
int leastTurnMostStaminaCity(int visited[], int stamina[], int citiesTraversedSize[], int numCities) {
    int next = -1;
    int leastTurns;
    int mostStamina;
    for (int i = 0; i < numCities; i++) {
        // if it isn't visited
        if (visited[i] == 0) {
            // if we havent chosen a reference city yet
            if (next == -1 && citiesTraversedSize[i] > 0) {
                next = i;
                leastTurns = citiesTraversedSize[i];
                mostStamina = stamina[i];
            }
            else {
                // if there is a city with less turns than the current lowest
                if (citiesTraversedSize[i] < leastTurns && citiesTraversedSize[i] > 0) {
                    next = i;
                    leastTurns = citiesTraversedSize[i];
                    mostStamina = stamina[i];
                }
                // if they are the same, compare the stamina of the agents at the city
                else if (citiesTraversedSize[i] == leastTurns && citiesTraversedSize[i] > 0){
                    if (stamina[i] > mostStamina) {
                        next = i;
                        leastTurns = citiesTraversedSize[i];
                        mostStamina = stamina[i];
                    }
                }
            }
        }
    }
    if (next == -1) {
        return 0;
    }
    return next;
}

////////////////////////////////////////////////////////////////////////
// Displaying state

/**
 * Prints information about the agent (for debugging purposes)
 */
void AgentShow(Agent agent) {
    // TODO: You can implement this function however you want
    //       You can leave this function blank if you want
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Node Functions //////////////////////////////////////////

struct node *newNode(int id) {
    struct node *new = malloc(sizeof(struct node));
    new->id = id;
    new->next = NULL;
    return new;
}

////////////////////////////////////////// Node Functions //////////////////////////////////////////


////////////////////////////////////////// DFS Functions //////////////////////////////////////////

// deletes the first node in the dfs, frees it and moves dfs to the next node

struct linkedList *newLinkedList(void) {
    struct linkedList *new = malloc(sizeof(struct linkedList));
    new->head = NULL;
    new->tail = NULL;
    return new;
}
void linkedListInsert(struct linkedList *list, int id) {
    struct node *node = newNode(id);
    // if the path is empty
    if (list->head == NULL) {
        list->head = node;
        list->tail = node;
    }
    else {
        list->tail->next = node;
        list->tail = node;
    }
}

void linkedListPop(struct linkedList *list) {
    struct node *nodeToFree = list->head;
    list->head = list->head->next;
    free(nodeToFree);
}

// frees the given linked list
void freeLinkedList(struct linkedList *list) {
    struct node *curr = list->head;
    while (curr != NULL) {
        struct node *temp = curr;
        curr = curr->next;
        free(temp);
    }
}

////////////////////////////////////////// DFS Functions //////////////////////////////////////////



///////////////////////////////////////// Stack Functions /////////////////////////////////////////

// creates a new stack
struct stack *newStack(void) {
    struct stack *stack = malloc(sizeof(struct stack));
    stack->head = NULL;
    return stack;
}

// inserts a new node at the start of the stack and removes any other nodes with the same value within the stack
void stackInsert(struct stack *stack, int next) {
    struct node *node = newNode(next);
    node->next = stack->head;
    stack->head = node;

    // we dont need to check if there duplicates in the very first node since we just inserted
    struct node *curr = stack->head->next;
    struct node *prev = stack->head;

    // removing duplicates
    int hasDuplicate = 0;
    while (curr != NULL && !hasDuplicate) {
        if (curr->id == next) {
            struct node *temp = curr;
            prev->next = curr->next;
            hasDuplicate = 1;
            free(temp);
        }
        else {
            prev = curr;
            curr = curr->next;
        }
    }
}

// returns the id of the head of the stack and then pops it
int stackPop(struct stack *stack) {
    struct node *node = stack->head;
    stack->head = node->next;
    int id = node->id;
    free(node);
    return id;
}

// Returns 1 if the stack is empty and 0 if it isnt
int stackIsEmpty(struct stack *stack) {
    if (stack->head == NULL) {
        return 1;
    }
    return 0;
}

void freeStack(struct stack *stack) {
    struct node *curr = stack->head;
    while (curr != NULL) {
        struct node *temp = curr;
        curr = curr->next;
        free(temp);
    }
    free(stack);
}

///////////////////////////////////////// Stack Functions /////////////////////////////////////////
