#include <iostream>
#include <vector>
#include <algorithm>
#include "../lib/utimer.hpp"
#include <ff/ff.hpp>


using namespace std;
using namespace ff;

// struct modeling a specimen for Genetic TSP
struct popMember{
    int fitness;
    vector<int> specimen;
};

vector<vector<int>> graph; // graph on which the algorithm will iterate
vector<popMember> population; // the population is represented by a vector of specimens
int n; // graph nodes number
int w; // number of workers
int pop_size; // population size auxiliary variable 
int mut_rate; // population size auxiliary variable 
int maxit; // number of times the algorithm must iterate
// auxiliary vector where the new generation will be stored before merging with previous population
vector<popMember> nextGen;


/*    
    @brief Metod which calculates the fitness function of a specimen in Genetic TSP problem by counting the total cost of the path.
    @param specimen -> input vector modeling the path.
    @return int -> fitness value of the given path.
*/
int fitnessFunction(vector<int> specimen){
    int w_sum = 0;
    int i = 0;

    for (int i = 0; i < n; i++)
        w_sum += graph[specimen[i]][specimen[(i + 1) % (n - 1)]];
    return w_sum;
}

/*    
    @brief Metod which selects a random specimen (Genetic Algorithm Selection Phase) in order to be a parent of the next generation.
            A given parent is selected by calculating a probability based on its fitness evaluation.
    @param seed_external -> thread local seed for generating pseudorandom numbers.
    @return int -> index of the selected parent respect to the population vector.
*/
int parentSelection(unsigned int seed_external){
    thread_local unsigned int seed = seed_external;
    bool done = false;
    int index;
    int rnd;
    float probability;

    while (true){
        index = rand_r(&seed) % pop_size;
        rnd = rand_r(&seed) % 101;

        probability = (1 / population[index].fitness) / population[pop_size - 1].fitness;

        if (rnd / 100 <= probability)
            return index;
    }
    return 0;
}

struct Emitter   : ff_monode{

    int generations;
    int parents;
    bool just_started = true;
    int w_counter;
    int w_num;
    vector<popMember>* nextGen;
    vector<popMember>* population;

    Emitter(int g, int w, vector<popMember>* ng, vector<popMember>* pop){
        generations = g;
        nextGen = ng;
        population = pop;
        parents = ng->size();
        w_counter = w;
        w_num = w;
    }

    void* svc(void*) {

        // The genetic algorithm has ended, begin termination sending End Of Stream
        if (generations == 0) {
            return EOS;
        }

        // The algorithm is just started, so there is no merge to perform, the workers
        // are now able to begin their work.
        if(just_started){
            just_started = false;
            broadcast_task(GO_ON);
            return GO_ON;
        }

        // Wait for all the workers to end the computation for the new generation
        if (!just_started && w_counter != 1){
            w_counter --;
            return GO_ON;
        }

        // The new generation has been computed by workers, the emitter merges the new specimens
        // with the population, substituing the less valuable specimens of the current population
        swap_ranges(population->end() - parents, population->end(), nextGen->begin());
        sort(population->begin(), population->end(), [](popMember &a, popMember &b)
            { return a.fitness < b.fitness; });
        generations--;
        w_counter = w_num;

        // The task is broadcasted to all workers, time to compute the next generation
        broadcast_task(GO_ON);
        return GO_ON;
    }
};

struct Worker   : ff_node {

    pair<int,int>* indexes;
    unsigned int s;
    vector<popMember>* nextGen;
    vector<popMember>* population;
    int count;
    int parents;
    int crossoverIndex = 0;
    int i1; // mutation swap index 1
    int i2; // mutation swap index 2

    popMember g1;
    popMember g2;
    popMember s1;
    popMember s2;

    // Vectors used for paths correction after crossovers and mutations, the occurences of each nodes is counted
    // in the corresponding index of the vector
    vector<int> founds1;
    bool c1 = false;
    vector<int> founds2;
    bool c2 = false;

    /*    
    @brief Worker ff_node constructior method.
    @param p -> pointer to a pair of indexes representing the start and end (excluded) of the partition of the nextGen vector
                        where the worker will operate (unique and independent form every other worker)
    @param ng -> pointer to the auxiliary vector where the new generation will be stored before merging with previous population
    @param pop -> pointer to the population vector modeling the population of specimens of the current computation
    @param seed -> thread local seed for generating pseudorandom numbers.
*/
    Worker(pair<int,int>* p, unsigned int seed, vector<popMember>* ng, vector<popMember>* pop){
        indexes = p;
        nextGen = ng;
        population = pop;
        s = seed;

        count = indexes->first;
        parents = indexes->second - indexes->first;
    }

    void* svc(void*) {

        thread_local unsigned int seed = s;
        count = indexes->first;

        while (count < indexes->second){
            c1 = false;
            c2 = false;
            g1 = (*population)[parentSelection(seed)];

            g2 = (*population)[parentSelection(seed)];

            // *** CROSSOVER ***
            crossoverIndex = rand_r(&seed) % n;
            if ((n - crossoverIndex) < 2) crossoverIndex = n - 2;
            else if (crossoverIndex < 2) crossoverIndex = 2;
            s1.specimen = vector<int>(n);
            s2.specimen = vector<int>(n);

            for (int j = 0; j < n; j++)
                {
                    s1.specimen[j] = (j < crossoverIndex) * g1.specimen[j] + (j >= crossoverIndex) * g2.specimen[j];
                    s2.specimen[j] = (j < crossoverIndex) * g2.specimen[j] + (j >= crossoverIndex) * g1.specimen[j];
                }

            // Correction of paths -> duplicate nodes are substituted with missing nodes
            founds1 = vector<int>(n, 0);
            founds2 = vector<int>(n, 0);

            for (int j = 0; j < n; j++){
                c1 = c1 || ++founds1[s1.specimen[j]] > 1;
                c2 = c2 || ++founds2[s2.specimen[j]] > 1;
            }

            if (c1)
                for (int j = 0; j < n; j++){
                    if (founds1[s1.specimen[j]] == 3)
                        s1.specimen[j] = distance(founds1.begin(), find(founds1.begin(), founds1.end(), 0));
                    founds1[s1.specimen[j]]++;
                }

            if (c2)
                for (int j = 0; j < n; j++){
                    if (founds2[s2.specimen[j]] == 3)
                        s2.specimen[j] = distance(founds2.begin(), find(founds2.begin(), founds2.end(), 0));
                    founds2[s2.specimen[j]]++;
                }

            // *** MUTATION ***
            //  With a mut_rate% probability 2 nodes in the paths are swapped
            if (rand_r(&seed) % 100 <= mut_rate){
                i1 = rand_r(&seed) % (n - 1) + 1;
                i2 = rand_r(&seed) % (n - 1) + 1;
                if (i1 == i2){
                    if (i2 != 1)
                        i2 = (i2 + rand_r(&seed) % (i2 - 1) + 1) % (n - 1) + 1;
                    else
                        i2 = rand_r(&seed) % (n - 2) + 2;
                }
                swap(s1.specimen[i1], s1.specimen[i2]);
            }

            if (rand_r(&seed) % 100 <= mut_rate){
                i1 = rand_r(&seed) % (n - 1) + 1;
                i2 = rand_r(&seed) % (n - 1) + 1;
                if (i1 == i2){
                    if (i2 != 1)
                        i2 = (i2 + rand_r(&seed) % (i2 - 1) + 1) % (n - 1) + 1;
                    else
                        i2 = rand_r(&seed) % (n - 2) + 2;
                }
                swap(s2.specimen[i1], s2.specimen[i2]);
            }

            // Fitness is evalued on children for the next generation
            s1.fitness = fitnessFunction(s1.specimen);
            s2.fitness = fitnessFunction(s2.specimen);

            (*nextGen)[count] = s1;
            count++;
            (*nextGen)[count] = s2;
            count++;
        }

        ff_send_out(GO_ON);
        return(GO_ON);
    }
};

int main(int argc, char *argv[]){
/* 
    INPUTS for Execution (in order):
    n -> number of nodes for generating a random graph.
    start -> starting node for TSP's paths.
    pop_size -> population size.
    ng_percentage -> percentage of specimen generating new generations (1 - 100)
    mut_rate -> mutation rate (0 - 100)
    maxit -> max iterations number (number of generations)
    w -> the number of workers for parallelisation
*/

    if (argc < 8)
    {
        cerr << "Not Enough Arguments\n";
        return -1;
    }
    n = atoi(argv[1]);
    int start = atoi(argv[2]);
    pop_size = atoi(argv[3]);
    int ng_percentage = atoi(argv[4]);
    mut_rate = atoi(argv[5]);
    if (start >= n || start < 0)
    {
        cerr << "Invalid Start \n";
        return -1;
    }
    if (ng_percentage <= 0 || ng_percentage > 100)
    {
        cerr << "Invalid ng_percentage: 1-100 \n";
        return -1;
    }
    if (mut_rate < 0 || mut_rate > 100)
    {
        cerr << "Invalid mt_rate: 0-100 \n";
        return -1;
    }
    maxit = atoi(argv[6]);
    w = atoi(argv[7]);
    unsigned int seed = 5400;

    // random graph generation (the seed is fixed for performance analysis reasons)
    graph = vector<vector<int>>(n);
    int weigth = 0;

    for (int i = 0; i < n; i++)
        graph[i] = vector<int>(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            weigth = rand_r(&seed) % 10; // generated weigths ranges from 1 to 10
            weigth++;
            graph[i][j] = weigth;
            graph[j][i] = weigth;
        }
    }

    // *** FIRST GENERATION - INITIALISATION ***

    utimer t(" - Whole Computation: ");
    population = vector<popMember>(pop_size);

    for (int i = 0; i < pop_size; i++){
        popMember aux;

        // 1 0 2 3 4 ...
        aux.specimen = vector<int>(n);
        aux.specimen[0] = start;
        for (int j = 0; j < n - 1; j++){
            aux.specimen[j + 1] = (j < start) * j + (j >= start) * (j + 1);
        }
        random_shuffle(aux.specimen.begin() + 1, aux.specimen.end());

        aux.fitness = fitnessFunction(aux.specimen);

        if (aux.specimen[0] != start){
            cout << "Error: Init";
            return -1;
        }
        population[i] = aux;
    }

    sort(population.begin(), population.end(), [](popMember &a, popMember &b)
         { return a.fitness < b.fitness; });

    int parents = pop_size * ng_percentage / 100;
    if (parents % 2 != 0)
        parents--;

    nextGen = vector<popMember>(parents);

    /* Indexes distibution for accessing the children specimens vector (nextGen)
     Workers know the partition of the nextGen vector that has been uniquely assigned to them, which is why,
     respecting this invariant, no concurrent accesses will be made to the same memory locations */

    if(parents/2 < w){
       cout << "\n Worker number is too high for the problem";
       cout << "\n Worker reduction: " << w <<" -> "<< parents/2;
       w = parents/2;
    }

    vector<pair<int, int>> w_indexes(w);
    int load = parents / w;
    if (load % 2 == 1)
        load--; // base value of load per worker
    int last_index = -1;

    for (int i = 0; i < w; i++)
    {
        w_indexes[i].first = ++last_index;
        w_indexes[i].second = w_indexes[i].first + load;
        last_index = w_indexes[i].second - 1;
    }

    // Distribution of the remaining load between workers
    int loads_to_add = (parents - last_index - 1) / 2;
    int count = 1;

    for (int i = w - loads_to_add; i < w; i++)
    {
        w_indexes[i].first = w_indexes[i - 1].second;
        w_indexes[i].second = w_indexes[i].second + 2 * count;
        count++;
    }

    // Workers Initialisation
    vector<unique_ptr<ff_node>> w_pool;
    for(int i = 0; i<w; i++) w_pool.push_back(make_unique<Worker>(&w_indexes[i], seed, &nextGen, &population));

    // Map Initialisation and Running
    ff_Farm<> mapGenerations(move(w_pool));
    Emitter e(maxit, w, &nextGen, &population);
    mapGenerations.add_emitter(e);
    mapGenerations.remove_collector();
    mapGenerations.wrap_around();
    mapGenerations.run_and_wait_end();

    cout << "\n\n TSP Solution: " << (*population.begin()).fitness << " ";
    cout << "\n\n FastFlow Parallel Version ";
}