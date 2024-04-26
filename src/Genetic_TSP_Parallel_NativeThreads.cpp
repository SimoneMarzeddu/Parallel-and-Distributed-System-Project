#include <iostream>
#include <vector>
#include <algorithm>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <atomic>
#include <barrier>
#include "../lib/utimer.hpp"

using namespace std;

// struct modeling a specimen for Genetic TSP
struct popMember{
    int index;
    int fitness;
    vector<int> specimen;
};

vector<vector<int>> graph; // graph on which the algorithm will iterate
vector<popMember> population; // the population is represented by a vector of specimens
atomic<int> id(0);
int n; // graph nodes number
int w; // number of workers
int pop_size; // population size auxiliary variable 
int mut_rate; // population size auxiliary variable 
int maxit; // number of times the algorithm must iterate
// auxiliary vector where the new generation will be stored before merging with previous population
vector<popMember> nextGen;

vector<vector<int>> random_values(w);

bool gen_ready = false; // condition waited for workers to start an iteration

mutex m1; // mutex associated with the next condition variable
condition_variable cv1; // condition variable where workers wait for the main thread to compute its merging operations


/*    
    @brief Metod called after barrier break, the "gen_ready" condition is set to false,
        the workers will wait to be awaken from the main thread after the new population has been formalised.
*/
void onCompletion ()noexcept{
    gen_ready = false;
}

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
int parentSelection(int id, int &dummy, int size){
    bool done = false;
    int index;
    int rnd;
    float probability;

    while (true)
    {
        index = random_values[id][dummy] % pop_size;
        dummy = (dummy + 1)%size;
        rnd = random_values[id][dummy] % 101;
        dummy = (dummy + 1)%size;
        probability = (1 / population[index].fitness) / population[pop_size - 1].fitness;

        if (rnd / 100 <= probability)
            return index;
    }
    return 0;
}


/*    
    @brief Core metod for workers execution which performs Selection Crossover Mutation end Evaluation phases of the
     Genetic TSP Algorithm. At the functional level it is the parallel transposition of the "while loop" within the main
      "for loop" of the sequential implementation of the algorithm.
    @param indexes -> a pair of indexes representing the start and end (excluded) of the partition of the nextGen vector
                        where the worker will operate (unique and independent form every other worker)
    @param b -> barrier for implementing the synchronization of all the algorithm iterations.
    @param seed_external -> thread local seed for generating pseudorandom numbers.
*/
void task(pair<int,int> indexes, barrier<void(*)(void) noexcept>&b, int id){

    int dummy = 0;
    int size = indexes.second - indexes.first;
    int count = indexes.first;
    int parents = indexes.second - indexes.first;
    int crossoverIndex = 0;
    int i1; // mutation swap index 1
    int i2; // mutation swap index 2

    popMember g1;
    popMember g2;
    popMember s1;
    popMember s2;
    int index_dummy = -1;

    // Vectors used for paths correction after crossovers and mutations, the occurences of each nodes is counted
    // in the corresponding index of the vector
    vector<int> founds1;
    bool c1 = false;
    vector<int> founds2;
    bool c2 = false;

    for (int i = 0; i < maxit; i++){

        count = indexes.first;

        // Workers wait for the main thread to formalize the current generation before calculating the next
        {
            unique_lock<mutex> l(m1);
            cv1.wait(l, [](){return gen_ready;});  
        }

        while (count < indexes.second){
            c1 = false;
            c2 = false;

            g1 = population[parentSelection(id, ref(dummy), size)];
            index_dummy = g1.index;

            while (index_dummy == g1.index){
                g2 = population[parentSelection(id, ref(dummy), size)];
                index_dummy = g2.index;
            }

            // *** CROSSOVER ***
            int r = random_values[id][dummy];
            dummy = (dummy + 1)%size;
            cout << "\n T:" << indexes.first << " " << r << "\n";
            crossoverIndex = r % n;
            if ((n - crossoverIndex) < 2) crossoverIndex = n - 2;
            else if (crossoverIndex < 2) crossoverIndex = 2;
            s1.specimen = vector<int>(n);
            s1.index = id;
            id++;
            s2.specimen = vector<int>(n);
            s2.index = id;
            id++;

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
            if (random_values[id][dummy] % 100 <= mut_rate){
                dummy = (dummy + 1)%size;
                i1 = random_values[id][dummy] % (n - 1) + 1;
                dummy = (dummy + 1)%size;
                i2 = random_values[id][dummy] % (n - 1) + 1;
                dummy = (dummy + 1)%size;
                if (i1 == i2){
                    if (i2 != 1)
                        i2 = (i2 + random_values[id][dummy] % (i2 - 1) + 1) % (n - 1) + 1;
                    else
                        i2 = random_values[id][dummy] % (n - 2) + 2;
                    dummy = (dummy + 1)%size;
                }
                swap(s1.specimen[i1], s1.specimen[i2]);
            }
            else dummy = (dummy + 1)%size;

            if (random_values[id][dummy] % 100 <= mut_rate){
                dummy = (dummy + 1)%size;
                i1 = random_values[id][dummy] % (n - 1) + 1;
                dummy = (dummy + 1)%size;
                i2 = random_values[id][dummy] % (n - 1) + 1;
                dummy = (dummy + 1)%size;
                if (i1 == i2){
                    if (i2 != 1)
                        i2 = (i2 + random_values[id][dummy] % (i2 - 1) + 1) % (n - 1) + 1;
                    else
                        i2 = random_values[id][dummy] % (n - 2) + 2;
                    dummy = (dummy + 1)%size;
                }
                swap(s2.specimen[i1], s2.specimen[i2]);
            }
            else dummy = (dummy + 1)%size;
            // Fitness is evalued on children for the next generation
            s1.fitness = fitnessFunction(s1.specimen);
            s2.fitness = fitnessFunction(s2.specimen);

            nextGen[count] = s1;
            count++;
            nextGen[count] = s2;
            count++;
        }
        b.arrive_and_wait();
    }

    return;
}

int main(int argc, char *argv[])
{
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

    // random graph generation (given seed)
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

    // set<popMember, decltype(&compareSpecimen)> population(&compareSpecimen);
    population = vector<popMember>(pop_size);

    { // Initialise timer block
        //utimer t(" - Initialisation: ");
        for (int i = 0; i < pop_size; i++)
        {
            popMember aux;
            aux.index = id;
            id++;

            // 1 0 2 3 4 ...
            aux.specimen = vector<int>(n);
            aux.specimen[0] = start;
            for (int j = 0; j < n - 1; j++)
            {
                aux.specimen[j + 1] = (j < start) * j + (j >= start) * (j + 1);
            }
            random_shuffle(aux.specimen.begin() + 1, aux.specimen.end());

            aux.fitness = fitnessFunction(aux.specimen);

            if (aux.specimen[0] != start)
            {
                cout << "Error: Init";
                return -1;
            }
            population[i] = aux;
        }
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

    barrier<void(*)(void) noexcept> b1(w + 1,onCompletion);

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

    //Random Number Sequential Generation
    srand(5400);

    //Workers Initialisation
    vector<thread> w_pool;
    for(int i = 0; i<w; i++){
        random_values[i] = vector<int>(w_indexes[i].second - w_indexes[i].first);
        for(int j = 0; j < w_indexes[i].second - w_indexes[i].first; j++)
            random_values[i][j] = rand();
        w_pool.push_back(thread(task, w_indexes[i], ref(b1), i)); 
    } 

    /*
        Main core loop interacting with workers:
        - Notify workers to start calculating the next generation
        - wait for ALL workers iteration execution (the main thread will arrive at a global barrier shared with all the workers)
        - Merge workers output (next generation formalisation)
    */
    for (int i = 0; i < maxit; i++){

        {
            unique_lock l(m1);
            gen_ready = true;
            cv1.notify_all();
        }

        b1.arrive_and_wait();  
        
        swap_ranges(population.end() - parents, population.end(), nextGen.begin());
        sort(population.begin(), population.end(), [](popMember &a, popMember &b)
            { return a.fitness < b.fitness; });
    }

    cout << "\n\n TSP Solution: " << (*population.begin()).fitness << " ";
    cout << "\n\n Native C++ Thread Parallel Version ";

    //Workers Join
    for(int i = 0; i<w; i++){
        w_pool[i].join(); 
    } 

}