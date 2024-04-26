#include <iostream>
#include <vector>
#include <algorithm>
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
int n; // graph nodes number
int pop_size; // population size auxiliary variable 

/*    
    @brief Metod which calculates the fitness function of a specimen in Genetic TSP problem by counting the total cost of the path.
    @param specimen -> input vector modeling the path.
    @return int -> fitness value of the given path.
*/
int fitnessFunction(vector<int> specimen){
    int w_sum = 0;
    int i = 0;

    for(int i = 0; i < n; i++) w_sum += graph[specimen[i]][specimen[(i+1)%(n-1)]];

    return w_sum;
}

/*    
    @brief Metod which selects a random specimen (Genetic Algorithm Selection Phase) in order to be a parent of the next generation.
            A given parent is selected by calculating a probability based on its fitness evaluation.
    @return int -> index of the selected parent respect to the population vector.
*/
int parentSelection(){
    bool done = false;
    int index;
    int rnd;
    float probability;

    while(true){
        index = rand()%pop_size;
        rnd = rand()%101;
        probability = (1/population[index].fitness)/population[pop_size-1].fitness;

        if(rnd/100 <= probability) return index; 
    }
    return 0;
}

int main(int argc, char* argv[]){
/* 
    INPUTS for Execution (in order):
    n -> number of nodes for generating a random graph.
    start -> starting node for TSP's paths.
    pop_size -> population size.
    ng_percentage -> percentage of specimen generating new generations (1 - 100)
    mut_rate -> mutation rate (0 - 100)
    maxit -> max iterations number (number of generations)
*/

    if(argc < 7){
        cerr << "Not Enough Arguments\n";
        return -1;
    }
    n = atoi(argv[1]);
    int start = atoi(argv[2]);
    pop_size = atoi(argv[3]);
    int ng_percentage = atoi(argv[4]);
    int mut_rate = atoi(argv[5]);
    if(start >= n || start < 0){
        cerr << "Invalid Start \n";
        return -1;
    }
    if(ng_percentage <= 0 || ng_percentage > 100 ){
        cerr << "Invalid ng_percentage: 1-100 \n";
        return -1;
    }
    if(mut_rate < 0 || mut_rate > 100 ){
        cerr << "Invalid mt_rate: 0-100 \n";
        return -1;
    }
    int maxit = atoi(argv[6]);
    srand(5400);

    // random graph generation (the seed is fixed for performance analysis reasons)
    graph = vector<vector<int>>(n);
    int weigth = 0;

    for(int i = 0; i<n; i++)
        graph[i] = vector<int>(n);

    for(int i = 0; i<n; i++){
        for(int j = 0; j<=i; j++){
            weigth = rand()%10; // generated weigths ranges from 1 to 10
            weigth++;
            graph[i][j] = weigth;
            graph[j][i] = weigth;
        }
    }

    // *** FIRST GENERATION - INITIALISATION ***

    utimer t(" - Whole Computation: ");

    int index = 0;
    population = vector<popMember>(pop_size);

    for(int i = 0; i<pop_size; i++){
        popMember aux;
        aux.index = index;
        index++;

        // 1 0 2 3 4 ...
        aux.specimen = vector<int>(n);
        aux.specimen[0] = start;
        for(int j = 0; j<n-1; j++){
            aux.specimen[j+1] = (j < start)*j + (j >= start)*(j+1);
        }
        random_shuffle(aux.specimen.begin() + 1, aux.specimen.end());

        aux.fitness = fitnessFunction(aux.specimen);

        if(aux.specimen[0] != start ){
            cout << "Error: Init";
            return -1;
        }
        population[i] = aux;
    }

    sort(population.begin(), population.end(), [](popMember& a, popMember& b) {return a.fitness < b.fitness; });

    int parents = pop_size * ng_percentage / 100;
    if(parents%2 != 0) parents--;

    int count = 0;
    int crossoverIndex = 0;
    int i1; // mutation swap index 1
    int i2; // mutation swap index 2

    popMember g1;
    popMember g2;
    popMember s1;
    popMember s2;
    int index_dummy = -1;

    // Vectors used for paths correction after crossovers and mutations, the occurences of each nodes is counted
    // in the location of corresponding index of the vector
    vector<int> founds1;
    bool c1 = false;
    vector<int> founds2;
    bool c2 = false;

    // auxiliary vector where the new generation will be stored before merging with previous population
    vector<popMember> nextGen(parents);

    /* 
        Sequential core loop, divided in 3 phases per iteration:
        - Parents selection
        - Next Generation (crossover and mutations)
        - Worst Specimens Elimination and Substitution
    */
    for(int i = 0; i < maxit; i++){

        count = 0;

        while(count < parents){

            c1 = false;
            c2 = false;

            g1 = population[parentSelection()];
            index_dummy = g1.index;

            while(index_dummy == g1.index){
                g2 = population[parentSelection()];
                index_dummy = g2.index;
            }    

            // *** CROSSOVER ***
            crossoverIndex = rand()%n;
            if((n - crossoverIndex) < 2) crossoverIndex = n - 2;
            else if(crossoverIndex < 2) crossoverIndex = 2;
            s1.specimen = vector<int>(n);
            s1.index = index;
            index++;
            s2.specimen = vector<int>(n);
            s2.index = index;
            index++;

            for(int j = 0; j<n; j++){
                s1.specimen[j] = (j < crossoverIndex)*g1.specimen[j] + (j >= crossoverIndex)*g2.specimen[j]; 
                s2.specimen[j] = (j < crossoverIndex)*g2.specimen[j] + (j >= crossoverIndex)*g1.specimen[j];
            }   
            
            // Correction of paths -> duplicate nodes are substituted with missing nodes
            founds1 = vector<int>(n,0);
            founds2 = vector<int>(n,0);

            for(int j = 0; j<n; j++){
                c1 = c1 || ++founds1[s1.specimen[j]] > 1;
                c2 = c2 || ++founds2[s2.specimen[j]] > 1;
            }

            if(c1)
                for(int j = 0; j<n; j++){
                    if(founds1[s1.specimen[j]] == 3) s1.specimen[j] = distance(founds1.begin(), find(founds1.begin(), founds1.end(), 0));
                    founds1[s1.specimen[j]] ++;
                }  
                    
            if(c2)
                for(int j = 0; j<n; j++){
                    if(founds2[s2.specimen[j]] == 3) s2.specimen[j] = distance(founds2.begin(), find(founds2.begin(), founds2.end(), 0));
                    founds2[s2.specimen[j]] ++;
                }

            // *** MUTATION ***
            // With a mut_rate% probability 2 nodes in the paths are swapped
            if(rand()%100 <= mut_rate){
                i1 = rand()%(n-1) + 1;
                i2 = rand()%(n-1) + 1;
                if(i1 == i2){
                    if(i2 != 1) i2 = (i2 + rand()%(i2-1) + 1)%(n-1) + 1;
                    else i2 = rand()%(n-2) + 2;
                }
                swap(s1.specimen[i1], s1.specimen[i2]);
            }
                    
            if(rand()%100 <= mut_rate){
                i1 = rand()%(n-1) + 1;
                i2 = rand()%(n-1) + 1;
                if(i1 == i2){
                    if(i2 != 1) i2 = (i2 + rand()%(i2-1) + 1)%(n-1) + 1;
                    else i2 = rand()%(n-2) + 2;
                }
                swap(s2.specimen[i1], s2.specimen[i2]);
            }

            // Fitness evaluation for generated children specimens
            s1.fitness = fitnessFunction(s1.specimen);
            s2.fitness = fitnessFunction(s2.specimen);
            
            nextGen[count] = s1;
            count++;
            nextGen[count] = s2;
            count++;
        }
        
        swap_ranges(population.end() - parents, population.end(), nextGen.begin());
        sort(population.begin(), population.end(), [](popMember& a, popMember& b) {return a.fitness < b.fitness; });
    }


    cout << "\n\n TSP Solution: " << (*population.begin()).fitness << " ";
    cout << "\n\n Sequential Version ";
}