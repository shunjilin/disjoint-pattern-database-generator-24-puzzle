//
// Created by shunji on 16/07/19.
// Disjoint Pattern Database for 24-Puzzle (Admissible, Inconsistent)
// According to Paper by Richard E. Korf and Ariel Felner
// "Disjoint pattern Database Heuristics"
// Artificial Intelligence 134 (2002) 9-22
//

#include "DPDB.h"

int main() {
    printf("start\n");
    try {
        DPDB search{{1, 2, 5, 6, 7, 12}};
        // DPDB search{{3, 4, 8, 9, 13, 14}};
        printf("search initialized\n");
        search.search();
        printf("end of search\n");
        search.output_db("pattern.tab");
    } catch(...) {
        throw;
    }
}
