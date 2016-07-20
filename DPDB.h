//
// Created by shunji on 16/07/19.
//
#ifndef DPDB_DPDB_H
#define DPDB_DPDB_H

#include <cstdlib>
#include <vector>
#include <string>
#include <array>

class DPDB {
public:
    static constexpr int kPatternSize = 6;
    struct Node {
        // inverse positions of pattern tiles,
        // takes 7 entries, where the last entry is the blank tile position
        std::array<unsigned char, kPatternSize + 1> pat_inv;
        unsigned char h_value;
        Node(std::vector<unsigned char> &pat_tiles);
        Node();
    };
    DPDB(const std::vector<unsigned char> &tiles);
    void search(); // perform retrograde search
    void output_db(std::string filename); // output into file
    std::vector<unsigned char> pat_tiles;

private:
    std::vector<unsigned char> h_table; // heuristic values
    std::vector<bool> visited;
    std::vector< std::vector<unsigned char> > move_table;

    // queues for BFS
    std::vector<Node> cur_buf;
    std::vector<Node> next_buf;
    int queue1_fd;
    int queue2_fd;
    int h_value;

    void init_node(Node &n);
    Node pop();
    void push(Node n);
    bool fill_buffer();
    void flush_cur_buffer();
    void flush_next_buffer();
    void switch_queues();
    size_t visited_hash(const Node &n);
    size_t dpdb_hash(const Node &n);

    // for debugging
    int manhattan_dist(unsigned char tile, unsigned char position);
    void debug(const Node &n);
};




#endif //DPDB_RETROGRADESEARCH_H
