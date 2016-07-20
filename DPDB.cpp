//
// Created by shunji on 16/07/19.
//
#include <cstdlib>
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>
#include <sstream>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include "DPDB.h"

//conexpr int kPatternSize = 6; // number of pattern tiles
constexpr size_t kTableSize = 244140625; // disjoint pattern database size (sparse) = 25^6
constexpr size_t kVisitedSize = 6103515625; // visited bit vector size (sparse) = 25^7
constexpr int kBoardSize = 25; // size of board
constexpr int kBoardWidth = 5; // width of board
constexpr int kMaxMoves = 125;
constexpr unsigned char kSentinel = std::numeric_limits<unsigned char>::max(); // for unused entries
constexpr size_t kBlockSize = 512; // block size for flushes, kBlockSize * sizeof(Node) = 16kB


// constructor takes as input pattern tiles
DPDB::DPDB(const std::vector<unsigned char> &tiles) :
        h_table(kTableSize, std::numeric_limits<unsigned char>::max()),
        visited(kVisitedSize, false),
        pat_tiles(tiles),
        move_table(kBoardSize),
        h_value(0)
{
    printf("kvisited size is %lu",kVisitedSize);
    printf("size of node is %lu\n", sizeof(Node));
    if (pat_tiles.size() != kPatternSize) {
        throw "Pattern tiles size mismatch.";
    }


    // initialize move table,
    // new blank positions reachable by current blank position
    for (int i = 0; i < kBoardSize; ++i) {
        //up move
        if (i > 4) {
            move_table[i].emplace_back(i - 5);
        }
        // down move
        if (i < 20) {
            move_table[i].emplace_back(i + 5);
        }
        // left move
        if (i % 5 != 0) {
            move_table[i].emplace_back(i - 1);
        }
        // right move
        if (i % 5 != 4) {
            move_table[i].emplace_back(i + 1);
        }
    }

    // initialize queues, 1 for current h-depth, 1 for next h-depth
    std::ostringstream oss1, oss2;
    oss1 << "q1";
    oss2 << "q2";
    std::string filename1 = oss1.str(), filename2 = oss2.str();
    // create files
    queue1_fd = open(filename1.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (queue1_fd == -1) {
        throw "Open queue1 file failed.";
    }
    queue2_fd = open(filename2.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (queue2_fd == -1) {
        throw "Open queue2 file failed.";
    }

}

// create initial node from pattern tiles
DPDB::Node::Node(std::vector<unsigned char> &pat_tiles) {
    h_value = 0;
    for (int i = 0; i < kPatternSize; ++i) {
        printf("pat tile %d\n", pat_tiles[i]);
        pat_inv[i] = pat_tiles[i];
    }
    pat_inv[kPatternSize] = 0; // blank position at last index
}

// empty constructor, sentinel node
DPDB::Node::Node() : h_value(kSentinel) {};

// pop a node from current queue
DPDB::Node DPDB::pop() {
    if (h_value == kMaxMoves) { // end of search
        Node terminate = Node();
        return terminate;
    }
    if (cur_buf.empty()) {
        if (!fill_buffer()) { // end of queue
            switch_queues();
            return pop();
        }
    }
    assert(!cur_buf.empty());
    Node n = cur_buf.back();
    cur_buf.pop_back();
    return n;
}

// push node into current or next queue
void DPDB::push(Node n)  {
    if (n.h_value == h_value) { // push to current queue
        if (cur_buf.size() == kBlockSize) {
            flush_cur_buffer();
        }
        cur_buf.emplace_back(n);
    } else { //push to next queue
        if (next_buf.size() == kBlockSize) {
            flush_next_buffer();
        }
        next_buf.emplace_back(n);
    }
}

// fill current buffer, return false if empty queue
bool DPDB::fill_buffer() {
    int &cur_fd = (h_value % 2 == 0) ? queue1_fd : queue2_fd;
    assert(cur_buf.empty());
    cur_buf.resize(kBlockSize);
    ssize_t read_sz = read(cur_fd, &cur_buf[0], kBlockSize * sizeof(Node));
    if (read_sz < 0) throw "Fill buffer failed.";
    cur_buf.resize(read_sz / sizeof(Node));
    return (read_sz == 0) ? false : true;
}

// flush current buffer to file
void DPDB::flush_cur_buffer() {
    int &cur_fd = (h_value % 2 == 0) ? queue1_fd : queue2_fd;
    if (cur_buf.size() == 0) return;
    ssize_t orig_offset = lseek(cur_fd, 0, SEEK_CUR); // save current offset (for reading)
    if (orig_offset < 0) throw ("Saving current offset failed.");
    if (lseek(cur_fd, 0, SEEK_END) < 0) throw ("Seek to end of file failed.");
    ssize_t write_sz = cur_buf.size() * sizeof(Node);
    if (write(cur_fd, &cur_buf[0], write_sz) != write_sz) throw "Flush current buffer failed.";
    if (lseek(cur_fd, orig_offset, SEEK_SET) < 0) throw "Restoring original offset failed."; // restore offset
    cur_buf.clear();
}

// flush next buffer to file
void DPDB::flush_next_buffer() {
    int &next_fd = (h_value % 2 == 0) ? queue2_fd : queue1_fd;
    if (next_buf.size() == 0) return;
    ssize_t write_sz = next_buf.size() * sizeof(Node);
    if (write(next_fd, &next_buf[0], write_sz) != write_sz) throw "Flush next buffer failed.";
    next_buf.clear();
}

// cleanup at the end of current h-depth
// clear current queue and reset offsets
// clear buffers
void DPDB::switch_queues() {
    printf("finished depth %d\n", h_value);
    flush_next_buffer();

    if (h_value % 2 == 0) {
        if (ftruncate(queue1_fd, 0) < 0) throw "Truncating queue1 failed.";
    } else {
        if (ftruncate(queue2_fd, 0) < 0) throw "Truncating queue2 failed.";
    }

    h_value++;

    if (lseek(queue1_fd, 0, SEEK_SET) < 0) throw "lseek on queue1 failed.";
    if (lseek(queue2_fd, 0, SEEK_SET) < 0) throw "lseek on queue2 failed.";

    cur_buf.clear();
    next_buf.clear();
}

// hash function for visited states
// base 25, indexed on pattern tile positions and blank tile position
size_t DPDB::visited_hash(const Node &n) {
    size_t hash = 0;
    for (int i = 0; i < kPatternSize + 1; ++i) {
        hash *= kBoardSize;
        hash += n.pat_inv[i];
    }
    return hash;
}

// hash function for disjoint pattern database
// base 25, indexed on pattern tiles
size_t DPDB::dpdb_hash(const Node &n) {
    size_t hash = 0;
    for (int i = 0; i < kPatternSize; ++i) {
        hash *= kBoardSize;
        hash += n.pat_inv[i];
    }
    return hash;

}

// output database to file
void DPDB::output_db(std::string filename) {
    FILE * outfile = fopen(filename.c_str(), "wb");
    if (!outfile) throw "Opening dpdb file failed.";
    for (size_t i = 0; i < kTableSize; ++i) {
        putc(h_table[i], outfile);
    }
    fclose(outfile);
}

/* For Debugging */
int DPDB::manhattan_dist(unsigned char tile, unsigned char pos) {
    return (abs((tile % kBoardWidth) - (pos % kBoardWidth)) + abs((tile / kBoardWidth) - (pos / kBoardWidth)));
}

void DPDB::debug(const Node &n) {
    int manhat = 0;
    for (int i = 0; i < kPatternSize; ++i) {
        manhat += manhattan_dist(pat_tiles[i], n.pat_inv[i]);
    }
    if (h_table[dpdb_hash(n)] < manhat) printf("underestimate!\n");
}

// Perform retrograde search from initial position of pattern tiles and blank space
// Perform BFS in increasing h-values
// h-values only incremented when pattern tile is moved
// visited state bit vector prevents cycles (visited)
// store minimum h-value in disjoint pattern database (h_table)
void DPDB::search() {
    Node initial_n(pat_tiles);
    cur_buf.emplace_back(initial_n);
    visited[visited_hash(initial_n)] = true;
    h_table[dpdb_hash(initial_n)] = 0;

    while (true) {
        Node popped_n = pop();
        if (popped_n.h_value == kSentinel) {
            return;
        }
        unsigned char old_blank_pos = popped_n.pat_inv[kPatternSize];
        for (auto new_blank_pos : move_table[old_blank_pos]) {
            Node child = popped_n;
            for (int i = 0; i < kPatternSize; ++i) {
                if (child.pat_inv[i] == new_blank_pos) {
                    child.pat_inv[i] = old_blank_pos;
                    child.h_value++;
                    break;
                }
            }
            child.pat_inv[kPatternSize] = new_blank_pos;

            if (!visited[visited_hash(child)]) {
                visited[visited_hash(child)] = true;
                if (child.h_value < h_table[dpdb_hash(child)]) {
                    h_table[dpdb_hash(child)] = child.h_value;
                    //debug(child);
                }
                push(child);
            }
        }
    }

}

