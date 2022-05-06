#pragma once
#include <assert.h>
#include <cstring>
#include <set>
#include <vector>
#include <cstdio>
#include <algorithm>

#ifndef INDEX
#define INDEX(x, y, n) ((x) * (n) + (y))
#endif

class Pattern {
public:
    Pattern(int _size);
    Pattern(int _size, std::string buffer);
    ~Pattern();
    Pattern(const Pattern &pat);
    void add_edge(int x, int y);
    void del_edge(int x, int y);
    inline void add_ordered_edge(int x, int y) { adj_mat[INDEX(x, y, size)] = '1'; }
    inline int get_size() const { return size; }
    inline std::string get_adj_mat() const { return adj_mat; }
    bool check_connected() const;
    void count_all_isomorphism(std::set<std::set<int>> &s) const;
    void print() const;
    bool is_dag() const;

private:
    Pattern &operator=(const Pattern &);
    void get_full_permutation(std::vector<std::vector<int>> &vec, bool use[], std::vector<int> tmp_vec, int depth) const;
    std::string adj_mat;
    int size;
};

Pattern::Pattern(int _size) {
    size = _size;
    adj_mat = std::string(size*size, '0');
}

Pattern::Pattern(int _size, std::string buffer) {
    size = _size;
    adj_mat = buffer;
}

Pattern::~Pattern() = default;

Pattern::Pattern(const Pattern &pat) {
    size = pat.get_size();
    adj_mat = pat.get_adj_mat();
}

void Pattern::add_edge(int x, int y) {
    adj_mat[INDEX(x, y, size)] = '1';
    adj_mat[INDEX(y, x, size)] = '1';
}

void Pattern::del_edge(int x, int y) {
    adj_mat[INDEX(x, y, size)] = '0';
    adj_mat[INDEX(y, x, size)] = '0';
}

void Pattern::print() const {
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (adj_mat[INDEX(i, j, size)] != '0')
                printf("(%d,%d) ", i, j);
    printf("\n");
}

bool Pattern::is_dag() const {
    int degree[size];
    memset(degree, 0, size * sizeof(int));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (adj_mat[INDEX(i, j, size)] > '0')
                ++degree[j];
    int que[size];
    int head = 0;
    int tail = -1;
    for (int i = 0; i < size; ++i)
        if (degree[i] == 0)
            que[++tail] = i;
    while (head <= tail) {
        int x = que[head++];
        for (int j = 0; j < size; ++j) {
            if (adj_mat[INDEX(x, j, size)] > '0') {
                --degree[j];
                if (degree[j] == 0)
                    que[++tail] = j;
            }
        }
    }
    if (tail == size - 1)
        return true;
    else
        return false;
}

void Pattern::get_full_permutation(std::vector<std::vector<int>> &vec, bool use[], std::vector<int> tmp_vec, int depth) const {
    if (depth == size) {
        vec.push_back(tmp_vec);
        return;
    }
    for (int i = 0; i < size; ++i) {
        if (use[i] == false) {
            use[i] = true;
            tmp_vec.push_back(i);
            get_full_permutation(vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }
    }
}
