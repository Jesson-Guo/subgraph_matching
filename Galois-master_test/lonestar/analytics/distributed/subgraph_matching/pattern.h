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
    Pattern(int _size, bool clique = false);
    Pattern(int _size, char* buffer);
    ~Pattern();
    Pattern(const Pattern& p);
    void add_edge(int x, int y);
    void del_edge(int x, int y);
    inline void add_ordered_edge(int x, int y) { adj_mat[INDEX(x, y, size)] = 1;}
    inline int get_size() const {return size;}
    inline const int* get_adj_mat_ptr() const {return adj_mat;}
    bool check_connected() const;
    void count_all_isomorphism(std::set< std::set<int> >& s) const;
    void print() const;
    bool is_dag() const;

private:
    Pattern& operator =(const Pattern&);
    void get_full_permutation(std::vector< std::vector<int> >& vec, bool use[], std::vector<int> tmp_vec, int depth) const;
    int* adj_mat;
    int size;
};

Pattern::Pattern(int _size, bool clique)
{
    size = _size;
    adj_mat = new int[size * size];
    memset(adj_mat, 0, size * size * sizeof(int));

    if( clique ) {
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < i; ++j)
                add_edge(i,j);
    }
}

Pattern::Pattern(int _size, char *buffer) {
    size = _size;
    adj_mat = new int[size * size];
    memset(adj_mat, 0, size * size * sizeof(int));

    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            if(buffer[INDEX(i,j,size)] == '1')
                add_edge(i,j);
}

Pattern::~Pattern()
{
    delete[] adj_mat;
}

Pattern::Pattern(const Pattern& p)
{
    size = p.get_size();
    adj_mat = new int[size * size];
    memcpy(adj_mat, p.get_adj_mat_ptr(), size * size * sizeof(int));
}

void Pattern::add_edge(int x, int y)
{
    adj_mat[INDEX(x, y, size)] = 1;
    adj_mat[INDEX(y, x, size)] = 1;
}

void Pattern::del_edge(int x, int y)
{
    adj_mat[INDEX(x, y, size)] = 0;
    adj_mat[INDEX(y, x, size)] = 0;
}

bool Pattern::check_connected() const
{
    bool vis[size];
    for (int i = 0; i < size; ++i)
        vis[i] = false;
    int que[size];
    int head = 0;
    int tail = 0;
    que[0] = 0;
    vis[0] = true;
    //bfs
    while (head <= tail)
    {
        int x = que[head++];
        for (int i = 0; i < size; ++i)
            if (adj_mat[INDEX(x, i, size)] != 0 && vis[i] == false)
            {
                vis[i] = true;
                que[++tail] = i;
            }
    }
    for (int i = 0; i < size; ++i)
        if (vis[i] == false)
            return false;
    return true;
}


void Pattern::count_all_isomorphism(std::set< std::set<int> >& s) const
{
    unsigned int pow = 1;
    for (int i = 2; i <= size; ++i)
        pow *= i;
    std::vector< std::vector<int> > vec;
    vec.clear();
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<int> tmp_vec;
    get_full_permutation(vec, use, tmp_vec, 0);
    assert(vec.size() == pow);
    std::set<int> edge_set;
    for (const std::vector<int>& v : vec)
    {
        edge_set.clear();
        for (int i = 0; i < size; ++i)
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != 0)
                {
                    if (v[i] < v[j])
                        edge_set.insert(v[i] * size + v[j]);
                    else
                        edge_set.insert(v[j] * size + v[i]);
                }
                    
        if (s.count(edge_set) == 0)
            s.insert(edge_set);
    }
}

void Pattern::print() const
{
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (adj_mat[INDEX(i, j, size)] != 0)
                printf("(%d,%d) ", i, j);
    printf("\n");
}

bool Pattern::is_dag() const
{
    int degree[size];
    memset(degree, 0, size * sizeof(int));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (adj_mat[INDEX(i, j, size)] > 0)
                ++degree[j];
    int que[size];
    int head = 0;
    int tail = -1;
    for (int i = 0; i < size; ++i)
        if (degree[i] == 0)
            que[++tail] = i;
    while (head <= tail)
    {
        int x = que[head++];
        for (int j = 0; j < size; ++j)
            if (adj_mat[INDEX(x, j, size)] > 0)
            {
                --degree[j];
                if (degree[j] == 0)
                    que[++tail] = j;
            }
    }
    if (tail == size - 1)
        return true;
    else
        return false;
}

void Pattern::get_full_permutation(std::vector< std::vector<int> >& vec, bool use[], std::vector<int> tmp_vec, int depth) const
{
    if (depth == size)
    {
        vec.push_back(tmp_vec);
        return;
    }
    for (int i = 0; i < size; ++i)
        if (use[i] == false)
        {
            use[i] = true;
            tmp_vec.push_back(i);
            get_full_permutation(vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }
}
