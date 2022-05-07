#pragma once
#include "dataloader.h"
#include "classes.h"

#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/Reduction.h"
#include "galois/PriorityQueue.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "Lonestar/BoilerPlate.h"
#include "Lonestar/Utils.h"

bool DataLoader::load_data(Graph* &g, const char* path, int oriented_type, long long tri_cnt) {
    return general_load_data(g, path, oriented_type, tri_cnt);
}

bool DataLoader::general_load_data(Graph *&g, const char* path, int oriented_type, long long tri_cnt) {
    if (freopen(path, "r", stdin) == NULL)
    {
        printf("File not found. %s\n", path);
        return false;
    }
    printf("Load begin in %s\n",path);
    g = new Graph();

    //load triangle counting information
    g->tri_cnt = tri_cnt;

    scanf("%d%u",&g->v_cnt,&g->e_cnt);
    int* degree = new int[g->v_cnt];
    memset(degree, 0, g->v_cnt * sizeof(int));
    g->e_cnt *= 2;
    std::pair<int,int> *e = new std::pair<int,int>[g->e_cnt];
    id.clear();
    int x,y;
    int tmp_v;
    unsigned int tmp_e;
    tmp_v = 0;
    tmp_e = 0;
    while(scanf("%d%d",&x,&y)!=EOF) {
        if(x == y) {
            printf("find self circle\n");
            g->e_cnt -=2;
            continue;
            //return false;
        }
        if(!id.count(x)) id[x] = tmp_v ++;
        if(!id.count(y)) id[y] = tmp_v ++;
        x = id[x];
        y = id[y];
        e[tmp_e++] = std::make_pair(x,y);
        e[tmp_e++] = std::make_pair(y,x);
        ++degree[x];
        ++degree[y];
        //if(tmp_e % 1000000u == 0u) {
        //    printf("load %u edges\n",tmp_e);
        //    fflush(stdout);
        //}
    }

    // oriented_type == 0 do nothing
    //               == 1 high degree first
    //               == 2 low degree first
    if ( oriented_type != 0 ) {
        std::pair<int,int> *rank = new std::pair<int,int>[g->v_cnt];
        int *new_id = new int[g->v_cnt];
        for(int i = 0; i < g->v_cnt; ++i) rank[i] = std::make_pair(i,degree[i]);
        if( oriented_type == 1) std::sort(rank, rank + g->v_cnt, cmp_degree_gt);
        if( oriented_type == 2) std::sort(rank, rank + g->v_cnt, cmp_degree_lt);
        for(int i = 0; i < g->v_cnt; ++i) new_id[rank[i].first] = i;
        for(unsigned int i = 0; i < g->e_cnt; ++i) {
            e[i].first = new_id[e[i].first];
            e[i].second = new_id[e[i].second];
        }
        delete[] rank;
        delete[] new_id;
    }
    std::sort(degree, degree + g->v_cnt);

    // The max size of intersections is the second largest degree.
    //TODO VertexSet::max_intersection_size has different value with different dataset, but we use a static variable now.
    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, degree[g->v_cnt - 2]);
    g->max_degree = degree[g->v_cnt - 1];
    delete[] degree;
    if(tmp_v != g->v_cnt) {
        printf("vertex number error!\n");
    }
    if(tmp_e != g->e_cnt) {
        printf("edge number error!\n");
    }
    if(tmp_v != g->v_cnt || tmp_e != g->e_cnt) {
        fclose(stdin);
        delete g;
        delete[] e;
        return false;
    }
    std::sort(e,e+tmp_e,cmp_pair);
    g->e_cnt = unique(e,e+tmp_e) - e;
    for(unsigned int i = 0; i < g->e_cnt - 1; ++i)
        if(e[i] == e[i+1]) {
            printf("have same edge\n");
            fclose(stdin);
            delete g;
            delete[] e;
            return false;
        }
    g->edge = new int[g->e_cnt];
    g->vertex = new unsigned int[g->v_cnt + 1];
    bool* have_edge = new bool[g->v_cnt];
    int lst_v = -1;
    for(int i = 0; i < g->v_cnt; ++i) have_edge[i] = false;
    for(unsigned int i = 0; i < g->e_cnt; ++i) {
        if(e[i].first != lst_v) {
            have_edge[e[i].first] = true;
            g->vertex[e[i].first] = i;
        }
        lst_v = e[i].first;
        g->edge[i] = e[i].second;
    }
    delete[] e;
    printf("Success! There are %d nodes and %u edges.\n",g->v_cnt,g->e_cnt);
    fflush(stdout);
    g->vertex[g->v_cnt] = g->e_cnt;
    for(int i = g->v_cnt - 1; i >= 0; --i)
        if(!have_edge[i]) {
            g->vertex[i] = g->vertex[i+1];
        }
    delete[] have_edge;

    return true;
}

bool DataLoader::load_complete(Graph* &g, int clique_size) {
    g = new Graph();

    g->v_cnt = clique_size;
    g->e_cnt = clique_size * (clique_size - 1) / 2;

    int* degree = new int[g->v_cnt];
    memset(degree, 0, g->v_cnt * sizeof(int));
    g->e_cnt *= 2;
    std::pair<int,int> *e = new std::pair<int,int>[g->e_cnt];
    id.clear();
    int tmp_v;
    unsigned int tmp_e;
    tmp_v = 0;
    tmp_e = 0;
    for(int i = 0; i < clique_size; ++i)
        for(int j = 0; j < i; ++j) {
            int x = i, y = j;
            if(!id.count(x)) id[x] = tmp_v ++;
            if(!id.count(y)) id[y] = tmp_v ++;
            x = id[x];
            y = id[y];
            e[tmp_e++] = std::make_pair(x,y);
            e[tmp_e++] = std::make_pair(y,x);
            ++degree[x];
            ++degree[y];
        }

    std::sort(degree, degree + g->v_cnt);

    // The max size of intersections is the second largest degree.
    //TODO VertexSet::max_intersection_size has different value with different dataset, but we use a static variable now.
    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, degree[g->v_cnt - 2]);
    g->max_degree = degree[g->v_cnt - 1];
    delete[] degree;
    if(tmp_v != g->v_cnt) {
        printf("vertex number error!\n");
    }
    if(tmp_e != g->e_cnt) {
        printf("edge number error!\n");
    }
    if(tmp_v != g->v_cnt || tmp_e != g->e_cnt) {
        fclose(stdin);
        delete g;
        delete[] e;
        return false;
    }
    std::sort(e,e+tmp_e,cmp_pair);
    g->e_cnt = unique(e,e+tmp_e) - e;
    g->edge = new int[g->e_cnt];
    g->vertex = new unsigned int[g->v_cnt + 1];
    bool* have_edge = new bool[g->v_cnt];
    int lst_v = -1;
    for(int i = 0; i < g->v_cnt; ++i) have_edge[i] = false;
    for(unsigned int i = 0; i < g->e_cnt; ++i) {
        if(e[i].first != lst_v) {
            have_edge[e[i].first] = true;
            g->vertex[e[i].first] = i;
        }
        lst_v = e[i].first;
        g->edge[i] = e[i].second;
    }
    delete[] e;
    g->vertex[g->v_cnt] = g->e_cnt;
    for(int i = g->v_cnt - 1; i >= 0; --i)
        if(!have_edge[i]) {
            g->vertex[i] = g->vertex[i+1];
        }
    delete[] have_edge;
    return true;
}

bool DataLoader::cmp_pair(std::pair<int,int>a, std::pair<int,int>b) {
    return a.first < b.first || (a.first == b.first && a.second < b.second);
}

bool DataLoader::cmp_degree_gt(std::pair<int,int> a,std::pair<int,int> b) {
    return a.second > b.second;
}

bool DataLoader::cmp_degree_lt(std::pair<int,int> a,std::pair<int,int> b) {
    return a.second < b.second;
}

long long DataLoader::comb(int n, int k) {
    long long ans = 1;
    for(int i = n; i > n - k; --i)
        ans = ans * i;
    for(int i = 1; i <= k; ++i)
        ans = ans / k;
    return ans;
}


int VertexSet::max_intersection_size = -1;

VertexSet::VertexSet()
:data(nullptr), size(0), allocate(false)
{}

void VertexSet::init()
{
    if (allocate == true && data != nullptr)
        size = 0; // do not reallocate
    else
    {
        size = 0;
        allocate = true;
        data = new int[max_intersection_size];
    }
}

//this function is only used for tmp_set in graph.cpp (i.e., init(Graph.max_degree))
void VertexSet::init(int input_size)
{
    if (allocate == true && data != nullptr)
        size = 0;
    else
    {
        size = 0;
        allocate = true;
        data = new int[input_size];
    }
}

void VertexSet::init(int input_size, int* input_data)
{
    if (allocate == true && data != nullptr)
        delete[] data;
    size = input_size;
    data = input_data;
    allocate = false;
}

void VertexSet::copy(int input_size, const int* input_data)
{
    size = input_size;
    for(int i = 0; i < input_size; ++i) data[i] = input_data[i];
}

VertexSet::~VertexSet()
{
    if (allocate== true && data != nullptr)
        delete[] data;
}

void VertexSet::intersection(const VertexSet& set0, const VertexSet& set1, int min_vertex, bool clique)
{
    if (&set0 == &set1)
    {
        copy(set0.get_size(), set0.get_data_ptr());
        return;
    }
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    // TODO : Try more kinds of calculation.
    // Like
    // while (true)
    //     ..., if (++i == size0) break;
    //     ..., if (++j == size1) break;
    //     ......
    // Maybe we can also use binary search if one set is very small and another is large.
    if (clique == true)
        if (set0.get_data(0) >= min_vertex || set1.get_data(0) >= min_vertex)
            return;
    int data0 = set0.get_data(0);
    int data1 = set1.get_data(0);
    if (clique)
        // TODO : Try more kinds of calculation.
        // For example, we can use binary search find the last element which is smaller than min_vertex, and set its index as loop_size.
        while (i < size0 && j < size1)
        {
            if (data0 < data1)
            {
                if ((data0 = set0.get_data(++i)) >= min_vertex)
                    break;
            }
            else if (data0 > data1)
            {
                if ((data1 = set1.get_data(++j)) >= min_vertex)
                    break;
            }
            else
            {
                push_back(data0);
                if ((data0 = set0.get_data(++i)) >= min_vertex)
                    break;
                if ((data1 = set1.get_data(++j)) >= min_vertex)
                    break;
            }
        }
    else
        while (i < size0 && j < size1)
        {
            data0 = set0.get_data(i);
            data1 = set1.get_data(j);
            if (data0 < data1)
                ++i;
            else if (data0 > data1)
                ++j;
            else
            {
                push_back(data0);
                ++i;
                ++j;
            }
        }
}

void VertexSet::intersection_with(const VertexSet& set1) {
    if (this == &set1)
        return;
    const VertexSet& set0 = *this;
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    // TODO : Try more kinds of calculation.
    // Like
    // while (true)
    //     ..., if (++i == size0) break;
    //     ..., if (++j == size1) break;
    //     ......
    // Maybe we can also use binary search if one set is very small and another is large.
    int data0 = set0.get_data(0);
    int data1 = set1.get_data(0);
    size = 0;
    while (i < size0 && j < size1)
    {
        data0 = set0.get_data(i);
        data1 = set1.get_data(j);
        if (data0 < data1)
            ++i;
        else if (data0 > data1)
            ++j;
        else
        {
            push_back(data0);
            ++i;
            ++j;
        }
    }
}

void VertexSet::build_vertex_set(const Schedule& schedule, const VertexSet* vertex_set, int* input_data, int input_size, int prefix_id, int min_vertex, bool clique)
{
    int father_id = schedule.get_father_prefix_id(prefix_id);
    if (father_id == -1)
        init(input_size, input_data);
    else
    {
        init();
        VertexSet tmp_vset;
        tmp_vset.init(input_size, input_data);
        intersection(vertex_set[father_id], tmp_vset, min_vertex, clique);
    }
}

void VertexSet::insert_ans_sort(int val)
{
    int i;
    for (i = size - 1; i >= 0; --i)
        if (data[i] >= val)
            data[i + 1] = data[i];
        else
        {
            data[i + 1] = val;
            break;
        }
    if (i == -1)
        data[0] = val;
    ++size;
}

bool VertexSet::has_data(int val)
{
    for (int i = 0; i < size; ++i)
        if (data[i] == val)
            return true;
    return false;
}

int VertexSet::unorderd_subtraction_size(const VertexSet& set0, const VertexSet& set1, int size_after_restrict)
{
    int size0 = set0.get_size();
    int size1 = set1.get_size();
    if (size_after_restrict != -1)
        size0 = size_after_restrict;

    int ret = size0;
    const int* set0_ptr = set0.get_data_ptr();
    for (int j = 0; j < size1; ++j)
        if (std::binary_search(set0_ptr, set0_ptr + size0, set1.get_data(j)))
            --ret;
    return ret;
}


Schedule::Schedule(const Pattern& pattern, bool &is_pattern_valid, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize ,int v_cnt, unsigned int e_cnt, long long tri_cnt)
{
    if( performance_modeling_type != 0 && tri_cnt == -1) {
        printf("Fatal: Can not use performance modeling if not have triangle number of this dataset.\n");
        fflush(stdout);
        assert(0);
    }

    is_pattern_valid = true;
    size = pattern.get_size();
    adj_mat = new int[size * size];
    
    // not use performance_modeling, simply copy the adj_mat from pattern
    memcpy(adj_mat, pattern.get_adj_mat_ptr(), size * size * sizeof(int));

    std::vector< std::pair<int,int> > best_pairs;
    best_pairs.clear();
    //Initialize adj_mat
    //If we use performance_modeling, we may change the order of vertex,
    //the best order produced by performance_modeling(...) is saved in best_order[]
    //Finally, we use best_order[] to relocate adj_mat
    if( performance_modeling_type != 0) { 
        unsigned int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;
        
        std::vector< std::vector<int> > candidate_permutations;
        candidate_permutations.clear();
        
        bool use[size];
        for (int i = 0; i < size; ++i) use[i] = false;
        std::vector<int> tmp_vec;
        get_full_permutation(candidate_permutations, use, tmp_vec, 0);
        assert(candidate_permutations.size() == pow);

        remove_invalid_permutation(candidate_permutations);

        if(performance_modeling_type == 1) {
            //reduce candidates
            int max_val = 0;
            for(const auto &vec : candidate_permutations) {
                max_val = std::max(max_val, get_vec_optimize_num(vec));
            }
            std::vector< std::vector<int> > tmp;
            tmp.clear();
            for(const auto &vec : candidate_permutations) 
                if( get_vec_optimize_num(vec) == max_val) {
                    tmp.push_back(vec);
                }
            candidate_permutations = tmp;
        }

        int *best_order = new int[size];
        double min_val;
        bool have_best = false;
        

        for(const auto &vec : candidate_permutations) {
            int rank[size];
            for(int i = 0; i < size; ++i) rank[vec[i]] = i;
        
            int* cur_adj_mat;
            cur_adj_mat = new int[size*size];
            for(int i = 0; i < size; ++i)
                for(int j = 0; j < size; ++j)
                    cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

            std::vector< std::vector< std::pair<int,int> > > restricts_vector;

            restricts_vector.clear();

            if(restricts_type == 1) {
                restricts_generate(cur_adj_mat, restricts_vector);
            }
            else {
                Schedule schedule(cur_adj_mat, size);

                std::vector< std::pair<int,int> > pairs;
                schedule.GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }

            if( restricts_vector.size() == 0) {
                std::vector< std::pair<int,int> > Empty;
                Empty.clear();

                double val;
                if(performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if(performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                    else {
                        val = Naive_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                }

                if(have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for(int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = Empty;
                }
            }


            for(const auto& pairs : restricts_vector) {
                double val;
                if(performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if(performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                    else {
                        val = Naive_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                }

                if(have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for(int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = pairs;
                }
            }

        }

        int rank[size];
        for(int i = 0; i < size; ++i) rank[best_order[i]] = i;

        const int* pattern_adj_mat = pattern.get_adj_mat_ptr();
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                adj_mat[INDEX(rank[i], rank[j], size)] = pattern_adj_mat[INDEX(i, j, size)]; 
        delete[] best_order;
    }
    else {
        std::vector< int > I;
        I.clear();
        for(int i = 0; i < size; ++i) I.push_back(i);

        std::vector< std::vector< std::pair<int,int> > > restricts_vector;
        restricts_vector.clear();

        if(restricts_type != 0) {

            if(restricts_type == 1) {
                restricts_generate(adj_mat, restricts_vector);
            }
            else {
                std::vector< std::pair<int,int> > pairs;
                GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }
        }

        bool have_best = false;
        double min_val;

        for(const auto& pairs : restricts_vector) {
            double val;
            if(restricts_type == 1) {
                val = our_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt, tri_cnt);
            }
            else {
                val = GraphZero_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt);
            }
            if(have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                best_pairs = pairs;
            }
        }

    }

    if( use_in_exclusion_optimize) {
        std::vector<int> I;
        I.clear();
        for(int i = 0; i < size; ++i) I.push_back(i);
        in_exclusion_optimize_num = get_vec_optimize_num(I);
        if( in_exclusion_optimize_num <= 1) {
            printf("Can not use in_exclusion_optimize with this schedule\n");
            in_exclusion_optimize_num = 0;
        }
        else {
            printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
            init_in_exclusion_optimize();
        }
    }
    else {
            in_exclusion_optimize_num = 0;
    }

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;

    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i)
    {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
            {
                valid = true;
                break;
            }
        if (valid == false)
        {
            //Invalid Schedule
            is_pattern_valid = false;
            return;
        }
    }

    build_loop_invariant();
    if( restricts_type != 0) add_restrict(best_pairs);
    
    set_in_exclusion_optimize_redundancy();
}

Schedule::Schedule(const int* _adj_mat, int _size)
{
    size = _size;
    adj_mat = new int[size * size];

    memcpy(adj_mat, _adj_mat, size * size * sizeof(int));

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;
    in_exclusion_optimize_num = 0;

    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i)
    {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
            {
                valid = true;
                break;
            }
        if (valid == false)
        {
            printf("invalid schedule!\n");
            assert(0);
        }
    }

    build_loop_invariant();

    set_in_exclusion_optimize_redundancy();
}

Schedule::~Schedule()
{
    delete[] adj_mat;
    delete[] father_prefix_id;
    delete[] last;
    delete[] next;
    delete[] loop_set_prefix_id;
    delete[] prefix;
    delete[] restrict_last;
    delete[] restrict_next;
    delete[] restrict_index;
}

int Schedule::get_max_degree() const{
    int mx = 0;
    for(int i = 0; i < size; ++i) {
        int cnt = 0;
        for(int j = 0; j < size; ++j)
            cnt += adj_mat[INDEX(i,j,size)];
        if(cnt > mx) mx = cnt;
    }
    return mx;
}

void Schedule::build_loop_invariant()
{
    int* tmp_data = new int[size];
    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size; ++i)
    {
        int data_size = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
                tmp_data[data_size++] = j;
        loop_set_prefix_id[i] = find_father_prefix(data_size, tmp_data);
    }
    assert(total_prefix_num <= size * (size - 1) / 2);
    delete[] tmp_data;
}

int Schedule::find_father_prefix(int data_size, const int* data)
{
    if (data_size == 0)
        return -1;
    int num = data[data_size - 1];
    for (int prefix_id = last[num]; prefix_id != -1; prefix_id = next[prefix_id])
        if (prefix[prefix_id].equal(data_size, data))
            return prefix_id;
    
    // not found, create new prefix and find its father prefix id recursively
    int father = find_father_prefix(data_size - 1, data);
    father_prefix_id[total_prefix_num] = father;
    next[total_prefix_num] = last[num];
    last[num] = total_prefix_num;
    prefix[total_prefix_num].init(data_size, data);
    ++total_prefix_num;
    return total_prefix_num - 1;
}

void Schedule::add_restrict(const std::vector< std::pair<int, int> >& restricts)
{
    restrict_pair = restricts;
    for(unsigned int i = 0; i < restrict_pair.size(); ) {
        bool tag = true;
        for(unsigned int j = 0; j < restrict_pair.size(); ++j) {
            if(i != j && restrict_pair[j].first == restrict_pair[i].first) 
                for(unsigned int k = 0; k < restrict_pair.size(); ++k)
                    if( i != k && j != k && restrict_pair[k].second == restrict_pair[i].second && restrict_pair[j].second == restrict_pair[k].first ) {
                        tag = false;
                        break;
                    }
            if(tag == false) break;
        }
        if(tag == false) {
            restrict_pair.erase(restrict_pair.begin() + i);
        }
        else ++i;
    }



    int max_prefix_num = size * (size - 1) / 2;
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    total_restrict_num = 0;
    for (const auto& p : restrict_pair)
    {
        // p.first must be greater than p.second
        restrict_index[total_restrict_num] = p.first;
        restrict_next[total_restrict_num] = restrict_last[p.second];
        restrict_last[p.second] = total_restrict_num;
        ++total_restrict_num;
    }
}

int Schedule::get_multiplicity() const{
    std::vector< std::vector<int> > isomorphism_vec = get_isomorphism_vec();
    return isomorphism_vec.size();
}

void Schedule::aggressive_optimize(std::vector< std::pair<int, int> >& ordered_pairs) const
{
    std::vector< std::vector<int> > isomorphism_vec = get_isomorphism_vec();

    std::vector< std::vector< std::vector<int> > > permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int>& v : isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int>& v : permutation_groups[i])
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
            }
            else if (v.size() != 1)
            {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1)
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        }
        else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int>& pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    bool changed = true;
    while (changed && isomorphism_vec.size() != 1)
    {
        // use restrictions to delete other isomophism
        for (unsigned int i = 0; i < isomorphism_vec.size(); )
        {
            Pattern test_dag(base_dag);
            const std::vector<int>& iso = isomorphism_vec[i];
            for (const std::pair<int, int>& pair : ordered_pairs)
                test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
            if (test_dag.is_dag() == false) // is not dag means conflict
            {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
            }
            else
                ++i;
        }

        changed = false;
        std::pair<int, int> found_pair;
        for (unsigned int i = 0; i < permutation_groups.size(); )
        {
            int two_element_number = 0;
            for (const std::vector<int>& v : permutation_groups[i])
                if (v.size() == 2)
                {
                    ++two_element_number;
                    found_pair = std::pair<int ,int>(v[0], v[1]);
                    break;
                }
            if (two_element_number >= 1)
            {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                ordered_pairs.push_back(found_pair);
                base_dag.add_ordered_edge(found_pair.first, found_pair.second);
                changed = true;
                break;
            }
            else
                ++i;
        }
    }
    assert(isomorphism_vec.size() == 1);
}

// Schedule::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT
void Schedule::aggressive_optimize_get_all_pairs(std::vector< std::vector< std::pair<int, int> > >& ordered_pairs_vector) 
{
    std::vector< std::vector<int> > isomorphism_vec = get_isomorphism_vec();

    std::vector< std::vector< std::vector<int> > > permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int>& v : isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs_vector.clear();

    std::vector< std::pair<int,int> > ordered_pairs;
    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int>& v : permutation_groups[i])
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
            }
            else if (v.size() != 1)
            {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1)
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        }
        else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int>& pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    aggressive_optimize_dfs(base_dag, isomorphism_vec, permutation_groups, ordered_pairs, ordered_pairs_vector);

}

void Schedule::aggressive_optimize_dfs(Pattern base_dag, std::vector< std::vector<int> > isomorphism_vec, std::vector< std::vector< std::vector<int> > > permutation_groups, std::vector< std::pair<int,int> > ordered_pairs, std::vector< std::vector< std::pair<int,int> > >& ordered_pairs_vector) {

    for (unsigned int i = 0; i < isomorphism_vec.size(); )
    {
        Pattern test_dag(base_dag);
        const std::vector<int>& iso = isomorphism_vec[i];
        for (const std::pair<int, int>& pair : ordered_pairs)
            test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
        if (test_dag.is_dag() == false) // is not dag means conflict
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
        }
        else
            ++i;
    }
    
    if(isomorphism_vec.size() == 1) {
        ordered_pairs_vector.push_back(ordered_pairs);
        return;
    }


    std::pair<int, int> found_pair;
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        for (const std::vector<int>& v : permutation_groups[i])
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
                std::vector< std::vector< std::vector<int> > > next_permutation_groups = permutation_groups;
                std::vector< std::vector<int> > next_isomorphism_vec = isomorphism_vec;
                std::vector< std::pair<int,int> > next_ordered_pairs = ordered_pairs;
                Pattern next_base_dag = base_dag;
                
                next_permutation_groups.erase(next_permutation_groups.begin() + i);
                next_isomorphism_vec.erase(next_isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                next_ordered_pairs.push_back(found_pair);
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);
                
                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec, next_permutation_groups, next_ordered_pairs, ordered_pairs_vector);
            }
        if( two_element_number >= 1) {
            break;
        }
        else {
           ++i;
        }
    }

}

void Schedule::GraphZero_aggressive_optimize(std::vector< std::pair<int, int> >& ordered_pairs) const { 
    std::vector< std::vector<int> > Aut;
    GraphZero_get_automorphisms(Aut);

    std::vector< std::pair<int,int> > L;
    L.clear();

    for(int v = 0; v < size; ++v) { // iterate all elements in schedule
        std::vector< std::vector<int> > stabilized_aut;
        stabilized_aut.clear();

        for(int i = 0; i < (int)Aut.size(); ++i) {
            std::vector<int>& x = Aut[i];
            if( x[v] == v) {
                stabilized_aut.push_back(x);
            }
            else {
                int x1  = v, x2 = x[v];
                if( x1 > x2) {
                    int tmp = x1;
                    x1 = x2;
                    x2 = tmp;
                }
                bool tag = true;
                std::pair<int,int> cur = std::make_pair(x1, x2);
                for(int j = 0; j < (int)L.size(); ++j)
                    if( L[j]  == cur) {
                        tag = false;
                        break;
                    }
                if(tag) {
                    L.push_back(cur);
                }
            }
        }
        Aut = stabilized_aut;
    }
    
    ordered_pairs.clear(); // In GraphZero paper, this vector's name is 'L'

    for(int i = 0; i < (int)L.size(); ++i) {
        bool tag = true;
        for(int j = 0; j < (int)ordered_pairs.size(); ++j)
            if( L[i].second == ordered_pairs[j].second) {
                tag = false;
                if( L[i].first > ordered_pairs[j].first) ordered_pairs[j].first = L[i].first;
                break;
            }
        if(tag) ordered_pairs.push_back(L[i]);
    }
}

void Schedule::GraphZero_get_automorphisms(std::vector< std::vector<int> > &Aut) const {
    int p[size];
    Aut.clear();
    for(int i = 0; i < size; ++i) p[i] = i;
    do{
        bool tag = true;
        for(int i = 0; i < size; ++i) {
            for(int j = 0; j < size; ++j)
                if( adj_mat[INDEX(i, j, size)] != adj_mat[INDEX(p[i], p[j], size)]) {
                    tag = false;
                    break;
                }
            if( !tag ) break;
        }
        if(tag) {
            std::vector<int> tmp;
            tmp.clear();
            for(int i = 0; i < size; ++i) tmp.push_back(p[i]);
            Aut.push_back(tmp);
        }
    } while( std::next_permutation(p, p + size) );

}

std::vector< std::vector<int> > Schedule::get_isomorphism_vec() const
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
    std::vector< std::vector<int> > isomorphism_vec;
    isomorphism_vec.clear();
    for (const std::vector<int>& v : vec)
    {
        bool flag = true;
        for (int i = 0; i < size; ++i)
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != 0)
                    if (adj_mat[INDEX(v[i], v[j], size)] == 0) // not isomorphism
                    {
                        flag = false;
                        break;
                    }
        if (flag == true)
            isomorphism_vec.push_back(v);
    }
    return isomorphism_vec;
}

void Schedule::get_full_permutation(std::vector< std::vector<int> >& vec, bool use[], std::vector<int> tmp_vec, int depth) const
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

std::vector< std::vector<int> > Schedule::calc_permutation_group(const std::vector<int> vec, int size)
{
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector< std::vector<int> > res;
    res.clear();
    for (unsigned int i = 0; i < vec.size(); ++i)
        if (use[i] == false)
        {
            std::vector<int> tmp_vec;
            tmp_vec.clear();
            tmp_vec.push_back(i);
            use[i] = true;
            int x = vec[i];
            while (use[x] == false)
            {
                use[x] = true;
                tmp_vec.push_back(x);
                x = vec[x];
            }
            res.push_back(tmp_vec);
        }
    return res;
}

void Schedule::performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt) {
    int* order;
    int* rank;

    double* p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    order = new int[size];
    rank = new int[size];
    
    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for(const std::vector<int>& vec : candidates) {
        for(int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;
        
        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::pair<int,int> > restricts;
        //TODO BUG!!!!!
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double* sum;
        sum = new double[restricts_size];
        for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
        int* tmp;
        tmp = new int[size];
        for(int i = 0; i < size; ++i) tmp[i] = i;
        do {
            for(int i = 0; i < restricts_size; ++i)
                if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                }
                else break;
        } while( std::next_permutation(tmp, tmp + size));
        double total = 1;
        for(int i = 2; i <= size; ++i) total *= i;
        for(int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] /total;
        for(int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for(int i = 0; i < size; ++i) invariant_size[i].clear();
        for(int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for(int j = 0; j < i; ++j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for(int j = i + 1; j < size; ++j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for(int j = i - 1; j >= 0; --j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for(int j = 0; j < (int)invariant_size[i].size(); ++j)
                if(invariant_size[i][j] > 1) 
                    val += p_size[invariant_size[i][j] - 1] + p_size[1];
            for(int j = 0; j < restricts_size; ++j)
                if(restricts[j].second == i)
                    val *=  sum[j];
            val *= p_size[cnt_forward];

        }
        if( have_best == false || val < min_val) {
            have_best = true;
            for(int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
        }
        delete[] sum;
        delete[] tmp;
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule::bug_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt) {
    int* order;
    int* rank;

    double* p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for(const std::vector<int>& vec : candidates) {
        for(int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;
        
        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::vector< std::pair<int,int> > > restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for(int restricts_rank = 0; restricts_rank < (int)restricts_vector.size(); ++restricts_rank) {
            std::vector< std::pair<int,int> >& restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double* sum;
            sum = new double[restricts_size];
            for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
            int* tmp;
            tmp = new int[size];
            for(int i = 0; i < size; ++i) tmp[i] = i;
            do {
                for(int i = 0; i < restricts_size; ++i)
                    if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    }
                    else break;
            } while( std::next_permutation(tmp, tmp + size));
            double total = 1;
            for(int i = 2; i <= size; ++i) total *= i;
            for(int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] /total;
            for(int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for(int i = 0; i < size; ++i) invariant_size[i].clear();
            for(int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for(int j = 0; j < i; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for(int j = i + 1; j < size; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for(int j = i - 1; j >= 0; --j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for(int j = 0; j < (int)invariant_size[i].size(); ++j)
                    if(invariant_size[i][j] > 1) 
                        val += p_size[invariant_size[i][j] - 1] + p_size[1];
                for(int j = 0; j < restricts_size; ++j)
                    if(restricts[j].second == i)
                        val *=  sum[j];
                val *= p_size[cnt_forward];

            }
            if( have_best == false || val < min_val) {
                have_best = true;
                for(int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule::new_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    int* order;
    int* rank;

    double* p_size;
    double* pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt; 
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    order = new int[size];
    rank = new int[size];
    
    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for(const std::vector<int>& vec : candidates) {
        for(int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;
        
        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::vector< std::pair<int,int> > > restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for(int restricts_rank = 0; restricts_rank < (int)restricts_vector.size(); ++restricts_rank) {
            std::vector< std::pair<int,int> >& restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double* sum;
            sum = new double[restricts_size];
            for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
            int* tmp;
            tmp = new int[size];
            for(int i = 0; i < size; ++i) tmp[i] = i;
            do {
                for(int i = 0; i < restricts_size; ++i)
                    if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    }
                    else break;
            } while( std::next_permutation(tmp, tmp + size));
            double total = 1;
            for(int i = 2; i <= size; ++i) total *= i;
            for(int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] /total;
            for(int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for(int i = 0; i < size; ++i) invariant_size[i].clear();
            for(int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for(int j = 0; j < i; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for(int j = i + 1; j < size; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for(int j = i - 1; j >= 0; --j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for(int j = 0; j < (int)invariant_size[i].size(); ++j)
                    if(invariant_size[i][j] > 1) 
                        val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
                val += 1;
                for(int j = 0; j < restricts_size; ++j)
                    if(restricts[j].second == i)
                        val *=  sum[j];
                val *= p_size[1] * pp_size[ cnt_forward - 1 ];

            }
            if( have_best == false || val < min_val) {
                have_best = true;
                for(int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] pp_size;
}

void Schedule::init_in_exclusion_optimize() {
    int optimize_num = in_exclusion_optimize_num;
    
    assert( in_exclusion_optimize_num > 1);

    int* id;
    id = new int[ optimize_num ];

    int* in_exclusion_val;
    in_exclusion_val = new int[ optimize_num * 2];

    for(int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;

        in_exclusion_val[ 2 * n - 2 ] = 0;
        in_exclusion_val[ 2 * n - 1 ] = 0;

        if( n == 1) {
            ++in_exclusion_val[0];
            continue;
        }

        std::pair<int,int> edge[m];
        int e_cnt = 0;
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < i; ++j)
                edge[e_cnt++] = std::make_pair(i,j);

        for(int s = 0; s < (1<<m); ++s) {
            dsu.init();
            int bit_cnt = 0;
            for(int i = 0; i < m; ++i) 
                if( s & (1<<i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);
                }
            if( dsu.get_set_size() == 1) {
                if( bit_cnt & 1) ++in_exclusion_val[2 * n -1];
                else ++in_exclusion_val[ 2 * n - 2];
            }
        }
    }        

    in_exclusion_optimize_group.clear();
    in_exclusion_optimize_val.clear();

    get_in_exclusion_optimize_group(0, id, 0, in_exclusion_val);

    delete[] id;
    delete[] in_exclusion_val;
}

void Schedule::get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val) {
    if( depth == in_exclusion_optimize_num) {
        int* size = new int[id_cnt];
        for(int i = 0; i < id_cnt; ++i)
            size[i] = 0;
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            size[ id[i] ] ++;
        int val[2];
        val[0] = in_exclusion_val[ size[0] * 2 - 2 ];
        val[1] = in_exclusion_val[ size[0] * 2 - 1 ];
        for(int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];

            val[0] = tmp0 * in_exclusion_val[ size[i] * 2 - 2] + tmp1 * in_exclusion_val[ size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[ size[i] * 2 - 1] + tmp1 * in_exclusion_val[ size[i] * 2 - 2];
        }

        std::vector< std::vector<int> > group;
        group.clear();
        for(int i = 0; i < id_cnt; ++i) {
            std::vector<int> cur;
            cur.clear();
            for(int j = 0; j < in_exclusion_optimize_num; ++j)
                if( id[j] == i) cur.push_back(j);
            group.push_back(cur);
        }

        in_exclusion_optimize_group.push_back(group);
        in_exclusion_optimize_val.push_back( val[0] - val[1] );

        delete[] size;
        return;
    }
    
    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);
    
    for(int i = 0; i < id_cnt; ++i) {
        id[depth] = i;
        get_in_exclusion_optimize_group(depth + 1, id, id_cnt, in_exclusion_val);
    }
}

void Schedule::print_schedule() const{
    printf("Schedule:\n");
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j)
            printf("%d", adj_mat[INDEX(i,j,size)]);
        puts("");
    }
}

void Schedule::GraphZero_performance_modeling(int* best_order, int v_cnt, unsigned int e_cnt) {
    int* order;
    int* rank;

    double* p_size;
    double* anti_p;
    p_size = new double[size];
    anti_p = new double[size];

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    printf("fuck p %.6lf\n",p);
    p_size[0] = v_cnt;
    for(int i = 1; i < size; ++i) {
        p_size[i] = p_size[i-1] * p;
        printf("p %d %.6lf\n", i, p_size[i]);
    }
    anti_p[0] = 1;
    for(int i = 1; i < size; ++i) {
        anti_p[i] = anti_p[i-1] * (1-p);
        printf("anti %d %.6lf\n", i, anti_p[i]);
    }

    order = new int[size];
    rank = new int[size];
    
    for(int i = 0; i < size; ++i) order[i] = i;
    double min_val;
    bool have_best = false;
    do {
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;
        
        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::pair<int,int> > restricts;
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double* sum;
        sum = new double[restricts_size];
        for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
        int* tmp;
        tmp = new int[size];
        for(int i = 0; i < size; ++i) tmp[i] = i;
        do {
            for(int i = 0; i < restricts_size; ++i)
                if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                }
                else break;
        } while( std::next_permutation(tmp, tmp + size));
        double total = 1;
        for(int i = 2; i <= size; ++i) total *= i;
        for(int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] /total;
        for(int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for(int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            for(int j = 0; j < i; ++j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;

            for(int j = 0; j < restricts_size; ++j)
                if(restricts[j].second == i)
                    val *=  sum[j];
      //      val *= p_size[cnt_forward + 1] * anti_p[i - cnt_forward];
            val *= v_cnt;
            for(int j = 0; j < i - cnt_forward; ++j)
                val *= (1-p);
            for(int j = 0; j < cnt_forward; ++j)
                val *= p;
        }
        if( have_best == false || val <= min_val) {
            have_best = true;
            for(int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
            printf("gz upd %.10lf\n", val);
        }
        
        delete[] cur_adj_mat;
        delete[] sum;
        delete[] tmp;

    } while( std::next_permutation(order, order + size) );

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] anti_p;
}

void Schedule::restrict_selection(int v_cnt, unsigned int e_cnt, long long tri_cnt, std::vector< std::vector< std::pair<int,int> > > ordered_pairs_vector, std::vector< std::pair<int,int> > &best_restricts) const{
    
    double* p_size;
    double* pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];

    for(int cur_restricts = 0; cur_restricts < (int)ordered_pairs_vector.size(); ++cur_restricts) {
        std::vector< std::pair<int,int> > &restricts = ordered_pairs_vector[cur_restricts];
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());

        double* sum;
        sum = new double[restricts_size];
        for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
        
        int* tmp;
        tmp = new int[size];
        for(int i = 0; i < size; ++i) tmp[i] = i;
        
        do {
            for(int i = 0; i < restricts_size; ++i)
                if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                }
                else break;
        } while( std::next_permutation(tmp, tmp + size));
        
        double total = 1;
        for(int i = 2; i <= size; ++i) total *= i;
        
        for(int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] /total;
        for(int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for(int i = 0; i < size; ++i) invariant_size[i].clear();
        
        for(int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for(int j = 0; j < i; ++j)
                if(adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for(int j = i + 1; j < size; ++j)
                if(adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for(int j = i - 1; j >= 0; --j)
                if(adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for(int j = 0; j < (int)invariant_size[i].size(); ++j)
                if(invariant_size[i][j] > 1) 
                    val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
            val += 1;
            for(int j = 0; j < restricts_size; ++j)
                if(restricts[j].second == i)
                    val *=  sum[j];
            val *= p_size[1] * pp_size[cnt_forward-1];
        
        }
        if( have_best == false || val < min_val) {
            have_best = true;
            best_restricts = restricts;
            min_val = val;
        }
        
        delete[] sum;
        delete[] tmp;

    }

    delete[] p_size;
    delete[] pp_size;
    assert(have_best);
}

void Schedule::restricts_generate(const int* cur_adj_mat, std::vector< std::vector< std::pair<int,int> > > &restricts) {
    Schedule schedule(cur_adj_mat, get_size());
    schedule.aggressive_optimize_get_all_pairs(restricts);
    int size = schedule.get_size();
    Graph* complete;
    DataLoader* D = new DataLoader();
    assert(D->load_complete(complete, size + 1));
    long long ans = complete->pattern_matching(schedule) / schedule.get_multiplicity();
    for(int i = 0; i < (int)restricts.size(); ) {
        Schedule cur_schedule(schedule.get_adj_mat_ptr(), schedule.get_size());
        cur_schedule.add_restrict(restricts[i]);
        long long cur_ans = complete->pattern_matching( cur_schedule);
        if( cur_ans != ans) {
            restricts.erase(restricts.begin() + i);
        }
        else {
            ++i;
        }
    }

    delete complete;
    delete D;
}

int Schedule::get_vec_optimize_num(const std::vector<int> &vec) {
    bool is_valid = true;
    for(int i = 1; i < size; ++i) {
        bool have_edge = false;
        for(int j = 0; j < i; ++j)
            if( adj_mat[INDEX(vec[i], vec[j], size)]) {
                have_edge = true;
                break;
            }
        if( have_edge == false) {
            is_valid = false;
            break;
        }
    }
    if( !is_valid) return -1;

    for(int k = 2; k <= size; ++k) {
        bool flag = true;
        for(int i = size - k + 1; i < size; ++i)
            if(adj_mat[INDEX(vec[size - k], vec[i], size)]) {
                flag = false;
                break;
            }
        if(flag == false) return k - 1;
    }
    // assert(0);
    return -1;
}

double Schedule::our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double pp_size[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt; 
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    
    int tmp[size];
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for(int i = 0; i < size; ++i) invariant_size[i].clear();
    
    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for(int j = i + 1; j < size; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for(int j = i - 1; j >= 0; --j)
            if(cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for(int j = 0; j < (int)invariant_size[i].size(); ++j)
            if(invariant_size[i][j] > 1) 
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        val += 1;
        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        if( i ) {
            val *= p_size[1] * pp_size[ cnt_forward - 1 ];
        }
        else {
            val *= p_size[0];
        }
    }
    delete[] cur_adj_mat;

    return val;
}

double Schedule::GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt) {
    int max_degree = get_max_degree();
    
    double p_size[max_degree];
    double p = e_cnt * 1.0 / v_cnt / v_cnt;
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());
    
    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    
    int tmp[size];
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for(int i = 0; i < size; ++i) invariant_size[i].clear();
    
    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for(int j = i + 1; j < size; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for(int j = i - 1; j >= 0; --j)
            if(cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for(int j = 0; j < (int)invariant_size[i].size(); ++j)
            if(invariant_size[i][j] > 1) 
                val += p_size[invariant_size[i][j] - 1] + p_size[1];
        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        val *= p_size[cnt_forward];

    }
    
    delete[] cur_adj_mat;

    return val;
}

double Schedule::Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt) {

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];
    
    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());
    
    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    int tmp[size];
    
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        val *= v_cnt;
        for(int j = 0; j < i - cnt_forward; ++j)
            val *= (1-p);
        for(int j = 0; j < cnt_forward; ++j)
            val *= p;
    }

    delete[] cur_adj_mat;

    return val;
}

// the i's vertex must connect to at lease one of i-1 vertices at front
void Schedule::remove_invalid_permutation(std::vector< std::vector<int> > &candidate_permutations) {
    for(unsigned int i = 0; i < candidate_permutations.size(); ) {
        const auto& vec = candidate_permutations[i];
        bool tag = true;
        for(int x = 1; x < size; ++x) {
            bool have_edge = false;
            for(int y = 0; y < x; ++y)
                if(adj_mat[INDEX(vec[x],vec[y],size)]) {
                    have_edge = true;
                    break;
                }
            if(!have_edge) {
                tag = false;
                break;
            }
        }
        if(tag) {
            ++i;
        }
        else {
            candidate_permutations.erase(candidate_permutations.begin() + i);
        }
    }
}

int Schedule::get_in_exclusion_optimize_num_when_not_optimize() {
    std::vector<int> I;
    for(int i = 0; i < size; ++i) I.push_back(i);
    return get_vec_optimize_num(I);
}

void Schedule::set_in_exclusion_optimize_redundancy() {
    int tmp = get_in_exclusion_optimize_num();
    if(tmp <= 1) {
        in_exclusion_optimize_redundancy = 1;
    }
    else {
        Graph* complete;
        DataLoader* D = new DataLoader();
        assert(D->load_complete(complete, get_size()));
        delete D;
        in_exclusion_optimize_redundancy = 1;
        long long ans = complete->pattern_matching( *this);
        set_in_exclusion_optimize_num(0);
        long long true_ans = complete->pattern_matching( *this);
        set_in_exclusion_optimize_num(tmp);
        delete complete;
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
}


constexpr static const unsigned CHUNK_SIZE = 128U;

int Graph::intersection_size(int v1,int v2) {
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    unsigned int l2, r2;
    get_edge_index(v2, l2, r2);
    int ans = 0;
    while(l1 < r1 && l2 < r2) {
        if(edge[l1] < edge[l2]) {
            ++l1;
        }
        else {
            if(edge[l2] < edge[l1]) {
                ++l2;
            }
            else {
                // v1 and v2 have the same neighbor
                ++l1;
                ++l2;
                ++ans;
            }
        }
    }
    return ans;
}

int Graph::intersection_size_clique(int v1,int v2) {
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    unsigned int l2, r2;
    get_edge_index(v2, l2, r2);
    int min_vertex = v2;
    int ans = 0;
    if (edge[l1] >= min_vertex || edge[l2] >= min_vertex)
        return 0;
    while(l1 < r1 && l2 < r2) {
        if(edge[l1] < edge[l2]) {
            if (edge[++l1] >= min_vertex)
                break;
        }
        else {
            if(edge[l2] < edge[l1]) {
                if (edge[++l2] >= min_vertex)
                    break;
            }
            else {
                ++ans;
                if (edge[++l1] >= min_vertex)
                    break;
                if (edge[++l2] >= min_vertex)
                    break;
            }
        }
    }
    return ans;
}

void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = vertex[v];
    r = vertex[v + 1];
}

long long Graph::pattern_matching(const Schedule& schedule) {
    long long global_ans = 0;
    std::vector<int> allNodes;
    for (int vertex = 0; vertex < v_cnt; ++vertex) allNodes.push_back(vertex);
    VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
    VertexSet subtraction_set;
    VertexSet tmp_set;
    subtraction_set.init();
    long long local_ans = 0;
    galois::do_all(
        galois::iterate(allNodes),
        [&](uint64_t vertex) {
            unsigned int l, r;
            get_edge_index(vertex, l, r);
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id);
            }
            subtraction_set.push_back(vertex);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
            subtraction_set.pop_back();
        },
        galois::loopname("matching"),
        galois::chunk_size<CHUNK_SIZE>(),
        galois::steal(),
        galois::no_stats()
    );
    delete[] vertex_set;
    // TODO : Computing multiplicty for a pattern
    global_ans += local_ans;

    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

void Graph::pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth) // 3 same # or @ in comment are useful in code generation ###
{
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);// @@@
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;

    int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    //Case: in_exclusion_optimize_num > 1
    if( depth == schedule.get_size() - schedule.get_in_exclusion_optimize_num() ) {
        int in_exclusion_optimize_num = schedule.get_in_exclusion_optimize_num();// @@@
        int loop_set_prefix_ids[ in_exclusion_optimize_num ];
        loop_set_prefix_ids[0] = loop_set_prefix_id;
        for(int i = 1; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id( depth + i );
        for(int optimize_rank = 0; optimize_rank < (int)schedule.in_exclusion_optimize_group.size(); ++optimize_rank) {
            const std::vector< std::vector<int> >& cur_graph = schedule.in_exclusion_optimize_group[optimize_rank];
            long long val = schedule.in_exclusion_optimize_val[optimize_rank];
            for(int cur_graph_rank = 0; cur_graph_rank < (int)cur_graph.size(); ++ cur_graph_rank) {
                //                VertexSet tmp_set;
                
                //if size == 1 , we will not call intersection(...)
                //so we will not allocate memory for data
                //otherwise, we need to copy the data to do intersection(...)
                if(cur_graph[cur_graph_rank].size() == 1) {
                    int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    val = val * VertexSet::unorderd_subtraction_size(vertex_set[id], subtraction_set);
                }
                else {
                    int id0 = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    int id1 = loop_set_prefix_ids[cur_graph[cur_graph_rank][1]];
                    tmp_set.init(this->max_degree);
                    tmp_set.intersection(vertex_set[id0], vertex_set[id1]);

                    for(int i = 2; i < (int)cur_graph[cur_graph_rank].size(); ++i) {
                        int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][i]];
                        tmp_set.intersection_with(vertex_set[id]);
                    }
                    val = val * VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
                }
                if( val == 0 ) break;

            }
            local_ans += val;
        }
        return;// @@@
            
    }
    //Case: in_exclusion_optimize_num <= 1
    if (depth == schedule.get_size() - 1)
    {
        // TODO : try more kinds of calculation. @@@
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0)
        {
            int min_vertex = v_cnt;
            for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
                if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
            const VertexSet& vset = vertex_set[loop_set_prefix_id];
            int size_after_restrict = std::lower_bound(vset.get_data_ptr(), vset.get_data_ptr() + vset.get_size(), min_vertex) - vset.get_data_ptr();
            if (size_after_restrict > 0)
                local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set, size_after_restrict);
        }
        else
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set); 
        return;// @@@
    }
  
    // TODO : min_vertex is also a loop invariant @@@
    int min_vertex = v_cnt;
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
    int ii = 0;
    for (int &i = ii; i < loop_size; ++i)
    {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int vertex = loop_data_ptr[i];
        if (subtraction_set.has_data(vertex))
            continue;
        unsigned int l, r;
        get_edge_index(vertex, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, vertex);
            if( vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex);
        pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
        subtraction_set.pop_back(); // @@@
    }
}
