#include "subgraph_matching_copy.h"
#include "llvm/Support/CommandLine.h"

#include <assert.h>
#include <iostream>
#include <string>
#include <algorithm>

namespace cll = llvm::cl;

static cll::opt<std::string> dataset(
    "dataset",
    cll::desc("dataset path"), 
    cll::Required);

static cll::opt<int> size(
    "size", 
    cll::desc("pattern size"), 
    cll::init(0)
);

static cll::opt<std::string> adj_mat(
    "mat", 
    cll::desc("pattern adjacent matrix"), 
    cll::Required
);

static cll::opt<long long> tri_cnt(
    "tricnt", 
    cll::desc("number of triangles of dataset"), 
    cll::init(0)
);

void test_pattern(Graph* g, const Pattern &pattern, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize = false) {
    bool is_pattern_valid;
    Schedule schedule(pattern, is_pattern_valid, performance_modeling_type, restricts_type, use_in_exclusion_optimize, g->v_cnt, g->e_cnt, g->tri_cnt);
    assert(is_pattern_valid);

    long long ans = g->pattern_matching(schedule);

    std::cout << "ans:\t" << ans << std::endl;
    schedule.print_schedule();
    for (auto p : schedule.restrict_pair) {
        std::cout << "(" << p.first << ", " << p.second << ")  ";
    }
    fflush(stdout);

}

int main(int argc,char *argv[]) {
    if(argc < 6) {
        printf("Usage: %s dataset_name graph_file pattern_size pattern_adjacency_matrix\n", argv[0]);
        printf("Example(Triangle counting on dataset WikiVote) : \n");
        printf("%s Wiki-Vote ../../dataset/wiki-vote_input 3 011101110 4\n", argv[0]);
        return 0;
    }

    Graph *g;
    DataLoader D;

    // comments in include/schedule.h explain the meaning of these parameters.
    int test_type = 1; // performance_modeling_type = restricts_type = use_in_exclusion_optimize = 1

    D.load_data(g, &dataset[0], 0, tri_cnt); 

    printf("Load data success!\n");
    fflush(stdout);

    Pattern p(size, &adj_mat[0]);
    test_pattern(g, p, test_type, test_type, test_type);
    
    delete g;
}
