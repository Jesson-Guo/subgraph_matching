#include "subgraph_matching.h"

#include "llvm/Support/CommandLine.h"

namespace cll = llvm::cl;

static const char* name = "subgraph matching";
static const char* desc = "Computes the number of subgraphs matching the pattern";
static const char* url = "subgraph_matching";

enum Algo {
    Base = 0,
};

const char* const ALGO_NAMES[] = {"base"};

static cll::opt<std::string> inputFile(
    "dataset",
    cll::desc("dataset path"), 
    cll::Required);

static cll::opt<std::string> completeFile(
    "complete", 
    cll::desc("complete graph"), 
    cll::Required
);

static cll::opt<Algo> algo(
    "algo", 
    cll::desc("Choose an algorithm (default value base):"),
    cll::values(clEnumVal(Base, "base")),
    cll::init(Base)
);

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

static int max_degree(0);

void baseMatching(Graph& graph, const Pattern p, Graph& complete, const long long tri_cnt) {
    int max_deg(0), sec_max_deg(0);
    galois::do_all(
        galois::iterate(complete),
        [&](uint64_t i) {
            complete.getData(i).degree = complete.size() - 1;
            for (int v=0; v<(int)complete.size(); v++) {
                if (v != (int)i) {
                    complete.getData(i).adj_node.push_back(v);
                }
            }
        },
        galois::loopname("complete"),
        galois::chunk_size<CHUNK_SIZE>(),
        galois::steal(),
        galois::no_stats()
    );
    galois::do_all(
        galois::iterate(graph),
        [&](uint64_t i) {
            int deg(0);
            for (auto e : graph.edges(i)) {
                GNode dst = graph.getEdgeDst(e);
                graph.getData(i).adj_node.push_back(dst);
                deg++;
            }
            std::sort(graph.getData(i).adj_node.begin(), graph.getData(i).adj_node.end());
            graph.getData(i).degree = deg;
            if (deg > max_deg) {
                sec_max_deg = max_deg;
                max_deg = deg;
            }
            else if (deg > sec_max_deg)
                sec_max_deg = deg;
        },
        galois::loopname("max intersection size calculate"),
        galois::chunk_size<CHUNK_SIZE>(),
        galois::steal(),
        galois::no_stats()
    );
    max_degree = max_deg;
    VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, sec_max_deg);
    bool is_pattern_valid;
    int perf_type(1), res_type(1);
    bool opt(true);
    Schedule schedule(p, is_pattern_valid, perf_type, res_type, opt, complete, max_deg, graph.size(), graph.sizeEdges(), tri_cnt);
    assert(is_pattern_valid);

    galois::StatTimer matchingTime("Timer_matching");
    matchingTime.start();
    uint64_t ans = pattern_matching(graph, schedule, max_deg);
    matchingTime.stop();
    std::cout << "ans:\t" << ans << std::endl;
    std::cout << "time:\t" << matchingTime.get_usec() << std::endl;
    std::cout << "restricts:\n";
    // for (auto p : schedule.restrict_pair) {
    for (int i=0; i<schedule.rest_size; i++) {
        std::cout << "(" << schedule.rest_first[i] << ", " << schedule.rest_second[i] << ")  ";
    }
    std::cout << std::endl;
}

int main(int argc, char** argv) {
    galois::SharedMemSys G;
    LonestarStart(argc, argv, name, desc, url, &inputFile);

    galois::StatTimer totalTime("TimerTotal");
    totalTime.start();

    Graph graph, complete;

    std::cout << "Reading from file: " << inputFile << "\n";
    galois::graphs::readGraph(graph, inputFile);
    std::cout << "Read " << graph.size() << " nodes, " << graph.sizeEdges() << " edges\n";

    galois::graphs::readGraph(complete, completeFile);

    size_t approxNodeData = graph.size() * 64;
    galois::preAlloc(numThreads + approxNodeData / galois::runtime::pagePoolSize());
    galois::reportPageAlloc("MeminfoPre");

    std::cout << "Running " << ALGO_NAMES[algo] << " algorithm\n";

    galois::StatTimer execTime("Timer_0");
    execTime.start();

    Pattern p(size, &adj_mat[0]);

    switch (algo) {
    case Base:
        baseMatching(graph, p, complete, tri_cnt);
        break;
    default:
        std::abort();
    }

    execTime.stop();

    galois::reportPageAlloc("MeminfoPost");

    totalTime.stop();

  return 0;
}
