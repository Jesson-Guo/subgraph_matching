#include "subgraph_matching_base.h"
#include "llvm/Support/CommandLine.h"
namespace cll = llvm::cl;

static const char* name = "subgraph matching";
static const char* desc = "Computes the number of subgraphs matching the pattern";
static const char* url = "subgraph_matching";

enum Algo {
    Base = 0,
};

const char* const ALGO_NAMES[] = {"base"};

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
    const auto& allNodes_c = complete.allNodesRange();
    galois::do_all(
        galois::iterate(allNodes_c),
        [&](uint64_t i) {
            complete.getData(i).degree = complete.globalSize() - 1;
            for (int v=0; v<(int)complete.globalSize(); v++) {
                if (v != (int)complete.getGID(i)) {
                    complete.getData(i).adj_node.push_back(v);
                }
            }
        },
        galois::loopname("complete"),
        galois::chunk_size<CHUNK_SIZE>(),
        galois::steal(),
        galois::no_stats()
    );
    const auto& allNodes = graph.allNodesRange();
    galois::do_all(
        galois::iterate(allNodes),
        [&](uint64_t i) {
            int deg(0);
            for (auto e : graph.edges(i)) {
                GNode dst = graph.getEdgeDst(e);
                graph.getData(i).adj_node.push_back(graph.getGID(dst));
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
    syncSubstrate->sync<writeSource, readAny, Reduce_set_degree>("InitializeGraph");
    syncSubstrate->sync<writeSource, readAny, Reduce_set_adj_node>("InitializeGraph");
    galois::runtime::getHostBarrier().wait();

    max_degree = max_deg;
    VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, sec_max_deg);
    bool is_pattern_valid;
    int perf_type(1), res_type(1);
    bool opt(true);
    Schedule schedule(p, is_pattern_valid, perf_type, res_type, opt, complete, max_deg, graph.size(), graph.sizeEdges(), tri_cnt);
    assert(is_pattern_valid);

    galois::StatTimer matchingTime("Timer_matching");
    matchingTime.start();
    long long ans = pattern_matching(graph, schedule, max_deg);
    matchingTime.stop();

    std::cout << "ans:\t" << ans << std::endl;
    std::cout << "time:\t" << matchingTime.get_usec() << std::endl;
    std::cout << "restricts:\n";
    for (auto p : schedule.restrict_pair) {
        std::cout << "(" << p.first << ", " << p.second << ")  ";
    }
    std::cout << std::endl;
}

int main(int argc, char** argv) {
    // sleep(20);
    galois::DistMemSys G;
    DistBenchStart(argc, argv, name, desc, url);

    auto& net = galois::runtime::getSystemNetworkInterface();

    if (net.ID == 0) {
        galois::runtime::reportParam("subgraph_matching", "pattern size", size);
        galois::runtime::reportParam("subgraph_matching", "pattern matrix", adj_mat);
    }

    galois::StatTimer totalTime("TimerTotal");
    totalTime.start();

    std::unique_ptr<Graph> hg, hg_c;

    std::tie(hg, syncSubstrate) = distGraphInitialization<NodeData, void>();
    // OEC default symmetric graph
    hg_c = galois::cuspPartitionGraph<NoCommunication, NodeData, void>(
        completeFile, galois::CUSP_CSR, galois::CUSP_CSR, false, inputFileTranspose, mastersFile);
    // Graph complete;
    // galois::graphs::readGraph(complete, completeFile);

    Pattern p(size, &adj_mat[0]);

    galois::StatTimer execTime("Timer_0");
    execTime.start();
    baseMatching(*hg, p, *hg_c, tri_cnt);
    execTime.stop();

    totalTime.stop();

    return 0;
}
