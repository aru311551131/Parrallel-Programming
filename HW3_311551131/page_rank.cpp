#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double *solution, double damping, double convergence)
{

  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
#pragma omp parallel for
  for (int i = 0; i < numNodes; ++i)
  {
    solution[i] = equal_prob;
  }

  // write here !
  double *score_old;
  score_old = (double*)malloc(sizeof(double) * g->num_nodes);

  while (true) {
    // copy new to old
#pragma omp parallel for
    for (int i = 0; i < numNodes; i++) {
      score_old[i] = solution[i];
    }
    // for task3
    double no_outgoing_num = 0.0;
#pragma omp parallel for reduction(+:no_outgoing_num)
    for (int i = 0; i < numNodes; i++) 
      if (outgoing_size(g, i) == 0)   no_outgoing_num += damping * score_old[i] / numNodes;
    // for all task
#pragma omp parallel for
    for (int i = 0; i < numNodes; i++) {
      const Vertex* start = incoming_begin(g, i);
      const Vertex* end = incoming_end(g, i);
      double sum = 0.0;
      for (const Vertex* v = start; v != end; v++)
        sum += score_old[*v] / (double)outgoing_size(g, *v);
      solution[i] = (damping * sum) + (1.0-damping) / numNodes + no_outgoing_num;
    }
    // for test
    double global_diff = 0.0;
#pragma omp parallel for reduction(+:global_diff)
    for (int i = 0; i < numNodes; i++) {
      global_diff += fabs(solution[i] - score_old[i]);
    }
    if(global_diff < convergence)
      break;
  }
  delete score_old;
  /*
     For PP students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes v in graph with no outgoing edges
                          { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */
}
