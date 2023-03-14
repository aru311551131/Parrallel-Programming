#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set *list)
{
    list->count = 0;
}

void vertex_set_init(vertex_set *list, int count)
{
    list->max_vertices = count;
    list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    int *distances)
{
    int current;
#pragma omp parallel for private(current) schedule(dynamic, 1024)
    for (int i = 0; i < frontier->count; i++)
    {

        int node = frontier->vertices[i];
        int dis = distances[node];
        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
        {
            int outgoing = g->outgoing_edges[neighbor];
            
            if (distances[outgoing] == NOT_VISITED_MARKER)
            {
                distances[outgoing] = dis + 1;
                do
                {
                    current = new_frontier->count;
                } while (!__sync_bool_compare_and_swap(&new_frontier->count, current, current + 1));
                new_frontier->vertices[current] = outgoing;
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);
    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}


void bottom_up_step(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    int *distances,
    int last_dis
    )
{
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < g->num_nodes; i++)
    {
        if (distances[i] == NOT_VISITED_MARKER) 
        {
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1) ? g->num_edges : g->incoming_starts[i + 1];
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int incoming = g->incoming_edges[neighbor];
                if (frontier->vertices[incoming])
                {
                    new_frontier->vertices[i] = 1;
                    distances[0] = 1;
                    distances[i] = last_dis + 1;
                    break;
                }
            }
        }
    }
}
void bfs_bottom_up(Graph graph, solution *sol)
{
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);
    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes; i++) 
        sol->distances[i] = -1;
    frontier->vertices[0] = 1;
    sol->distances[0] = 1;
    int last_dis = 0;
    
    while (sol->distances[ROOT_NODE_ID])
    {
        sol->distances[ROOT_NODE_ID] = 0;
        bottom_up_step(graph, frontier, new_frontier, sol->distances, last_dis);
        last_dis++;
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}
void bottom_up_step_h(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    int *distances,
    int last_dis
    )
{
    int current; //hybrid
#pragma omp parallel for schedule(dynamic, 1024) private(current)
    for (int i = 0; i < g->num_nodes; i++)
    {
        if (distances[i] == NOT_VISITED_MARKER) 
        {
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1) ? g->num_edges : g->incoming_starts[i + 1];
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int incoming = g->incoming_edges[neighbor];
                if (distances[i] == last_dis) //hybrid
                {
                    do
                    {
                        current = new_frontier->count;
                    } while (!__sync_bool_compare_and_swap(&new_frontier->count, current, current + 1));
                    new_frontier->vertices[current] = incoming;
                    new_frontier->vertices[i] = 1;
                    distances[i] = last_dis + 1;
                    break;
                }
            }
        }
    }
}
void bfs_hybrid(Graph graph, solution *sol)
{
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);
    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes; i++) 
        sol->distances[i] = -1;
    frontier->vertices[frontier->count++] = 0;
    sol->distances[0] = 0;

    int last_dis = 0;
    int some_node = graph->num_nodes / 3;
    int cnt_sum = 1;
    
    while (some_node > cnt_sum)
    {
        vertex_set_clear(new_frontier);
        top_down_step(graph, frontier, new_frontier, sol->distances);
        last_dis++;
        cnt_sum += new_frontier->count;
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
    vertex_set list3;
    vertex_set_init(&list3, graph->num_nodes);
    vertex_set *n_frontier = &list3;
    sol->distances[0] = 1;

#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes; i++) {
        n_frontier->vertices[i] = 0;
        new_frontier->vertices[i] = 0;
    }
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < frontier->count; i++) {
        n_frontier->vertices[frontier->vertices[i]] = 1;
        new_frontier->vertices[frontier->vertices[i]] = 1;
    }
    
    while (sol->distances[ROOT_NODE_ID])
    {
        sol->distances[ROOT_NODE_ID] = 0;
        bottom_up_step(graph, n_frontier, new_frontier, sol->distances, last_dis);
        last_dis++;
        vertex_set *tmp = n_frontier;
        n_frontier = new_frontier;
        new_frontier = tmp;
    }
}


/* bad version
void bottom_up_step(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    int *distances,
    int last_dis)
{
    int current;
#pragma omp parallel for private(current) schedule(dynamic, 1024)
    for (int i = 0; i < frontier->count; i++)
    {
        int node = frontier->vertices[i];
        int start_edge = g->incoming_starts[node];
        int end_edge = (node == g->num_nodes - 1) ? g->num_edges : g->incoming_starts[node + 1];
        bool found = false;
        for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
        {
            int incoming = g->incoming_edges[neighbor];
            if (distances[incoming] == last_dis)
            {
                distances[node] = last_dis + 1;
                found = true;
                break;
            }
        }
        if (!found) {
            do
            {
                current = new_frontier->count;
            } while (!__sync_bool_compare_and_swap(&new_frontier->count, current, current + 1));
            new_frontier->vertices[current] = node;
        }
    }


}
void bfs_bottom_up(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
//#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier without the root node
//#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes - 1; i++) 
        frontier->vertices[i] = i + 1;
    frontier->count = graph->num_nodes - 1;
    sol->distances[ROOT_NODE_ID] = 0;
    int last_dis = 0;
    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step(graph, frontier, new_frontier, sol->distances, last_dis);
        
#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        if (frontier->count == new_frontier->count) break;
        last_dis++;
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}
*/


void bottom_up_step_v2_slowQQ(
    Graph g,
    int *distances,
    int last_dis,
    int *flag    
    )
{
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < g->num_nodes; i++)
    {
        if (distances[i] == NOT_VISITED_MARKER) 
        {
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1) ? g->num_edges : g->incoming_starts[i + 1];
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int incoming = g->incoming_edges[neighbor];
                if (distances[incoming] == last_dis)
                {
                    *flag = 1;
                    distances[i] = last_dis + 1;
                    break;
                }
            }
        }
    }


}
void bfs_bottom_up_v2_slowQQ(Graph graph, solution *sol)
{
#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < graph->num_nodes; i++) 
        sol->distances[i] = NOT_VISITED_MARKER;
    sol->distances[ROOT_NODE_ID] = 0;
    int last_dis = 0;
    int flag = 1;
    
    while (flag)
    {
        flag = 0;
        bottom_up_step_v2_slowQQ(graph, sol->distances, last_dis, &flag);
        last_dis++;
    }
}