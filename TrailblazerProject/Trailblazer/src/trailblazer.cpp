/** Author: David Flatow
 *  Date: 11/21/2014
 *  File: trailblazer.cpp
 * -----------------------------------
 * This file implements various graph search algorithims and
 * Kruskal's algorithm for finding a minimum spanning tree.
 */

#include "costs.h"
#include "trailblazer.h"
#include "queue.h"
#include "pqueue.h"
#include "stack.h"

Vector<Vertex*> retraceSteps(Vertex* vert);
bool visitUsingDFS(BasicGraph& graph, Vertex* start, Vertex* end, Vertex* prev);
Vector<Vertex*> aStarWrapper(BasicGraph& graph, Vertex* start, Vertex* end, bool useHeuristic);
void mergeClusters(Vector<Set<Vertex*>>& clusters, Edge* e, Set<Edge*>& mst);

using namespace std;

/*
 * breadthFirstSearch implements a breadth first search on the graph from start to
 * end. If no path exisists an empty vector is returned. If start = end then a
 * vector containing only start is returned. Otherwise, a vector with the path from
 * start to end is returned.
 */
Vector<Vertex*> breadthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end) {

    graph.resetData();

    // initialize queue that will hold vertexes
    Queue<Vertex*> verts;

    start->setColor(YELLOW);
    start->visited = true;
    verts.enqueue(start);

    Vertex* currVertex;
    while (!verts.isEmpty()){

        currVertex = verts.dequeue();
        currVertex->setColor(GREEN);

        // we've found the end, search done
        if (currVertex == end) break;

        for (Edge* edge : currVertex->edges){

            Vertex* neighbor = edge->finish;

            // want only outgoing edges
            if (neighbor == currVertex) continue;

            // vertex already visited
            if (neighbor->visited) continue;

            // set neighbor's parent to be curr_vertex
            neighbor->previous = currVertex;

            neighbor->setColor(YELLOW);
            neighbor->visited = true;
            verts.enqueue(neighbor);

        }
    }

    // retrace steps and return path
    return retraceSteps(end);
}

/*
 * retraceSteps returns a vector representing the path from some starting vertex to
 * the vertex vert. The path is determined, via the Hansel and Gretel method, by retracing
 * over the 'previous' pointers at each vertex until such pointer is NULL. Since the path
 * is discovered in reverse order, we store an intermediary copy in stack, then pop all
 * vertexes into a vector and return it.
 */
Vector<Vertex*> retraceSteps(Vertex* vert){

    Vector<Vertex*> path;

    // no path leading to vert
    if (vert->previous == NULL) return path;

    // retrace steps (collect breadcrumbs)
    Stack<Vertex*> pathStack; // intermediary stack, used to reverse order of vertexes
    while (vert != NULL){
        pathStack.push(vert);
        vert = vert->previous;
    }

    // pop into Vector, populate Vector path in correct forwards order
    while (!pathStack.isEmpty()){
        path.add(pathStack.pop());
    }

    return path;
}

/*
 * depthFirstSearch carries out a depth first search on the graph to find a
 * path from start to end. this implementation is recursive.
 */
Vector<Vertex*> depthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end) {
    graph.resetData();
    visitUsingDFS(graph, start, end, NULL);
    return retraceSteps(end);
}

bool visitUsingDFS(BasicGraph& graph, Vertex* start, Vertex* end, Vertex* prev) {

    // base case: already visisted, search fail
    if (start->visited) return false;

    // TRY
    start->previous = prev;
    start->visited = true;
    start->setColor(GREEN);

    // base case: end found, search success
    if (start->name == end->name) return true;

    // recursive step
    for (Edge* edge : start->edges){
        if (visitUsingDFS(graph, edge->finish, end, start)) return true;
    }

    /* UNDO
     * note, we don't set visited to false since, if any path was visited
     * at this point was involved in a search that failed to find end, and
     * so any search that gets to this vertex will also fail.
     */
    start->previous = NULL; // probably don't need this, just feel a little safer with it though
    start->setColor(GRAY);
    return false;
}

/*
 * dijkstrasAlgorithm search algorithm is a special case of aStar using a constant
 * heuristicFunction (or none at all). see aStar for more details
 */
Vector<Vertex*> dijkstrasAlgorithm(BasicGraph& graph, Vertex* start, Vertex* end) {
    return aStarWrapper(graph, start, end, false);
}

/*
 * aStar calls aStarWrapper specifying whether or not to use heuristic function. we
 * do this so we can implement dijkstrasAlgorithm without copying code.
 */
Vector<Vertex*> aStar(BasicGraph& graph, Vertex* start, Vertex* end) {
    return aStarWrapper(graph, start, end, true);
}

/*
 * aStar search algorith is a sort of breadth first search; however, paths with lower cost and
 * lower heuristic cost are explored first.
 */
Vector<Vertex*> aStarWrapper(BasicGraph& graph, Vertex* start, Vertex* end, bool useHeuristic) {
    graph.resetData();

    start->setColor(YELLOW);
    start->visited = true;
    start->cost = 0.0;

    // initialize queue that will hold vertexes
    PriorityQueue<Vertex*> verts;
    double priority = (useHeuristic)? heuristicFunction(start, end) : 0.0;
    verts.enqueue(start, priority);

    Vertex* currVertex;
    while (!verts.isEmpty()){

        currVertex = verts.dequeue();
        currVertex->setColor(GREEN);
        currVertex->visited = true;

        // we've found the end, search done
        if (currVertex == end) break;

        for (Edge* edge : currVertex->edges){

            Vertex* neighbor = edge->finish;

            // vertex already visited
            if (neighbor->visited) continue;

            // new proposed total cost to get to neighbor
            double newCost = currVertex->cost + edge->cost;

            // new cost is higher than existing, ignore
            if ((neighbor->previous != NULL) && (newCost >= neighbor->cost)) continue;

            // queue or update in queue
            priority = (useHeuristic)? newCost + heuristicFunction(neighbor, end) : newCost;
            if (neighbor->previous == NULL){
                verts.enqueue(neighbor, priority);
                neighbor->setColor(YELLOW);
            } else {
                verts.changePriority(neighbor, priority);
            }

            // update values and links
            neighbor->cost = newCost;
            neighbor->previous = currVertex;
        }
    }

    // retrace steps and return path
    return retraceSteps(end);
}

/*
 * kruskal returns a minimum spanning tree for the graph by implementing
 * kruskal's algorithm. the algorithm is achieved by itteratively merging
 * clusters of nodes so long as there's an edge connecting any member in
 * one cluster to any member in the other cluster.
 */
Set<Edge*> kruskal(BasicGraph& graph) {

    Set<Edge*> mst;

    // put all edges into pqueue using weights as priorities
    PriorityQueue<Edge*> edges;
    for (Edge* e : graph.getEdgeSet()){
        edges.enqueue(e, e->cost);
    }

    // place each vertex of graph as single element sets into a vector
    Vector<Set<Vertex*>> clusters;
    for (Vertex* v: graph.getVertexSet()){
        Set<Vertex*> s;
        s.add(v);
        clusters.add(s);
    }

    // while there are two or more separate clusters remaining
    while (!edges.isEmpty() && clusters.size() >= 2){
        Edge* e = edges.dequeue();
        mergeClusters(clusters, e, mst);
    }

    return mst;
}

void mergeClusters(Vector<Set<Vertex*>>& clusters, Edge* e, Set<Edge*>& mst){
    // vector to hold indexes of clusters to merge
    Vector<int> indexesToMerge;

    // find indexes of clusters containing start and finish of edge
    for (int i = 0; i < clusters.size(); i++){
        if (clusters[i].contains(e->start)) indexesToMerge.add(i);
        if (clusters[i].contains(e->finish)) indexesToMerge.add(i);
        if (indexesToMerge.size() == 2) break;
    }

    // found fewer than two clusters containing start and finish. no merging.
    if (indexesToMerge.size() < 2) return;

    // start and finish in same cluster. no merging needed.
    if (indexesToMerge[0] == indexesToMerge[1]) return;

    // merge clusters containing start and finish
    clusters[indexesToMerge[0]].addAll(clusters[indexesToMerge[1]]);
    clusters.remove(indexesToMerge[1]);
    mst.add(e);
}

/*
 * prim returns a minimum spanning tree for the graph by implementing
 * prim's algorithm. the algorithm is achieved by iteratively adding the
 * minimum weight edge connecting any node in the tree of edges to a node
 * not in the tree. this is repeated until all nodes are in the tree.
 */
Set<Edge*> prim(BasicGraph& graph) {

    Set<Edge*> mst;

    Set<Vertex*> vertexes;
    Vertex* initialVertex = graph.getVertexSet().first();
    vertexes.add(initialVertex);

    while (true){
        for (Edge* e: graph.getEdgeSet()){

            if ((vertexes.contains(e->start) && !vertexes.contains(e->finish)) ||
                (!vertexes.contains(e->start) && vertexes.contains(e->finish))){
                vertexes.add(e->start);
                vertexes.add(e->finish);
                mst.add(e);
                break;
            }
        }
        // done
        if (vertexes == graph.getVertexSet()) break;
    }
    return mst;
}
