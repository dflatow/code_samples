/*
 * CS 106X Trailblazer
 * This file declares code to serve as an "adapter" to convert between various
 * data types used by legacy support code and new code for the current version
 * of this assignment.
 * Specifically, past quarters represented a graph as a Grid<double> and each
 * of its vertices/edges as a Loc (TBLoc) and Edge (TBEdge), while the
 * current version of the assignment represents the graph as a BasicGraph
 * (an extension of Stanford's Graph class) and each vertex/edge as a
 * Vertex / Edge object.
 * This code converts between the two formats to help localize changes made
 * to the assignment starting in 13au by Marty Stepp.
 *
 * Please do not modify this provided file. Your turned-in files should work
 * with an unmodified version of all provided code files.
 *
 * Author: Marty Stepp
 * Version: 2014/11/10 (initial version for 14aut)
 */

#ifndef _adapter_h
#define _adapter_h

#include "basicgraph.h"
#include "grid.h"
#include "set.h"
#include "types.h"

/* Type: AlgorithmType
 *
 * An enumerated type representing the various algorithms to choose from.
 */
enum AlgorithmType {
    AUTODETECT,
    BFS,
    DFS,
    DIJKSTRA,
    A_STAR,
    BIDIRECTIONAL
};

/*
 * Finds the shortest path between the locations given by start and end in the
 * specified world.  The cost of moving from one edge to the next is specified
 * by the given cost function.  The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.    If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this function prototype and the
 * implementation inside of Trailblazer.cpp.
 */
Vector<TBLoc>
shortestPath(TBLoc start,
             TBLoc end,
             const Grid<double>& world,
             double costFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world),
             double heuristicFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world),
             AlgorithmType algorithm = AUTODETECT);

// Support functions called by the GUI to improve loading times for large graphs.

/*
 * Converts the given graph into a 2D grid structure, which is returned;
 * essentially the opposite of gridToGraph.
 * Used to convert randomly generated mazes into grids for legacy code support.
 */
Grid<double>* graphToGrid(BasicGraph* graph);

/*
 * Converts the given 2D grid structure into a BasicGraph object,
 * using the Grid to grab the vertices and the given cost function to grab
 * the edge weights between neighbors.
 * Returns a pointer to the graph that was created.
 */
BasicGraph* gridToGraph(const Grid<double>& world,
                        double costFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world));

/*
 * Makes sure that the given world has already been converted into an equivalent
 * BasicGraph.  If it has not, does so (by calling gridToGraph) and caches the
 * result to be used on future calls.
 * This is done to improve runtime when very large/huge mazes and terrains are
 * loaded and then searched multiple times by the user.
 */
void ensureWorldCache(const Grid<double>& world,
                      double costFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world));

/*
 * Removes all entries from the internal cache of BasicGraphs and frees
 * any memory associated with them.
 */
void flushWorldCache();

/*
 * Returns a grid representing a maze of the given size with randomly
 * filled walls.  The walls are generated through Kruskal's algorithm.
 */
Grid<double> createRandomMaze(int size);


/*
 * Returns a name for the given row/column location in the given world,
 * such as "r08c17".
 * The reason you have to pass the world is so that I know how many rows/cols
 * it has so I can zero-pad the numbers in the string.
 */
string vertexName(int r, int c, const Grid<double>& world);
string vertexName(int r, int c, int digits);

#endif
