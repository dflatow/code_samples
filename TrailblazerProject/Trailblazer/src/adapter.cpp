/*
 * CS 106X Trailblazer
 * This file implements code to serve as an "adapter" to convert between various
 * data types used by legacy support code and new code for the current version
 * of this assignment.
 * See adapter.h for more details.
 *
 * Please do not modify this provided file. Your turned-in files should work
 * with an unmodified version of all provided code files.
 *
 * Author: Marty Stepp
 * Version: 2014/11/10 (initial version for 14aut)
 */

#include "adapter.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include "basicgraph.h"
#include "gmath.h"
#include "graph.h"
#include "map.h"
#include "random.h"
#include "costs.h"
#include "trailblazer.h"
#include "trailblazergui.h"
#include "types.h"

// global variables
static Map<Grid<double>*, BasicGraph*> WORLD_CACHE;
static double (*s_heuristicFunction)(const TBLoc& from, const TBLoc& to, const Grid<double>& world) = NULL;
static Grid<double>* s_world = NULL;

// function prototype declarations
static int getRow(Vertex* v);
static int getCol(Vertex* v);
static int getGridValue(Vertex* v);

double heuristicFunction(Vertex* from, Vertex* to) {
    if (s_world == NULL || s_heuristicFunction == NULL) {
        return 0.0;
    } else {
        TBLoc* fromLoc = (TBLoc*) from->extraData;
        TBLoc* toLoc = (TBLoc*) to->extraData;
        return s_heuristicFunction(*fromLoc, *toLoc, *s_world);
    }
}


// function implementations

class VertexObserver : public Observer {
public:
    const Grid<double>* world;

    void update(Observable* obs, void* /*arg*/) {
        Vertex* v = (Vertex*) obs;
        TBLoc* loc = (TBLoc*) v->extraData;
        Color color = v->getColor();
        colorCell(*world, *loc, color);
    }
};

BasicGraph* gridToGraph(const Grid<double>& world,
                        double costFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world)) {
    BasicGraph* graph = new BasicGraph();
    VertexObserver* obs = new VertexObserver();
    obs->world = &world;

    // add vertices
    int rows = world.numRows();
    int cols = world.numCols();
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            string name = vertexName(r, c, world);
            Vertex* v = new Vertex(name);
            v->extraData = new TBLoc(r, c, world.get(r, c));
            v->addObserver(obs);
            graph->addVertex(v);
        }
    }

    // add edges
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            Vertex* v = graph->getVertex(vertexName(r, c, world));
            for (int dr = -1; dr <= 1; dr++) {
                for (int dc = -1; dc <= 1; dc++) {
                    int nr = r + dr;
                    int nc = c + dc;
                    if ((dr == 0 && dc == 0) || !world.inBounds(nr, nc)) {
                        continue;
                    }
                    Vertex* neighbor = graph->getVertex(vertexName(nr, nc, world));
                    double cost = costFn(TBLoc(r, c), TBLoc(nr, nc), world);
                    if (cost != POSITIVE_INFINITY) {
                        Edge* e = new Edge(v, neighbor, cost);
                        e->extraData = new TBEdge(TBLoc(r, c), TBLoc(nr, nc));
                        graph->addEdge(e);
                    }
                }
            }
        }
    }

    return graph;
}

Grid<double> createRandomMaze(int size) {
    Grid<double> maze(size, size, kMazeFloor);
    BasicGraph* graph = gridToGraph(maze, mazeCost);
    
    // assign random weights to the edges
    // give each edge a 'random' weight;
    // put all edges into a priority queue, sorted by weight
    Set<Edge*> edgeSet = graph->getEdgeSet();
    int edgeCount = edgeSet.size();
    for (Edge* edge : edgeSet) {
        int weight = randomInteger(1, edgeCount * 1000);
        edge->cost = weight;
    }
    
    // run the student's Kruskal algorithm to get a minimum spanning tree (MST)
    //Set<Edge*> mst = kruskal(*graph);
    Set<Edge*> mst = prim(*graph);
    
    // convert the MST/graph back into a maze grid
    // (insert a 'wall' between any neighbors that do not have a connecting edge)
    graph->clearEdges();
    for (Edge* edge : mst) {
        graph->addEdge(edge->start, edge->finish);
        graph->addEdge(edge->finish, edge->start);
    }

    // physical <-> logical size; a maze of size MxN has 2M-1 x 2N-1 grid cells.
    // cells in row/col 0, 2, 4, ... are open squares (floors), and cells in 
    // row/col 1, 3, 5, ... are blocked (walls).
    int digits = countDigits(size);
    int worldSize = size * 2 - 1;
    Grid<double> world(worldSize, worldSize, kMazeWall);
    for (int row = 0; row < worldSize; row++) {
        for (int col = 0; col < worldSize; col++) {
            if (row % 2 == 0 && col % 2 == 0) {
                world[row][col] = kMazeFloor;
            }
        }
    }
    
    for (int row = 0; row < size; row++) {
        int gridRow = row * 2;
        for (int col = 0; col < size; col++) {
            int gridCol = col * 2;
            string name = vertexName(row, col, digits);
            
            // decide whether to put open floor between neighbors
            // (if there is an edge between them)
            for (int dr = -1; dr <= 1; dr++) {
                int nr = row + dr;
                int gridNr = gridRow + dr;
                for (int dc = -1; dc <= 1; dc++) {
                    int nc = col + dc;
                    int gridNc = gridCol + dc;
                    if ((nr != row && nc != col)
                            || (nr == row && nc == col)
                            || !world.inBounds(gridNr, gridNc)) {
                        continue;
                    }
                    string neighborName = vertexName(nr, nc, digits);
                    if (graph->containsEdge(name, neighborName)) {
                        world[gridNr][gridNc] = kMazeFloor;
                    }
                }
            }
        }
    }
    
    delete graph;
    return world;
}

Grid<double>* graphToGrid(BasicGraph* graph) {
    int maxRow = -1;
    int maxCol = -1;
    for (Vertex* v : graph->getVertexSet()) {
        maxRow = max(maxRow, getRow(v));
        maxCol = max(maxCol, getCol(v));
    }
    if (maxRow < 0 || maxCol < 0) {
        return NULL;
    }

    // initially set all to be walls
    Grid<double>* grid = new Grid<double>(maxRow + 1, maxCol + 1);
    for (int r = 0; r <= maxRow; r++) {
        for (int c = 0; c <= maxCol; c++) {
            grid->set(r, c, kMazeWall);
        }
    }
    for (Vertex* v : graph->getVertexSet()) {
        grid->set(getRow(v), getCol(v), getGridValue(v));
    }

    return grid;
}

void ensureWorldCache(const Grid<double>& world,
                      double costFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world)) {
    Grid<double>* const pWorld = const_cast<Grid<double>*>(&world);
    if (!WORLD_CACHE.containsKey(pWorld)) {
        cout << "Preparing world model ..." << endl;
        BasicGraph* graph = gridToGraph(world, costFn);
        cout << "World model completed." << endl;
        WORLD_CACHE[pWorld] = graph;
    }
}

void flushWorldCache() {
    for (Grid<double>* grid : WORLD_CACHE) {
        BasicGraph* graph = WORLD_CACHE[grid];
        delete graph;
    }
    WORLD_CACHE.clear();
}

Vector<TBLoc>
shortestPath(TBLoc start,
             TBLoc end,
             const Grid<double>& world,
             double costFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world),
             double heuristicFn(const TBLoc& from, const TBLoc& to, const Grid<double>& world),
             AlgorithmType algorithm) {
    // modified by Marty to use an actual Graph object
    ensureWorldCache(world, costFn);
    cout << endl;

    Grid<double>* const pWorld = const_cast<Grid<double>*>(&world);
    BasicGraph* graph = WORLD_CACHE[pWorld];
    // graph->resetData();   // make the student worry about this

    s_world = pWorld;
    s_heuristicFunction = heuristicFn;

    // convert start/end from Loc to Vertex
    Vertex* startVertex = graph->getVertex(vertexName(start.row, start.col, world));
    Vertex* endVertex   = graph->getVertex(vertexName(end.row, end.col, world));
    if (startVertex == NULL) {
        error(string("Graph can not find start vertex with name \"") + vertexName(start.row, start.col, world) + "\"");
        for (Vertex* v : graph->getVertexSet()) {
            cout << v->name << " ";
        }
        cout << endl;
    }
    if (endVertex == NULL) {
        error(string("Graph can not find end vertex with name \"") + vertexName(start.row, start.col, world) + "\"");
        for (Vertex* v : graph->getVertexSet()) {
            cout << v->name << " ";
        }
        cout << endl;
    }

    cout << "Looking for a path from " << startVertex->name
         << " to " << endVertex->name << "." << endl;

    Vector<Vertex*> result;
    switch (algorithm) {
    case BFS:
        cout << "Executing breadth-first search algorithm ..." << endl;
        result = breadthFirstSearch(*graph, startVertex, endVertex);
        break;
    case DIJKSTRA:
        cout << "Executing Dijkstra's algorithm ..." << endl;
        result = dijkstrasAlgorithm(*graph, startVertex, endVertex);
        break;
    case A_STAR:
        cout << "Executing A* algorithm ..." << endl;
        result = aStar(*graph, startVertex, endVertex);
        break;
#ifdef BIDIRECTIONAL_SEARCH_ALGORITHM_ENABLED
    case BIDIRECTIONAL:
        cout << "Executing Bidirectional Search algorithm ..." << endl;
        extern Vector<Vertex*> bidirectionalSearch(BasicGraph& graph, Vertex* start, Vertex* end);
        result = bidirectionalSearch(*graph, startVertex, endVertex);
        break;
#endif // BIDIRECTIONAL_SEARCH_ALGORITHM_ENABLED
    case DFS:
    default:
        cout << "Executing depth-first search algorithm ..." << endl;
        result = depthFirstSearch(*graph, startVertex, endVertex);
        break;
    }

    cout << "Algorithm complete." << endl;

    // convert Vector<Vertex*> to Vector<Loc>
    Vector<TBLoc> locResult;
    for (Vertex* v : result) {
        locResult.add(TBLoc(getRow(v), getCol(v)));
    }
    return locResult;
}

Set<TBEdge> createMaze(int /* numRows */, int /* numCols */) {
    Set<TBEdge> set;
    return set;
}

static int getRow(Vertex* v) {
    TBLoc* loc = (TBLoc*) v->extraData;
    return loc->row;
}

static int getCol(Vertex* v) {
    TBLoc* loc = (TBLoc*) v->extraData;
    return loc->col;
}

static int getGridValue(Vertex* v) {
    TBLoc* loc = (TBLoc*) v->extraData;
    return loc->gridValue;
}

string vertexName(int r, int c, const Grid<double>& world) {
    // zero-pad the number of rows/cols for better alphabetic sorting
    int rowCols = max(world.numRows(), world.numCols());
    return vertexName(r, c, countDigits(rowCols));
}

string vertexName(int r, int c, int digits) {
    ostringstream out;
    out << 'r' << setw(digits) << setfill('0') << r
        << 'c' << setw(digits) << setfill('0') << c;
    return out.str();
}
