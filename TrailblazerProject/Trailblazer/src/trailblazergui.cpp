/*
 * CS 106X Trailblazer
 * This file implements functions to perform drawing in the graphical user
 * interface (GUI).
 * This file also contains the main function to launch the program.
 * See trailblazergui.h for documentation of each public function.
 *
 * Please do not modify this provided file. Your turned-in files should work
 * with an unmodified version of all provided code files.
 *
 * Author: Marty Stepp, Keith Schwarz, et al
 * Version: 2014/11/10 (initial version for 14aut)
 */

#include "trailblazergui.h"
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include "error.h"
#include "filelib.h"
#include "gevents.h"
#include "gfilechooser.h"
#include "ginteractors.h"
#include "gobjects.h"
#include "gwindow.h"
#include "set.h"
#include "simpio.h"
#include "strlib.h"
#include "vector.h"
#include "adapter.h"
#include "costs.h"
#include "trailblazer.h"
using namespace std;

// internal constants
const int MIN_ANIMATION_DELAY = 0;
const int MAX_ANIMATION_DELAY = 2000;
const int DEFAULT_ANIMATION_DELAY = 1000;
static int animationDelay = 0;
static bool pathSearchInProgress = false;

// various UI strings and settings
const string kSelectedLocationColor("RED");
const string kPathColor("RED");
const string kBackgroundColor("Black");
const string GUI_STATE_FILE("trailblazer-gui-state.sav");
const string OTHER_FILE_LABEL("Other file ...");
const bool SHOULD_SAVE_GUI_STATE = true;
const bool RANDOM_MAZE_OPTION_ENABLED = true;
const bool RANDOM_USE_SEED = false;   // true to get predictable random mazes

/*
 * The size, in pixels, of the displayed world.
 * Number of padding pixels between the border of the window and the
 * start of the content.
 * Constants controlling the amount we should adjust the size of the window
 * to reflect the extra space used up by the border and controls.
 */
const int WINDOW_WIDTH = 500;
const int WINDOW_HEIGHT = 500;
const int kMargin = 5;
int kWidthAdjustment  = 0;
int kHeightAdjustment = 0;

/*
 * Maximum number of rows or columns we allow in a world.    This is mostly a
 * safety feature to prevent an OOM on a malformed input file.
 */
const int kMaxRows = 400;
const int kMaxCols = 400;

/*
 * Constants controlling the number of rows/cols in the different sizes of
 * worlds.  The indices into these arrays are controlled by the enumerated
 * type WorldSize.
 */
const int kTerrainNumRows[] = { 10, 33, 65, 129, 257 };
const int kTerrainNumCols[] = { 10, 33, 65, 129, 257 };
const int kMazeNumRows[]    = {  6, 11, 31,  81, 161 };
const int kMazeNumCols[]    = {  6, 11, 31,  81, 161 };

/*
 * Constants reflecting the colors that are used to highlight cells different
 * colors.  This is an array of arrays, where each array contains RGB triplets
 * for the collors.  The ordering here should be consistent with the ordering
 * of the Color type.
 */
const int kColorMultipliers[6][3] = {
    {   0,   0,   0 },   // UNCOLORED
    { 255, 255, 255 },   // WHITE
    { 192, 192, 192 },   // GRAY
    { 255, 255,   0 },   // YELLOW
    {   0, 255,   0 },   // GREEN
    { 255,   0,   0 }    // RED
};


// Internal global variables

static GWindow* gWindow = NULL;
static GChooser* gWorldChooser = NULL;
static GChooser* gAlgorithmList = NULL;
static GSlider* gDelaySlider = NULL;
static GRect* gFirstSelected = NULL;      // which cells were seletected so far
static GRect* gSecondSelected = NULL;
static GTextField* gPositionField = NULL;
static string gPositionFieldText = " r???c???";
static Grid<Color> gMarked;               // which cells have been colored by the user
static Grid<double> gMarkedValues;        // the values we've marked them with
static Vector<GLine*> gHighlightedPath;   // all highlighting lines
static TBLoc gStartLocation;              // their corresponding locations
static TBLoc gEndLocation;
static double gPixelsPerWidth;            // width of each cell
static double gPixelsPerHeight;           // height of each cell
static State state;


// internal function prototypes (see implementations for documentation)

static void colorLocation(TBLoc loc, double value, Color locColor);
static TBLoc coordToLoc(double x, double y);
static void drawWorld(Grid<double>& world);
static void fillRect(int x, int y, int width, int height, string color);
static void findMidpoint(TBLoc loc, double& xc, double& yc);
static Set<string> getFiles(string substr);
static WorldSize getWorldSize(string worldText);
static void intro();
static bool loadGUIState();
static void locToCoord(TBLoc loc, double& x, double& y);
static bool readWorldFile(ifstream& input, Grid<double>& world, WorldType& worldType);
static bool regenerateWorld(Grid<double>& world, WorldType& worldType);
static bool registerClick(Grid<double>& world, double x, double y, WorldType worldType);
static void removeAndDelete(GObject* object);
static void removeMarkedSquares();
static void removeOverlays();
static void restoreWorldDisplay();
static void runSearch(State& state);
static void runShortestPath(const Grid<double>& world, WorldType worldType,
                            const TBLoc& start, const TBLoc& end, Vector<TBLoc>& path, double& pathCost);
static bool saveGUIState();
static void setAnimationDelay(int delay);
static void shutdownProgram();
static void uncolorSquares();
static void updateAnimationDelayFromSlider(bool forbidZero = false);
static string valueToColor(double value, Color locColor);
static void worldUpdated(Grid<double>& world, WorldType worldType);


// function implementations

/*
 * Main program.
 */
int main() {
    if (RANDOM_USE_SEED) {
        setRandomSeed(106);
    }
    
    setConsoleLocation(WINDOW_WIDTH + 6 * kMargin, 30);   // temporary
    setConsoleSize(450, 400);
    setConsoleEventOnClose(true);
    initializeGUI();
    Point guiLoc = gWindow->getLocation();
    GDimension guiSize = gWindow->getSize();
    setConsoleLocation(
            guiLoc.getX() + guiSize.getWidth() + kMargin,
            guiLoc.getY());

    // Main event loop to process events as they happen.
    while (true) {
        GEvent e = waitForEvent(ACTION_EVENT | MOUSE_EVENT | WINDOW_EVENT);
        if (e.getEventType() == MOUSE_CLICKED || e.getEventType() == MOUSE_MOVED) {
            processMouseEvent(GMouseEvent(e));
        } else if (e.getEventClass() == ACTION_EVENT) {
            processActionEvent(GActionEvent(e));
        } else if (e.getEventClass() == WINDOW_EVENT) {
            processWindowEvent(GWindowEvent(e));
        }
    }
    return 0;
}

/*
 * (public function)
 * Colors the specified cell in the world the indicated color.
 */
void colorCell(const Grid<double>& world, TBLoc loc, Color locColor) {
    colorLocation(loc, world[loc.row][loc.col], locColor);
    gMarked[loc.row][loc.col] = locColor;

    // possibly delay and then repaint the square, for animation
    if (animationDelay > 0) {
        gWindow->repaint();

        // while I'm here, check if delay updated
        updateAnimationDelayFromSlider(/* forbidZero */ false);

        pause(animationDelay);
    }
}

/*
 * (public function)
 * Initializes state of the GUI subsystem.
 */
void initializeGUI() {
    // Calculate the intended width and height of the window based on the content
    // area size, the margin size, and the adjustment amount.
    int windowWidth = WINDOW_WIDTH + 2 * kMargin + kWidthAdjustment;
    int windowHeight = WINDOW_HEIGHT + 2 * kMargin + kHeightAdjustment;

    gWindow = new GWindow(windowWidth, windowHeight);
    gWindow->setWindowTitle("CS 106X Trailblazer");
    // gWindow->setExitOnClose(true);
    // gWindow->setLocation(1, gWindow->getLocation().getY());

    // Add the algorithms list.
    gAlgorithmList = new GChooser();
    gAlgorithmList->addItem("DFS");
    gAlgorithmList->addItem("BFS");
    gAlgorithmList->addItem("Dijkstra");
    gAlgorithmList->addItem("A*");
#ifdef BIDIRECTIONAL_SEARCH_ALGORITHM_ENABLED
    gAlgorithmList->addItem("Bidirectional");
#endif // BIDIRECTIONAL_SEARCH_ALGORITHM_ENABLED

    gDelaySlider = new GSlider(MIN_ANIMATION_DELAY, MAX_ANIMATION_DELAY, MIN_ANIMATION_DELAY);

    gPositionField = new GTextField(7);
    gPositionField->setText(gPositionFieldText);
    gPositionField->setEditable(false);

    // Add in the list of existing world files.
    gWorldChooser = new GChooser();
    Set<string> worldFiles = getFiles("maze") + getFiles("terrain");
    for (string worldFile : worldFiles) {
        gWorldChooser->addItem(worldFile);
    }
    gWorldChooser->addItem(OTHER_FILE_LABEL);
    if (RANDOM_MAZE_OPTION_ENABLED) {
        // gWorldChooser->addItem("Random Maze (tiny)");
        gWorldChooser->addItem("Random Maze (small)");
        gWorldChooser->addItem("Random Maze (medium)");
        gWorldChooser->addItem("Random Maze (large)");
        gWorldChooser->addItem("Random Maze (huge)");
        // gWorldChooser->addItem("Random Terrain (tiny)");
        // gWorldChooser->addItem("Random Terrain (small)");
        // gWorldChooser->addItem("Random Terrain (medium)");
        // gWorldChooser->addItem("Random Terrain (large)");
        // gWorldChooser->addItem("Random Terrain (huge)");
    }

    // north layout
    gWindow->addToRegion(gAlgorithmList, "NORTH");
    gWindow->addToRegion(new GLabel("Delay:"), "NORTH");
    gWindow->addToRegion(gDelaySlider, "NORTH");
    gWindow->addToRegion(new GButton("Run"), "NORTH");
    gWindow->addToRegion(new GButton("Clear"), "NORTH");
    gWindow->addToRegion(gPositionField, "NORTH");

    // south layout
    gWindow->addToRegion(new GLabel("World:"), "SOUTH");
    gWindow->addToRegion(gWorldChooser, "SOUTH");
    gWindow->addToRegion(new GButton("Load"), "SOUTH");
    // gWindow->addToRegion(new GButton("Exit"), "SOUTH");

    intro();
    setAnimationDelay(DEFAULT_ANIMATION_DELAY);
    gDelaySlider->setValue(animationDelay);
    if (SHOULD_SAVE_GUI_STATE) {
        loadGUIState();
    }
    
    gWindow->pack();   // correct the window size
    
    if (!regenerateWorld(state.world, state.worldType)) {
        error("Cannot set up initial world properly!");
    }
    state.uiState = FRESH;
    drawWorld(state.world);
}

/*
 * (public function)
 * Reacts to an action event in the window.
 */
void processActionEvent(GActionEvent e) {
    string cmd = e.getActionCommand();
    if (cmd == "Load") {
        // Want to load a new world?  Try to do so and update the UI accordingly.
        if (regenerateWorld(state.world, state.worldType)) {
            drawWorld(state.world);
            state.uiState = FRESH;
        }
    } else if (cmd == "Run") {
        // Rerunning the search is only possible if we already did a search.
        if (state.uiState == DRAWN) {
            uncolorSquares();
            removeOverlays();
            runSearch(state);
        } else {
            cout << "Cannot rerun a search; no search has been done." << endl;
        }
    } else if (cmd == "Clear") {
        // Clearing the display just sets us back to the fresh state.
        restoreWorldDisplay();
        state.uiState = FRESH;
    } else if (cmd == "Exit") {
    }
}

/*
 * (public function)
 * Reacts to a mouse event in the window.
 */
void processMouseEvent(GMouseEvent e) {
    if (e.getEventType() == MOUSE_CLICKED) {
        switch (state.uiState) {
        case DRAWN:
            // Already have something drawn?  Clear it and pretend we're fresh.
            restoreWorldDisplay();
            state.uiState = FRESH;
            // deliberate: don't break.

        case FRESH:
            // In a fresh state?  Try to select what the user clicked on.
            if (registerClick(state.world, e.getX(), e.getY(), state.worldType)) {
                state.uiState = MARKED;
            }
            break;

        case MARKED:
            // Already marked?  Try to select what the user clicked on, then
            // try to find a path if it worked.
            if (registerClick(state.world, e.getX(), e.getY(),
                              state.worldType)) {
                runSearch(state);
                state.uiState = DRAWN;
            }
        }
    } else if (e.getEventType() == MOUSE_MOVED) {
        // update display of current mouse row/col position to aid testing
        TBLoc loc = coordToLoc(e.getX(), e.getY());
        if (state.world.inBounds(loc.row, loc.col)) {
            string newPositionText = vertexName(loc.row, loc.col, state.world);
            if (newPositionText != gPositionFieldText) {
                gPositionFieldText = newPositionText;
                gPositionField->setText(newPositionText);
            }
        }
    }
}

/*
 * Reacts to a window event in the window.
 */
void processWindowEvent(GWindowEvent e) {
    if (e.getEventType() == WINDOW_CLOSED
            || e.getEventType() == CONSOLE_CLOSED) {
        shutdownProgram();
    }
}

/*
 * Given a physical screen location, maps it to a logical grid location.
 */
static TBLoc coordToLoc(double x, double y) {
    TBLoc loc;
    loc.row = (int) ((y - kMargin) / gPixelsPerHeight);
    loc.col = (int) ((x - kMargin) / gPixelsPerWidth);
    return loc;
}

/*
 * Colors the given location, which has the specified intensity, the indicated
 * color.
 */
static void colorLocation(TBLoc loc, double value, Color locColor) {
    double x, y;
    locToCoord(loc, x, y);

    fillRect(x, y, gPixelsPerWidth, gPixelsPerHeight,
             valueToColor(value, locColor));
}

/*
 * Given a path, returns the cost of that path.
 */
static double costOf(Vector<TBLoc>& path,
                     const Grid<double>& world,
                     double costFn(const TBLoc&, const TBLoc&, const Grid<double>&)) {
    double result = 0.0;
    for (int i = 1; i < path.size(); i++) {
        result += costFn(path[i - 1], path[i], world);
    }
    return result;
}

/*
 * Given a path, draws that path on the screen.
 */
static void drawPath(Vector<TBLoc>& path) {
    for (int i = 1; i < path.size(); i++) {
        // highlight connection between path[i - 1] and path[i]
        double srcX, srcY, dstX, dstY;
        findMidpoint(path[i - 1], srcX, srcY);
        findMidpoint(path[i], dstX, dstY);

        GLine* connection = new GLine(srcX, srcY, dstX, dstY);
        connection->setColor(kPathColor);
        connection->setLineWidth(3.0);
        gWindow->add(connection);
        gHighlightedPath.add(connection);
    }
}

/*
 * Draws the specified grid in the world.
 */
static void drawWorld(Grid<double>& world) {
    if (gWindow == NULL) error("Cannot draw world before window is set up.");

    // This should act as if we are repainting the top-level display, so we
    // will remove all overlays.
    removeMarkedSquares();
    removeOverlays();

    // Recompute the size of a cell in the display.
    gPixelsPerWidth = (double) WINDOW_WIDTH / world.numCols();
    gPixelsPerHeight = (double) WINDOW_HEIGHT / world.numRows();

    // Draw the background.
    fillRect(0, 0,
             WINDOW_WIDTH + 2 * kMargin,
             WINDOW_HEIGHT + 2 * kMargin,
             kBackgroundColor);

    // With the redraw, no locations are marked anymore.
    gMarked.resize(world.numRows(), world.numCols());
    gMarkedValues = world;

    // Draw each cell.
    for (int row = 0; row < world.numRows(); row++) {
        for (int col = 0; col < world.numCols(); col++) {
            TBLoc loc = { row, col };
            double x, y;
            locToCoord(loc, x, y);

            fillRect(x, y, gPixelsPerWidth, gPixelsPerHeight,
                     valueToColor(world[row][col], WHITE));
        }
    }

    // Issue a repaint so that we see the changes.
    gWindow->repaint();
}

/*
 * A convenience wrapper that sets the given color and then fills a rectangle
 * of the given dimensions.  Saves a line of code on each call.
 */
static void fillRect(int x, int y, int width, int height, string color) {
    gWindow->setColor(color);
    gWindow->fillRect(x, y, width, height);
}

/*
 * Given a logical grid location, returns the physical screen coordinates of
 * its midpoint.
 */
static void findMidpoint(TBLoc loc, double& xc, double& yc) {
    locToCoord(loc, xc, yc);
    xc += gPixelsPerWidth / 2;
    yc += gPixelsPerHeight / 2;
}

/*
 * Returns which type of algorithm is currently selected in the drop-down
 * menu.
 */
static AlgorithmType getAlgorithmType() {
    if (gWindow == NULL) error("Window has not yet been initialized.");
    string algorithmLabel = gAlgorithmList->getSelectedItem();

    if (algorithmLabel == "DFS") {
        return DFS;
    } else if (algorithmLabel == "BFS") {
        return BFS;
    } else if (algorithmLabel == "Dijkstra") {
        return DIJKSTRA;
    } else if (algorithmLabel == "A*") {
        return A_STAR;
    } else if (algorithmLabel == "Bidirectional") {
        return BIDIRECTIONAL;
    } else {
        error("Invalid algorithm provided.");
        return DIJKSTRA;
    }
}

/*
 * Returns all files in the current directory that start with the given
 * substring prefix.  Used to find all maze and/or terrain files to display.
 */
static Set<string> getFiles(string substr) {
    substr = toLowerCase(substr);
    Vector<string> files;
    listDirectory(".", files);
    Set<string> result;
    for (string file : files) {
        string fileLC = toLowerCase(file);
        if (startsWith(fileLC, substr) && endsWith(fileLC, ".txt")) {
            result.add(file);
        }
    }
    return result;
}

/*
 * Given the contents of the world size selector, returns a WorldSize encoding
 * the desired world size.
 */
static WorldSize getWorldSize(string worldText) {
    string worldLC = toLowerCase(worldText);
    if (worldLC.find("tiny") != string::npos) {
        return TINY_WORLD;
    } else if (worldLC.find("small") != string::npos) {
        return SMALL_WORLD;
    } else if (worldLC.find("medium") != string::npos) {
        return MEDIUM_WORLD;
    } else if (worldLC.find("large") != string::npos) {
        return LARGE_WORLD;
    } else if (worldLC.find("huge") != string::npos) {
        return HUGE_WORLD;
    } else {
        error("Invalid world size provided.");
        return SMALL_WORLD;
    }
}

/*
 * Prints an introductory text message on the screen.
 */
static void intro() {
    cout << "Welcome to CS 106X Trailblazer!" << endl;
    cout << "This program searches for paths through graphs" << endl;
    cout << "representing mazes and rocky terrains." << endl;
    cout << "It demonstrates several graph algorithms for" << endl;
    cout << "finding paths, such as depth-first search (DFS)," << endl;
    cout << "breadth-first search (BFS), Dijkstra's Algorithm," << endl;
    cout << "and A* search. You can also generate random mazes" << endl;
    cout << "using Kruskal's algorithm." << endl;
}

/*
 * Restores the previously saved GUI state, including which algorithm is
 * currently selected, the animation delay, and the world file chosen.
 * If the saved state does not exist or is corrupt, returns false and
 * uses a default initial state.
 */
static bool loadGUIState() {
    ifstream input;
    input.open(GUI_STATE_FILE.c_str());
    if (input.fail()) {
        return false;
    }
    string algorithm;
    getline(input, algorithm);
    if (input.fail()) {
        return false;
    }

    string line;
    getline(input, line);
    if (input.fail()) {
        return false;
    }
    setAnimationDelay(stringToInteger(line));

    string worldFile;
    getline(input, worldFile);
    if (input.fail()) {
        return false;
    }
    input.close();

    // delete the save state file in case there is a crash loading a world
    deleteFile(GUI_STATE_FILE);

    gAlgorithmList->setSelectedItem(algorithm);
    gDelaySlider->setValue(animationDelay);
    if (worldFile != OTHER_FILE_LABEL) {
        gWorldChooser->setSelectedItem(worldFile);
    }
    return true;
}

/*
 * Given a logical grid location, maps it to a physical screen location.
 */
static void locToCoord(TBLoc loc, double& x, double& y) {
    x = loc.col * gPixelsPerWidth + kMargin;
    y = loc.row * gPixelsPerHeight + kMargin;
}

/*
 * Tries to read a world file from the specified stream.  On success, returns
 * true and updates the input parameters to mark the type of the world and
 * the world contents.  On failure, returns false, but may still modify the
 * input parameters.
 */
static bool readWorldFile(ifstream& input, Grid<double>& world, WorldType& worldType) {
    try {
        // Enable exceptions on the stream so that we can handle errors using try-
        // catch rather than continuously testing everything.
        input.exceptions(ios::failbit | ios::badbit);

        // The file line of the file identifies the type, which should be either
        // "terrain" or "maze."
        string type;
        input >> type;

        if (type == "terrain") {
            worldType = TERRAIN_WORLD;
        } else if (type == "maze") {
            worldType = MAZE_WORLD;
        } else {
            cerr << "world file does not contain type (terrain/maze) as first line." << endl;
            return false;
        }

        // Read the size of the world.
        int numRows, numCols;
        input >> numRows >> numCols;

        if (numRows <= 0 || numCols <= 0 ||
                numRows >= kMaxRows || numRows >= kMaxCols) {
            cerr << "world file contains invalid number of rows/cols: "
                 << numRows << "," << numCols << endl;
            return false;
        }

        world.resize(numRows, numCols);

        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                double value;
                input >> value;
                if (input.fail()) {
                    cerr << "Illegal input file format; row #" << (row+1)
                         << "does not contain " << numCols << " valid numbers" << endl;
                    return false;
                }

                // Validate the input based on the type of world.
                if (worldType == MAZE_WORLD) {
                    if (value != kMazeWall && value != kMazeFloor) {
                        cerr << "world file contains invalid square value of " << value
                             << ", must be " << kMazeFloor << " or " << kMazeWall << endl;
                        return false;
                    }
                } else {  // worldType == TERRAIN_WORLD
                    if (value < 0.0 || value > 1.0) {
                        cerr << "world file contains invalid terrain value of " << value
                             << ", must be 0.0 - 1.0" << endl;
                        return false;
                    }
                }
                world[row][col] = value;
            }
        }

        return true;
    } catch (...) {
        // Something went wrong, so report an error.
        cerr << "exception thrown while reading world file" << endl;
        return false;
    }
}

/*
 * Generates a new world based on the user's preferences.
 */
static bool regenerateWorld(Grid<double>& world, WorldType& worldType) {
    string worldFile = gWorldChooser->getSelectedItem();
    if (startsWith(worldFile, "Random")) {
        // generate a random maze
        cout << endl;
        cout << "Generating a random maze ..." << endl;
        int size = kMazeNumCols[getWorldSize(worldFile)];
        
        // a maze of size M x N has a logical size of (M/2+1 x N/2+1)
        // for the purposes of Kruskal's algorithm
        world = createRandomMaze(size / 2 + 1);
        worldType = MAZE_WORLD;
        worldUpdated(world, worldType);
        return true;
    } else {
        // not random (most common); load a pre-defined maze or terrain input file
        cout << endl;
        if (worldFile == OTHER_FILE_LABEL) {
            // prompt for file name
            worldFile = GFileChooser::showOpenDialog(getCurrentDirectory());
        }

        if (worldFile.empty() || !fileExists(worldFile)) {
            cout << "File not found; aborting." << endl;
            return false;
        }

        cout << "Loading world from " << getTail(worldFile) << " ..." << endl;
        ifstream input;
        input.open(worldFile.c_str());
        if (input.fail()) {
            error(string("Unable to open input file ") + worldFile);
            return false;
        }

        // Try reading in the world file.  If we can't, report an error.
        Grid<double> newWorld;
        WorldType newWorldType;
        if (!readWorldFile(input, newWorld, newWorldType)) {
            cerr << worldFile << " is not a valid world file." << endl;
            return false;
        }

        world = newWorld;
        worldType = newWorldType;
        worldUpdated(world, worldType);
        return true;
    }
}

/*
 * Registers that a click occurred, storing the location appropriately and
 * updating the graphics.  If the click was on an invalid location, returns
 * false.  If the click was at a valid location, returns true.
 */
static bool registerClick(Grid<double>& world, double x, double y,
                          WorldType worldType) {
    if (gWindow == NULL) {
        error("Window has not yet been initialized.");
    }
    if (gFirstSelected != NULL && gSecondSelected != NULL) {
        error("Two tiles have already been selected.");
    }

    // Determine where we clicked and make sure it's in-bounds.
    TBLoc loc = coordToLoc(x, y);
    if (!world.inBounds(loc.row, loc.col) ||
            (worldType == MAZE_WORLD && world[loc.row][loc.col] == kMazeWall)) {
        return false;
    }

    // Add the selection to the display.
    double selectionX, selectionY;
    locToCoord(loc, selectionX, selectionY);
    GRect* selection = new GRect(selectionX, selectionY, gPixelsPerWidth, gPixelsPerHeight);
    selection->setLineWidth(2.0);
    selection->setColor(kSelectedLocationColor);
    gWindow->add(selection);

    // Store the selection appropriately based on whether this is the first or
    // second click.
    if (gFirstSelected == NULL) {
        gFirstSelected = selection;
        gStartLocation = loc;
    } else {
        gSecondSelected = selection;
        gEndLocation = loc;
    }

    return true;
}

/*
 * Utility function to remove a GObject from the display window and deallocate
 * it.
 */
static void removeAndDelete(GObject* object) {
    if (object != NULL) {
        gWindow->remove(object);
        delete object;
    }
}

/*
 * Removes the highlighted squares from the display.
 */
static void removeMarkedSquares() {
    removeAndDelete(gFirstSelected);
    removeAndDelete(gSecondSelected);
    gFirstSelected = gSecondSelected = NULL;
}

/*
 * Removes the line overlays.
 */
static void removeOverlays() {
    // Clear the lines drawn from the path.
    for (int i = 0; i < gHighlightedPath.size(); i++) {
        removeAndDelete(gHighlightedPath[i]);
    }
    gHighlightedPath.clear();
    if (!gWindow->isRepaintImmediately()) {
        gWindow->repaint();
    }
}

/*
 * Clears all cells that are currently selected and deselects any currently-
 * selected squares.
 */
static void restoreWorldDisplay() {
    if (gWindow == NULL) {
        error("Window has not yet been initialized.");
    }

    uncolorSquares();
    removeMarkedSquares();
    removeOverlays();

    // Repaint the window so the changes show.
    gWindow->repaint();
}

/*
 * Runs the currently selected graph path-searching algorithm on the current
 * world, then displays information about the path that was found.
 */
static void runSearch(State& state) {
    try {
        updateAnimationDelayFromSlider();
        Vector<TBLoc> path;
        double pathCost;
        runShortestPath(state.world, state.worldType, gStartLocation,
                        gEndLocation, path, pathCost);
        
        if (animationDelay == 0) {
            // GUI will not have repainted itself to show the path being drawn;
            // manually repaint it
            gWindow->repaint();
        }
        
        cout << "Path length: " << path.size() << endl;
        cout << "Path cost: " << pathCost << endl;
        int greenGray = 0;
        int yellow = 0;
        for (int r = 0; r < gMarked.numRows(); r++) {
            for (int c = 0; c < gMarked.numCols(); c++) {
                Color color = gMarked.get(r, c);
                if (color == GREEN || color == GRAY) {
                    greenGray++;
                } else if (color == YELLOW) {
                    yellow++;
                }
            }
        }
        cout << "Locations visited: " << greenGray << endl;

    } catch (const ErrorException& e) {
        cout << e.what() << endl;
    }
}

/*
 * Computes the shortest path between the start and end locations, displaying
 * it on the screen and returning its length.
 */
static void runShortestPath(const Grid<double>& world, WorldType worldType,
                            const TBLoc& start, const TBLoc& end,
                            Vector<TBLoc>& path, double& pathCost) {
    AlgorithmType algType = getAlgorithmType();

    // Determine which cost/heuristic functions to use.
    double (*costFn)(const TBLoc&, const TBLoc&, const Grid<double>&);
    double (*hFn)(const TBLoc&, const TBLoc&, const Grid<double>&);

    if (worldType == TERRAIN_WORLD) {
        costFn = terrainCost;
        hFn = terrainHeuristic;
    } else if (worldType == MAZE_WORLD) {
        costFn = mazeCost;
        hFn = mazeHeuristic;
    } else {
        error("Unknown world type.");
    }

    // Invoke the student's shortestPath function to find out the cost of the path.
    pathSearchInProgress = true;
    path = shortestPath(start, end, world, costFn, hFn, algType);
    pathSearchInProgress = false;

    if (path.isEmpty()) {
        cout << "No path was found. (The returned path is empty.)" << endl;
    } else if (path[0] != start) {
        cout << "Warning: Start of path is not the start location." << endl;
    } else if (path[path.size() - 1] != end) {
        cout << "Warning: End of path is not the end location." << endl;
    }

    drawPath(path);
    pathCost = costOf(path, world, costFn);
}

/*
 * Saves the GUI's current state, including which algorithm is
 * currently selected, the animation delay, and the world file chosen.
 * If the saved state does not exist or is corrupt, returns false and
 * uses a default initial state.
 */
static bool saveGUIState() {
    string algorithm = gAlgorithmList->getSelectedItem();
    int delay = gDelaySlider->getValue();
    string worldFile = gWorldChooser->getSelectedItem();
    ofstream output;
    output.open(GUI_STATE_FILE.c_str());
    if (output.fail()) {
        return false;
    }
    output << algorithm << endl;
    output << delay << endl;
    output << worldFile << endl;
    if (output.fail()) {
        return false;
    }
    output.flush();
    output.close();
    return true;
}

/*
 * Sets the delay between animation frames to the given value.
 */
static void setAnimationDelay(int delay) {
    int oldDelay = animationDelay;
    animationDelay = delay;

    // tell the window whether or not to repaint on every square colored
    if ((animationDelay == 0) != (oldDelay == 0)) {
        if (gWindow) {
            gWindow->setRepaintImmediately(animationDelay != 0);
        }
    }
}

/*
 * Reads the delay slider and uses its value to decide on an animation delay.
 */
static void updateAnimationDelayFromSlider(bool forbidZero) {
    int delay = gDelaySlider->getValue();
    double percent = 100.0 * delay / MAX_ANIMATION_DELAY;
    // convert scale so delays don't increase so rapidly
    if (percent == 0.0) {
        delay = forbidZero ? 1 : 0;
    } else if (percent <= 10) {
        delay = MAX_ANIMATION_DELAY / 1000;
    } else if (percent <= 20) {
        delay = MAX_ANIMATION_DELAY / 500;
    } else if (percent <= 30) {
        delay = MAX_ANIMATION_DELAY / 200;
    } else if (percent <= 40) {
        delay = MAX_ANIMATION_DELAY / 100;
    } else if (percent <= 50) {
        delay = MAX_ANIMATION_DELAY / 50;
    } else if (percent <= 60) {
        delay = MAX_ANIMATION_DELAY / 25;
    } else if (percent <= 70) {
        delay = MAX_ANIMATION_DELAY / 10;
    } else if (percent <= 80) {
        delay = MAX_ANIMATION_DELAY / 5;
    } else if (percent <= 90) {
        delay = MAX_ANIMATION_DELAY / 2;
    } else {  // percent > 90
        delay = MAX_ANIMATION_DELAY;
    }

    setAnimationDelay(delay);
}

/*
 * Given an elevation and a color, determines what #rrggbb format color to use
 * to draw it.
 */
static string valueToColor(double value, Color locColor) {
    // To make non-gray colors show up better, we artificially cap the lowest
    // possible intensity at 0.2, rather than 0.0.  This next step is a linear
    // remapping of the range [0, 1] to [0.2, 1].
    if (locColor != WHITE) {
        value = 0.8 * value + 0.2;
    }

    stringstream hexValue;
    hexValue << "#" << hex << setfill('0');
    for (int i = 0; i < 3; i++) {
        int intensity = int(value * kColorMultipliers[locColor][i]);
        hexValue << setw(2) << intensity;
    }

    return hexValue.str();
}

/*
 * Restores all squares to their original colors.
 */
static void uncolorSquares() {
    // Restore all colored cells to their original condition.
    for (int row = 0; row < gMarked.numRows(); row++) {
        for (int col = 0; col < gMarked.numCols(); col++) {
            if (gMarked[row][col] != UNCOLORED && gMarked[row][col] != WHITE) {
                TBLoc loc = { row, col };
                colorLocation(loc, gMarkedValues[row][col], WHITE);

                // Unmark this cell; it's no longer colored.
                gMarked[row][col] = UNCOLORED;
            }
        }
    }
}

/*
 * Called anytime the current world is changed, so that we can update the
 * cache of Grid -> Graph.
 */
static void worldUpdated(Grid<double>& world, WorldType worldType) {
    flushWorldCache();
    if (worldType == TERRAIN_WORLD) {
        ensureWorldCache(world, terrainCost);
    } else if (worldType == MAZE_WORLD) {
        ensureWorldCache(world, mazeCost);
    } else {
        error("Unknown world type.");
    }
}

/*
 * Closes the console and GUI, saves current world state, and exits
 * the program.
 */
static void shutdownProgram() {
    // Trying to quit exits the entire program.
    cout << endl;
    cout << "Exiting." << endl;
    if (SHOULD_SAVE_GUI_STATE) {
        saveGUIState();
    }
    exitGraphics();
}
