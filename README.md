This project implementes the [Affine Particle-In-Cell](https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf) method.

The 2D simulator starter code and renderer is provided by Doug James [here](http://graphics.stanford.edu/courses/cs348c/PA3_APIC2017/index.html). The starter code itself is based on [Robert Bridson's 2D simulator](https://www.cs.ubc.ca/%7Erbridson/).

# Dependency
Requires gcc and OpenGL.

# How to Build
1. Compile the simulator
    1. Run `make depend` at root folder.
    2. Run `make` to build. Or use `make_release` to build the release version.
2. Compile the renderer
    1. `cd` to `viewer` folder.
    2. Run `make depend`.
    3. Run `make` to build

# How to Run
1. Simulator
    1. Use `./simulate <output-folder> <simulation-type>` (or `./simulate_release <output-folder> <simulation-type>`) to run. The simulator will generate one output file per frame into the output folder. Simulation type should be one of `pic`, `apic` or `flip`. If no argument is provided, output will be saved at root folder, and simulatrion type is deafult to `apic`.
2. Viewer
    1. From the `viewer` folder, use `./view <input-folder>` to run the renderer. The input folder should be the same as the output folder of the simulator.
    2. When the viewer is running, use the left and right arrow keys to move between frames.

Note: See the original starter code's compile instruction [here](https://docs.google.com/document/d/1cdG9zB3fslVxtV5L-pFG9xcujWVvxdeyQRKzzXbp4OI/edit).
