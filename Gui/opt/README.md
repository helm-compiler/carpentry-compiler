# Memory-efficient E-Graph with Multi-objective Optimization in C++

In this project, we developed a C++ memory-efficient E-Graph with multi-objective optimization supports for SIGGRAPH Asia 2019 submission.

## Usage

* Use FreeCAD - Carpentry to populate e-graphs to `D:\e-graph`
* Modify `.\scripts\parse.py` `Line 408: taskName = "your task name"`, run this script by `python parse.py`
* After step-2, you can see the creation of `.\benchmarks\your task name`. Then execute compiled `EGraph.exe`, input your task name.
* Results will be generated in `.\benchmarks\your task name\collapse\result.xml`.
* You can explore other sciripts in `.\scripts\` such as `plot.py` to visualize the generated `result.xml`.

## Known problems

- [ ] All `References` are created since we assign concrete floating numbers to its `vector`
