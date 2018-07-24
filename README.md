# shortest path problem with circle obstacles

**ongoing project**

[Demo](https://lhttjdr.github.io/shortest-path/)

- [X] plot tools (HTML5 canvas)
- [X] adjacent map of obstacles (can not walk through "wall")
- [X] convex hull of circles (avoid enter convex hull unless start point or destination point are in them)
 - [An algorithm for constructing the convex hull of a set of spheres in dimension d (1996, Jean-DanielBoissonnat et al.)](https://doi.org/10.1016/0925-7721(95)00024-0)
 - incremental convex hull 3d
- [X] Delaunay triangulation (**current task:** replace it by constrained Delaunay triangulation)
- [ ] A* search for shortest path
