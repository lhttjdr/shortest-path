let circles = remove_duplicate_circles( get_circles( 30 ) );
let [ start, end ] = get_start_end( circles );

plot_scene_1( myGraph, circles, start, end );
// 将start和end看作是两个半径为0的圆，避免以后分情况讨论
circles.push( new Circle( start, 0 ) );
circles.push( new Circle( end, 0 ) );
// 计算邻接矩阵adjacent，相交或相切的圆为邻接，相离的圆为不邻接。
// 内含的圆会被剔除，保存在inner_circle
let [ adjacent, inner_circle ] = get_adjacent( circles );
// 根据邻接矩阵adjacent计算联通分量blocks
// 每个圆所在的分量索引circle_in_block
let [ blocks, circle_in_block ] = get_blocks( circles, adjacent, inner_circle );

let [ circle_convex_hulls, inside_convex_circle ] = get_circle_convex_hulls(
  circles, blocks );
let ploygons = get_polygons( circle_convex_hulls );
[ circle_convex_hulls, polygons ] = remove_inside_segments( circle_convex_hulls,
  ploygons );

// plot_circles( myGraph, circles, inner_circle, inside_convex_circle );
// plot_circle_convex_hulls( circles, circle_convex_hulls );
// plot_polygons( circle_convex_hulls, ploygons );
//  plot_tree(circles,blocks);
// plot_st( start, end );


let points = get_points( circle_convex_hulls )
  .concat( circles.map( c => c.center ) );
// circles=circles.filter((c,i)=>!inner_circle[i]);
// let centers=circles.map((c,i)=>c.center);
let [ triangles, baseTriangle ] = Geometry2D.delaunay( points );

//plot_triangles( points, triangles, baseTriangle );