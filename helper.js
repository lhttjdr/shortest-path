const plot_circles = ( canvas, circles, inner_circle, inside_convex_circle ) => {
  for(let i=0;i<circles.length;i++){
    if( inner_circle && inner_circle[ i ] ) {
      canvas.plot( circles[ i ], "rgba(255,255,153,0.5)", 1,
        "rgba(255,255,153,0.5)" );
    } else if( inside_convex_circle && inside_convex_circle[ i ] ) {
      canvas.plot( circles[ i ], "rgba(204,255,153,0.5)", 1,
        "rgba(204,255,153,0.5)" );
    } else {
      canvas.plot( circles[ i ], "rgba(155,255,255,0.5)", 1,
        "rgba(155,255,255,0.5)" );
    }
  }
}

const plot_start_end = ( canvas, start, end ) => {
  canvas.plot( start, 'blue', 2, 'red' );
  canvas.plot( end, 'blue', 2, 'red' );
}

const plot_tree = ( canvas, circles, blocks ) => {
  for( let block of blocks ) {
    if( block.vertices.length === 1 ) {
      continue;
    }
    for( let e of block.edges ) {
      canvas.plot( new Segment( circles[ e[ 0 ] ].center, circles[ e[ 1 ] ].center ),
        "green" );
    }
  }
}

const plot_circle_convex_hull = ( canvas, circles, convex_hull ) => {
  let points = convex_hull.points;
  for( let c of convex_hull.circles ) canvas.plot( circles[ c ] );
  for( let s of convex_hull.segments ) canvas.plot( new Segment( points[ s[
    2 ] ], points[ s[ 3 ] ] ) );
}

const plot_circle_convex_hulls = ( canvas, circles, circle_convex_hulls ) => {
  for( let convex_hull of circle_convex_hulls ) {
    plot_circle_convex_hull( canvas, circles, convex_hull );
  }
}

const plot_polygons = ( canvas, circle_convex_hulls, polygons ) => {
  for( let i = 0, n = polygons.length; i < n; i++ ) {
    let points = polygons[ i ].map( c => circle_convex_hulls[ i ].points[ c ] );
    if( ConvexPolygon.check_points( points ) ) {
      let poly = new ConvexPolygon( points );
      canvas.plot( poly, 'Chocolate', 2, 'rgba(189,183,107,0.3)' );
    }
  }
}

const plot_triangles = ( canvas, points, triangles, baseTriangle ) => {
  let base = baseTriangle.points;
  //triangles.filter(t=>t.some(p=>p<0)).forEach(triangle=>{
  //  myGraph.plot(new Triangle(triangle.map(p=>p<0? base[p+3]: points[p])),'rgba(8, 239, 4, 0.67)',0.5);
  //});
  triangles.filter( t => t.every( p => p >= 0 ) )
    .forEach( triangle => {
      canvas.plot( new Triangle( triangle.map( p => p < 0 ? base[ p + 3 ] :
        points[ p ] ) ), 'rgba(200, 0, 200, 0.67)', 0.5 );
    } );
}