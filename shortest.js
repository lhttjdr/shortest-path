/****tools****/
let getRandom = Math.random;
let getRandomArbitary = ( min, max ) => getRandom() * ( max - min ) + min;
let getRandomInt = ( min, max ) => Math.floor( getRandom() * ( max - min + 1 ) ) +
  min;

let myGraph = new Graph( {
  canvasId: 'myCanvas',
  minX: -900,
  minY: -600,
  maxX: 900,
  maxY: 600,
  unitsPerTick: 100
} );

/********data***********/
const get_circles = ( n ) => {
  let circles = [];
  for( let i = 0; i < n; i++ ) {
    let [ x, y, r ] = [ getRandomArbitary( -700, 700 ), getRandomArbitary( -
      400,
      400 ), getRandomArbitary( 50, 150 ) ];
    c = new Circle( x, y, r );
    circles.push( c );
  }
  return circles;
}
const remove_duplicate_circles = ( circles ) => circles.filter( ( c, i, cs ) =>
  cs.findIndex( ( t ) => t.radius ===
    c.radius && t.center.equals( c.center ) === i ) );
const get_start_end = ( circles ) => {
  let start, end;
  do {
    start = new Point( getRandomArbitary( -700, 700 ), getRandomArbitary( -
      400,
      400 ) );
  } while ( circles.some( c => c.includes( start ) ) );
  do {
    end = new Point( getRandomArbitary( -700, 700 ), getRandomArbitary( -400,
      400 ) );
  } while ( circles.some( c => c.includes( end ) ) );
  return [ start, end ];
}

const get_polygons = function ( circle_convex_hulls ) { // 广义的，包括线段和单独的点
  let polygons = [];
  for( let circle_convex_hull of circle_convex_hulls ) {
    let points = circle_convex_hull.points;
    polygons.push( Geometry2D.GrahamScan( points ) ); //返回的是点的下标
  }
  return polygons;
}

const remove_inside_segments = function ( circle_convex_hulls, polygons ) {
  for( let i = 0, n = circle_convex_hulls.length; i < n; i++ ) {
    let points = circle_convex_hulls[ i ].points;
    let segments = circle_convex_hulls[ i ].segments;
    let reverse_dict = Array( points.length )
      .fill( -1 );
    let valid_points = polygons[ i ];
    for( let i = 0; i < valid_points.length; i++ ) reverse_dict[ valid_points[
      i ] ] = i;
    // 注意，仅仅为了polygons和circle_convex_hulls要配套。点调整为凸包顺序后，polygons其实就没用了。
    polygons[ i ] = Array.from( valid_points.keys() );
    circle_convex_hulls[ i ].points = valid_points.map( x => points[ x ] ); //点变为凸包顺序
    valid_points = new Set( valid_points );
    circle_convex_hulls[ i ].segments = segments.filter( s => valid_points.has(
        s[ 2 ] ) && valid_points.has( s[ 3 ] ) )
      .map( s => [ s[ 0 ], s[ 1 ], reverse_dict[ s[ 2 ] ], reverse_dict[ s[ 3 ] ] ] );
  }
  return [ circle_convex_hulls, polygons ];
}

const get_points = function ( circle_convex_hulls ) {
  return circle_convex_hulls.reduce( ( set, hull ) => set.concat( hull.points ), [] );
}

const get_adjacent = function ( circles ) {
  let N = circles.length;
  let adjacent = Array( N )
    .fill()
    .map( ( r ) => [] );
  let inner_circle = Array( N )
    .fill( false );
  for( let i = 0; i < N; i++ ) {
    for( let j = i + 1; j < N; j++ ) {
      let d = Geometry2D.distance( circles[ i ].center, circles[ j ].center );
      let rdiff = circles[ i ].radius - circles[ j ].radius;
      if( d <= Math.abs( rdiff ) ) { // 内含或内切
        if( circles[ i ].radius <= circles[ j ].radius ) {
          adjacent[ j ].push( i );
          inner_circle[ i ] = true;
        } else {
          adjacent[ i ].push( j );
          inner_circle[ j ] = true;
        }
        continue;
      }
      let rsum = circles[ i ].radius + circles[ j ].radius;
      if( d <= rsum ) {
        adjacent[ i ].push( j );
        adjacent[ j ].push( i );
      }
    }
  }
  return [ adjacent, inner_circle ];
}

const get_blocks = function ( circles, adjacent, inner_circle ) {
  let N = circles.length; // number of circles
  let blocks = []; // spanning tree of connected component
  let circle_in_block = []; // Which connected component a circle in
  // DFS
  let vis = Array( N )
    .fill( false ); // visited flag
  for( let i = 0; i < N; i++ ) { // try a new circle
    // visited circle or invalid circle (=circle in another circle)
    if( vis[ i ] || inner_circle[ i ] ) continue;
    // DFS tree = Connected Component
    let tree = {};
    let blockId = blocks.length; // label of current connected component
    // add circle i to current connected component
    tree.vertices = [ i ];
    tree.edges = [];
    let stack = [ i ];
    circle_in_block[ i ] = blockId;
    vis[ i ] = true;
    // loop while DFS stack is non-empty
    while( stack.length > 0 ) {
      let now = stack[ stack.length - 1 ]; // top
      let expand = false; // suppose children not exist
      for( let j of adjacent[ now ] ) { // check adjacent circle j
        if( !vis[ j ] && !inner_circle[ j ] ) {
          expand = true;
          // add circle j to current connected component
          tree.edges.push( [ now, j ] );
          tree.vertices.push( j );
          stack.push( j );
          vis[ j ] = true;
          circle_in_block[ j ] = blockId;
        }
      }
      // if no children, backtrack
      if( !expand ) stack.pop();
    }
    blocks.push( tree );
  }
  return [ blocks, circle_in_block ];
}

const get_circle_convex_hulls = function ( circles, blocks ) {
  let circle_convex_hulls = [];
  let inside_convex_circle = Array( circles.length )
    .fill( true ); // 第二类无效圆, 起始点所在的block，若可见则有效，否则无效
  for( let block of blocks ) {
    let circle_convex_hull = new Object();
    if( block.vertices.length === 1 ) { // 单独的圆
      inside_convex_circle[ block.vertices[ 0 ] ] = false;
      circle_convex_hull.circles = [ block.vertices[ 0 ] ];
      circle_convex_hull.segments = [];
      circle_convex_hull.points = [];
      circle_convex_hulls.push( circle_convex_hull );
      continue;
    }
    let block_circles_idx = [ ...new Set( block.edges.reduce( ( u, e ) => u.concat(
      e ), [] ) ) ];
    let block_circles = block_circles_idx.map( ( x ) => circles[ x ] );
    let convex;
    convex = Geometry2D.circle_convex( block_circles );
    if( typeof convex === "number" ) { // 仅有一个圆，collinear or concurrent
      inside_convex_circle[ block_circles_idx[ convex ] ] = false;
      circle_convex_hull.arcs = [ block_circles_idx[ convex ] ];
      circle_convex_hull.segments = [];
    } else if( Array.isArray( convex[ 0 ] ) ) { // 圆的对
      circle_convex_hull.circles = convex.reduce( ( a, b ) => a.concat( b ), [] )
        .filter( ( c, i, a ) => a.indexOf( c ) === i )
        .map( bi => block_circles_idx[ bi ] );
      let convex_segments = [];
      let points = []; // 这里最适合用一个二叉树，查找为O(log n)。目前是O(n)
      for( let i = 0, n = convex.length; i < n; i++ ) {
        let a = block_circles_idx[ convex[ i ][ 0 ] ],
          b = block_circles_idx[ convex[ i ][ 1 ] ];
        inside_convex_circle[ a ] = inside_convex_circle[ b ] = false;
        let segments = Geometry2D.common_tangent_segments( circles[ a ],
          circles[ b ] );
        if( segments.length > 0 ) {
          segments = segments.slice( 0, 2 );
          for( let s of segments ) {
            let id = points.findIndex( p => p.equals( s.points[ 0 ] ) );
            if( id === -1 ) {
              id = points.length;
              points.push( s.points[ 0 ] );
            }
            let id_nxt = points.findIndex( p => p.equals( s.points[ 1 ] ) );
            if( id_nxt === -1 ) {
              id_nxt = points.length;
              points.push( s.points[ 1 ] );
            }
            convex_segments.push( [ a, b, id, id_nxt ] ); // 从圆a到圆b，点的编号id, id+1
          }
        }
      }
      circle_convex_hull.points = points;
      circle_convex_hull.segments = convex_segments;
    }
    circle_convex_hulls.push( circle_convex_hull );
  }
  return [ circle_convex_hulls, inside_convex_circle ];
}

const intersect_with_convex = function ( circle_convex_hull, segment ) {
  let points = circle_convex_hull.points;
  circle_convex_hull.segments.forEach( pair => {
    let edge = new Segment( points[ 2 ], points[ 3 ] );
    if( Geometry2D.intersect( segment, edge ) ) return true;
  } );
  return false;
}

const get_search_map = function ( circles, triangles ) {
  let search_map = Array( circles.length )
    .fill()
    .map( r => new Array() );
  triangles.forEach( ( triangle ) => {
    let pair = triangle.map( ( p, i, t ) => t.map( p1 => [ p, p1 ] ) )
      .reduce( ( a, b ) => a.concat( b ) );
    pair.forEach( pr => {
      if( pr[ 0 ] > 0 && pr[ 1 ] > 0 && pr[ 0 ] !== pr[ 1 ] )
        search_map[ pr[ 0 ] ].push( pr[ 1 ] );
    } );
  } );
  search_map.forEach( ( adj, cur ) => {
    let base = circles[ cur ];
    adj.sort( ( a, b ) => {
      return Geometry2D.distance( base, circles[ a ] ) - Geometry2D.distance(
        base, circles[ b ] );
    } );
  } );
  return search_map;
}

const get_visable = function ( circles, i, search_map ) {
  let n = circles.length;
  let base = circles[ i ];
  let visable = [];
  let used = Array( n )
    .fill( false );
  let Q = [ i ];
  used[ i ] = true;
  let head = 0;
  while( head < Q.length ) {
    let cur = Q[ head++ ];
    for( let adj of search_map[ cur ] ) {
      if( used[ adj ] ) continue;
      let segments = Geometry2D.common_tangent_segments( base, circles[ adj ] );
      let visable_adj = [];
      console.log( segments );
      for( let segment of segments ) {
        let not_blocked = true;
        for( let k = 1; k < Q.length; k++ ) { // 已搜过的圆
          if( Geometry2D.intersect( segment, circles[ Q[ k ] ] ) ) { // 被圆挡住
            not_blocked = false;
            break;
          }
        }
        if( not_blocked ) {
          visable_adj.push( [ adj, segment ] );
        }
      }
      if( visable_adj.length > 0 ) {
        visable.push( visable_adj );
        Q.push( adj );
        used[ adj ] = true;
      }
    }
  }
  return visable;
}

let search_map = get_search_map( circles, triangles );
let visable = get_visable( circles, 0, search_map );
// 按距离排序后，总是近的遮住远的。O(n log n)

const plot_visable = function ( visable ) {
  visable.forEach( vis => {
    vis.forEach( segment => {
      myGraph.plot( segment[ 1 ], "red", 2 );
    } );
  } );
}

plot_visable( visable );

const neighbors = function ( cur ) {
  if( cur === start ) {

  } else {
    let c1 = circles[ cur ];
    for( let i = 0; i < circles.length; i++ ) {
      if( inner_circle[ i ] ) continue;

    }
  }
}

const A_star = function ( start ) {
  let frontier = new PriorityQueue( ( a, b ) => a[ 1 ] < b[ 1 ] );
  frontier.init( [ start, 0 ] );
  let came_from = [];
  let cost_so_far = [];
  came_from[ start ] = start;
  cost_so_far[ start ] = 0;

  while( !frontier.empty() ) {
    let [ current, dist ] = frontier.remove();
    if( current === goal ) {
      break;
    }
    let visable = graph.neighbors( current );
    for( let next of visable ) {
      let new_cost = cost_so_far[ current ] + graph.cost( current, next );
      if( !cost_so_far.count( next ) || new_cost < cost_so_far[ next ] ) {
        cost_so_far[ next ] = new_cost;
        let priority = new_cost + heuristic( next, goal );
        frontier.insert( next, priority );
        came_from[ next ] = current;
      }
    }
  }
}

const PriorityQueue = function ( compare ) {
  var compare = compare || function ( a, b ) {
      return a < b;
    },
    heap = [],
    slice = Array.prototype.slice;

  // Initialize the priority queue with an array of items.
  this.init = function ( items ) {
    if( heap.length ) heap = [];
    this.insert.apply( null, items );
    return this;
  };

  // Insert one or more items into the priority queue.
  this.insert = function () {
    slice.call( arguments, 0 )
      .forEach( function ( item ) {
        heap.push( item );
        float( heap.length - 1 );
      } );
    return this;
  };

  // Remove and return the next queue element.
  this.remove = function () {
    var head = null,
      last;
    if( heap.length ) {
      head = heap[ 0 ];
      last = heap.length - 1;
      swap( 0, last );
      heap.splice( last, 1 );
      sink( 0 );
    }
    return head;
  };

  // Clear all items from the queue.
  this.clear = function () {
    heap = [];
    return this;
  };

  // Return the top item without modifying the queue.
  this.peek = function () {
    return heap[ 0 ];
  };

  // Return the number of items in the queue.
  this.size = function () {
    return heap.length;
  };

  // Return true if the queue is empty.
  this.empty = function () {
    return !heap.length;
  };

  // Sink an item from the top down to heap order.
  function sink( h ) {
    var left = 2 * h + 1,
      right = 2 * h + 2,
      child = left;
    if( left < heap.length ) {
      right < heap.length && compare( heap[ right ], heap[ left ] ) && child++;
      if( compare( heap[ child ], heap[ h ] ) ) {
        swap( h, child );
        sink( child );
      }
    }
  }

  // Float an item from the bottom up to heap order.
  function float( h ) {
    var parent = Math.floor( ( h - 1 ) / 2 );
    if( parent >= 0 && compare( heap[ h ], heap[ parent ] ) ) {
      swap( h, parent );
      float( parent );
    }
  }

  // Swap two heap elements by index.
  function swap( i, j ) {
    var tmp = heap[ i ];
    heap[ i ] = heap[ j ];
    heap[ j ] = tmp;
  }
}