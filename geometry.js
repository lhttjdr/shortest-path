class Point {
  // n coordinates can be given by
  // 1. an array type argument with n numbers
  // 2. n arguments that all is number type
  constructor( x ) {
    this.x = Array.isArray( x ) ? x : Array.from( arguments );
    this._type_guard();
  }
  coordinates() {
    return this.x;
  }
  dimension() {
    return this.x.length;
  }
  move( v ) {
    let vx = v.components();
    return new Point( this.x.map( ( x, i ) => x + vx[ i ] ) );
  }
  static ORIGIN( n ) {
    return new Point( Array( n )
      .fill( 0 ) );
  }
  equals( p ) {
    let px = p.coordinates();
    for( let i = 0, n = px.length; i < n; i++ ) {
      if( this.x[ i ] !== px[ i ] ) return false;
    }
    return true;
  }
  toString() {
    return "(" + this.x.join( "," ) + ")";
  }
  static isPoint( p ) {
    return p instanceof Point;
  }
  // private
  _type_guard() {
    if( !this.x.every( x => typeof x === "number" ) ) {
      throw new TypeError( "Point must be an array of numbers!" );
    }
  }
};

class Vector {
  // n components can be given by
  // 1. an array type argument with n numbers
  // 2. n arguments that all is number type
  // 3. two points in n-dimensional space
  constructor( x ) {
    this.x = Array.isArray( x ) ? x : Array.from( arguments );
    if( this.x.length === 2 && this.x.every( x => x instanceof Point ) ) {
      let [ a, b ] = this.x;
      if( a.dimension() !== b.dimension() ) {
        throw new Error( "Two points must have same dimension!" );
      }
      let [ ax, bx ] = [ a.coordinates(), b.coordinates() ];
      this.x = bx.map( ( x, i ) => x - ax[ i ] );
    }
    this._type_guard();
  }
  components() {
    return this.x;
  }
  add( v ) {
    this._dimension_check( v );
    return new Vector( this.x.map( ( x, i ) => x + v.x[ i ] ) );
  }
  mult( a ) {
    return new Vector( this.x.map( ( x ) => a * x ) );
  }
  negated() {
    return this.mult( -1 );
  }
  minus( v ) {
    this._dimension_check( v );
    return this.add( v.negated() );
  }
  divide( a ) {
    return this.mult( 1.0 / a );
  }
  dot( v ) {
    this._dimension_check( v );
    return this.x.map( ( x, i ) => x * v.x[ i ] )
      .reduce( ( x, sum ) => x + sum, 0 );
  }
  static dot( u, v ) {
    return u.dot( v );
  }
  cross( v ) {
    this._dimension_check( v );
    let u = this;
    if( this.x.length !== 3 ) {
      if( this.x.length === 2 ) {
        // WARNNING: only return the vertical component
        return u.x[ 0 ] * v.x[ 1 ] - u.x[ 1 ] * v.x[ 0 ];
      }
      throw new Error(
        "Cross product is only defined in 3 dimensional space!" );
    }
    return new Vector( [
      u.x[ 1 ] * v.x[ 2 ] - u.x[ 2 ] * v.x[ 1 ],
      u.x[ 2 ] * v.x[ 0 ] - u.x[ 0 ] * v.x[ 2 ],
      u.x[ 0 ] * v.x[ 1 ] - u.x[ 1 ] * v.x[ 0 ],
    ] );
  }
  static cross( u, v ) {
    return u.cross( v );
  }
  norm() {
    return Math.sqrt( this.dot( this ) );
  }
  normlize() {
    return this.divide( this.norm() );
  }
  dimension() {
    return this.x.length;
  }
  equals( v ) {
    return this.x.minus( v )
      .components()
      .every( x => x === 0 );
  }
  isZero( v ) {
    return this.x.every( x => x === 0 );
  }
  toString() {
    return "(" + this.x.join( "," ) + ")";
  }
  static isVector( v ) {
    return v instanceof Vector;
  }
  // private
  _dimension_check( v ) {
    if( v.x.length !== this.x.length ) {
      throw new Error( "Two vectors must have same dimension!" );
    }
    return v;
  }
  _type_guard() {
    if( !this.x.every( x => typeof x === "number" ) ) {
      throw new TypeError( "Vector must be an array of numbers!" );
    }
  }
};

class Geometry {
  constructor() {

  }
  static distance( c1, c2 ) {
    if( Point.isPoint( c1 ) && Point.isPoint( c2 ) ) {
      return Geometry._distance_point_point( c1, c2 );
    }
  }
  // private
  static _distance_point_point( p1, p2 ) {
    return( new Vector( p1, p2 ) )
      .norm();
  }
}

class Geometry2D extends Geometry {
  constructor() {

  }
  static distance( c1, c2 ) {
    if( c1.dimension() === 2 && c2.dimension() === 2 ) {
      if( Point.isPoint( c1 ) && Line.isLine( c2 ) ) {
        return Geometry2D._distance_point_line( c1, c2 );
      } else if( Point.isPoint( c2 ) && Line.isLine( c1 ) ) {
        return Geometry2D._distance_point_line( c2, c1 );
      } else if( Point.isPoint( c1 ) && Circle.isCircle( c2 ) ) {
        return Geometry2D._distance_point_circle( c1, c2 );
      } else if( Point.isPoint( c2 ) && Circle.isCircle( c1 ) ) {
        return Geometry2D._distance_point_circle( c2, c1 );
      } else if( Circle.isCircle( c1 ) && Circle.isCircle( c2 ) ) {
        return Geometry2D._distance_circle_circle( c1, c2 );
      } else {
        return super.distance( c1, c2 );
      }
    }
  }
  static intersect( c1, c2 ) {
    if( Segment.isSegment( c1 ) && Circle.isCircle( c2 ) ) {
      return Geometry2D._intersect_circle_segment( c2, c1 );
    } else if( Circle.isCircle( c1 ) && Segment.isSegment( c2 ) ) {
      return Geometry2D._intersect_circle_segment( c1, c2 );
    } else if( Segment.isSegment( c1 ) && Segment.isSegment( c2 ) ) {
      return Geometry2D._intersect_segment_segment( c1, c2 );
    }
  }
  static intersection( c1, c2 ) {
    if( Line.isLine( c1 ) && Line.isLine( c2 ) ) {
      let v1 = new Vector( c1.general_form() );
      let v2 = new Vector( c2.general_form() );
      let d = v1.cross( v2 )
        .components();
      if( d[ 2 ] !== 0 ) return new Point( d[ 0 ] / d[ 2 ], d[ 1 ] / d[ 2 ] );
      return null;
    }
  }
  static inConvexPolygon( point, polygon ) {
    let points = polygon.points;
    let n = points.length;
    let v = new Vector( points[ 0 ], point );
    let v1 = new Vector( points[ 0 ], points[ 1 ] );
    let vn_1 = new Vector( points[ 0 ], points[ n - 1 ] );
    if( v.cross( v1 ) > 0 || v.cross( vn_1 ) < 0 ) {
      return false; //不在夹角p[1]p[0]p[n-1]内
    }
    // 二分查找落在哪两条对角线夹角内: 角p[low-1], p[0], p[low]
    let low = 2,
      high = n - 1;
    while( low <= high ) {
      let mid = ( low + high ) / 2;
      let vm = new Vector( points[ 0 ], vector( points[ mid ] ) );
      if( v.cross( vm ) > 0 ) high = mid - 1;
      else low = mid + 1;
    }
    /*
    let low = 2;
    for (; low < n; low++) {
        var vm = points[0].vector(points[low]);
        if (v.cross(vm).isPositive()) {
            break;
        }
    }*/
    // p0,p[low-1],p[low] 构成的夹角内
    if( ( new Vector( points[ low - 1 ], points[ low ] ) )
      .cross( new Vector( points[ low - 1 ], point ) ) > 0 ) {
      return true; // 三角形内
    }
    return false;
  }
  static tangent_segments( point, circle ) {
    if( Point.isPoint( circle ) && Circle.isCircle( point ) ) {
      [ point, circle ] = [ circle, point ];
    }
    //注意避免三角函数运算
    let segs = new Array();
    let d = Geometry2D.distance( point, circle.center );
    let r = circle.radius;
    if( d <= r ) return segs; // 点在圆内
    // 切线长 L, d为到圆心距离，r为圆半径，弦长一半l，点到弦的距离为d1。弦心距为d2。d=d1+d2
    // d/L=L/d1=r/l, d/r=L/l=r/d2 （相似三角形）
    // L=sqrt(d^2-r^2) (勾股定理)
    // d1=L^2/d，d2=r^2/d, l=rL/d
    let L = Math.sqrt( d * d - r * r );
    let d1 = L * L / d,
      l = r * L / d;
    let v = new Vector( point, circle.center )
      .normlize();
    let a = v.mult( d1 ),
      b = Geometry2D.orthogonal( v )
      .mult( l );
    let p1 = point.move( a.add( b ) ),
      p2 = point.move( a.minus( b ) );
    segs.push( new Segment( point, p1 ) );
    segs.push( new Segment( point, p2 ) );
    return segs;
  }
  static common_tangent_segments( c1, c2 ) {
    if( c1.radius === 0 && c2.radius === 0 ) { // 退化为两点
      return [ new Segment( c1.center, c2.center ) ];
    } else if( c1.radius === 0 ) { // 点与圆，求两条切线即可
      return Geometry2D.tangent_segments( c1.center, c2 );
    } else if( c2.radius === 0 ) {
      return Geometry2D.tangent_segments( c2.center, c1 )
        .map( s => s.reverse() );
    }
    // 求两圆的公切线
    let segs = new Array(); // there will be 0,2,4 segments
    let swap = false;
    let [ A, B ] = [ c1, c2 ];
    if( A.radius < B.radius ) {
      swap = true;
      [ A, B ] = [ B, A ];
    }
    let AB = new Vector( A.center, B.center );
    let [ R, r ] = [ A.radius, B.radius ];
    let d = Geometry2D.distance( A.center, B.center ); // 圆心距
    let [ rdiff, rsum ] = [ R - r, R + r ]; // 半径差，半径和
    if( d === 0 || d < rdiff ) return segs; // 1、内含 2: 重合
    if( d === rdiff ) { //3.内切
      // A.point(base); 存在公切线，但不存在线段
      return segs;
    }
    //4.相交（外切、外离的外公切线也在此求出）
    // 过圆B的圆心向圆A的半径做垂线，得到直角三角形。切线长L=sqrt(d^2-(R-r)^2)
    // 过圆A的切点向圆心连线做垂线，即半条弦，得到直角三角形。相似。
    // l/R=L/d => l=LR/d，即半条弦长。(R-r)/d=d1/R => d1=R(R-r)/d，即弦心距。
    // 小圆切点向圆连心线做垂线，由相似比 r/R。
    let L = Math.sqrt( d * d - rdiff * rdiff );
    let l = L * R / d,
      d1 = R * rdiff / d;
    let v = AB.normlize();
    let a = v.mult( d1 ),
      b = Geometry2D.orthogonal( v )
      .mult( l );
    let [ v1, v2 ] = [ a.add( b ), a.minus( b ) ];
    let r_R = r / R;
    let [ p1, p2 ] = [ A.center.move( v1 ), B.center.move( v1.mult( r_R ) ) ];
    let [ p3, p4 ] = [ A.center.move( v2 ), B.center.move( v2.mult( r_R ) ) ];
    if( !swap ) {
      segs.push( new Segment( p1, p2 ) );
      segs.push( new Segment( p3, p4 ) );
    } else {
      segs.push( new Segment( p2, p1 ) );
      segs.push( new Segment( p4, p3 ) );
    }
    if( d === rsum ) { //5.外切
      // A.point(base); 存在公切线，但不存在线段
      return segs;
    }
    if( d > rsum ) { //6.外离
      // 反向延长小圆B的半径，过大圆A的圆心垂线。
      // 切线长 L=sqrt(d^2-(R+r)^2)
      // l/R=L/d => l=RL/d, d1/R=(R+r)/d =>d1=R(R+r)/d
      let L = Math.sqrt( d * d - rsum * rsum );
      let l = L * R / d,
        d1 = R * rsum / d;
      let v = AB.normlize();
      let [ a, b ] = [ v.mult( d1 ), Geometry2D.orthogonal( v )
        .mult( l )
      ];
      let [ v1, v2 ] = [ a.add( b ), a.minus( b ) ];
      let r_R = -r / R; //大圆和小圆的径向是反的
      let [ p1, p2 ] = [ A.center.move( v1 ), B.center.move( v1.mult( r_R ) ) ];
      let [ p3, p4 ] = [ A.center.move( v2 ), B.center.move( v2.mult( r_R ) ) ];
      if( !swap ) {
        segs.push( new Segment( p1, p2 ) );
        segs.push( new Segment( p3, p4 ) );
      } else {
        segs.push( new Segment( p2, p1 ) );
        segs.push( new Segment( p4, p3 ) );
      }
    }
    return segs;
  }
  static GrahamScan( points ) {
    if( points.length < 3 ) return [ ...points.keys() ];
    points = points.map( ( p, i ) => {
      p.idx = i;
      return p;
    } );
    //选定基点
    let base = points.reduce( ( p0, p ) => {
      let xp0 = p0.coordinates();
      let xp = p.coordinates();
      for( let [ i, x ] of xp.entries() ) {
        if( x < xp0[ i ] ) return p;
        if( x > xp0[ i ] ) return p0;
      }
      return p0;
    }, points[ 0 ] );
    let ps = points.filter( ( p ) => p !== base )
      .sort( function ( p1, p2 ) { // 极角排序
        let v1 = new Vector( base, p1 );
        let v2 = new Vector( base, p2 );
        let d = v2.cross( v1 ); // d>0 v1 在 v2 逆时针方向
        if( d === 0 ) { //共线时，按距离，即模比较
          return v1.norm() - v2.norm();
        }
        return d;
      } )
      .filter( function ( p, i, ary ) {
        return i === 0 || !p.equals( ary[ i - 1 ] );
      } );
    let vector = Array( ps.length - 1 )
      .fill()
      .map( ( x, i ) => [ i, i + 1 ] );
    //依次删除不在凸包上的向量
    for( let i = 1, N = vector.length; i < N; i++ ) {
      //回溯删除旋转方向相反的向量，使用外积判断旋转方向
      for( let j = i - 1; j >= 0; j-- ) {
        if( vector[ j ] == null ) continue;
        let vi = new Vector( vector[ i ].map( ( x ) => ps[ x ] ) ); //i0-->i1
        let vj = new Vector( vector[ j ].map( ( x ) => ps[ x ] ) ); //j0-->j1
        var d = vi.cross( vj );
        /*
                         i1
                        / vi
                vj    /
           j0-----j1(i0)
        */
        if( d < 0 ) break; //无逆向旋转
        if( d === 0 ) { // 共线
          if( vi.dot( vj ) > 0 ) break; // 同向
        }
        //删除前一个向量后，需更新当前向量，与前面的向量首尾相连
        //向量三角形计算公式
        /*
            j0-----j1(i0)
                       \
                        \
                        i1
        */
        if( vector[ j ][ 1 ] !== vector[ i ][ 0 ] ) console.log( "ERROR!" );
        vector[ i ][ 0 ] = vector[ j ][ 0 ];
        vector[ j ] = null;
      }
    }
    vector = vector.filter( ( v ) => v !== null );
    //计算凸包点
    return [ base.idx, ps[ vector[ 0 ][ 0 ] ].idx ].concat( vector.map( ( v ) =>
      ps[ v[ 1 ] ].idx ) );
  }
  static circle_convex( circles ) {
    const eps = 1e-12;
    if( circles.length === 1 ) return 0;
    // 变换到三维空间
    let points = circles.map( c => {
      let [ x, y ] = c.center.coordinates();
      return new Point( x, y, c.radius );
    } );
    let convex = Geometry3D.quickhull( points );
    if( Array.isArray( convex ) ) { //多面体
      console.log( "polyhedron" );
      //求重心, 任何一个位于凸包内的点都可以，比较保险的就是重心了
      let O = new Point( convex.map( ( f ) =>
          f.reduce( ( o, p ) => points[ p ].coordinates()
            .map( ( x, i ) => x / 3.0 + o[ i ] ), [ 0, 0, 0 ] ) // 每个面的重心
        )
        .reduce( ( o, p, j, c ) => p.map( ( x, i ) => x * 1.0 / c.length +
          o[ i ] ), [ 0, 0, 0 ] ) ); // 整个凸包的重心
      let v = new Vector( O, Point.ORIGIN( 3 ) );
      let dual = convex.map( ( f ) => {
        // 原点指向面的单位法向量
        let n = ( new Vector( points[ f[ 0 ] ], points[ f[ 1 ] ] ) )
          .cross( new Vector( points[ f[ 0 ] ], points[ f[ 2 ] ] ) )
          .normlize();
        let d = n.dot( new Vector( O, points[ f[ 1 ] ] ) );
        if( d < 0 ) {
          n = n.negated();
          d = -d; // 原点到面的距离
        }
        if( Math.abs( d ) < eps ) console.log( "possible error" );
        return new Point( n.mult( 1.0 / d )
          .components() );
      } );
      // 锥面方程 x^2+y^2=z^2， z>=0
      // P.(0,0,1)=|P|sqrt(2) ==> z^2=sqrt(x^2+y^2+z^2)/sqrt(2)
      // 直线： P=P0+vt ===> x=x0+vxt, y=y0+vyt, z=z0+vzt
      // (x0+vxt)^2+(y0+vyt)^2=(z0+vzt)^2
      // (vx^2+vy^2-vz^2)t^2 + 2(x0vx+y0vy-z0vz)t+(x0^2+y0^2-z0^2)=0
      // delta=B^2-4AC
      // t1, t2 ==> z>=0
      let circle_pair = [];
      // 两个面有且只有一条交线，因此交线等价于面的对
      let edges = Array.from( dual.keys() )
        .map( ( f, i, a ) => a.map( ( f1 ) => [ f, f1 ] ) )
        .reduce( ( a, b ) => a.concat( b ) )
        .filter( ( p ) => p[ 0 ] < p[ 1 ] )
        .map( ( e ) => e.concat( [ convex[ e[ 0 ] ].filter( ( p ) => convex[
          e[ 1 ] ].includes( p ) ) ] ) ); // [面1, 面2, [交线的点（两个）]]
      for( let e of edges ) {
        if( e[ 2 ].length === 2 ) {
          let [ x0, y0, z0 ] = dual[ e[ 0 ] ].coordinates();
          let [ vx, vy, vz ] = ( new Vector( dual[ e[ 0 ] ], dual[ e[ 1 ] ] ) )
          .components();
          if( z0 < 0 && z0 + vz < 0 ) continue; // [z0,z0+vz]
          let A = vx * vx + vy * vy - vz * vz;
          let B = 2 * ( x0 * vx + y0 * vy - z0 * vz );
          let C = x0 * x0 + y0 * y0 - z0 * z0;
          if( A !== 0 ) { // 二次方程
            let d = B * B - 4 * A * C;
            if( d < 0 ) continue; // 无解
            let t0 = [ ( -B + Math.sqrt( d ) ) / ( 2 * A ), ( -B - Math.sqrt(
              d ) ) / ( 2 * A ) ];
            if( t0.some( ( t ) => z0 + vz * t > 0 && 0 <= t && t <= 1 ) ) {
              circle_pair.push( e[ 2 ] );
            }
          } else if( B !== 0 ) { //一次方程
            let t0 = -C / B;
            if( z0 + vz * t0 > 0 && 0 <= t0 && t0 <= 1 ) {
              circle_pair.push( e[ 2 ] );
            }
          } else if( C !== 0 ) { // 无解
            continue;
          } else { // 恒等式，无数解
            circle_pair.push( e[ 2 ] );
          }
        }
      }
      return circle_pair; // 并不能表示成一个环
    } else if( convex.error === "coplanar" ) {
      console.log( "coplanar" );
      let normal = convex.normal;
      let convex2d = Geometry3D.convex2d( points, normal ); // 必须内有内含和内切的圆
      return convex2d.map( ( p, i, a ) => [ p, a[ ( i + 1 ) % a.length ] ] );
    } else if( convex.error === "collinear" ) {
      console.log( "collinear" );
      let u = convex.direction.normlize();
      let cos_angle = Math.abs( u.dot( new Vector( 0, 0, 1 ) ) );
      if( cos_angle >= 0.5 * Math.sqrt( 2 ) ) { //只有最大一个圆
        let [ max_r, i ] = circles.reduce( ( res, cur, idx ) => ( res[ 0 ] <
          cur.radius ? [ cur, idx ] : res ), [ circles[ 0 ].radius, 0 ] );
        return i;
      } else { // 全部圆，或者最大和最小两个圆即可
        let [ max_r, max_ri, min_r, min_ri ] = circles.map( c => c.radius )
          .reduce( ( res, cur, idx ) => {
            if( res[ 0 ] < cur ) res[ 0 ] = cur, res[ 1 ] = idx;
            if( res[ 0 ] > cur ) res[ 2 ] = cur, res[ 3 ] = idx;
            return res;
          }, [ circles[ 0 ].radius, 0, circles[ 0 ].radius, 0 ] );
        if( max_r !== min_r ) return [
          [ max_ri, min_ri ]
        ];
        // 最远的两个
        let signed_distance = circles.map( ( c, i, a ) => u.dot( new Vector(
          a[ 0 ].center, c.center ) ) );
        let [ max_d, max_di, min_d, min_di ] = signed_distance.reduce( ( res,
          cur, idx ) => {
          if( res[ 0 ] < cur ) res[ 0 ] = cur, res[ 1 ] = idx;
          if( res[ 0 ] > cur ) res[ 2 ] = cur, res[ 3 ] = idx;
          return res;
        }, [ signed_distance[ 0 ], 0, signed_distance[ 0 ], 0 ] );
        return [
          [ max_di, min_di ]
        ];
      }
    } else if( convex.error === "concurrent" ) {
      return 0; // 第一个圆，其实任一个圆都可以
    }
  }
  // delaunay 分割
  static delaunay( points ) {
    let [ x0, y0 ] = points[ 0 ].coordinates();
    let [ max_x, min_x, max_y, min_y ] = points.reduce( ( res, cur ) => {
      let [ x, y ] = cur.coordinates();
      if( res[ 0 ] < x ) res[ 0 ] = x;
      if( res[ 1 ] > x ) res[ 1 ] = x;
      if( res[ 2 ] < y ) res[ 2 ] = y;
      if( res[ 3 ] > y ) res[ 3 ] = y;
      return res;
    }, [ x0, x0, y0, y0 ] );
    //矩形框的外接圆
    let [ w, h ] = [ max_x - min_x, max_y - min_y ];
    let r = 0.5 * Math.sqrt( w * w + h * h );
    let [ xc, yc ] = [ ( max_x + min_x ) / 2, ( max_y + min_y ) / 2 ];
    //let circumscribed_circle=new Circle(new Point(xc,yc),r);
    let a = new Point( xc - Math.sqrt( 3.0 ) * r, yc - r );
    let b = new Point( xc + Math.sqrt( 3.0 ) * r, yc - r );
    let c = new Point( xc, yc + 2 * r );
    a.idx = -3, b.idx = -2, c.idx = -1;
    //points=points.concat([a,b,c]);
    points = points.map( ( p, i ) => {
      let p1 = new Point( p.coordinates() );
      p1.idx = i;
      return p1;
    } );
    let baseTriangle = new Triangle( a, b, c );
    //外部三角形を格納
    let triangles = [ baseTriangle ];
    points.forEach( function ( point ) {
      let edges = [];
      let unique_edges = [];
      //ポイントを外接円に含む三角形を抽出、辺に分解
      triangles.forEach( function ( targetTriangle, i ) {
        if( targetTriangle.circumscribed_circle()
          .includes( point ) ) {
          edges = edges.concat( targetTriangle.edges() );
          delete triangles[ i ];
        }
      } );
      //分解した辺リストから重複する辺を削除
      edges.forEach( function ( edge0, i ) {
        let unique_flag = true;
        edges.forEach( function ( edge1, j ) {
          //重複する辺がある場合
          if( i != j && edge0.equals( edge1 ) ) {
            unique_flag = false;
          }
        } );
        if( unique_flag ) {
          unique_edges.push( edge0 )
        }
      } );
      edges = unique_edges;
      //重複しない辺リストから三角形を生成
      edges.forEach( function ( edge, i ) {
        triangles.push( new Triangle( edge.points[ 0 ], edge.points[
          1 ], point ) );
      } );
    } );
    return [ triangles.filter( t => t instanceof Triangle )
      .map( t => t.points.map( p => p.idx ) ), baseTriangle
    ];
  }
  // private
  static _distance_point_line( p, l ) {
    let v = new Vector( p, l.p );
    let u = Geometry2D.orthogonal( l.direction() )
      .normlize();
    return Math.abs( v.dot( u ) ); // project v on u;
  }
  static _distance_point_circle( p, c ) {
    return Geometry2D.distance( p, c.center ) - c.radius;
  }
  static _distance_circle_circle( c1, c2 ) {
    return Geometry2D.distance( c1.center, c2.center ) - ( c1.radius + c2.radius );
  }
  static orthogonal( v ) {
    let [ x, y ] = v.components();
    return new Vector( -y, x );
  }
  static _intersect_circle_segment( circle, segment ) {
    // 线段的端点落在圆内或圆上
    if( segment.points.some( ( p ) => Geometry2D.distance( p, circle.center ) <=
        circle.radius ) ) {
      return true;
    }
    // 垂线的方向向量=>该线段法向量
    let L = segment.toLine();
    let normal = Geometry2D.orthogonal( L.direction() );
    let sign = segment.points.map( ( p ) => Math.sign( ( new Vector( circle.center,
        p ) )
      .cross( normal ) ) );
    return sign[ 0 ] !== sign[ 1 ] && Geometry2D.distance( circle.center, L ) <=
      circle.radius; //p1p2在垂线两侧，距离小于半径
  }
  static _intersect_segment_segment( segment1, segment2 ) {
    const [ p1, p2 ] = segment1.points;
    const [ p3, p4 ] = segment2.points;
    // 快速排斥实验，线段所在的矩形区域不相交===>线段不相交
    if( Math.max( p1.x[ 0 ], p2.x[ 0 ] ) < Math.min( p3.x[ 0 ], p4.x[ 0 ] ) || // p1p2 位于 p3p4之左
      Math.max( p1.x[ 1 ], p2.x[ 1 ] ) < Math.min( p3.x[ 1 ], p4.x[ 1 ] ) || // p1p2 位于 p3p4之下
      Math.max( p3.x[ 0 ], p4.x[ 0 ] ) < Math.min( p1.x[ 0 ], p2.x[ 0 ] ) || // p3p4 位于 p1p2之左
      Math.max( p3.x[ 1 ], p4.x[ 1 ] ) < Math.min( p1.x[ 1 ], p2.x[ 1 ] ) ) { // p3p4 位于 p1p2之左
      return false;
    }
    // 快速跨立实验
    const v12 = new Vector( p1, p2 ),
      v23 = new Vector( p2, p3 ),
      v24 = new Vector( p2, p4 );
    const s1 = Math.sign( v12.cross( v23 ) ); // p3在p1p2左
    const s2 = Math.sign( v12.cross( v24 ) ); // p4在p1p2左
    const v34 = new Vector( p3, p4 ),
      v41 = new Vector( p4, p1 ),
      v42 = new Vector( p4, p2 );
    const s3 = Math.sign( v34.cross( v41 ) ); // p1在p3p4左
    const s4 = Math.sign( v34.cross( v42 ) ); // p2在p3p4左
    // p1p2同侧，或p3p4同侧
    if( s1 === s2 && s3 === s4 ) return false;
    //相交
    return true;
  }
  static _intersect_segment_convex_polygon( segment, polygon ) {
    if( segment.points.some( ( p ) => Geometry2D.inConvexPolygon( p, polygon ) ) ) {
      return false;
    }
    var v = segment.toLine()
      .direction(); // 方向向量
    // 用两条平行于该直线的直线把多边形夹住，可以找到两个顶点。
    // TODO::二分O(log n)查找峰值
    let peaks = [];
    let points = polygon.points;
    let n = points.length;
    let det_sign = polygon.points.map( ( p, i, a ) => Math.sign( v.cross( new Vector(
      p, a[ ( i + 1 ) % n ] ) ) ) );
    // 相邻
    det_sign.push( det_sign[ 0 ] ); // loop
    for( let i = 0; i < n; i++ ) {
      if( det_sign[ i ] !== det_sign[ i + 1 ] || det_sign[ i ] !== 0 &&
        det_sign[ i + 1 ] == 0 ) {
        peaks.push( ( i + 1 ) % n );
      }
    }
    if( peaks.length != 2 ) {
      throw new RangeError( "Impossile" );
    }
    let diameter = new Segment( points[ peaks[ 0 ] ], points[ peaks[ 1 ] ] );
    return Geometry2D.intersect( segment, diameter );
  }
};

class Line {
  // 1. 2n numbers (x1,y1,z1,..., x2,y2,z2,...)
  // 2. 2 points (p1,p2)
  // 3. a point and a direction vector (p, v)
  constructor( p, v ) {
    let illegal = false;
    if( Point.isPoint( p ) ) {
      if( Vector.isVector( v ) ) {
        if( p.dimension() !== v.dimension() ) {
          throw new Error(
            "The point and direction vector must have same dimension!" );
        }
        this.p = p;
        this.v = v;
      } else if( Point.isPoint( v ) ) {
        if( p.dimension() !== v.dimension() ) {
          throw new Error( "Two points must have same dimension!" );
        }
        this.p = p;
        this.v = new Vector( p, v );
      } else {
        illegal = true;
      }
    } else {
      let p = Array.from( arguments );
      if( !( p.length % 2 ) && p.every( ( x ) => typeof x === "number" ) ) {
        this.p = new Point( p.slice( 0, p.length / 2 ) );
        this.v = new Vector( this.p, new Point( p.slice( p.length / 2, p.length ) ) );
      } else {
        illegal = true;
      }
    }
    if( illegal ) {
      throw new Error( "Illegal arguments!" );
    }
  }
  direction() {
    return this.v;
  }
  general_form() {
    let [ x1, y1 ] = this.p.coordinates();
    let [ x2, y2 ] = this.p.move( this.v )
      .coordinates();
    return [ y2 - y1, x1 - x2, x2 * y1 - x1 * y2 ];
  }
  dimension() {
    return this.p.dimension();
  }
  static isLine( l ) {
    return l instanceof Line;
  }
};

class Segment {
  // 2 points
  // 2n numbers
  constructor( p1, p2 ) {
    if( Point.isPoint( p1 ) && Point.isPoint( p2 ) && p1.dimension() === p2.dimension() ) {
      this.points = [ p1, p2 ];
    } else {
      let p = Array.from( arguments );
      if( !( p.length % 2 ) && p.every( ( x ) => typeof x === "number" ) ) {
        this.points = [ new Point( p.slice( 0, p.length / 2 ) ), new Point( p
          .slice( p.length / 2, p.length ) ) ];
      } else if( p.length === 2 && p.every( ( x ) => Point.isPoint( x ) ) &&
        p[ 0 ].dimension() === p[ 1 ].dimension() ) {
        this.points = p;
      }
    }
  }
  length() {
    return Geometry.distance( this.points[ 0 ], this.points[ 1 ] );
  }
  dimension() {
    return this.points[ 0 ].dimension();
  }
  toLine() {
    return new Line( this.points[ 0 ], this.points[ 1 ] );
  }
  equals( l ) {
    return l.points[ 0 ].equals( this.points[ 0 ] ) && l.points[ 1 ].equals(
        this.points[ 1 ] ) ||
      l.points[ 0 ].equals( this.points[ 1 ] ) && l.points[ 1 ].equals( this.points[
        0 ] );
  }
  reverse() {
    return new Segment( this.points.reverse() );
  }
  static isSegment( s ) {
    return s instanceof Segment;
  }
};

class Circle { // only in 2d space
  // 1. center and radius
  // 2. coordinates of center and radius
  constructor( c, r ) {
    if( Point.isPoint( c ) && c.dimension() === 2 && typeof r === "number" ) {
      this.center = c;
      this.radius = r;
    } else {
      let p = Array.from( arguments );
      if( p.length === 3 && p.every( ( x ) => typeof x === "number" ) ) {
        this.center = new Point( p[ 0 ], p[ 1 ] );
        this.radius = p[ 2 ];
      }
    }
  }
  dimension() {
    return 2; // must be in 2d space
  }
  includes( p ) {
    return Geometry2D.distance( p, this.center ) < this.radius;
  }
  static isCircle( c ) {
    return c instanceof Circle;
  }
};

class Quaternion {
  constructor( r, v ) {
    if( Vector.isVector( v ) && v.dimension() === 3 && typeof r === "number" ) {
      this.real = a;
      this.imag = v;
    } else {
      let x = Array.from( arguments );
      if( x.length === 4 ) {
        this.real = x[ 0 ];
        this.imag = new Vector( x.slice( 1 ) );
      }
    }
  }
  add( q ) {
    return new Quaternion( this.real + q.real, this.imag.add( q.imag ) );
  }
  conjugate( q ) {
    return new Quaternion( this.real, this.imag.negated() );
  }
  mult( q ) {
    if( q instanceof Quaternion ) {
      return new Quaternion( this.real * q.real - this.imag.dot( q.imag ),
        q.imag.mult( this.real )
        .add( this.imag.mult( q.real ) )
        .add( this.imag.cross( q.imag ) ) );
    } else if( typeof q === "number" ) {
      return new Quaternion( this.real * q, this.imag.mult( q ) );
    }
  }
  norm() {
    return Math.sqrt( this.real * this.real + this.imag.dot( this.imag ) );
  }
  reciprocal() {
    let norm = this.norm();
    return this.conjugate()
      .mult( 1.0 / ( norm * norm ) );
  }
}

class Geometry3D extends Geometry {
  constructor() {

  }
  static rotate( point, q ) {

  }
  static distance( c1, c2 ) {
    if( Point.isPoint( c1 ) && Point.isPoint( c2 ) ) {
      return super.distance( c1, c2 );
    } else if( Point.isPoint( c1 ) && Line.isLine( c2 ) ) {
      return Geometry3D._distance_point_line( c1, c2 );
    } else if( Point.isPoint( c2 ) && Line.isLine( c1 ) ) {
      return Geometry3D._distance_point_line( c2, c1 );
    }
  }
  static _distance_point_line( point, line ) {
    let s = new Vector( point, line.p );
    return s.cross( line.v )
      .norm() / line.v.norm();
  }
  static convex2d( points, normal ) {
    let base = points[ 0 ];
    let points_project = points.map( ( p ) => {
      let r = new Vector( base, p );
      let p1 = r.minus( normal.mult( r.dot( normal ) ) );
      return new Point( p1.components()
        .slice( 0, 2 ) );
    } );
    return Geometry2D.GrahamScan( points_project );
  }
  /*
  {error:"concurrent", point:xxx} 共点则返回点
  {error:"collinear", direction:xxx} 共线则返回线的方向向量
  {error:"coplanar",normal"xxx} 共面则返回法向量
  [face1,face2,...] 返回多面体的面
  */
  static quickhull( points ) {
    const eps = 1e-12;
    points = points.map( ( p, i ) => { // 保存下标后...
        p.idx = i;
        return p;
      } )
      .filter( ( p, i, ps ) => ps.findIndex( p1 => p.equals( p1 ) === i ) ); // ...去重
    if( points.length === 1 ) return {
      error: "concurrent",
      point: points[ 0 ]
    };
    if( points.length === 2 ) {
      return {
        error: "collinear",
        direction: new Vector( points[ 0 ], points[ 1 ] )
      };
    }
    if( points.length === 3 ) {
      let [ p1, p2, p3 ] = points;
      let v12 = new Vector( p1, p2 ),
        v13 = new Vector( p1, p3 );
      let normal = v12.cross( v13 );
      if( normal.isZero() ) {
        return {
          error: "collinear",
          direction: new Vector( p1, p2 )
        };
      }
      return {
        error: "coplanar",
        normal: normal
      };
    }
    let faces = [];
    let faceStack = [];
    const origin = Point.ORIGIN( 3 );
    const process = function ( points ) {
      // Iterate through all the faces and remove
      while( faceStack.length > 0 ) { //待定面集，队列
        cull( faceStack.shift(), points ); //取出第一个面
      }
    }
    const norm = function ( a, b, c ) { // 单位法向量
      return( new Vector( a, c ) )
        .cross( new Vector( a, b ) )
        .normlize();
    }
    //                n    p[1]
    //                 |  /    \
    //                 |/       \
    //               / |         \
    //            p[0]-----------p[2]
    //
    const getNormal = function ( face, points ) {
      if( face.normal !== undefined ) return face.normal;
      let p = face.map( ( idx ) => points[ idx ] );
      return face.normal = norm( p[ 0 ], p[ 1 ], p[ 2 ] );
    }
    const assignPoints = function ( face, pointset, points ) {
      let p0 = points[ face[ 0 ] ];
      let norm = getNormal( face, points );
      let dis = pointset.map( ( p ) => [ p, norm.dot( new Vector( p0,
        points[ p ] ) ) ] );
      face.visiblePoints = dis.filter( ( p ) => p[ 1 ] > eps ) // 法向量指向凸包外
        .sort( function ( a, b ) {
          return a[ 1 ] - b[ 1 ];
        } )
        .map( ( p ) => p[ 0 ] );
    };

    let vIndex = [];
    // 待定面集的每个面
    const cull = function ( face, points ) {
      // 外部点集中最远的点
      let apex = face.visiblePoints.pop();
      // 两个三角面等价 ==> f1的三个点都在f2上 ==> f1.every((p)=>f2.includes(p))
      // 计算可见面 =========> 法向量与视线夹角小于90度
      let visibleFaces = [ face ].concat( faces.filter( ( f ) => !f.every(
          ( p ) => face.includes( p ) ) && // 不是同一个面
        getNormal( f, points )
        .dot( new Vector( points[ f[ 0 ] ], points[ apex ] ) ) > eps ) ); //可见
      // 计算临界边
      let allPoints = [];
      let perimeter = []; // 临界边
      for( let currentFace of visibleFaces ) {
        //可见面必然不在凸包上
        faceStack = faceStack.filter( ( f ) => !f.every( ( p ) =>
          currentFace.includes( p ) ) ); // 相同的面就不需要再拓展了
        if( currentFace.visiblePoints !== undefined && currentFace.visiblePoints
          .length > 0 )
          allPoints = allPoints.concat( currentFace.visiblePoints ); //汇总外部点集
        faces = faces.filter( ( f ) => !f.every( ( p ) => currentFace.includes(
          p ) ) ); // 从凸包移除
        // 边
        let edges = currentFace.map( ( p, i, f ) => [ p, f[ ( i + 1 ) % 3 ] ] ); // 3 edges of triangle
        for( let edge of edges ) {
          // 如果当前边不是被两个可见面共享的, 或者当前只有一个面
          if( !visibleFaces.some( ( f ) => !f.every( ( p ) => currentFace.includes(
              p ) ) && edge.every( ( p ) => f.includes( p ) ) ) ) {
            perimeter.push( edge );
          }
        }
      }
      allPoints = allPoints.filter( ( p, i, a ) => a.indexOf( p ) === i ); // unique
      // 构造新的面
      for( let edge of perimeter ) {
        let f = [ edge[ 1 ], apex, edge[ 0 ] ]; //顺序很重要, 保证从凸包外面看, 顺时针.这样法向量才会向外
        assignPoints( f, allPoints, points );
        faces.push( f );
        if( f.visiblePoints !== undefined && f.visiblePoints.length > 0 )
          faceStack.push( f );
      }
    }

    const initTetrahedron = function () {
      // 左右前后上下，找到六个最远的位置
      let extremes = [];
      let px = points.reduce( ( xs, p ) => p.coordinates()
        .map( ( x, i ) => xs[ i ].concat( [ x ] ) ), [
          [],
          [],
          []
        ] );
      for( let i = 0; i < 3; i++ ) {
        extremes[ i + i ] = px[ i ].indexOf( Math.min( ...px[ i ] ) );
        extremes[ i + i + 1 ] = px[ i ].indexOf( Math.max( ...px[ i ] ) );
      }
      // 最远的一对点, 这并非所有点对的最远点对, 尽量远即可
      let pair = Array.from( extremes.keys() )
        .map( ( x, i, a ) => a.map( ( y ) => [ x, y ] ) )
        .reduce( ( a, b ) => a.concat( b ) ) //下标的笛卡尔集
        .filter( ( p ) => p[ 0 ] < p[ 1 ] )
        .map( ( p ) => p.map( ( i ) => extremes[ i ] ) ); // 仅保留不相等的组合
      let dis = pair.map( p => Geometry3D.distance( points[ p[ 0 ] ],
        points[ p[ 1 ] ] ) );
      let [ max_dis, v ] = dis.reduce( ( res, cur, idx ) => ( res[ 0 ] <
        cur ? [ cur, pair[ idx ] ] : res ), [ 0, [ 0, 0 ] ] );
      if( max_dis === 0 ) return {
        error: "concurrent",
        point: points[ 0 ]
      }; //仅当所有的点相同时,最远点对距离为0
      // 到该线段最远的点，于是v[0],v[1],v[2]构成一个面
      let L = new Line( points[ v[ 0 ] ], points[ v[ 1 ] ] );
      let dis_line = points.map( p => Geometry3D.distance( L, p ) );
      let [ max_dis_line, v2 ] = dis_line.reduce( ( res, cur, idx ) => (
        res[ 0 ] < cur ? [ cur, idx ] : res ), [ 0, v[ 0 ] ] );
      if( max_dis_line === 0 ) return {
        error: "collinear",
        direction: L.direction()
      }; //所有的点共线
      v[ 2 ] = v2;
      // 到该面最远的点
      let N = norm( points[ v[ 0 ] ], points[ v[ 1 ] ], points[ v[ 2 ] ] ); // 法向量
      let D = N.dot( new Vector( origin, points[ v[ 0 ] ] ) );
      // 平面方程 N*P+D=0
      let dis_plane = points.map( ( p ) => Math.abs( new Vector( origin, p )
        .dot( N ) - D ) );
      let [ max_dis_plane, v3 ] = dis_plane.reduce( ( res, cur, idx ) => (
        res[ 0 ] < cur ? [ cur, idx ] : res ), [ 0, v[ 2 ] ] );
      if( max_dis_plane === 0 ) return {
        error: "coplanar",
        normal: N
      }; //所有的点共面
      v[ 3 ] = v3;
      //至此，v[0..3]构成四面体
      let tetrahedron = [ // 从凸包外看去, 顺时针方向存取每个点
        [ v[ 2 ], v[ 1 ], v[ 0 ] ],
        [ v[ 1 ], v[ 3 ], v[ 0 ] ],
        [ v[ 2 ], v[ 3 ], v[ 1 ] ],
        [ v[ 0 ], v[ 3 ], v[ 2 ] ],
      ];
      // 对于面v0-v1-v2,法向量固定为(v0,v2)x(v0,v1). 因此顺序要保证法向量指向凸包外
      let sign = ( new Vector( points[ v[ 0 ] ], points[ v[ 3 ] ] ) )
        .dot( new Vector( points[ v[ 0 ] ], points[ v[ 2 ] ] )
          .cross( new Vector( points[ v[ 0 ] ], points[ v[ 1 ] ] ) ) );
      if( sign < 0 ) tetrahedron = tetrahedron.map( ( x ) => x.reverse() );
      return tetrahedron;
    }
    //  初始化一个由四个点构成的四面体；
    let tetrahedron = initTetrahedron();
    if( !Array.isArray( tetrahedron ) ) return tetrahedron; // 无法构成四面体
    // 未被分配的点
    let v = Array.from( new Set( tetrahedron.reduce( ( a, b ) => a.concat( b ) ) ) );
    let pointsCloned = Array.from( points.keys() )
      .filter( ( p ) => !v.includes( p ) );
    // 四面体的每个面
    for( let f of tetrahedron ) {
      // p上方的点分配到外部点集中
      assignPoints( f, pointsCloned, points );
      if( f.visiblePoints !== undefined && f.visiblePoints.length > 0 ) {
        faceStack.push( f ); // 外部点集非空的面保存到待定面集中
      }
      faces.push( f );
    }
    process( points );
    for( let f of faces ) {
      f = f.map( ( x ) => points[ x ].idx );
    }
    return faces;
  }
}

class ConvexPolygon {
  constructor( points ) {
    if( !Array.isArray( points ) ) points = Array.from( arguments );
    let error_info = ConvexPolygon._check_points( points );
    if( error_info !== "success" ) throw new TypeError( error_info );
    this.points = points;
  }
  static check_points( points ) {
    return ConvexPolygon._check_points( points ) === "success";
  }
  static _check_points( points ) {
    if( points.length >= 3 && points.every( ( p ) => Point.isPoint( p ) ) ) {
      if( points.every( ( p ) => p.dimension() === 2 ) ) {
        let collinear = true;
        for( let i = 0, n = points.length; i < n; i++ ) {
          let sign = Math.sign( ( new Vector( points[ i ], points[ ( i + 1 ) %
              n ] ) )
            .cross( new Vector( points[ i ], points[ ( i + 2 ) % n ] ) ) );
          if( sign < 0 ) return "The polygon should be convex, and points are given in counterclockwise order!";
          if( sign != 0 ) collinear = false;
        }
        if( collinear ) return "All points are collinear!";
        return "success";
      }
      return "All points must in 2d plane";
    }
    return "At least 3 points are required! Please check the number and the type of elements.";
  }
  dimension() {
    return 2;
  }
};

class Triangle extends ConvexPolygon {
  constructor( points ) {
    if( !Array.isArray( points ) ) points = Array.from( arguments );
    if( points.length === 3 ) {
      if( ( new Vector( points[ 0 ], points[ 1 ] ) )
        .cross( new Vector( points[ 1 ], points[ 2 ] ) ) < 0 ) {
        points = points.reverse();
      }
      super( points );
    }
  }
  circumscribed_circle() {
    let [ p1, p2, p3 ] = this.points;
    let a = new Vector( p1, p2 ),
      b = new Vector( p1, p3 );
    let L1 = new Line( p1.move( a.mult( 0.5 ) ), Geometry2D.orthogonal( a ) ); // 两条中垂线方程
    let L2 = new Line( p1.move( b.mult( 0.5 ) ), Geometry2D.orthogonal( b ) );
    let O = Geometry2D.intersection( L1, L2 );
    let r = Geometry2D.distance( O, p1 );
    return new Circle( O, r );
  }
  edges() {
    return this.points.map( ( p, i, a ) => new Segment( p, a[ ( i + 1 ) % 3 ] ) );
  }
  static isTriangle( t ) {
    return t instanceof Triangle;
  }
};