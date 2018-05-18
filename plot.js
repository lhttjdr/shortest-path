function Graph(config) {
    // user defined properties
    this.canvas = document.getElementById(config.canvasId);
    this.minX = config.minX;
    this.minY = config.minY;
    this.maxX = config.maxX;
    this.maxY = config.maxY;
    this.unitsPerTick = config.unitsPerTick;

    // constants
    this.axisColor = '#aaa';
    this.font = '8pt Calibri';
    this.tickSize = 20;

    // relationships
    this.context = this.canvas.getContext('2d');
    this.rangeX = this.maxX - this.minX;
    this.rangeY = this.maxY - this.minY;
    this.unitX = this.canvas.width / this.rangeX;
    this.unitY = this.canvas.height / this.rangeY;
    this.centerY = Math.round(Math.abs(this.minY / this.rangeY) * this.canvas.height);
    this.centerX = Math.round(Math.abs(this.minX / this.rangeX) * this.canvas.width);
    this.iteration = (this.maxX - this.minX) / 1000;
    this.scaleX = this.canvas.width / this.rangeX;
    this.scaleY = this.canvas.height / this.rangeY;

    // draw x and y axis
    this.drawXAxis();
    this.drawYAxis();
}

Graph.prototype = {
    prepare: function() {

    },
    drawXAxis: function() {
        var context = this.context;
        context.save();
        context.beginPath();
        context.moveTo(0, this.centerY);
        context.lineTo(this.canvas.width, this.centerY);
        context.strokeStyle = this.axisColor;
        context.lineWidth = 2;
        context.stroke();

        // draw tick marks
        var xPosIncrement = this.unitsPerTick * this.unitX;
        var xPos, unit;
        context.font = this.font;
        context.textAlign = 'center';
        context.textBaseline = 'top';

        // draw left tick marks
        xPos = this.centerX - xPosIncrement;
        unit = -1 * this.unitsPerTick;
        while (xPos > 0) {
            context.moveTo(xPos, this.centerY - this.tickSize / 2);
            context.lineTo(xPos, this.centerY + this.tickSize / 2);
            context.stroke();
            context.fillText(unit, xPos, this.centerY + this.tickSize / 2 + 3);
            unit -= this.unitsPerTick;
            xPos = Math.round(xPos - xPosIncrement);
        }

        // draw right tick marks
        xPos = this.centerX + xPosIncrement;
        unit = this.unitsPerTick;
        while (xPos < this.canvas.width) {
            context.moveTo(xPos, this.centerY - this.tickSize / 2);
            context.lineTo(xPos, this.centerY + this.tickSize / 2);
            context.stroke();
            context.fillText(unit, xPos, this.centerY + this.tickSize / 2 + 3);
            unit += this.unitsPerTick;
            xPos = Math.round(xPos + xPosIncrement);
        }
        context.restore();
    },
    drawYAxis: function() {
        var context = this.context;
        context.save();
        context.beginPath();
        context.moveTo(this.centerX, 0);
        context.lineTo(this.centerX, this.canvas.height);
        context.strokeStyle = this.axisColor;
        context.lineWidth = 2;
        context.stroke();

        // draw tick marks
        var yPosIncrement = this.unitsPerTick * this.unitY;
        var yPos, unit;
        context.font = this.font;
        context.textAlign = 'right';
        context.textBaseline = 'middle';

        // draw top tick marks
        yPos = this.centerY - yPosIncrement;
        unit = this.unitsPerTick;
        while (yPos > 0) {
            context.moveTo(this.centerX - this.tickSize / 2, yPos);
            context.lineTo(this.centerX + this.tickSize / 2, yPos);
            context.stroke();
            context.fillText(unit, this.centerX - this.tickSize / 2 - 3, yPos);
            unit += this.unitsPerTick;
            yPos = Math.round(yPos - yPosIncrement);
        }

        // draw bottom tick marks
        yPos = this.centerY + yPosIncrement;
        unit = -1 * this.unitsPerTick;
        while (yPos < this.canvas.height) {
            context.moveTo(this.centerX - this.tickSize / 2, yPos);
            context.lineTo(this.centerX + this.tickSize / 2, yPos);
            context.stroke();
            context.fillText(unit, this.centerX - this.tickSize / 2 - 3, yPos);
            unit -= this.unitsPerTick;
            yPos = Math.round(yPos + yPosIncrement);
        }
        context.restore();
    },
    drawEquation: function(equation, color, thickness) {
        var context = this.context;
        context.save();
        context.save();
        this.transformContext();

        context.beginPath();
        context.moveTo(this.minX, equation(this.minX));

        for (var x = this.minX + this.iteration; x <= this.maxX; x += this.iteration) {
            context.lineTo(x, equation(x));
        }

        context.restore();
        context.lineJoin = 'round';
        context.lineWidth = thickness;
        context.strokeStyle = color;
        context.stroke();
        context.restore();
    },
    plot: function(curve, color, thickness, fillcolor, alpha){
        if(curve.dimension()!==2) return;

        var context = this.context;
        context.save();
        context.save();
        this.transformContext();

        context.beginPath();
        if(curve instanceof Segment){
          let [x1,y1]=curve.points[0].coordinates();
          context.moveTo(x1, y1);
          let [x2,y2]=curve.points[1].coordinates();
          context.lineTo(x2, y2);
        }if(curve instanceof Line){
          let [vx,vy]=curve.v.components();
          let [x0,y0]=curve.p.coordinates();
          if(vx===0){
            context.moveTo(x0, this.minY);
            context.lineTo(x0, this.maxY);
          }else if(vy==0){
            context.moveTo(this.minX, y0);
            context.lineTo(this.minY, y0);
          }else{
            let k=vy*1.0/vx;
            let xmin=this.minX, xmax=this.maxX;
            let ymin=k*(xmin-x0)+y0, ymax=k*(xmax-x0)+y0;
            let k1=vx*1.0/vy;
            if(ymin< this.minY){
              ymin=this.minY;
              xmin=k1*(ymin-y0)+x0;
            }
            if(ymax> this.maxY){
              ymax=this.maxY;
              xmax=k1*(ymax-y0)+x0;
            }
            context.moveTo(xmin, ymin);
            context.lineTo(xmax, ymax);
          }
        }else if(curve instanceof Circle){
          let [x,y]=curve.center.coordinates();
          let r=curve.radius;
          context.arc(x, y, r, 0, 2 * Math.PI, false);
        }else if(curve instanceof Point){
          let [x,y]=curve.coordinates();
          let r=5;
          context.arc(x, y, r, 0, 2 * Math.PI, false);
        }else if(curve instanceof ConvexPolygon || curve instanceof Triangle){
          let points=curve.points.slice();
          let [x0,y0]=points[0].coordinates();
          points.push(points[0]);
          points.unshift();
          context.moveTo(x0,y0);
          for(let p of points){
            let [x,y]=p.coordinates();
            context.lineTo(x,y);
          }
        }
        context.restore();
        context.lineJoin = 'round';
        context.lineWidth = thickness || 1;
        context.strokeStyle = color || 'grey';
        context.fillStyle = fillcolor || 'rgba(0,0,0,0)';
        context.stroke();
        context.fill();
        context.restore();
    },
    transformContext: function() {
        var context = this.context;

        // move context to center of canvas
        this.context.translate(this.centerX, this.centerY);

        /*
         * stretch grid to fit the canvas window, and
         * invert the y scale so that that increments
         * as you move upwards
         */
        context.scale(this.scaleX, -this.scaleY);
    },
    clear: function() {
        var context = this.context;
        context.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.drawXAxis();
        this.drawYAxis();
    }
};
