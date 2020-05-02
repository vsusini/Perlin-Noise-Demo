/*
 * A speed-improved perlin and simplex noise algorithms for 2D.
 *
 * Based on example code by Stefan Gustavson (stegu@itn.liu.se).
 * Optimisations by Peter Eastman (peastman@drizzle.stanford.edu).
 * Better rank ordering method by Stefan Gustavson in 2012.
 * Converted to Javascript by Joseph Gentle.
 *
 * Version 2012-03-09
 *
 * This code was placed in the public domain by its original author,
 * Stefan Gustavson. You may use it as you see fit, but
 * attribution is appreciated.
 *
 */

(function(global){
    var module = global.noise = {};

    function Grad(x, y, z) {
        this.x = x; this.y = y; this.z = z;
    }
    
    Grad.prototype.dot2 = function(x, y) {
        return this.x*x + this.y*y;
    };

    Grad.prototype.dot3 = function(x, y, z) {
        return this.x*x + this.y*y + this.z*z;
    };

    var grad3 = [new Grad(1,1,0),new Grad(-1,1,0),new Grad(1,-1,0),new Grad(-1,-1,0),
                    new Grad(1,0,1),new Grad(-1,0,1),new Grad(1,0,-1),new Grad(-1,0,-1),
                    new Grad(0,1,1),new Grad(0,-1,1),new Grad(0,1,-1),new Grad(0,-1,-1)];
    var p = [151,160,137,91,90,15,
        131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
        190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
        88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
        77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
        102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
        135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
        5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
        223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
        129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
        251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
        49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
        138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180];
    // To remove the need for index wrapping, double the permutation table length
    var perm = new Array(512);
    var gradP = new Array(512);

    // This isn't a very good seeding function, but it works ok. It supports 2^16
    // different seed values. Write something better if you need more seeds.
    module.seed = function(seed) {
      if(seed > 0 && seed < 1) {
        // Scale the seed out
        seed *= 65536;
      }
  
      seed = Math.floor(seed);
      if(seed < 256) {
        seed |= seed << 8;
      }
  
      for(var i = 0; i < 256; i++) {
        var v;
        if (i & 1) {
          v = p[i] ^ (seed & 255);
        } else {
          v = p[i] ^ ((seed>>8) & 255);
        }
  
        perm[i] = perm[i + 256] = v;
        gradP[i] = gradP[i + 256] = grad3[v % 12];
      }
    };

    module.seed(0);

    // Skewing and unskewing factors for 2, 3, and 4 dimensions
    var G2 = (3-Math.sqrt(3))/6;

    // 2D simplex noise
    module.simplex2 = function(xin, yin) {
      var n0, n1, n2; // Noise contributions from the three corners
      // Skew the input space to determine which simplex cell we're in
      var s = (xin+yin)*0.5*(Math.sqrt(3)-1);; // Hairy factor for 2D
      var i = Math.floor(xin+s);
      var j = Math.floor(yin+s);
      var t = (i+j)*G2;
      var x0 = xin-i+t; // The x,y distances from the cell origin, unskewed.
      var y0 = yin-j+t;
      // For the 2D case, the simplex shape is an equilateral triangle.
      // Determine which simplex we are in.
      var i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
      if(x0>y0) { // lower triangle, XY order: (0,0)->(1,0)->(1,1)
        i1=1; j1=0;
      } else {    // upper triangle, YX order: (0,0)->(0,1)->(1,1)
        i1=0; j1=1;
      }
      // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
      // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
      // c = (3-sqrt(3))/6
      var x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
      var y1 = y0 - j1 + G2;
      var x2 = x0 - 1 + 2 * G2; // Offsets for last corner in (x,y) unskewed coords
      var y2 = y0 - 1 + 2 * G2;
      // Work out the hashed gradient indices of the three simplex corners
      i &= 255;
      j &= 255;
      var gi0 = gradP[i+perm[j]];
      var gi1 = gradP[i+i1+perm[j+j1]];
      var gi2 = gradP[i+1+perm[j+1]];
      // Calculate the contribution from the three corners
      var t0 = 0.5 - x0*x0-y0*y0;
      if(t0<0) {
        n0 = 0;
      } else {
        t0 *= t0;
        n0 = t0 * t0 * gi0.dot2(x0, y0);  // (x,y) of grad3 used for 2D gradient
      }
      var t1 = 0.5 - x1*x1-y1*y1;
      if(t1<0) {
        n1 = 0;
      } else {
        t1 *= t1;
        n1 = t1 * t1 * gi1.dot2(x1, y1);
      }
      var t2 = 0.5 - x2*x2-y2*y2;
      if(t2<0) {
        n2 = 0;
      } else {
        t2 *= t2;
        n2 = t2 * t2 * gi2.dot2(x2, y2);
      }
      // Add contributions from each corner to get the final noise value.
      // The result is scaled to return values in the interval [-1,1].
      return 70 * (n0 + n1 + n2);
    };

    // ##### Perlin noise stuff

    function fade(t) {
      return t*t*t*(t*(t*6-15)+10);
    }

    function lerp(a, b, t) {
      return (1-t)*a + t*b;
    }

    // 2D Perlin Noise
    module.perlin2 = function(x, y) {
      // Find unit grid cell containing point
      var X = Math.floor(x), Y = Math.floor(y);
      // Get relative xy coordinates of point within that cell
      x = x - X; y = y - Y;
      // Wrap the integer cells at 255 (smaller integer period can be introduced here)
      X = X & 255; Y = Y & 255;
  
      // Calculate noise contributions from each of the four corners
      var n00 = gradP[X+perm[Y]].dot2(x, y);
      var n01 = gradP[X+perm[Y+1]].dot2(x, y-1);
      var n10 = gradP[X+1+perm[Y]].dot2(x-1, y);
      var n11 = gradP[X+1+perm[Y+1]].dot2(x-1, y-1);
  
      // Compute the fade curve value for x
      var u = fade(x);
  
      // Interpolate the four results
      return lerp(lerp(n00, n10, u), lerp(n01, n11, u), fade(y));
    };

})(this);

//effective animation code
var wWidth = window.innerWidth;
var wHeight = window.innerHeight;
var scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera(75, wWidth / wHeight, 0.01, 1000);

  //Camera Positioning
camera.position.x = 1.1;
camera.position.y = 150;
camera.position.z = 87;
camera.lookAt(0,60,-50);

var renderer = new THREE.WebGLRenderer({
    alpha: true
});
renderer.setClearColor(new THREE.Color(0x3d3d3d));
document.getElementById('canvas').appendChild(renderer.domElement);

var noiseTypeList = ["Perlin 2D","Simplex 2D"];

var controls = new function () {
    this.noiseType = "Perlin 2D"
    this.noiseFlag = 0;
    //Wave Information
    this.floorColor = "rgb(83,131,248)"
    this.middleColor = "rgb(39,175,86)"
    this.topMiddleColor = "rgb(0,0,150)"
    this.topColor = "rgb(255,255,255)"
    this.rows = 125;
    this.cols = 200;
    this.perlinScale = 0.01
    this.waveSpeed = 0.2
    this.waveHeight = 30;
    this.dotSize = 0.01;
    this.fpAngleFlag = 0;

    //Buttons
    this.angle = 0;

    //Controler Functions
    this.getColor = function(color){
        var colorsRGB = color.split("(")[1].split(",");
        colorsRGB[0] = colorsRGB[0]/255;
        colorsRGB[1] = colorsRGB[1]/255;
        colorsRGB[2] = colorsRGB[2].substring(0,colorsRGB[2].length-1)/255;
        colorsRGB = new THREE.Color(colorsRGB[0], colorsRGB[1], colorsRGB[2]);
        return colorsRGB
    }
    this.getFloorColor = function(){
        return this.getColor(this.floorColor)
    }
    this.getMiddleColor = function (){
        return this.getColor(this.middleColor)
    }
    this.getTopColor = function (){
        return this.getColor(this.topColor)
    }
    this.getTopMiddleColor = function (){
        return this.getColor(this.topMiddleColor)
    }
    this.setFloorColor = function(value){
        this.FloorColor = value;
        updateGeo(controls.getFloorColor());
    }
    this.setMiddleColor = function(value){
        this.middleColor = value;
        updateGeo(controls.getMiddleColor());
    }
    this.setTopColor = function(value){
        this.topColor = value;
        updateGeo(controls.getTopColor());
    }
    this.setTopMiddleColor = function(value){
        this.topMiddleColor = value;
        updateGeo(controls.getTopMiddleColor());
    }
    this.updateSeed = function(){
        updateSeed();
    }
    this.updateRows = function(value){
        this.rows = value
        updateGeo(controls.getFloorColor());
    }
    this.updateCols = function(value){
        this.cols = value
        updateGeo(controls.getFloorColor());
    }
    // this.updateSpeed = function(value){
    //     this.waveSpeed = value;
    //     updateGeo(controls.getFloorColor());
    // }
    this.updateWaveHeight = function(value){
        this.waveHeight = value;
        updateGeo(controls.colorArray)
    }
    this.updatePerlinScale = function(value){
        this.perlinScale = value;
        updateGeo(controls.colorArray)
    }
    this.updateDotSize = function (value){
      this.dotSize = value;
      updateGeo(controls.colorArray);
    }
    this.skyAngle = function(){
      this.fpAngleFlag = 0;
      camera.position.y = 150;
      camera.position.x = 1.1;
      camera.position.z = 87; 
      camera.lookAt(0,60,-50);
    }
    this.normAngle = function(){
      this.fpAngleFlag = 0;
      camera.position.y = controls.waveHeight;
      camera.position.x = 1.1;
      camera.position.z = 87; 
      camera.lookAt(0,60,-50);
    }
    this.fpAngle = function(){
      this.fpAngleFlag = 1;
      camera.position.z = -200;
      camera.position.x = -200;
      
      camera.lookAt(-60,60,-50);
    }
    this.changeNoise = function(){
      if (this.noiseFlag == 0){
        this.noiseType = noiseTypeList[1]
        this.noiseFlag = 1;
      }else{
        this.noiseType = noiseTypeList[0];
        this.noiseFlag = 0;
        
      }
    }
}

  //Animation parameters
var startTime = new Date().getTime();

function updateSeed(){
    noise.seed(Math.random());
}
updateSeed();

function createGeometry(color) {
    var geometry = new THREE.Geometry();

    for (var y = 0; y < controls.rows; y++) {
        for (var x = 0; x < controls.cols; x++) {
            geometry.vertices.push(
                new THREE.Vector3(x, 0, y)
            );
            geometry.colors.push(color)
        }
    }
    geometry.dynamic = true;
    //Place the geo into the correct location.
    geometry.translate(-200, 0, -200);
    return geometry;
}

  //Set to check if pointGroup exists or not
var pointGroup = 0;
function updateGeo(colorString){
    //Check if it exists, then remove it
    if(pointGroup != 0){
        scene.remove(pointGroup);
    }
    var geo = createGeometry(colorString);
    pointGroup = new THREE.Points(geo, new THREE.PointsMaterial({
        size: controls.dotSize,
        vertexColors: THREE.VertexColors
    }));
    scene.add(pointGroup);
}
  //Create the org Geo
updateGeo(controls.getFloorColor());

function animatePoints(colorArray,func) {
    var currentTime = new Date().getTime();
    var i = 0;
    for (var y = 0; y < controls.rows; y++) {
        for (var x = 0; x < controls.cols; x++) {
            var pX = (x * controls.perlinScale) + ((currentTime - startTime) / 1000) * controls.waveSpeed;
            var pY = (y * controls.perlinScale) + ((currentTime - startTime) / 1000) * controls.waveSpeed;
            pointGroup.geometry.vertices[i].y = (func(pX, pY)) * controls.waveHeight;

            //Floor 
            if (pointGroup.geometry.vertices[i].y <= -controls.waveHeight/3){ //-10 to Start
                pointGroup.geometry.colors[i] = colorArray[0];
            }else if (pointGroup.geometry.vertices[i].y <= -controls.waveHeight/30){
                pointGroup.geometry.colors[i] = colorArray[1];
            }else if (pointGroup.geometry.vertices[i].y <= controls.waveHeight/3){
                pointGroup.geometry.colors[i] = colorArray[2];
            }else{
                pointGroup.geometry.colors[i] = colorArray[3];
            }
            i += 1;
        }
    }
    pointGroup.geometry.verticesNeedUpdate = true;
    pointGroup.geometry.colorsNeedUpdate = true;
    return pointGroup.geometry.vertices[0].y
}

function animate() {
  var camHeight;
    if (controls.noiseFlag == 0){
      camHeight = animatePoints([controls.getFloorColor(),controls.getMiddleColor(), controls.getTopMiddleColor(), controls.getTopColor()],noise.perlin2);
    }else{
      camHeight = animatePoints([controls.getFloorColor(),controls.getMiddleColor(), controls.getTopMiddleColor(), controls.getTopColor()],noise.simplex2)
    }
    if (controls.fpAngleFlag == 1){
      camera.position.y = camHeight + 3;
    }
    renderer.render(scene, camera);
    requestAnimationFrame(animate);
}

  //GUI Setup
var gui = new dat.GUI();
  var f1 = gui.addFolder("Wave Settings")
    f1.open();
    f1.add(controls,"rows",10,300,1).name("Row Number").onChange(controls.updateRows)
    f1.add(controls,"cols",10,400,1).name("Col Number").onChange(controls.updateCols)
    // f1.add(controls,"waveSpeed",0.001,0.5,0.01).name("Wave Speed").onChange(controls.updateSpeed)
    f1.add(controls,"waveHeight",5,50,1).name("Wave Height").onChange(controls.updateWaveHeight)
    f1.add(controls,"perlinScale",0.01,0.1,0.01).name("Perlin Scale").onChange(controls.updatePerlinScale)
    f1.add(controls,"dotSize",0.0001,3,0.01).name("Dot Size").onChange(controls.updateDotSize)
    f1.add(controls,"updateSeed").name("Generate Seed");
    f1.add(controls,"changeNoise").name("Change Noise")
    f1.add(controls, 'noiseType').name("Noise Algorithm").listen(); 

  var f2 = gui.addFolder("Color Settings")
    f2.addColor(controls, 'topColor').name("Top Color").onChange(controls.setTopColor);
    f2.addColor(controls, 'topMiddleColor').name("Top-Middle Color").onChange(controls.setTopMiddleColor);
    f2.addColor(controls, 'middleColor').name("Middle Color").onChange(controls.setMiddleColor);
    f2.addColor(controls, 'floorColor').name("Floor Color").onChange(controls.setFloorColor);

  var f3 = gui.addFolder("Camera Options");
    f3.add(controls,"skyAngle").name("Sky Camera");
    f3.add(controls,"normAngle").name("Wave-Height Camera")
    f3.add(controls,"fpAngle").name("First Person");


  //Events needed to create animation and any listeners such as Canas
window.addEventListener('resize', refreshCanvasState, false);
animate();
refreshCanvasState();

function refreshCanvasState() {
    wWidth = window.innerWidth;
    wHeight = window.innerHeight*0.95;
    camera.aspect = wWidth / wHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(wWidth, wHeight);
}
