

public class Grid{
  public int n;
  public float[] u, v, u_prev, v_prev;
  public float[] dens;
  public float[] dens_prev;

  public Grid(int n){
    this.n = n;
    int size = (n+2) * (n+2);
    u = new float[size];
    v = new float[size];
    u_prev = new float[size];
    v_prev = new float[size];
    dens = new float[size];
    dens_prev = new float[size];
  }

  int IX(int i, int j, int N) {
    return i + (N + 2) * j;
  }

  void swap(float[] a, float[] b) {
    float[] temp = a.clone();
    for (int i = 0; i < a.length; i++) {
      a[i] = b[i];
      b[i] = temp[i];
    }
  }

  void add_source(int N, float[] x, float[] s, float dt) {
    for (int i = 0; i < N * N; i++) {
      x[i] += dt * s[i];
    }
  }

  void diffuse(int N, int b, float[] x, float[] x0, float diff, float dt) {
    int i, j, k;
    float a = dt * diff * N * N;
    for (k = 0; k < 20; k++) {
      for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
          x[IX(i, j, N)] = (x0[IX(i, j, N)] + a * (x[IX(i - 1, j, N)] + x[IX(i + 1, j, N)] +
                          x[IX(i, j - 1, N)] + x[IX(i, j + 1, N)])) / (1 + 4 * a);
        }
      }
      set_bnd(N, b, x);
    }
  }

  void advect(int N, int b, float[] d, float[] d0, float[] u, float[] v, float dt) {
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt * N;
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        x = i - dt0 * u[IX(i, j, N)];
        y = j - dt0 * v[IX(i, j, N)];
        if (x < 0.5) x = 0.5; if (x > N + 0.5) x = N + 0.5; i0 = (int) x; i1 = i0 + 1;
        if (y < 0.5) y = 0.5; if (y > N + 0.5) y = N + 0.5; j0 = (int) y; j1 = j0 + 1;
        s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
        d[IX(i, j, N)] = s0 * (t0 * d0[IX(i0, j0, N)] + t1 * d0[IX(i0, j1, N)]) +
                        s1 * (t0 * d0[IX(i1, j0, N)] + t1 * d0[IX(i1, j1, N)]);
      }
    }
    set_bnd(N, b, d);
  }

  void dens_step(int N, float[] x, float[] x0, float[] u, float[] v, float diff, float dt) {
  add_source(N, x, x0, dt);
  swap(x0, x); 
  diffuse(N, 0, x, x0, diff, dt);
  swap(x0, x); 
  advect(N, 0, x, x0, u, v, dt);
}

  void vel_step(int N, float[] u, float[] v, float[] u0, float[] v0, float visc, float dt) {
    add_source(N, u, u0, dt);
    add_source(N, v, v0, dt);
    swap(u0, u);
    diffuse(N, 1, u, u0, visc, dt);
    swap(v0, v);
    diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    swap(u0, u);
    swap(v0, v);
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
  }

  void project(int N, float[] u, float[] v, float[] p, float[] div) {
    int i, j, k;
    float h = 1.0 / N;
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        div[IX(i, j, N)] = -0.5 * h * (u[IX(i + 1, j, N)] - u[IX(i - 1, j, N)] +
                                      v[IX(i, j + 1, N)] - v[IX(i, j - 1, N)]);
        p[IX(i, j, N)] = 0;
      }
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);
    for (k = 0; k < 20; k++) {
      for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
          p[IX(i, j, N)] = (div[IX(i, j, N)] + p[IX(i - 1, j, N)] + p[IX(i + 1, j, N)] +
                            p[IX(i, j - 1, N)] + p[IX(i, j + 1, N)]) / 4;
        }
      }
      set_bnd(N, 0, p);
    }
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        u[IX(i, j, N)] -= 0.5 * (p[IX(i + 1, j, N)] - p[IX(i - 1, j, N)]) / h;
        v[IX(i, j, N)] -= 0.5 * (p[IX(i, j + 1, N)] - p[IX(i, j - 1, N)]) / h;
      }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
  }

  void set_bnd(int N, int b, float[] x) {
    int i;
    for (i = 1; i <= N; i++) {
      x[IX(0, i, N)] = b == 1 ? -x[IX(1, i, N)] : x[IX(1, i, N)];
      x[IX(N + 1, i, N)] = b == 1 ? -x[IX(N, i, N)] : x[IX(N, i, N)];
      x[IX(i, 0, N)] = b == 2 ? -x[IX(i, 1, N)] : x[IX(i, 1, N)];
      x[IX(i, N + 1, N)] = b == 2 ? -x[IX(i, N, N)] : x[IX(i, N, N)];
    }
    x[IX(0, 0, N)] = 0.5 * (x[IX(1, 0, N)] + x[IX(0, 1, N)]);
    x[IX(0, N + 1, N)] = 0.5 * (x[IX(1, N + 1, N)] + x[IX(0, N, N)]);
    x[IX(N + 1, 0, N)] = 0.5 * (x[IX(N, 0, N)] + x[IX(N + 1, 1, N)]);
    x[IX(N + 1, N + 1, N)] = 0.5 * (x[IX(N, N + 1, N)] + x[IX(N + 1, N, N)]);
  }
}

float g = 9.81;








// physics and shape
public class Line{
 public Vec2 pt1, pt2;
  
 public Line(Vec2 pt1, Vec2 pt2){
   this.pt1 = pt1;
   this.pt2 = pt2;
 }

 public Vec2 vec(){
   return pt2.minus(pt1);
 }

 public float length(){
   return pt1.minus(pt2).length();
 }

 public void display(){
   line(pt1.x, pt1.y, pt2.x, pt2.y);
 }

 public String toString(){
   return pt1.toString();
 }
}


public class Circle{
 public Vec2 center;
 public float radius;
  
 public Circle(Vec2 center, float radius){
   this.center = center;
   this.radius = radius;
 }

 public void display(){
   circle(center.x, center.y, radius*2);
 }

 public String toString(){
   return center.toString();
 }
}


public boolean colliding(Line l, Circle c){ // check
 Vec2 toCircle1 = c.center.minus(l.pt1);
 Vec2 toCircle2 = c.center.minus(l.pt2);
 if (toCircle1.length() <= c.radius || toCircle2.length() <= c.radius) return true;

 float a = 1;  //Lenght of l_dir (we noramlized it)
 float b = -2*dot((l.vec()).normalized(),toCircle1); //-2*dot(l_dir,toCircle)
 float c_val = toCircle1.lengthSqr() - (c.radius)*(c.radius); //different of squared distances
  
 float d = b*b - 4*a*c_val; //discriminant 
  
 if (d >=0){ 
   //If d is positive we know the line is colliding, but we need to check if the collision line within the line segment
   //  ... this means t will be between 0 and the lenth of the line segment
   float t1 = (-b - sqrt(d))/(2*a); //Optimization: we only take the first collision [is this safe?]
   if (t1 > 0 && t1 < l.length()){
     return true;
   } 
 }
 return false;
}

Line collisionResponseStatic(Circle ball, Line line){
 Vec2 v1 = ball.center.minus(line.pt1);
 Vec2 v2 = line.pt2.minus(line.pt1);
 float proj = dot(v2, v1) / v2.length();
 Vec2 closest = line.pt1.plus(v2.normalized().times(proj));
 Vec2 dist = ball.center.minus(closest);

 Vec2 normal = new Vec2(-v2.y, v2.x).normalized();

 float d = dot(normal, dist);
 if (d < 0){
   normal.mul(-1);
 }
 normal.mul(-1);

 Vec2 outside_point = ball.center.plus(normal.times(ball.radius));

 Vec2 new_pt2 = outside_point.minus(line.pt1).normalized().times(v2.length());
 return new Line(new_pt2, line.pt1);
 //float angle = atan2(new_pt2.y - line.pt1.y, new_pt2.x - line.pt1.x);
 // ball.center = closest.plus(normal.times(ball.radius));
 // Vec2 velNormal = normal.times(dot(ball.vel,normal));
 // ball.vel.subtract(velNormal.times(1 + cor));
 //return angle;
}




// Camera

class Camera
{
Camera()
{
  position      = new PVector( 0, 0, 0 ); // initial position
  theta         = 0; // rotation around Y axis. Starts with forward direction as ( 0, 0, -1 )
  phi           = 0; // rotation around X axis. Starts with up direction as ( 0, 1, 0 )
  moveSpeed     = 50;
  turnSpeed     = 1.57; // radians/sec
  boostSpeed    = 10;  // extra speed boost for when you press shift
    
  // dont need to change these
  shiftPressed = false;
  negativeMovement = new PVector( 0, 0, 0 );
  positiveMovement = new PVector( 0, 0, 0 );
  negativeTurn     = new PVector( 0, 0 ); // .x for theta, .y for phi
  positiveTurn     = new PVector( 0, 0 );
  fovy             = PI / 4;
  aspectRatio      = width / (float) height;
  nearPlane        = 0.1;
  farPlane         = 10000;
}
void Update(float dt)
{
  theta += turnSpeed * ( negativeTurn.x + positiveTurn.x)*dt;
    
  // cap the rotation about the X axis to be less than 90 degrees to avoid gimble lock
  float maxAngleInRadians = 85 * PI / 180;
  phi = min( maxAngleInRadians, max( -maxAngleInRadians, phi + turnSpeed * ( negativeTurn.y + positiveTurn.y ) * dt ) );
    
  // re-orienting the angles to match the wikipedia formulas: https://en.wikipedia.org/wiki/Spherical_coordinate_system
  // except that their theta and phi are named opposite
  float t = theta + PI / 2;
  float p = phi + PI / 2;
  PVector forwardDir = new PVector( sin( p ) * cos( t ),   cos( p ),   -sin( p ) * sin ( t ) );
  PVector upDir      = new PVector( sin( phi ) * cos( t ), cos( phi ), -sin( t ) * sin( phi ) );
  PVector rightDir   = new PVector( cos( theta ), 0, -sin( theta ) );
  PVector velocity   = new PVector( negativeMovement.x + positiveMovement.x, negativeMovement.y + positiveMovement.y, negativeMovement.z + positiveMovement.z );
  position.add( PVector.mult( forwardDir, moveSpeed * velocity.z * dt ) );
  position.add( PVector.mult( upDir,      moveSpeed * velocity.y * dt ) );
  position.add( PVector.mult( rightDir,   moveSpeed * velocity.x * dt ) );
    
  aspectRatio = width / (float) height;
  perspective( fovy, aspectRatio, nearPlane, farPlane );
  camera( position.x, position.y, position.z,
          position.x + forwardDir.x, position.y + forwardDir.y, position.z + forwardDir.z,
          upDir.x, upDir.y, upDir.z );
}
  
// only need to change if you want difrent keys for the controls
void HandleKeyPressed()
{
  if ( key == 'w' || key == 'W' ) positiveMovement.z = 1;
  if ( key == 's' || key == 'S' ) negativeMovement.z = -1;
  if ( key == 'a' || key == 'A' ) negativeMovement.x = -1;
  if ( key == 'd' || key == 'D' ) positiveMovement.x = 1;
  if ( key == 'q' || key == 'Q' ) positiveMovement.y = 1;
  if ( key == 'e' || key == 'E' ) negativeMovement.y = -1;
    
  if ( key == 'r' || key == 'R' ){
   //  Camera defaults = new Camera();
   //  position = defaults.position;
   //  theta = defaults.theta;
   //  phi = defaults.phi;
  //  camera.position = new PVector(719.04254, 117.900024, 545.1439);
  //  camera.theta = 0.131;
  //  camera.phi = -0.6536;
  }
    
  if ( keyCode == LEFT )  negativeTurn.x = 1;
  if ( keyCode == RIGHT ) positiveTurn.x = -0.5;
  if ( keyCode == UP )    positiveTurn.y = 0.5;
  if ( keyCode == DOWN )  negativeTurn.y = -1;
    
  if ( keyCode == SHIFT ) shiftPressed = true; 
  if (shiftPressed){
    positiveMovement.mult(boostSpeed);
    negativeMovement.mult(boostSpeed);
  }

  if ( key == 'p' || key == 'P'){
    println("position:", position.x, position.y, position.z);
    println("theta:", theta);
    println("phi:", phi);

  }
}

  void HandleKeyReleased()
{
  if ( key == 'w' || key == 'W' ) positiveMovement.z = 0;
  if ( key == 'q' || key == 'Q' ) positiveMovement.y = 0;
  if ( key == 'd' || key == 'D' ) positiveMovement.x = 0;
  if ( key == 'a' || key == 'A' ) negativeMovement.x = 0;
  if ( key == 's' || key == 'S' ) negativeMovement.z = 0;
  if ( key == 'e' || key == 'E' ) negativeMovement.y = 0;
    
  if ( keyCode == LEFT  ) negativeTurn.x = 0;
  if ( keyCode == RIGHT ) positiveTurn.x = 0;
  if ( keyCode == UP    ) positiveTurn.y = 0;
  if ( keyCode == DOWN  ) negativeTurn.y = 0;
    
  if ( keyCode == SHIFT ){
    shiftPressed = false;
    positiveMovement.mult(1.0/boostSpeed);
    negativeMovement.mult(1.0/boostSpeed);
  }
}
  
// only necessary to change if you want different start position, orientation, or speeds
PVector position;
float theta;
float phi;
float moveSpeed;
float turnSpeed;
float boostSpeed;
  
// probably don't need / want to change any of the below variables
float fovy;
float aspectRatio;
float nearPlane;
float farPlane;  
PVector negativeMovement;
PVector positiveMovement;
PVector negativeTurn;
PVector positiveTurn;
boolean shiftPressed;
};




public class Vec2 {
 public float x, y;
  
 public Vec2(float x, float y){
   this.x = x;
   this.y = y;
 }
  
 public String toString(){
   return "(" + x+ "," + y +")";
 }
  
 public float length(){
   return sqrt(x*x+y*y);
 }
  
 public float lengthSqr(){
   return x*x+y*y;
 }

 public Vec2 plus(Vec2 rhs){
   return new Vec2(x+rhs.x, y+rhs.y);
 }
  
 public void add(Vec2 rhs){
   x += rhs.x;
   y += rhs.y;
 }
  
 public Vec2 minus(Vec2 rhs){
   return new Vec2(x-rhs.x, y-rhs.y);
 }
  
 public void subtract(Vec2 rhs){
   x -= rhs.x;
   y -= rhs.y;
 }
  
 public Vec2 times(float rhs){
   return new Vec2(x*rhs, y*rhs);
 }
  
 public void mul(float rhs){
   x *= rhs;
   y *= rhs;
 }
  
 public void clampToLength(float maxL){
   float magnitude = sqrt(x*x + y*y);
   if (magnitude > maxL){
     x *= maxL/magnitude;
     y *= maxL/magnitude;
   }
 }
  
 public void setToLength(float newL){
   float magnitude = sqrt(x*x + y*y);
   x *= newL/magnitude;
   y *= newL/magnitude;
 }
  
 public void normalize(){
   float magnitude = sqrt(x*x + y*y);
   x /= magnitude;
   y /= magnitude;
 }
  
 public Vec2 normalized(){
   float magnitude = sqrt(x*x + y*y);
   return new Vec2(x/magnitude, y/magnitude);
 }
  
 public float distanceTo(Vec2 rhs){
   float dx = rhs.x - x;
   float dy = rhs.y - y;
   return sqrt(dx*dx + dy*dy);
 }
}

Vec2 interpolate(Vec2 a, Vec2 b, float t){
 return a.plus((b.minus(a)).times(t));
}

float interpolate(float a, float b, float t){
 return a + ((b-a)*t);
}

float dot(Vec2 a, Vec2 b){
 return a.x*b.x + a.y*b.y;
}

float cross(Vec2 a, Vec2 b){
 return a.x*b.y - a.y*b.x;
}

Vec2 mid_point(Vec2 a, Vec2 b){
 return new Vec2((a.x + b.x) / 2, (a.y + b.y) / 2);
}

Vec2 projAB(Vec2 a, Vec2 b){
 return b.times(a.x*b.x + a.y*b.y);
}

float clamp(float f, float min, float max){
 if (f < min) return min;
 if (f > max) return max;
 return f;
}
