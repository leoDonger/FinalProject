import java.util.Iterator;
import java.util.Random;
import beads.*;

AudioContext ac;
SamplePlayer player;
Gain masterGain;
Static rateUGen;
float dyingRate = 1;


Cigarette cig;
Grid grid;
Grid fire;
int N = 256; // Grid size
float dt = 1.0 / 120;
// float dt = 1.0 / frameRate;
float diff = 0.001; // Diffusion rate
float visc = 0.01; // Viscosity
boolean paused;

int arraySize = (N+2) * (N+2);
float[] antiGravity = new float[arraySize];
float[] chaoticEffects = new float[arraySize];

int IX(int i, int j) {
  return i + (N + 2) * j;
}

public class Grid{
  public int n;
  public float[] u, v, u_prev, v_prev;
  public float[] dens;
  public float[] dens_prev;
  private float[] temp_for_swap;

  // Fire
  public float[] oxy, fuel, heat, heat_prev, oxy_prev;
  public float energyBarrier, rateConstant, maxRate;
  public float exothermicness, convectivness;
  public float viscC, diffC, diffC_Heat, diffC_O2;
  

  // Professor Guy's code
  void updateFire(){
  // Compute the "burn equation" for each cell (O,F -> H)
    int nx = n+2;
    int ny = n+2;
    for (int i = 0; i < nx; i++){
      for (int j = 0; j < ny; j++){
        int index = IX(i, j);
        float O = oxy[index];
        float F = fuel[index];
        float H = heat[index] + 30;  //Ambient temp
              
        float reactionRate = (O*F*H - energyBarrier) * rateConstant;
              
        if (reactionRate < 0) reactionRate = 0;
        if (reactionRate > maxRate) reactionRate = maxRate;

        oxy[index] -= reactionRate * dt;
        fuel[index] -= reactionRate * dt;
        heat[index] += reactionRate * dt * exothermicness;
                  
        if (oxy[index] < 0) oxy[index] = 0;
        if (fuel[index] < 0) fuel[index] = 0;
      }
    }

    // Add an updraft for convection due to heat
    for (int i=1 ; i<nx ; i++){
      for (int j=1 ; j<ny ; j++){
        int index = IX(i, j);
        v[index] -= heat[index] * convectivness; //V is y velocity
      }
    }

    // Diffusion and advection
    vel_step(u,v,u_prev,v_prev,viscC,dt);
    dens_step(dens,dens_prev,u,v,diffC,dt);
    dens_step(heat,heat_prev,u,v,diffC_Heat,dt);
    dens_step(oxy,oxy_prev,u,v,diffC_O2,dt);
  }


  // Based on Stam's explicit code
  public Grid(int n){
    this.n = n;
    int size = (n+2) * (n+2);
    u = new float[size];
    v = new float[size];
    u_prev = new float[size];
    v_prev = new float[size];
    dens = new float[size];
    dens_prev = new float[size];

    oxy = new float[size];
    fuel = new float[size];
    heat = new float[size];

    oxy = new float[arraySize];
    oxy_prev = new float[arraySize];
    fuel = new float[arraySize];
    heat = new float[arraySize];
    heat_prev = new float[arraySize];

    energyBarrier = 1.0f;
    rateConstant = 0.1f;
    maxRate = 0.5f;
    exothermicness = 5.0f;
    convectivness = 0.001f;
    viscC = 0.001f;
    diffC = 0.001f;
    diffC_Heat = 0.001f;
    diffC_O2 = 0.001f;

    for (int i = 0; i < arraySize; i++) {
      oxy[i] = 1.0f; 
      fuel[i] = 0.5f;
      heat[i] = 0.0f;
    }
  }

  void addSmoke(int x, int y, float d_u, float d_v, float d_dens) {
    int index = IX(x, y);
    dens_prev[index] += d_dens;
    u_prev[index] += d_u;
    v_prev[index] += d_v;
  }

  int IX(int i, int j) {
    return i + (n + 2) * j;
  }

  void add_source(float[] x, float[] s, float dt) {
    for (int i = 0; i < n * n; i++) {
      x[i] += dt * s[i];
    }
  }

  void diffuse(int b, float[] x, float[] x0, float diff, float dt) {
    int i, j, k;
    float a = dt * diff * n * n;
    for (k = 0; k < 10; k++) {
      for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
          x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                          x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
        }
      }
      set_bnd(n, b, x);
    }
  }

  void advect(int b, float[] d, float[] d0, float[] u, float[] v, float dt) {
    // b: a boundary condition identifier
    // d: density array
    // d0: Previous density array
    // u: velocity array in x
    // v: velocity array in y
    // dt: time step

    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt * n;
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
        x = i - dt0 * u[IX(i, j)];
        y = j - dt0 * v[IX(i, j)];
        if (x < 0.5) x = 0.5; if (x > n + 0.5) x = n + 0.5; i0 = (int) x; i1 = i0 + 1;
        if (y < 0.5) y = 0.5; if (y > n + 0.5) y = n + 0.5; j0 = (int) y; j1 = j0 + 1;
        s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1; // bilinear intepolation
        d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                        s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
      }
    }

    set_bnd(n, b, d);
  }

  void dens_step(float[] x, float[] x0, float[] u, float[] v, float diff, float dt) {
    add_source(x, x0, dt);

    // swap(x0, x); 
    temp_for_swap = x;
    x = x0;
    x0 = temp_for_swap;

    diffuse(0, x, x0, diff, dt);

    // swap(x0, x); 
    temp_for_swap = x;
    x = x0;
    x0 = temp_for_swap;

    advect(0, x, x0, u, v, dt);
  }

  void vel_step(float[] u, float[] v, float[] u0, float[] v0, float visc, float dt) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);

    // swap(u0, u);
    temp_for_swap = u;
    u = u0;
    u0 = temp_for_swap;

    diffuse(1, u, u0, visc, dt);

    // swap(u0, u);
    temp_for_swap = v;
    v = v0;
    v0 = temp_for_swap;

    diffuse(2, v, v0, visc, dt);

    project(u, v, u0, v0);

    // swap(u0, u);
    temp_for_swap = u;
    u = u0;
    u0 = temp_for_swap;

    // swap(v0, v);
    temp_for_swap = v;
    v = v0;
    v0 = temp_for_swap;

    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);
    project(u, v, u0, v0);
  }

  void project(float[] u, float[] v, float[] p, float[] div) {
    // extract the gradient field from velocity field to obtain incompressible field
    int i, j, k;
    float h = 1.0 / n;
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
        div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] +
                                      v[IX(i, j + 1)] - v[IX(i, j - 1)]);
        p[IX(i, j)] = 0;
      }
    }

    set_bnd(n, 0, div);
    set_bnd(n, 0, p);

    for (k = 0; k < 20; k++) {
      for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
          p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                            p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
        }
      }
      set_bnd(n, 0, p);
    }
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
        u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
        v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
      }
    }
    set_bnd(n, 1, u);
    set_bnd(n, 2, v);
  }

  void set_bnd(int N, int b, float[] x) {
    int i;
    for (i = 1; i <= N; i++) {
      x[this.IX(0, i)] = b == 1 ? -x[this.IX(1, i)] : x[this.IX(1, i)];
      x[this.IX(N + 1, i)] = b == 1 ? -x[this.IX(N, i)] : x[this.IX(N, i)];
      x[this.IX(i, 0)] = b == 2 ? -x[this.IX(i, 1)] : x[this.IX(i, 1)];
      x[this.IX(i, N + 1)] = b == 2 ? -x[this.IX(i, N)] : x[this.IX(i, N)];
    }
    x[this.IX(0, 0)] = 0.5 * (x[this.IX(1, 0)] + x[this.IX(0, 1)]);
    x[this.IX(0, N + 1)] = 0.5 * (x[this.IX(1, N + 1)] + x[this.IX(0, N)]);
    x[this.IX(N + 1, 0)] = 0.5 * (x[this.IX(N, 0)] + x[this.IX(N + 1, 1)]);
    x[this.IX(N + 1, N + 1)] = 0.5 * (x[this.IX(N, N + 1)] + x[this.IX(N + 1, N)]);
  }
}

public class CigaretteTip {
  public float x, y;
  public float u, v;
  public float density;
  private Random random;

  public CigaretteTip(float x, float y) {
    this.x = x;
    this.y = y;
    this.u = -10;
    this.v = -2;
    this.density = 50;
    this.random = new Random();
  }

  public void updatePosition(float x, float y) {
    this.x = x;
    this.y = y;
  }

  public void randomize() {
    this.u = (random.nextFloat() - 0.5f) * 2 * 4;
    this.v = random.nextFloat() * -4;
  }
}

public class Cigarette{
  public Vec2 pos;
  public float burningRate;
  public CigaretteTip cTip;
  public ParticleSystem ps;
  public float frontPos;

  public float totalLength = 550;
  public float frontLength = 400;
  public float tipLength = 150;

  public Cigarette(Vec2 pos) {
    this.pos = pos;
    this.cTip = new CigaretteTip(pos.x - frontLength, pos.y);
    frontPos = 0;
    burningRate = 0.2;
    ps = new ParticleSystem(new Vec2(pos.x - frontLength -10, pos.y));

  }

  void update(){

    if (frontPos <= frontLength){
      frontPos += burningRate;
      cTip.x += burningRate;
      ps.origin.x += burningRate;
      cTip.density += 0.05;
    }else if(frontPos > frontLength){
      cTip.density = 0;
      ps.clear();
      player.kill();
    }
    ps.addParticle(10, 30);
    if(frontPos > 50){
      ac.start();
    }
  }

  void updatePosition(float x, float y){
    pos = new Vec2(x, y);
    cTip.updatePosition(x - frontLength + frontPos, y);
    ps.updatePosition(x - frontLength + frontPos  -15, y);
  }


  void draw(){
    float top = pos.y - 30;
    float bottom = pos.y + 30;
    fill(255);
    rect(pos.x - frontLength + frontPos, top, frontLength - frontPos, 60);
    fill(255, 0, 0);
    ellipse(pos.x - frontLength + frontPos, pos.y, 20, 60);

    fill(139, 69, 19);
    rect(pos.x, top, tipLength, 60);
    ellipse(pos.x, pos.y, 20, 60);
    ellipse(pos.x + tipLength, pos.y, 20, 60);

    ps.run(pos.x - frontLength + frontPos -15, pos.x - frontLength + frontPos+5);
  } 
}


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


public class Particle{
  Vec2 pos;
  Vec2 vel;
  float lifespan;
  color particleColor;

  public Particle(Vec2 p){
    pos = p;
    vel = new Vec2(random(0, 2), random(-1, 1));
    lifespan = 255.0;
    particleColor = color(128);
  }

  public void update() {
    pos.add(vel);
    lifespan -= 10.0 * dyingRate;
  }

  public void updateColor(float leftBound, float rightBound) {
    float t = map(pos.x, leftBound, rightBound, 0, 1);
    t = constrain(t, 0, 1);

    particleColor = lerpColor(color(50), color(255, 0, 0), t);
  }


  public void draw() {
    fill(particleColor, lifespan);
    noStroke();
    circle(pos.x, pos.y, 2.5);
  }

  public boolean isDead(){
    return lifespan < 0;
  }
}


class ParticleSystem {
  ArrayList<Particle> particles;
  Vec2 origin;

  ParticleSystem(Vec2 pos) {
    origin = pos;
    particles = new ArrayList<Particle>();
  }

  void addParticle(float width, float height) {
    for(int i = 0; i< 40; i++){
      float angle = random(TWO_PI);
      float r = sqrt(random(1));
      float x = r * cos(angle) * width / 2;
      float y = r * sin(angle) * height / 2;

      Vec2 randomPos = new Vec2(origin.x + x, origin.y + y);
      particles.add(new Particle(randomPos));
    }
  }

  void run() {
    Iterator<Particle> it = particles.iterator();
    while (it.hasNext()) {
      Particle p = it.next();
      p.update();
      p.draw();
      if (p.isDead()) {
        it.remove();
      }
    }
  }

  void run(float leftBound, float rightBound) {
    Iterator<Particle> it = particles.iterator();
    while (it.hasNext()) {
      Particle p = it.next();
      p.update();
      p.updateColor(leftBound, rightBound);
      p.draw();
      if (p.isDead()) {
        it.remove();
      }
    }
  }

  void updatePosition(float x, float y){
    origin = new Vec2(x, y);
  }
  
  void clear(){
    particles = new ArrayList<Particle>();
  }
}
