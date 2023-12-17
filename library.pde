import java.util.Random;

Grid grid;
Grid3D grid3d;

CigaretteTip cTip;
int N = 128; // Grid size
float dt = 1.0 / 120; // Time step
// float dt = 1.0 / frameRate; // Time step
float diff = 0.001; // Diffusion rate
float visc = 0.001; // Viscosity
boolean paused;

// int U_FIELD = 0;
// int V_FIELD = 1;
// int W_FIELD = 2;
// int D_FIELD = 3;
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
  private float[] age;
  private float max_den = 45.0;

  public Grid(int n){
    this.n = n;
    int size = (n+2) * (n+2);
    u = new float[size];
    v = new float[size];
    u_prev = new float[size];
    v_prev = new float[size];
    dens = new float[size];
    dens_prev = new float[size];
    age = new float[size];
  }

  void updateSmokeAge(float threshold) {
    for (int i = 0; i < age.length; i++) {
        if (dens[i] > 0) {
            age[i]++;  // Increment the age if there is smoke
        } else {
            age[i] = 0;  // Reset age if no smoke
        }

        if (age[i] > threshold) {
            dens[i] = 0;
            age[i] = 0;
        }
    }
  }

  void removeOldSmoke(float threshold) {
    for (int i = 0; i < age.length; i++) {
        if (age[i] > threshold) {
            dens[i] = 0;
            age[i] = 0;
        }
    }
  }

  // void integrate(float dt, float force){
  //   for (i = 1; i <= n; i++) {
  //     for (j = 1; j <= n; j++) {
        
  //       x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
  //                       x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
  //     }
  //   }
  // }


  // void addSmoke(int x, int y, float d_dens) {
  //   int index = IX(x, y);
  //   dens_prev[index] += d_dens;
  //   // u_prev[index] += d_u;
  //   // v_prev[index] += d_v;
  //   // age[index] = 0;  // Reset age for new smoke
  // }


  void addSmoke(int x, int y, float d_u, float d_v, float d_dens) {
    int index = IX(x, y);
    dens_prev[index] += d_dens;
    u_prev[index] += d_u;
    v_prev[index] += d_v;
    // age[index] = 0;  // Reset age for new smoke
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

    set_bnd(N, b, d);
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

    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    for (k = 0; k < 20; k++) {
      for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
          p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                            p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
        }
      }
      set_bnd(N, 0, p);
    }
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
        u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
        v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
      }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
  }



  void set_bnd(int N, int b, float[] x) {
    int i;
    for (i = 1; i <= N; i++) {
      x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
      x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
      x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
      x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
  }


}

public class CigaretteTip {
    public int x, y;
    public float u, v;
    public float density;
    private Random random;

    public CigaretteTip(int x, int y) {
        this.x = x;
        this.y = y;
        this.u = -10;
        this.v = -2;
        this.density = 100;
        this.random = new Random();
    }

    public void updatePosition(int x, int y) {
      this.x = x;
      this.y = y;
    }

    public void randomize() {
      this.u = (random.nextFloat() - 0.5f) * 2 * 4;
      // this.v += (random.nextFloat() - 0.5f) * 2 * 10;
      this.v = random.nextFloat() * -4;
      // this.density += (random.nextFloat() - 0.5f) * 2 * 100;
    }
}




float g = 9.81;



public class Grid3D {
  public int n;
  public float[] u, v, w, u_prev, v_prev, w_prev;
  public float[] dens;
  public float[] dens_prev;
  private float[] temp_for_swap;

  public Grid3D(int n) {
    this.n = n;
    int size = (n + 2) * (n + 2) * (n + 2);
    u = new float[size];
    v = new float[size];
    w = new float[size];
    u_prev = new float[size];
    v_prev = new float[size];
    w_prev = new float[size];
    dens = new float[size];
    dens_prev = new float[size];
  }

  int IX(int i, int j, int k) {
      return i + (n + 2) * (j + (n + 2) * k);
  }

  void add_source(float[] x, float[] source, float dt){
    for (int i = 0; i < n * n * n; i++){
      x[i] += dt * source[i];
    }
  }

  void diffuse(int b, float[] x, float[] x0, float diff, float dt) {
    float a = dt * diff * n * n * n; 

    for (int iter = 0; iter < 20; iter++) {
      for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
          for (int k = 1; k <= n; k++) {
            x[IX(i, j, k)] = (x0[IX(i, j, k)] +
                            a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] +
                                  x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] +
                                  x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / (1 + 6 * a);
          }
        }
      }
      set_bnd(b, x);
    }
  }


// void advect(int b, float[] d, float[] d0, float[] u, float[] v, float[] w, float dt) {
//     int i, j, k, i0, j0, k0, i1, j1, k1;
//     float x, y, z, s0, t0, r0, s1, t1, r1, dt0;

//     dt0 = dt * n;
//     for (i = 1; i <= n; i++) {
//       for (j = 1; j <= n; j++) {
//         for (k = 1; k <= n; k++) {
//           x = i - dt0 * u[IX(i, j, k)];
//           y = j - dt0 * v[IX(i, j, k)];
//           z = k - dt0 * w[IX(i, j, k)];

//           if (x < 0.5) x = 0.5;
//           if (x > n + 0.5) x = n + 0.5;
//           i0 = (int) x;
//           i1 = i0 + 1;

//           if (y < 0.5) y = 0.5;
//           if (y > n + 0.5) y = n + 0.5;
//           j0 = (int) y;
//           j1 = j0 + 1;

//           if (z < 0.5) z = 0.5;
//           if (z > n + 0.5) z = n + 0.5;
//           k0 = (int) z;
//           k1 = k0 + 1;

//           s1 = x - i0; s0 = 1 - s1;
//           t1 = y - j0; t0 = 1 - t1;
//           r1 = z - k0; r0 = 1 - r1;

//           d[IX(i, j, k)] =
//               s0 * (t0 * (r0 * d0[IX(i0, j0, k0)] + r1 * d0[IX(i0, j0, k1)]) +
//                     t1 * (r0 * d0[IX(i0, j1, k0)] + r1 * d0[IX(i0, j1, k1)])) +
//               s1 * (t0 * (r0 * d0[IX(i1, j0, k0)] + r1 * d0[IX(i1, j0, k1)]) +
//                     t1 * (r0 * d0[IX(i1, j1, k0)] + r1 * d0[IX(i1, j1, k1)]));
//         }
//       }
//     }
//     set_bnd(b, d);
//   }

  void advect(int b, float[] d, float[] d0, float[] u, float[] v, float[] w, float dt) {
    int i, j, k, i0, j0, k0, i1, j1, k1;
    float x, y, z, s0, t0, r0, s1, t1, r1, dt0;

    dt0 = dt * n;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            for (k = 1; k <= n; k++) {
                x = i - dt0 * u[IX(i, j, k)];
                y = j - dt0 * v[IX(i, j, k)];
                z = k - dt0 * w[IX(i, j, k)];
                if (x < 0.5) x = 0.5; if (x > n + 0.5) x = n + 0.5; i0 = (int) x; i1 = i0 + 1;
                if (y < 0.5) y = 0.5; if (y > n + 0.5) y = n + 0.5; j0 = (int) y; j1 = j0 + 1;
                if (z < 0.5) z = 0.5; if (z > n + 0.5) z = n + 0.5; k0 = (int) z; k1 = k0 + 1;

                s1 = x - i0; s0 = 1 - s1;
                t1 = y - j0; t0 = 1 - t1;
                r1 = z - k0; r0 = 1 - r1;

                d[IX(i, j, k)] = 
                    s0 * (t0 * (r0 * d0[IX(i0, j0, k0)] + r1 * d0[IX(i0, j0, k1)]) +
                          t1 * (r0 * d0[IX(i0, j1, k0)] + r1 * d0[IX(i0, j1, k1)])) +
                    s1 * (t0 * (r0 * d0[IX(i1, j0, k0)] + r1 * d0[IX(i1, j0, k1)]) +
                          t1 * (r0 * d0[IX(i1, j1, k0)] + r1 * d0[IX(i1, j1, k1)]));
            }
        }
    }
    set_bnd(b, d);
}


void dens_step(float[] x, float[] x0, float[] u, float[] v, float[] w, float diff, float dt) {
    add_source(x, x0, dt);

    // Swap x0 and x for diffusion step
    temp_for_swap = x;
    x = x0;
    x0 = temp_for_swap;

    diffuse(0, x, x0, diff, dt);

    // Swap x0 and x for advection step
    temp_for_swap = x;
    x = x0;
    x0 = temp_for_swap;

    advect(0, x, x0, u, v, w, dt);
}

void vel_step(float[] u, float[] v, float[] w, float[] u0, float[] v0, float[] w0, float visc, float dt) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);
    add_source(w, w0, dt);

    // Swap for diffusion step
    temp_for_swap = u;
    u = u0;
    u0 = temp_for_swap;

    diffuse(1, u, u0, visc, dt);

    temp_for_swap = v;
    v = v0;
    v0 = temp_for_swap;

    diffuse(2, v, v0, visc, dt);

    temp_for_swap = w;
    w = w0;
    w0 = temp_for_swap;

    diffuse(3, w, w0, visc, dt);

    project(u, v, w, u0, v0);

    // Swap for advection step
    temp_for_swap = u;
    u = u0;
    u0 = temp_for_swap;

    temp_for_swap = v;
    v = v0;
    v0 = temp_for_swap;

    temp_for_swap = w;
    w = w0;
    w0 = temp_for_swap;

    advect(1, u, u0, u0, v0, w0, dt);
    advect(2, v, v0, u0, v0, w0, dt);
    advect(3, w, w0, u0, v0, w0, dt);

    project(u, v, w, u0, v0);
}


void project(float[] u, float[] v, float[] w, float[] p, float[] div) {
    int i, j, k;
    float h = 1.0f / n;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            for (k = 1; k <= n; k++) {
                div[IX(i, j, k)] = -0.5f * h * (u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)] +
                                                v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)] +
                                                w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]);
                p[IX(i, j, k)] = 0;
            }
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);

    // Solve for pressure
    for (int iter = 0; iter < 20; iter++) {
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++) {
                for (k = 1; k <= n; k++) {
                    p[IX(i, j, k)] = (div[IX(i, j, k)] + p[IX(i - 1, j, k)] +
                                        p[IX(i + 1, j, k)] + p[IX(i, j - 1, k)] +
                                        p[IX(i, j + 1, k)] + p[IX(i, j, k - 1)] +
                                        p[IX(i, j, k + 1)]) / 6.0f;
                }
            }
        }
        set_bnd(0, p);
    }

    // Subtract pressure gradient
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            for (k = 1; k <= n; k++) {
                u[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]) / h;
                v[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]) / h;
                w[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]) / h;
            }
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
    set_bnd(3, w);
}

  void set_bnd(int b, float[] x) {
      int i, j, k;
      for (i = 1; i <= n; i++) {
          for (j = 1; j <= n; j++) {
              for (k = 1; k <= n; k++) {
                  int ijk = IX(i, j, k);
                  // Set boundaries for 3D case
              }
          }
      }
      // Handle corners and edges in the 3D case
  }

  // Add the third dimension to all other methods: diffuse, advect, project, etc.

  // You will also need to extend the logic of each method to consider the 3D neighbors.
}









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
