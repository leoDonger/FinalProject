Grid grid;
int N = 64; // Grid size
float dt = 1.0 / frameRate; // Time step
float diff = 0.0001; // Diffusion rate
float visc = 0.0001; // Viscosity
boolean paused;


void setup() {
  size(400, 400);
  grid = new Grid(N);
  // Initialize density or add initial conditions here if needed
}

void draw() {
  background(255);
  if(!paused){
    addForces();

    grid.vel_step(N, grid.u, grid.v, grid.u_prev, grid.v_prev, visc, dt);
    grid.dens_step(N, grid.dens, grid.dens_prev, grid.u, grid.v, diff, dt);
  }

  // Render the density
  renderDensity();
}

void addForces() {
    // Add forces to the grid based on user input or other criteria
    // For example, add density or velocity when the mouse is pressed
    if (mousePressed) {
        int x = mouseX * N / width;
        int y = mouseY * N / height;
        int index = grid.IX(x, y, N);
        grid.dens_prev[index] += 100; // Add density
        grid.u_prev[index] += (mouseX - pmouseX) * 5; // Add velocity
        grid.v_prev[index] += (mouseY - pmouseY) * 5;
    }
}

// void renderDensity() {
//     loadPixels();
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             int x = i * width / N;
//             int y = j * height / N;
//             int index = grid.IX(i, j, N);
//             float d = grid.dens[index];
//             float colorValue = 255 - constrain(d * 5, 0, 255); // Multiply by a factor to increase contrast
//             pixels[x + y * width] = color(colorValue);
//         }
//     }
//     updatePixels();
// }

void renderDensity() {
    noStroke(); // Disable drawing outlines for the rectangles
    float cellWidth = (float) width / N;
    float cellHeight = (float) height / N;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int x = i * width / N;
            int y = j * height / N;
            int index = grid.IX(i, j, N);
            float d = grid.dens[index];
            float colorValue = 255 - constrain(d * 5, 0, 255); // Adjusted for contrast as before
            fill(colorValue);
            rect(x, y, cellWidth, cellHeight);
        }
    }
}


void keyPressed()
{
  if (key == ' ') {
    paused = !paused;
  }
}
