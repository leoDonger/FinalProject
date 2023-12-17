void setup() {
  background(0);
  size(1280, 1280);
  cTip = new CigaretteTip(640, 768);
  grid = new Grid(N);

  for (int i = 0; i < arraySize; i++) {
    antiGravity[i] = -40;
    chaoticEffects[i] = 0;
  }
}

void draw() {
  // background(0);
  if(!paused){
    addForcesFromCigarette();
    addVirticalForces();
    grid.vel_step(grid.u, grid.v, grid.u_prev, grid.v_prev, visc, dt);
    grid.dens_step(grid.dens, grid.dens_prev, grid.u, grid.v, diff, dt);
  }

  renderDensity();
  // renderDensityInterpolated();
}


void addForcesFromCigarette() {
  int x = (cTip.x * N / width) - 5;
  int y = (cTip.y * N / height);
  for (int i = 0; i < 10; i++){
    grid.addSmoke(x+i, y, cTip.u, cTip.v, cTip.density);
    // grid.addSmoke(x+i, y, cTip.density);
  }
  // grid.addSmoke(x, y, cTip.u * 5, cTip.v * 5, cTip.density);
  cTip.randomize();
}

void addVirticalForces(){
  grid.add_source(grid.v_prev, antiGravity, dt);
}



void addForces() {
  // Add forces to the grid based on user input or other criteria
  // For example, add density or velocity when the mouse is pressed
  if (mousePressed) {
    int x = mouseX * N / width;
    int y = mouseY * N / height;
    int index = grid.IX(x, y);
    grid.dens_prev[index] += 100; // Add density
    grid.u_prev[index] += (mouseX - pmouseX) * 5; // Add velocity
    println("mouseX: ", mouseX, "pmouseX:", pmouseX, "u_dif: ", mouseX - pmouseX);
    grid.v_prev[index] += (mouseY - pmouseY) * 5;
    println("mouseY: ", mouseY, "pmouseY:", pmouseY, "v_dif: ", mouseY - pmouseY);
  }
}


// void renderDensityInterpolated() {
//     float cellWidth = (float) width / N;
//     float cellHeight = (float) height / N;
//     color[][] colors = new color[N][N]; // Array to store colors for interpolation

//     // First, calculate and store colors for each cell
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             int index = grid.IX(i, j);
//             float d = grid.dens[index];
//             colors[i][j] = color(255 - constrain(d * 5, 0, 255));
//         }
//     }

//     // Render using interpolated colors
//   noStroke();
//   for (int i = 0; i < N - 1; i++) {
//     for (int j = 0; j < N - 1; j++) {
//       float x = i * cellWidth;
//       float y = j * cellHeight;

//       // Use lerpColor for interpolation
//       color c1 = lerpColor(colors[i][j], colors[i+1][j], 0.5);
//       color c2 = lerpColor(colors[i][j+1], colors[i+1][j+1], 0.5);
//       color avgColor = lerpColor(c1, c2, 0.5);

//       fill(avgColor);
//       rect(x, y, cellWidth+0.7, cellHeight+0.7);
//     }
//   }
// }


void renderDensity() {
    float cellWidth = (float) width / N;
    float cellHeight = (float) height / N;
    noStroke();
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        int x = i * width / N;
        int y = j * height / N;
        int index = grid.IX(i, j);
        float d = grid.dens[index];
        // float colorValue = 255 - constrain(d * 10, 0, 255);
        float colorValue = 0 + constrain(d * 5, 0, 255);
        // colorValue = colorValue == 255 ? 0 : colorValue;
        fill(colorValue);
        rect(x, y, cellWidth+0.731, cellHeight+0.731);
      }
    }
}

void keyPressed()
{
  if (key == ' ') {
    paused = !paused;
  }
}
