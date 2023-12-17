void setup() {
 background(0);
 size(1280, 1280);
 cTip = new CigaretteTip(640, 768);
 grid = new Grid(N);

 for (int i = 0; i < arraySize; i++) {
   antiGravity[i] = -60;
   chaoticEffects[i] = 0;
 }
}

void draw() {
 if(!paused){
   addForcesFromCigarette();
   addVirticalForces();
   grid.vel_step(grid.u, grid.v, grid.u_prev, grid.v_prev, visc, dt);
   grid.dens_step(grid.dens, grid.dens_prev, grid.u, grid.v, diff, dt);
 }

 renderDensity();
}


void addForcesFromCigarette() {
 int x = (cTip.x * N / width) - 5;
 int y = (cTip.y * N / height);
 for (int i = 0; i < 10; i++){
   grid.addSmoke(x+i, y, cTip.u, cTip.v, cTip.density);
 }
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


void renderDensity() {
   float cellWidth = (float) width / N;
   float cellHeight = (float) height / N;
   noStroke();
   for (int i = 0; i < N; i++) {
     for (int j = 0; j < N; j++) {
       int x = i * width / N;
       int y = j * height / N;
       int index = grid.IX(i, j);
       float density = grid.dens[index];
       float colorValue = 0 + constrain(density * 5, 0, 255);
       // colorValue = colorValue == 255 ? 0 : colorValue;
      //  fill(255, 255, 255, density);
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
