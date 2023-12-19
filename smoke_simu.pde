// void setup() {
//   size(1280, 1280);
//   grid = new Grid(N);
//   cTip = new CigaretteTip(475, 768);
//   cig = new Cigarette(new Vec2(1000, 768), cTip);
//   minim = new Minim(this);
//   player = minim.loadFile("cough.wav");
//   paused = true;
//   for (int i = 0; i < arraySize; i++) {
//     antiGravity[i] = -30;
//     chaoticEffects[i] = random(-20, 20);
//   }
// }

// void draw() {
//  background(0);
//  if(!paused){
//     addForcesFromCigarette();
//     addVirticalForces();
//     grid.vel_step(grid.u, grid.v, grid.u_prev, grid.v_prev, visc, dt);
//     grid.dens_step(grid.dens, grid.dens_prev, grid.u, grid.v, diff, dt);
//     renderDensity();
//     cig.draw();
//  }

//   // renderDensity();
//   // cig.draw();
// }


// void addForcesFromCigarette() {
//  int x = int(cTip.x * N / width) - 2;
//  int y = int(cTip.y * N / height) - 5;
//  for (int i = 0; i < 2; i++){
//    for (int j = 0; j < 10; j++){
//     grid.addSmoke(x+i, y+j, cTip.u, cTip.v, cTip.density);
//    }
//  }
//  cTip.randomize();
//  updateWind();
// }

// void addVirticalForces(){
//  grid.add_source(grid.v_prev, antiGravity, dt);
//  grid.add_source(grid.u_prev, chaoticEffects, dt);
// }

// void updateWind() {
//     // Shift wind values
//   for (int i = 0; i < arraySize; i++) {
//     chaoticEffects[i] = random(-20, 20);
//   }

//     // Optional: Visualize the wind forces (for debugging)
//     // for (int i = 0; i < WIND_ARRAY_SIZE; i++) {
//     //     line(i * (width / WIND_ARRAY_SIZE), height / 2, i * (width / WIND_ARRAY_SIZE), height / 2 + windForces[i]);
//     // }
// }



// // void addForces() {
// //  // Add forces to the grid based on user input or other criteria
// //  // For example, add density or velocity when the mouse is pressed
// //  if (mousePressed) {
// //    int x = mouseX * N / width;
// //    int y = mouseY * N / height;
// //    int index = grid.IX(x, y);
// //    grid.dens_prev[index] += 100; // Add density
// //    grid.u_prev[index] += (mouseX - pmouseX) * 5; // Add velocity
// //    println("mouseX: ", mouseX, "pmouseX:", pmouseX, "u_dif: ", mouseX - pmouseX);
// //    grid.v_prev[index] += (mouseY - pmouseY) * 5;
// //    println("mouseY: ", mouseY, "pmouseY:", pmouseY, "v_dif: ", mouseY - pmouseY);
// //  }
// // }


// void renderDensity() {
//    float cellWidth = (float) width / N;
//    float cellHeight = (float) height / N;
//    noStroke();
//    for (int i = 0; i < N; i++) {
//      for (int j = 0; j < N; j++) {
//        int x = i * width / N;
//        int y = j * height / N;
//        int index = grid.IX(i, j);
//        float density = grid.dens[index];
//        float colorValue = 0 + constrain(density * 40, 0, 255);
//        // colorValue = colorValue == 255 ? 0 : colorValue;
//       //  fill(255, 255, 255, density);
//        fill(colorValue);
//        rect(x, y, cellWidth+0.731, cellHeight+0.731);
//      }
//    }
// }

// void keyPressed()
// {
//  if (key == ' ') {
//   paused = !paused;
//  }else if(key == 'd' || key == 'D'){
//   println("Increase diffusion Rate!");
//   diff *= 10;
//  }else if(key == 'f' || key == 'F'){
//   println("Decrease diffusion Rate!");
//   diff /= 10;
//  }else if(key == 'v' || key == 'V'){
//   println("Increase Viscosity Rate!");
//   visc *= 10;
//  }else if(key == 'b' || key == 'B'){
//   println("Decrease viscosity Rate!");
//   visc /= 10;
//  }
// }

// void mouseDragged(){
//   cig.pos = new Vec2(mouseX, mouseY);
//   Vec2 cTip_pos = cig.pos.plus(new Vec2(-525, 0));
//   cTip.updatePosition(cTip_pos.x + cig.frontPos, cTip_pos.y);
// }

