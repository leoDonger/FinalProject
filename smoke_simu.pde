// void setup() {
//   size(1280, 1280);
//   paused = true;
//   grid = new Grid(N);
//   fire = new Grid(128);
//   cig = new Cigarette(new Vec2(1000, 768));

//   ac = new AudioContext();
//   try {
//     Sample sample = new Sample("cough.wav");
//     player = new SamplePlayer(ac, sample);
//     player.setLoopType(SamplePlayer.LoopType.LOOP_FORWARDS);

//     rateUGen = new Static(ac, 1.0f);
//     player.setRate(rateUGen);

//     masterGain = new Gain(ac, 2, 0.5f);
//     masterGain.addInput(player);
//     ac.out.addInput(masterGain);

//     } catch (Exception e) {
//     println("Error loading audio file: " + e.getMessage());
//   }

//   for (int i = 0; i < arraySize; i++) {
//     antiGravity[i] = -50;
    
//     // chaoticEffects[i] = random(-1, 1) * 40 + 20;
//     chaoticEffects[i] = (random(1) < 0.5) ? -50 : 50;
//   }
// }

// void updateSmoke(){
//   addForcesFromCigarette();
//   addVirticalForces();
//   grid.vel_step(grid.u, grid.v, grid.u_prev, grid.v_prev, visc, dt);
//   grid.dens_step(grid.dens, grid.dens_prev, grid.u, grid.v, diff, dt);

//   cig.update();
// }


// // void updateFire(float dt){
// // // Compute the "burn equation" for each cell (O,F -> H)
// //  for (int i = 0; i < nx; i++){
// //    for (int j = 0; j < ny; j++){
// //      int index = IX_FIRE(i, j);
// //      float O = oxy[index];
// //      float F = fuel[index];
// //      float H = heat[index] + 30;  //Ambient temp
          
// //      float reactionRate = (O*F*H - energyBarrier) * rateConstant;
          
// //      if (reactionRate < 0) reactionRate = 0;
// //      if (reactionRate > maxRate) reactionRate = maxRate;

// //      oxy[index] -= reactionRate * dt;
// //      fuel[index] -= reactionRate * dt;
// //      heat[index] += reactionRate * dt * exothermicness;
              
// //      if (oxy[index] < 0) oxy[index] = 0;
// //      if (fuel[index] < 0) fuel[index] = 0;
// //    }
// //  }

// //  // Add an updraft for convection due to heat
// //  for (int i=1 ; i<=nx ; i++){
// //    for (int j=1 ; j<=ny ; j++){
// //      V[index] += heat[index] * convectivness; //V is y velocity
// //    }
// //  }

// //  // Diffusion and advection
// //  vel_step(U,V,U_prev,V_prev,viscC,dt);
// //  dens_step(dens,dens_prev,U,V,diffC,dt);
// //  dens_step(heat,heat_prev,U,V,diffC_Heat,dt);
// //  dens_step(oxy,oxy_prev,U,V,diffC_O2,dt);
// // }


// void draw() {
//  background(0);
//  if(!paused){
//   updateSmoke();
//   updateFire();
//  }

//   renderDensity();
//   cig.draw();
// }


// void addForcesFromCigarette() {
//  int x = int(cig.cTip.x * N / width) - 2;
//  int y = int(cig.cTip.y * N / height) - 5;
//  for (int i = 0; i < 2; i++){
//    for (int j = 0; j < 10; j++){
//     grid.addSmoke(x+i, y+j, cig.cTip.u, cig.cTip.v, cig.cTip.density);
//    }
//  }
//  cig.cTip.randomize();
//  updateWind();
// }

// void addVirticalForces(){
//  grid.add_source(grid.v_prev, antiGravity, dt);
//  grid.add_source(grid.u_prev, chaoticEffects, dt);
// }

// void updateWind() {
//   for (int i = 0; i < arraySize; i++) {
//     chaoticEffects[i] = (random(1) < 0.5) ? -50 : 50;
//   }
// }

// void renderDensity() {
//    float cellWidth = (float) width / N;
//    float cellHeight = (float) height / N;
//    noStroke();
//    for (int i = 0; i < N; i++) {
//      for (int j = 0; j < N; j++) {
//        int x = i * width / N;
//        int y = j * height / N;
//        int index = IX(i, j);
//        float density = grid.dens[index];
//        float colorValue = 0 + constrain(density * 60, 0, 255);
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
//  }else if(key == 'n' || key == 'N'){
//   println("Speed up the burning rate of the cigarette!");
//   cig.burningRate *= 2;
//   rateUGen.setValue(rateUGen.getValue() + 0.2f);
//   dyingRate *= 2;
//  }else if(key == 'm' || key == 'M'){
//   println("Slow down the burning rate of the cigarette!");
//   cig.burningRate *= 0.5;
//   rateUGen.setValue(rateUGen.getValue() - 0.2f);
//   dyingRate *= 0.5;
//  }
// }

// void mouseDragged(){
//   cig.updatePosition(mouseX, mouseY);
// }
