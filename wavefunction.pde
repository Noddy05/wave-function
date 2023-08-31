//Bohr radius
float a_0 = 0.01;
//Atomic number
int Z = 1;
int N = 4, L = 0, M = 0;

float[][] densities;
float[][] notSquaredDensities;
float totalProbability = 0;
float max = 0;

int location(int n, int l, int m){
  return n * (n - 1) * (2 * n - 1) / 6 + l * (l + 1) + m;
}

void setup(){
  size(800, 800, P2D);
  background(30);
  
  int maxN = 1;
  int totalWeights = (maxN * (maxN + 1) * (2 * maxN + 1)) / 6;
  String combined = "float weights[" + totalWeights + "] = float[" + totalWeights + "](";
  int index = 0;
  for(int n = 1; n <= maxN; n++){
    for(int l = 0; l < n; l++){
      for(int m = -l; m <= l; m++){
        index++;
        densities = new float[width][height];
        totalProbability = 0;
        for(int x = 0; x < width; x++){
          densities[x] = new float[height];
          for(int y = 0; y < height; y++){
            float plotX = x - width / 2;
            float plotY = y - height / 2;
            float r = sqrt(pow(plotX, 2) + pow(plotY, 2)) / 500.0;
            float theta = atan2(y - height / 2, x - width / 2) + HALF_PI;
            float phi = frameCount / 50.0;
            float waveFunction = waveFunction(n, l, m, r, theta, phi);
            float probabilityDensity = pow(waveFunction, 2);
            if(probabilityDensity != probabilityDensity)
              continue;
            
            max = max(max, min(probabilityDensity, 100));
            totalProbability += min(probabilityDensity, 100);
            densities[x][y] = probabilityDensity;
          }
        }
        if(n != 1)
          combined += ", ";
        combined += str(min(round(totalProbability * 100) / 100.0, 500000));
        //println("Max: ", max);
        //println("For n: ", n, ", l: ", l, " and m: ", m, " the total probability is ", totalProbability);
        //println("For index: ", index, " the location is ", location(n, l, m));
      }
    }
  }
  combined += ");";
  //println(combined);
  notSquaredDensities = new float[width][height];
  for(int x = 0; x < width; x++){
    notSquaredDensities[x] = new float[height];
    densities[x] = new float[height];
    for(int y = 0; y < height; y++){
      float plotX = x - width / 2;
      float plotY = y - height / 2;
      float r = sqrt(pow(plotX, 2) + pow(plotY, 2)) / 500.0;
      float theta = atan2(plotX, plotY);
      float phi = frameCount / 50.0;
      float waveFunction = waveFunction(N, L, M, r, theta, phi);
      float probabilityDensity = pow(waveFunction, 2);
      
      if(probabilityDensity != probabilityDensity)
        continue;
      
      max = max(max, min(probabilityDensity, 100));
      totalProbability += min(probabilityDensity, 100);
      densities[x][y] = probabilityDensity;
      notSquaredDensities[x][y] = waveFunction;
      if(SphericalHarmonics(L, M, theta, phi).x < 0)
        notSquaredDensities[x][y] *= -1;
    }
  }
  //println("total:" + totalProbability);
  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      densities[x][y] /= totalProbability;
    }
  }
  totalProbability = 1;
}

void draw(){
  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      if(notSquaredDensities[x][y] < 0)
        set(x, y, color(0, 0, densities[x][y] * 255 * 5000));
      else
        set(x, y, color(densities[x][y] * 255 * 5000, 0, 0));
    }
  }
  
  fill(255);
  textSize(64);
  textAlign(CENTER, CENTER);
  text("Ψ     (r, θ, φ)", width / 2, height - 50);
  textSize(32);
  text(N + ", " + L + ", " + M, width / 2 - 85, height - 20);
  strokeWeight(3);
  stroke(255, 255, 0);
  //line(width / 2, 50, mouseX, mouseY);
  PVector start = new PVector(width / 2, 50); 
  PVector line = new PVector(mouseX - width / 2, mouseY - 50);
  float stepSize = 5;
  float distance = 0;
  float absorptionCoefficient = 500;
  noStroke();
  for(int i = 1; i < line.mag() / stepSize; i++){
    PVector pos = start.copy().add(line.copy().normalize().mult(stepSize * i));
    PVector prevPos = start.copy().add(line.copy().normalize().mult(stepSize * (i - 1)));
    distance += (densities[(int)pos.x][(int)pos.y] + densities[(int)prevPos.x][(int)prevPos.y]) / 2.0 * stepSize;
    float localAbsorption = exp(-distance * absorptionCoefficient);
    fill(localAbsorption * 255);
    circle(pos.x, pos.y, stepSize);
  }
  float absorption = exp(-distance * absorptionCoefficient);
  
  // Wave function collapse
  /*
  noStroke();
  fill(255, 0, 0, 20);
  
  Boolean found = false;
  while(!found){
    int x = (int)random(0, width);
    int y = (int)random(0, height);
    float probability = densities[x][y] / totalProbability;
    float rand = random(0, 1);
    if(rand <= probability){
      circle(x, y, 5);
      found = true;
    }
  }*/
}

//https://en.wikipedia.org/wiki/File:Hydrogen_Density_Plots.png
//https://en.wikipedia.org/wiki/Atomic_orbital
//https://chem.libretexts.org/Courses/University_of_California_Davis/
//    UCD_Chem_107B%3A_Physical_Chemistry_for_Life_Scientists/Chapters/4%3A_Quantum_Theory/4.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom
//https://en.wikipedia.org/wiki/Bohr_radius
float waveFunction(int n, int l, int m, float r, float theta, float phi){
  float underSquareRoot = sqrt(pow(2 / (n * a_0), 3) * Factorial(n - l - 1) / (2 * n * Factorial(n + l)));
  float rho = Z * r / a_0;
  float exponentials = exp(-rho/2) * pow(rho, l);
  float laguerre = LaguerrePolynomial(n - l - 1, 2 * l + 1, rho);
  PVector harmonics = SphericalHarmonics(l, m, theta, phi);
  
  return harmonics.mag() * underSquareRoot * exponentials * laguerre;
}

//https://en.wikipedia.org/wiki/Laguerre_polynomials
float LaguerrePolynomial(int n, int alpha, float x){
  if(n == 0)
    return 1;
  
  float L_ksub1 = 1;
  float L_k = 1 + alpha - x;
  float L_kadd1 = 0;
  for(int k = 1; k < n; k++){
    L_kadd1 = ((2 * k + 1 + alpha - x) * L_k - (k + alpha) * L_ksub1)/(k+1);
    L_ksub1 = L_k;
    L_k = L_kadd1;
  }
  
  return L_k;
}

float Choose(int top, int bottom){
  return Factorial(top)/(Factorial(bottom)*Factorial(top-bottom));
}

//https://en.wikipedia.org/wiki/Legendre_polynomials
float LegendrePolynomial(int n, float x){
  float sum = 0;
  for(int k = 0; k <= n; k++){
    sum += pow(Choose(n, k), 2) * pow(x - 1, n - k) * pow(x + 1, k);
  }
  return 1.0 / pow(2, n) * sum;
}

//https://functions.wolfram.com/Polynomials/LegendreP2/02/0001/
PVector AssociatedLegendrePolynomial(int l, int m, float x){
  if(m == 0)
    return new PVector(LegendrePolynomial(l, x), 0);
    
  float sum = 0;
  for(int k = 0; k <= l; k++){
    if(k - m < 0)
      continue;
    sum += Pochhammer(-l, k) * Pochhammer(l + 1, k) / (Factorial(k - m) * Factorial(k))
      * pow((1 - x) / 2, k);
  }
  float imTop = 0;
  float reTop = 0;
  float imBottom = 0;
  float reBottom = 0;
  if(1 + x < 0){
    imTop = pow(-(1 + x), m / 2.0);
  } else{
    reTop = pow(1 + x, m / 2.0);
  }
  if(1 - x < 0){
    imBottom = pow(-(1 - x), m / 2.0);
  } else{
    reBottom = pow(1 - x, m / 2.0);
  }
  return ComplexDivision(new PVector(reTop, imTop), new PVector(reBottom, imBottom)).mult(sum);
}

// a/b
PVector ComplexDivision(PVector a, PVector b){
  return new PVector(a.x * b.x + a.y * b.y, a.y * b.x - a.x * b.y).div(b.x * b.x + b.y * b.y);
}
// a*b
PVector ComplexMultiplication(PVector a, PVector b){
  return new PVector(a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y);
}

//https://mathworld.wolfram.com/PochhammerSymbol.html
float Pochhammer(float x, float n){
  float product = 1;
  
  for(int i = 0; i < n; i++){
    product *= (x + i);
  }
  
  return product;
}

//https://youtu.be/yRJ2Xf2oH5w?t=517

PVector SphericalHarmonics(int l, int m, float theta, float phi){
  float underSquareRoot = sqrt((2.0 * l + 1.0) * Factorial(l - m)/(4.0 * PI * Factorial(l + m)));
  float exponentialReal = cos(m * phi);
  float exponentialImaginary = sin(m * phi);
  PVector polynomial = AssociatedLegendrePolynomial(l, m, cos(theta));

  return ComplexMultiplication(new PVector(underSquareRoot * exponentialReal, underSquareRoot * exponentialImaginary), polynomial);
}

int Factorial(int n)
{
  int product = 1;
  for(int i = 2; i <= n; i++){
    product *= i;
  }
  
  return product;
}
