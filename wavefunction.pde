//Bohr radius
float a_0 = 0.1;
//Atomic number
int Z = 1;

void setup(){
  size(800, 800);
  background(30);
  for(int x = 0; x < width; x++){
    fill(255, 120, 140);
    noStroke();
    circle(x, height / 2 - (int)(LaguerrePolynomial(25, 14, (x - width / 2) / 16.0) / 2500.0), 5);
  }
  println(LaguerrePolynomial(4, 5, 5));
}

void draw(){
  background(30);
  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      float r = sqrt(pow(x - width / 2, 2) + pow(y - height / 2, 2)) / 100.0;
      float waveFunction = waveFunction(2, 0, 0, r, 0, 0);
      set(x, y, color(pow(waveFunction, 2) * 255));
    }
  }
}

//https://en.wikipedia.org/wiki/File:Hydrogen_Density_Plots.png
//https://en.wikipedia.org/wiki/Atomic_orbital
//https://en.wikipedia.org/wiki/Laguerre_polynomials
//https://chem.libretexts.org/Courses/University_of_California_Davis/
//    UCD_Chem_107B%3A_Physical_Chemistry_for_Life_Scientists/Chapters/4%3A_Quantum_Theory/4.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom
//https://en.wikipedia.org/wiki/Bohr_radius

//add caching
float waveFunction(int n, int l, int m, float r, float theta, float phi){
  float underSquareRoot = sqrt(pow(2 / (n * a_0), 3) * factorial(n - l - 1) / (2 * n * factorial(n + l)));
  float rho = Z * r / a_0;
  float exponentials = exp(-rho/2) * pow(rho, l);
  float laguerre = LaguerrePolynomial(n - l - 1, 2 * l + 1, rho);
  float harmonics = SphericalHarmonics(theta, phi);
  
  return underSquareRoot * exponentials * laguerre * harmonics;
}

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

float SphericalHarmonics(float theta, float phi){
  return 1;
}

//Caching
HashMap<Integer, Integer> factorialsEvaluated = new HashMap<Integer, Integer>();
int factorial(int n)
{
  var output = factorialsEvaluated.get(n);
  if(output != null){
    return (int)output;
  }
  int product = 1;
  for(int i = 2; i <= n; i++){
    product *= i;
  }
  
  println("No cached data, calculated " + n + "! = " + product);
  factorialsEvaluated.put(n, product);
  return product;
}
