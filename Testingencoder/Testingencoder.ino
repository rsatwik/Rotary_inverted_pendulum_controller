#define sel1 38
#define sel2 39
#define rst_ 40
#define clk 12

int alpha, theta, alphadot, thetadot, prev_alpha=0, prev_theta=0;
byte alpha1, alpha2, theta1, theta2;

void setup()
 {
  pinMode(sel1, OUTPUT);
  pinMode(sel2, OUTPUT);
  pinMode(rst_, OUTPUT);
  pinMode(clk, OUTPUT);

  TCCR1B=0x01;
  analogWrite(clk, 127);

  digitalWrite(rst_, LOW);
  delay(1000);
  digitalWrite(rst_, HIGH);

 alpha = 0;
 DDRC=0b00000000;

  Serial.begin(9600);

}

void loop()
 {
  digitalWrite(sel1, HIGH);
  digitalWrite(sel2, LOW);
  alpha1=PINC;
  theta1=PINA;
  digitalWrite(sel1, LOW);
  digitalWrite(sel2, LOW);
  alpha2=PINC;
  theta2=PINA;
  alpha=word(alpha2, alpha1);
  theta=word(theta2,theta1);
  //Serial.print(alpha);
  //Serial.print("\t");
  Serial.println(theta);
}
