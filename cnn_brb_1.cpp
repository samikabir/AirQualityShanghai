#include <iostream>
#include <math.h> 
#include <fstream>
#include <string>  

using namespace std;  
  
float a1 = 0;
float a2 = 0; 
 
//float PM_H = 500; //220, 500.4;
float PM_H = 500.4; //Revised PM2.5
//float PM_M = 35.5;       
//float PM_M = 150.5;  
float PM_M = 55.5; //Revised PM2.5  
//float PM_L = 0.0; 
float PM_L = 0.0; //Revised PM2.5  

float CP_H = 1.00; 
float CP_M = 0.50; 
float CP_L = 0.00; //0.00, 0.24;              
 
float RH_H = 0.90;   
float RH_M = 0.75;        
float RH_L = 0.35;//0.55, 0;                

float H_Cl = 0.0;
float M_Cl = 0.0; 
float L_Cl = 0.0;  

float H_Hm = 0.0;
float M_Hm = 0.0;
float L_Hm = 0.0;  

float H1 = 0.0;
float M1 = 0.0;
float L1 = 0.0;  
  
float H2 = 0.0;   
float M2 = 0.0;
float L2 = 0.0;  
    

int numberOfAntAttributes = 2;
float matchingDegree[27];
float relativeWeight = 1.0; 
float totalWeight = 0;
float consequentBeliefDegree[81];  
float updatedConsequentBeliefDegree[81];
float beliefDegreeChangeLevel = 0; 
float activationWeight[27];
float ruleWiseBeliefDegreeSum[27];  
string line;
string cnn_mild;
string cnn_nominal;
string cnn_severe;
int counter = 0;
string temp;
string dew_point; 
string cloud;
string humidity;

float normalized_cnn_severe_degree = 1.0; //transformed value after image prediction
float normalized_cnn_mild_degree = 1.0;
float normalized_cnn_nominal_degree = 1.0;

float image_normalized_cnn_severe_degree = 1.0; //directly from image
float image_normalized_cnn_mild_degree = 1.0;
float image_normalized_cnn_nominal_degree = 1.0; 

float cnn_pm25 = 1.0;
float pm25_revised = 1.0;
float PO = 1; 
float aggregatedBeliefDegreeH = 1;
float aggregatedBeliefDegreeM = 1;
float aggregatedBeliefDegreeL = 1;   
float finalAggregatedBeliefDegreeH = 1.0;
float finalAggregatedBeliefDegreeM = 1.0;
float finalAggregatedBeliefDegreeL = 1.0;
float brbH = 0;
float brbM = 0;
float brbL = 0;
//int aqi = 1;
float aqi = 1.0;
float aqi1 = 1.0;
float aqi2 = 1.0;
float aqi3 = 1.0;
float aqi4 = 1.0;
float aqi5 = 1.0;
float aqi6 = 1.0;     
 
void transformInputCNN(float i);        
void transformInputCloud(float i);
void transformInputRH(float i);

   
void ruleBase()        
{        
    consequentBeliefDegree[0] = 0;   
    consequentBeliefDegree[1] = 0;
    consequentBeliefDegree[2] = 1;  
    consequentBeliefDegree[3] = 0; 
    consequentBeliefDegree[4] = 0.20;//0.20
    consequentBeliefDegree[5] = 0.80;//0.80 
    consequentBeliefDegree[6] = 1; 
    consequentBeliefDegree[7] = 0; 
    consequentBeliefDegree[8] = 0;//
    consequentBeliefDegree[9] = 0;
    consequentBeliefDegree[10] = 0;
    consequentBeliefDegree[11] = 1; //
    consequentBeliefDegree[12] = 0;
    consequentBeliefDegree[13] = 0.20;
    consequentBeliefDegree[14] = 0.80;//
    consequentBeliefDegree[15] = 1;
    consequentBeliefDegree[16] = 0;
    consequentBeliefDegree[17] = 0;//
    consequentBeliefDegree[18] = 1;
    consequentBeliefDegree[19] = 0;
    consequentBeliefDegree[20] = 0;//
    consequentBeliefDegree[21] = 1;
    consequentBeliefDegree[22] = 0;
    consequentBeliefDegree[23] = 0;//
    consequentBeliefDegree[24] = 1;
    consequentBeliefDegree[25] = 0;
    consequentBeliefDegree[26] = 0;//
    consequentBeliefDegree[27] = 0;  
    consequentBeliefDegree[28] = 0;
    consequentBeliefDegree[29] = 1;//
    consequentBeliefDegree[30] = 0;
    consequentBeliefDegree[31] = 0.40;
    consequentBeliefDegree[32] = 0.60;//
    consequentBeliefDegree[33] = 0;
    consequentBeliefDegree[34] = 1;
    consequentBeliefDegree[35] = 0;//
    consequentBeliefDegree[36] = 0;
    consequentBeliefDegree[37] = 0; 
    consequentBeliefDegree[38] = 1;//
    consequentBeliefDegree[39] = 0;
    consequentBeliefDegree[40] = 0;
    consequentBeliefDegree[41] = 1;//
    consequentBeliefDegree[42] = 0;
    consequentBeliefDegree[43] = 1;
    consequentBeliefDegree[44] = 0;//
    consequentBeliefDegree[45] = 0; 
    consequentBeliefDegree[46] = 1;
    consequentBeliefDegree[47] = 0;//
    consequentBeliefDegree[48] = 0;
    consequentBeliefDegree[49] = 1;
    consequentBeliefDegree[50] = 0;//  
    consequentBeliefDegree[51] = 0;
    consequentBeliefDegree[52] = 1;
    consequentBeliefDegree[53] = 0;//
    consequentBeliefDegree[54] = 0;
    consequentBeliefDegree[55] = 0;
    consequentBeliefDegree[56] = 1;//  
    consequentBeliefDegree[57] = 0;
    consequentBeliefDegree[58] = 0;
    consequentBeliefDegree[59] = 1;//
    consequentBeliefDegree[60] = 0;
    consequentBeliefDegree[61] = 0;
    consequentBeliefDegree[62] = 1;//
    consequentBeliefDegree[63] = 0;
    consequentBeliefDegree[64] = 0;
    consequentBeliefDegree[65] = 1; //
    consequentBeliefDegree[66] = 0;
    consequentBeliefDegree[67] = 0;
    consequentBeliefDegree[68] = 1;//
    consequentBeliefDegree[69] = 0;
    consequentBeliefDegree[70] = 0;
    consequentBeliefDegree[71] = 1;//
    consequentBeliefDegree[72] = 0;
    consequentBeliefDegree[73] = 0;
    consequentBeliefDegree[74] = 1;//
    consequentBeliefDegree[75] = 0;
    consequentBeliefDegree[76] = 0;
    consequentBeliefDegree[77] = 1;//
    consequentBeliefDegree[78] = 0;
    consequentBeliefDegree[79] = 0;   
    consequentBeliefDegree[80] = 1; //              
}        
     
void takeInput()
{
     
    cout<<"Insert value for Cloud percentage (between 0 and 1): ";  
    cin>>a1;    
     
    cout<<"Insert value for Relative Humidity percentage (between 0 and 1): ";
    cin>>a2;  
      
    transformInputCloud(a1);
    transformInputRH(a2);       
}  
   
void transformInputCloud(float i)
{
    if (i >= CP_H)
    {
        H_Cl = 1; 
        M_Cl = 0;
        L_Cl = 0; 
    }
    else if (i == CP_M)
    {
        H_Cl = 0;
        M_Cl = 1;
        L_Cl = 0; 
    } 
    else if (i <= CP_L) 
    {
        H_Cl = 0;
        M_Cl = 0;  
        L_Cl = 1; 
    }    
    else if ((i <= CP_H) && (i >= CP_M)) 
    {
        M_Cl = (CP_H-i)/(CP_H-CP_M);
        H_Cl = 1 - M_Cl; 
        L_Cl = 0.0; 
    } 
    else if ((i <= CP_M) && (i >= CP_L))    
    {
        L_Cl = (CP_M-i)/(CP_M-CP_L);
        M_Cl = 1 - L_Cl; 
        H_Cl = 0.0; 
    } 
} 
 
void transformInputRH(float i) //done till this 2020.7.8 Wed end-of-day  
{//RH_H, H_Hm
    if (i >= RH_H)
    {
        H_Hm = 1; 
        M_Hm = 0;
        L_Hm = 0; 
    }
    else if (i == RH_M)
    {
        H_Hm = 0;
        M_Hm = 1;
        L_Hm = 0; 
    }  
    else if (i <= RH_L) 
    {
        H_Hm = 0;
        M_Hm = 0;  
        L_Hm = 1; 
    }    
    else if ((i <= RH_H) && (i >= RH_M)) 
    {
        M_Hm = (RH_H-i)/(RH_H-RH_M);
        H_Hm = 1 - M_Hm; 
        L_Hm = 0.0; 
    } 
    else if ((i <= RH_M) && (i >= RH_L))     
    {
        L_Hm = (RH_M-i)/(RH_M-RH_L);
        M_Hm = 1 - L_Hm;  
        H_Hm = 0.0; 
    } 
}

void transformInputCNN(float i)   
{//normalized_cnn_severe_degree, normalized_cnn_mild_degree, normalized_cnn_nominal_degree 
    if (i >= PM_H)
    {
        normalized_cnn_severe_degree = 1;
        normalized_cnn_mild_degree = 0;
        normalized_cnn_nominal_degree = 0;
    }
    else if (i == PM_M)
    {
        normalized_cnn_severe_degree = 0;
        normalized_cnn_mild_degree = 1;
        normalized_cnn_nominal_degree = 0;
    }
    else if (i <= PM_L)
    {
        normalized_cnn_severe_degree = 0;
        normalized_cnn_mild_degree = 0;
        normalized_cnn_nominal_degree = 1;
    }    
    else if ((i <= PM_H) && (i >= PM_M))
    {
        normalized_cnn_mild_degree = (PM_H-i)/(PM_H-PM_M);
        normalized_cnn_severe_degree = 1 - normalized_cnn_mild_degree;
        normalized_cnn_nominal_degree = 0.0; 
    }
    else if ((i <= PM_M) && (i >= PM_L))
    {
        normalized_cnn_nominal_degree = (PM_M-i)/(PM_M-PM_L);
        normalized_cnn_mild_degree = 1 - normalized_cnn_nominal_degree; 
        normalized_cnn_severe_degree = 0.0; 
    } 
    cout<<" After transformInputCNN(), normalized_cnn_severe_degree: "<< normalized_cnn_severe_degree << ", normalized_cnn_mild_degree: " << normalized_cnn_mild_degree << ", normalized_cnn_nominal_degree: " << normalized_cnn_nominal_degree <<endl;          
} 

void showTransformedInput() 
{//H_Hm, H_Cl  
    cout<< endl << "Transformed Inputs are as follow." << endl;  
    cout<< "Revised PM2.5 due to fog = {(H, " << normalized_cnn_severe_degree << "); (M, " << normalized_cnn_mild_degree << "); (L, " << normalized_cnn_nominal_degree << ")}" << endl;
    cout<< "Cloud % = {(H, " << H_Cl << "); (M, " << M_Cl << "); (L, " << L_Cl << ")}" << endl;
    cout<< "RH % = {(H, " << H_Hm << "); (M, " << M_Hm << "); (L, " << L_Hm << ")}" << endl;   
} 
  
void calculateMatchingDegreeBrbCnn()    
{ 
    int increment = 0;
    float ti1[3] = {normalized_cnn_severe_degree, normalized_cnn_mild_degree, normalized_cnn_nominal_degree}; 
    float ti2[3] = {H_Cl, M_Cl, L_Cl};     
    float ti3[3] = {H_Hm, M_Hm, L_Hm};            
    
    
    for (int c = 0; c < 3; c++)
        for (int d = 0; d < 3; d++)
        	for (int e = 0; e < 3; e++){ 
                //weight[increment] = ti1[c] * ti2[d] * ti3[e];  
                matchingDegree[increment] = pow(ti1[c], relativeWeight) * pow(ti2[d], relativeWeight) * pow(ti3[e], relativeWeight);
                increment++; 
            } 
} 
  
void showMatchingDegree()
{
    int track = 1; 
    //cout << endl << "Matching degrees of the rules are as follow." << endl; 
    for (int counter = 0; counter < 27; counter++)
    { 
        //cout<< "Matching Degree of Rule " << track << " = " << matchingDegree[counter] << endl;
        track++; 
    } 
} 

void showActivationWeight()
{   
    int trace = 1; 
    for (int x = 0; x < 27; x++)
    {
        totalWeight += matchingDegree[x];
    } 
    
    //cout << endl << "Activation Weights of the rules are as follow."<< endl; 
     
    for (int counter = 0; counter < 27; counter++)
    {  
        activationWeight[counter] = matchingDegree[counter]/totalWeight; 
        cout<< "Activation weight of Rule " << trace << " = " << activationWeight[counter] << endl;
        trace++;  
    }    
}

void updateBeliefDegree()
{
    int update = 0;
    float sumAntAttr1 = 1;
    float sumAntAttr2 = 1;    
    
    if ((H1 + M1 + L1) < 1)
    {
        sumAntAttr1 = H1 + M1 + L1;
        update = 1;
    }
    
    if ((H2 + M2 + L2) < 1)
    {
        sumAntAttr2 = H2 + M2 + L2;
        update = 1;
    }
    
    if (update == 1)
    {
        beliefDegreeChangeLevel = (sumAntAttr1 + sumAntAttr2)/numberOfAntAttributes;
        //cout << "Belief Degree Change Level = " << beliefDegreeChangeLevel << endl;
        for (int go = 0; go < 27; go++)
        {
            consequentBeliefDegree[go] = beliefDegreeChangeLevel * consequentBeliefDegree[go];
            //cout << "Updated Consequent Belief Degree : " << consequentBeliefDegree[go] << endl;
        }
    }
    else
    {
        //cout << endl << "No upgradation of belief degree required." << endl;
    }
    
}


void takeCnnOutput() 
{ 
         
  ifstream myfile ("cnn_prediction.txt"); //cnn output 
  //ifstream myfile ("cnn_prediction1.txt"); //severe 408          
  //ifstream myfile ("cnn_prediction2.txt"); //nominal 36   
  //ifstream myfile ("cnn_prediction3.txt"); //mild 117        
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      //cout << line << '\n';
      if(counter == 0)
      {
         cnn_severe = line;        
      }
      else if (counter == 1)
      {
          cnn_nominal = line;        
      }
      else
      {
          cnn_mild = line;        
      } 
        
      counter++;     
     
    }
    myfile.close();
  } 
 
  else cout << "Unable to open file";  
  float a = stof(cnn_mild);
  float b = stof(cnn_nominal);
  float c = stof(cnn_severe);  
  
  //cout<<" cloud: "<<cloud<<", dew_point: "<<dew_point<<endl; 
  //cout<<" Cloud: "<<cloud<<", Humidity: "<< humidity << endl;
  //cout<<" cnn_mild: "<< a << ", cnn_nominal: " << b << ", cnn_severe: " << c <<endl; 


  float mild_degree = a/100;
  float nominal_degree = b/100;
  float severe_degree = c/100;
    
  float sum_degree = severe_degree + mild_degree + nominal_degree;
   
  normalized_cnn_severe_degree = severe_degree/sum_degree; 
  normalized_cnn_mild_degree = mild_degree/sum_degree;   
  normalized_cnn_nominal_degree = nominal_degree/sum_degree;

  image_normalized_cnn_severe_degree = normalized_cnn_severe_degree;
  image_normalized_cnn_mild_degree = normalized_cnn_mild_degree;
  image_normalized_cnn_nominal_degree = normalized_cnn_nominal_degree;

  cout<<" Before transformInputCNN(), normalized_cnn_severe_degree: "<< image_normalized_cnn_severe_degree << ", normalized_cnn_mild_degree: " << image_normalized_cnn_mild_degree << ", normalized_cnn_nominal_degree: " << image_normalized_cnn_nominal_degree <<endl;         
    
  if ((normalized_cnn_severe_degree > normalized_cnn_mild_degree) && (normalized_cnn_severe_degree > normalized_cnn_nominal_degree))
  {
      //cnn_pm25 = (150.5 + 349.9*normalized_cnn_severe_degree) + ((150.4*normalized_cnn_mild_degree)/2);  
      //cnn_pm25 = (150.5 + 275.9*normalized_cnn_severe_degree) + ((150.4*normalized_cnn_mild_degree)/2);   
      //cnn_pm25 = (150.5 + 69.5*normalized_cnn_severe_degree);    
      cnn_pm25 = (150.5 + 349.5*normalized_cnn_severe_degree);    
        
      cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;           
  }
  else if ((normalized_cnn_nominal_degree > normalized_cnn_mild_degree) && (normalized_cnn_nominal_degree > normalized_cnn_severe_degree))
  { 
      if (normalized_cnn_nominal_degree == 1){ 
        cnn_pm25 = 10;      
      }
      else if (normalized_cnn_nominal_degree >= 0.95){ 
      	cnn_pm25 = 30;
      }
      else if (normalized_cnn_nominal_degree >= 0.90){
      	cnn_pm25 = 40;  
      }    
      else{
        cnn_pm25 = (35.4*(1 - normalized_cnn_nominal_degree)) + ((150.4*normalized_cnn_mild_degree)/2);
      }            
      cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;          
  }
  else if ((normalized_cnn_mild_degree > normalized_cnn_severe_degree) && (normalized_cnn_mild_degree > normalized_cnn_nominal_degree))
  { 
      if (normalized_cnn_mild_degree == 1){  
      	cnn_pm25 = 92; 
      }

      else if (normalized_cnn_severe_degree > normalized_cnn_nominal_degree)
      { 
        cnn_pm25 = (35.5 + 114.9*normalized_cnn_mild_degree) + ((500.4*normalized_cnn_severe_degree)/2);
        cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;          
      }
      
      else if ((normalized_cnn_nominal_degree > normalized_cnn_severe_degree))
      {
        cnn_pm25 = (35.5 + 114.9*normalized_cnn_mild_degree) + ((35.4*normalized_cnn_nominal_degree)/2);    
        cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;           
      }  
  }     
  
/*  if((d >= 0.25) && (e > 0.55)){     
	  if(cnn_pm25 >= 55.50){ 
   		if ((d <= 0.50)){ 
       		PO = (cnn_pm25)/500.4;  
       		cnn_pm25 = 55.4 * PO * 0.90;           
   		}
   		else if ((d > 0.50) && (d <= 0.75)){
       		PO = (cnn_pm25)/500.4;  
       		cnn_pm25 = 55.4 * PO * 0.80;
   		}
   		else if ((d > 0.75) && (d <= 1.00)){ 
       		PO = (cnn_pm25)/500.4; 
       		cnn_pm25 = 55.4 * PO * 0.70; 
   		}
   	}  
 }*/ 


      /*if((d - e) < 2.5){
        PO = (cnn_pm25)/500.4; 
        cnn_pm25 = 55.4 * PO;
     }*/

  cout << endl << "Revised PM2.5 computed by CNN (due to Cloud): " << cnn_pm25 << " µg/m3 " << endl;  
  transformInputCNN(cnn_pm25);   
  //cout << endl << "CNN-generated Belief Degree for Severe Pollution: " << normalized_cnn_severe_degree << endl;        
  //cout << "CNN-generated Belief Degree for Mild Pollution: " << normalized_cnn_mild_degree << endl;   
  //cout << "CNN-generated Belief Degree for Nominal Pollution: " << normalized_cnn_nominal_degree << endl;       
} 
   
void aggregateER_BrbCnn()
{ 
    int parse = 0;
    int move1 = 0; 
    int move2 = 1; 
    int move3 = 2; 
    int action1 = 0;
    int action2 = 1;
    int action3 = 2;
    
    float part11 = 1;
    float part12 = 1;
    float part13 = 1;
    float part1 = 1;
    float part2 = 1;
    float value = 1;
    float meu = 1;
    
    float numeratorH1 = 1;
    float numeratorH2 = 1;
    float numeratorH = 1;
    float denominatorH1 = 1;
    float denominatorH = 1;
    
    float numeratorM1 = 1;
    float numeratorM = 1;
    
    float numeratorL1 = 1;
    float numeratorL = 1;
    
    float utilityScoreH = 1;
    float utilityScoreM = 0.5;
    float utilityScoreL = 0;
    float crispValue = 1;
    float degreeOfIncompleteness = 1;
    float utilityMax = 1;
    float utilityMin = 1;
    float utilityAvg = 1;       
    
    for (int t = 0; t < 27; t++)         
    {
        parse = t * 3;
        ruleWiseBeliefDegreeSum[t] = consequentBeliefDegree[parse] + consequentBeliefDegree[parse+1] + consequentBeliefDegree[parse+2];
    }   
    
    for (int rule = 0; rule < 27; rule++){
        part11 *= (activationWeight[rule] * consequentBeliefDegree[move1] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move1 += 3;
    }
     
    for (int rule = 0; rule < 27; rule++){
        part12 *= (activationWeight[rule] * consequentBeliefDegree[move2] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move2 += 3;
    }

    for (int rule = 0; rule < 27; rule++){
        part13 *= (activationWeight[rule] * consequentBeliefDegree[move3] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move3 += 3;
    }  
    
    part1 = (part11 + part12 + part13);
    
    for (int rule = 0; rule < 27; rule++){
        part2 *= (1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));           
    }    
    
    value = part1 - part2;
    
    meu = 1/value;

    for (int rule = 0; rule < 27; rule++){
        numeratorH1 *= (activationWeight[rule] * consequentBeliefDegree[action1] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action1 += 3;
    }
    
    for (int rule = 0; rule < 27; rule++){
        numeratorH2 *= (1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
    }    
    
    numeratorH = meu * (numeratorH1 - numeratorH2);
    
    for (int rule = 0; rule < 27; rule++){
        denominatorH1 *= (1 - activationWeight[rule]);        
    }
    
    denominatorH = 1 - (meu * denominatorH1);
    
    aggregatedBeliefDegreeH = (numeratorH/denominatorH);
    //cout << endl << "ER Aggregated Belief Degree for Severe Pollution: " << aggregatedBeliefDegreeH << endl;
    
    for (int rule = 0; rule < 27; rule++){
        numeratorM1 *= (activationWeight[rule] * consequentBeliefDegree[action2] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action2 += 3;     
    }
    
    numeratorM = meu * (numeratorM1 - numeratorH2); 
    aggregatedBeliefDegreeM = (numeratorM/denominatorH); 
    //cout << "ER Aggregated Belief Degree for Mild Pollution: " << aggregatedBeliefDegreeM << endl;
    
    for (int rule = 0; rule < 27; rule++){
        numeratorL1 *= (activationWeight[rule] * consequentBeliefDegree[action3] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action3 += 3;
    }
     
    numeratorL = meu * (numeratorL1 - numeratorH2);
    aggregatedBeliefDegreeL = (numeratorL/denominatorH); 
    //cout << "ER Aggregated Belief Degree for Nominal Pollution: " << aggregatedBeliefDegreeL << endl;    
    
    if ((aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL) == 1){
        crispValue = (aggregatedBeliefDegreeH * utilityScoreH) + (aggregatedBeliefDegreeM * utilityScoreM) + (aggregatedBeliefDegreeL * utilityScoreL);
        //cout << "Crisp or numerical value is: " << crispValue << endl;        
        brbH = aggregatedBeliefDegreeH;
        brbM = aggregatedBeliefDegreeM;
        brbL = aggregatedBeliefDegreeL;    
        
        cout << endl << "BRB-CNN integrated Belief Degree for Hazardous AQI: " << aggregatedBeliefDegreeH << endl;
        cout << "BRB-CNN integrated Belief Degree for Unhealthy AQI: " << aggregatedBeliefDegreeM << endl; 
        cout << "BRB-CNN integrated Belief Degree for Good AQI: " << aggregatedBeliefDegreeL << endl;
        //cout << "brbH: " << brbH << " brbM: " << brbM << " brbL: " << brbL <<endl; 
    }   
 
        
    else{
        
        degreeOfIncompleteness = 1 - (aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        //cout << "Usassigned Degree of Belief: " << degreeOfIncompleteness << endl; 
        
        utilityMax = ((aggregatedBeliefDegreeH + degreeOfIncompleteness) * utilityScoreH + (aggregatedBeliefDegreeM*utilityScoreM) + (aggregatedBeliefDegreeL*utilityScoreL));
        
        utilityMin = (aggregatedBeliefDegreeH*utilityScoreH) + (aggregatedBeliefDegreeM*utilityScoreM) + (aggregatedBeliefDegreeL + degreeOfIncompleteness) * utilityScoreL;
        
        utilityAvg = (utilityMax + utilityMin)/2;
        
        //cout << "Maximum expected utility: " << utilityMax << endl;
        //cout << "Minimum expected utility: " << utilityMin << endl; 
        //cout << "Average expected utility: " << utilityAvg << endl; 
        
        cout << endl << "BRB-CNN integrated Belief Degrees considering degree of Incompleteness:" << endl;   
        
        finalAggregatedBeliefDegreeH = aggregatedBeliefDegreeH/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
         
        finalAggregatedBeliefDegreeM = aggregatedBeliefDegreeM/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        
        finalAggregatedBeliefDegreeL = aggregatedBeliefDegreeL/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
                
        
        cout << endl << "BRB-CNN integrated Belief Degree for High PM2.5: " << finalAggregatedBeliefDegreeH << endl; 
        cout << "BRB-CNN integrated Belief Degree for Medium PM2.5: " << finalAggregatedBeliefDegreeM << endl; 
        cout << "BRB-CNN integrated Belief Degree for Low PM2.5: " << finalAggregatedBeliefDegreeL << endl << endl; 
         
        brbH = finalAggregatedBeliefDegreeH;
        brbM = finalAggregatedBeliefDegreeM;
        brbL = finalAggregatedBeliefDegreeL;  
          
        //cout << "brbH: " << brbH << " brbM: " << brbM << " brbL: " << brbL <<endl;
        if ((finalAggregatedBeliefDegreeH > finalAggregatedBeliefDegreeM) && (finalAggregatedBeliefDegreeH > finalAggregatedBeliefDegreeL))
        { 
            //aqi = (201 + 299*finalAggregatedBeliefDegreeH) + ((200*finalAggregatedBeliefDegreeM)/2);
            /*if(finalAggregatedBeliefDegreeH >= 0.80){  
                aqi = 301 + (199 * finalAggregatedBeliefDegreeH);
            } 
            else if(finalAggregatedBeliefDegreeH<0.80){  
                aqi = 201 + (99 * (finalAggregatedBeliefDegreeH));  
            }*/
            if(finalAggregatedBeliefDegreeH >= 0.90){  
                //aqi = 230 + (40 * finalAggregatedBeliefDegreeH); 
              pm25_revised = 180 + (40 * finalAggregatedBeliefDegreeH); 
            }             
            else if((finalAggregatedBeliefDegreeH >= 0.80) && (finalAggregatedBeliefDegreeH < 0.90)){  
                //aqi = 230 + (40 * finalAggregatedBeliefDegreeH); 
              pm25_revised = 120 + (25 * finalAggregatedBeliefDegreeH); 
            } 
            else if(finalAggregatedBeliefDegreeH<0.80){
                //aqi = 201 + (29 * (finalAggregatedBeliefDegreeH));  
              pm25_revised = 100.5 + (25 * (finalAggregatedBeliefDegreeH)); 
            }                                     
            cout << endl << "PM2.5 (revised) predicted by BRB-CNN: " << pm25_revised << endl;           
        }   
        else if ((finalAggregatedBeliefDegreeL > finalAggregatedBeliefDegreeM) && (finalAggregatedBeliefDegreeL > finalAggregatedBeliefDegreeH))   
        { 
            if(finalAggregatedBeliefDegreeL >= 0.98){
                //aqi = (50 * finalAggregatedBeliefDegreeL)/2;
                pm25_revised = (12 * finalAggregatedBeliefDegreeL)/2; 
                cout << " image_normalized_cnn_nominal_degree: " << image_normalized_cnn_nominal_degree << ", image_normalized_cnn_mild_degree: " << image_normalized_cnn_mild_degree << endl; 
            }   
            else if((finalAggregatedBeliefDegreeL < 0.98) && (finalAggregatedBeliefDegreeL >= 0.95)){
            	//aqi = (70 * finalAggregatedBeliefDegreeL);     
            	//pm25_revised = (22 * finalAggregatedBeliefDegreeL);      
            	pm25_revised = (10 * finalAggregatedBeliefDegreeL);      
            }     
            else if((finalAggregatedBeliefDegreeL<0.95) && (finalAggregatedBeliefDegreeL >= 0.90)){
                //aqi = 50 + (50 * (finalAggregatedBeliefDegreeL)); 
                //aqi = 40 + (50 * (1 - finalAggregatedBeliefDegreeL));  
                if((image_normalized_cnn_nominal_degree <= 0.85) && (image_normalized_cnn_mild_degree >= 0.15))
                { 
                	cout << " normalized_cnn_nominal_degree: " << image_normalized_cnn_nominal_degree << ", normalized_cnn_mild_degree: " << image_normalized_cnn_mild_degree << endl; 
                	//aqi = 190 + (49 * finalAggregatedBeliefDegreeM);
                	pm25_revised = 136 + (94.9 * finalAggregatedBeliefDegreeM);
                	/*float temp = 1;                 		
                	temp = finalAggregatedBeliefDegreeL;	
                	finalAggregatedBeliefDegreeL = finalAggregatedBeliefDegreeM; 
                	finalAggregatedBeliefDegreeM = temp;*/
                } 
                else if (finalAggregatedBeliefDegreeM >= 0.08)
                {	
                	cout << " normalized_cnn_nominal_degree: " << image_normalized_cnn_nominal_degree << ", normalized_cnn_mild_degree: " << image_normalized_cnn_mild_degree << endl; 
	                //aqi = 35 + (50 * (1 - finalAggregatedBeliefDegreeM));                 	
	                //pm25_revised = 8 + (23.3 * (1 - finalAggregatedBeliefDegreeM));                 	
	                pm25_revised = 8 + (23.3 * (finalAggregatedBeliefDegreeM));                 	
                }
                else 
                {
                	//aqi = 40 + (50 * (1 - finalAggregatedBeliefDegreeL));  
                	pm25_revised = 10 + (23.3 * (1 - finalAggregatedBeliefDegreeL)); //corrected till this from bottom to up 15:05 10.8.2021
                }


            }
            else if((finalAggregatedBeliefDegreeL<0.90) && (finalAggregatedBeliefDegreeL >= 0.85)){
               
                //aqi = 80 + (50 * (1 - finalAggregatedBeliefDegreeL)); 
                pm25_revised = 26 + (19.9 * (1 - finalAggregatedBeliefDegreeL));
            }
            else if((finalAggregatedBeliefDegreeL<0.85) && (finalAggregatedBeliefDegreeL >= 0.80)){
               
                //aqi = 85 + (50 * (1 - finalAggregatedBeliefDegreeL)); 
                pm25_revised = 28 + (19.9 * (1 - finalAggregatedBeliefDegreeL)); 
            }
            else if((finalAggregatedBeliefDegreeL<0.80) && (finalAggregatedBeliefDegreeL >= 0.75)){
               
                //aqi = 90 + (50 * (1 - finalAggregatedBeliefDegreeL)); 
                pm25_revised = 31 + (19.9 * (1 - finalAggregatedBeliefDegreeL)); 
            }
            else if((finalAggregatedBeliefDegreeL<0.75) && (finalAggregatedBeliefDegreeL >= 0.70)){
               
                //aqi = 95 + (50 * (1 - finalAggregatedBeliefDegreeL)); 
                pm25_revised = 33 + (19.9 * (1 - finalAggregatedBeliefDegreeL)); 
            }
            else if(finalAggregatedBeliefDegreeL<0.70){ 
                //aqi = 101 + (60 * (1-finalAggregatedBeliefDegreeL));
                pm25_revised = 35.5 + (24 * (1-finalAggregatedBeliefDegreeL));
            }            
            //aqi = (100*(1 - finalAggregatedBeliefDegreeL)) + ((200*finalAggregatedBeliefDegreeM)/2);             
            //cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;          
            cout << endl << "Revised PM2.5 predicted by BRB-CNN: " << pm25_revised << endl;          
        }
        else if ((finalAggregatedBeliefDegreeM > finalAggregatedBeliefDegreeH) && (finalAggregatedBeliefDegreeM > finalAggregatedBeliefDegreeL))
        {
            if (finalAggregatedBeliefDegreeH > finalAggregatedBeliefDegreeL)
            { 
                //aqi = (151 + 49*finalAggregatedBeliefDegreeM);  
                pm25_revised = (55.5 + 94.9*finalAggregatedBeliefDegreeM);  
                cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;            
            }
      
            else if ((finalAggregatedBeliefDegreeL > finalAggregatedBeliefDegreeH))
            { 	
            	if(finalAggregatedBeliefDegreeM >= 0.95)
            	{
            		//aqi = (51 + 49*(1 - finalAggregatedBeliefDegreeM));	
            		pm25_revised = (12.1 + 23.3*(1 - finalAggregatedBeliefDegreeM));	
            	}
            	else if ((finalAggregatedBeliefDegreeM >= 0.90) && (finalAggregatedBeliefDegreeM < 0.95))
            	{
            		//aqi = (51 + 49*(finalAggregatedBeliefDegreeM));		
            		pm25_revised = (12.1 + 23.3*(finalAggregatedBeliefDegreeM));		 
            	}
            	else if ((finalAggregatedBeliefDegreeM >= 0.85) && (finalAggregatedBeliefDegreeM < 0.90))
            	{
            		//aqi = (64 + 49*(finalAggregatedBeliefDegreeM)); 		
            		pm25_revised = (19 + 23.3*(finalAggregatedBeliefDegreeM));		 
            	}            	
            	else if ((finalAggregatedBeliefDegreeM >= 0.80) && (finalAggregatedBeliefDegreeM < 0.85))
            	{
            		//aqi = (76 + 49*(finalAggregatedBeliefDegreeM));		
            		pm25_revised = (24 + 23.3*(finalAggregatedBeliefDegreeM));		
            	}           
            	else if ((finalAggregatedBeliefDegreeM >= 0.75) && (finalAggregatedBeliefDegreeM < 0.80))
            	{
            		//aqi = (89 + 49*(finalAggregatedBeliefDegreeM));		 
            		pm25_revised = (31 + 19.9*(finalAggregatedBeliefDegreeM));		 
            	}
            	else
            	{
            		//aqi = (101 + 49*finalAggregatedBeliefDegreeM);  
            		pm25_revised = (35.5 + 19.9*finalAggregatedBeliefDegreeM); 
            	}            	 	            	
                	
                //aqi = (51 + 49*finalAggregatedBeliefDegreeM);
                //cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;           
                cout << endl << "Revised PM2.5 predicted by BRB-CNN: " << pm25_revised << endl;           
            }
        }
        
        /*if(aqi >= 301)
        {
            aqi6 = (aqi- 301)/199.0;  
        }
        else if((aqi >= 201)&& (aqi <= 300))
        {
            aqi6 = (aqi- 201)/99.0; 
        } 
        else if((aqi >= 151)&& (aqi <= 200)) 
        {
            aqi6 = (aqi- 151)/49.0;  
        }
        else if((aqi >= 101)&& (aqi <= 150))  
        {
            aqi6 = (aqi- 101)/49.0; 
        }
        else if((aqi >= 51)&& (aqi <= 100))  
        {
            aqi6 = (aqi- 51)/49.0; 
        }
        else if(aqi <= 50)    
        {
            aqi6 = (aqi/49.0);  
        }   
        //cout << "aqi6: " << aqi6 << endl;       
        Cloud 0.90, RH 0.54       
        cout << endl << "BRB-CNN integrated Belief Degree for Hazardous AQI: " << finalAggregatedBeliefDegreeH*aqi6 << endl;
        cout << "BRB-CNN integrated Belief Degree for Very Unhealthy AQI: " << finalAggregatedBeliefDegreeH*(1-aqi6) << endl;
        cout << "BRB-CNN integrated Belief Degree for Unhealthy AQI: " << finalAggregatedBeliefDegreeM*aqi6 << endl;    
        cout << "BRB-CNN integrated Belief Degree for Unhealthy (Sensitive Groups) AQI: " << finalAggregatedBeliefDegreeM*(1-aqi6) << endl; 
        cout << "BRB-CNN integrated Belief Degree for Moderate AQI: " << finalAggregatedBeliefDegreeL*aqi6 << endl;
        cout << "BRB-CNN integrated Belief Degree for Good AQI: " << finalAggregatedBeliefDegreeL*(1-aqi6) << endl << endl; */ 
         
    }     
}

int main()
{
    ruleBase();
    takeInput();   
    takeCnnOutput(); 
    showTransformedInput(); 
    calculateMatchingDegreeBrbCnn();
    showMatchingDegree();
    showActivationWeight();
    //updateBeliefDegree(); 
    aggregateER_BrbCnn(); 
      
    return 0; 
}
