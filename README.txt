1. Run demo1.m for a short and quick demonstration:
   (a) It applies three different FRAs to problems with at most 100 variables;
   (b) It applies FRA and then solve Shor's relaxation with/without 
   preprocessing for problems with at most 100 variables;

2. Run demo2.m to reproduce the results in the paper:
   (a) It applies three different FRAs to problems with at most 10000 variables;
   (b) It applies FRA and then solve Shor's relaxation with/without 
   preprocessing for problems with at most 350 variables;
   (c) It takes many hours to complete demo2.m.
       For your convinience, we run demo2.m at our local machine, and save
       the results in file "demo2SDPoutput.mat" and "demo2FRoutput.mat".
       The terminal output is saved in file "demo2output.txt".
   (d) The server does not have enough memory for solving SDP
   relaxation with 200~350 variables. It can require up to 128GB RAM.

3. How can I customized the tests?
   BIP_test(maxsize,solveSDP) apply FRAs to all problems at most "maxsize" 
   variables. If solveSDP is true, then we also solve the Shor's SDP 
   relaxation with and without FRA when there exists a positive reduction.

4. Structure of the folder/codes
-autoFRcode_MPC_rev1
-data 	             "the MILPLIB2017 data is saved here"
-YALMIP-master       "YALMIP is needed for SDP test"
demo1.m 	         "a short test"
demo2.m 	         "a complete test"
affineFR.m           "implementation of affine FR algorithm"
partialFR.m          "implementation of partial FR algorithm"
BIP_test.m	         "apply FRA/solve Shor's relaxation"
FRAformat.m          "convert gurobi format to proper format for FRA"
FRAcompute.m         "wrapper for calling one of the FRAs"
SDP_Shor.m           "solve the Shor's SDP relaxation"
SDP_Shor_FR.m        "solve the facially reduced Shor's SDP relaxation"
prob_list.mat        "the instance names sorted by #variables"
demo2FRoutput.mat    "output of demo2.m"
demo2SDPoutput.mat   "output of demo2.m"
demo2output.txt      "terminal output of demo2.m"