# sambinder
Prediction of SAM binding sites

SAMbinder is a Python-based tool for predicting SAM interacting residue in a protein chain. It comprises of number of machine learning models for example, SVC model, Random Forest model, Artificial Neural Network model, which is implemented using Scikit package. These models are developed using widely used features like (i) Binary Profile of patterns and (ii) Evolutionary Information in the form of PSSM matrix generated using PSI-BLAST.
Results are produced in the form of propensity score in between value 0-9. Residues showing the propensity score equal or above the selected threshold are said to be “Interacting” whereas residues showing value lesser than selected threshold are said to be “Non-Interacting”. Prediction model developed using evolutionary information where SVC machine learning technique was implemented performed best in our study. 

Installation

Command for downloading SAMbinder
```
git clone https://github.com/raghavagps/sambinder
```

SAMbinder is an open source Python based software, which operates depending on the Python environment (Python version 3.3 or above) and can be run on multi-OS systems (such as Windows, Linux and Mac operating systems). Before running SAMbinder, user should make sure of all the following packages are installed in their Python environment: sys, os, shutil, scipy, numpy(), pandas(), sklearn version 0.19.1, math and re. For convenience, we strongly recommend users to install the Anaconda Python 3.3 (or above) in their local machine. The software can be freely downloadable from https://www.anaconda.com/download/ .
User also needs to download the model directories and blastpr floder in order to run the prediction. Please run the commands given below to download the model folders and blastpr folder

**COMMANDS**
```
wget -c http://webs.iiitd.edu.in/gpsr2/blastpr.zip 
wget -c http://webs.iiitd.edu.in/gpsr2/sam_models.zip  
```

After downloading these folders, untar these folders in the present directory using the command
```
unzip blastpr.zip
unzip sam_models.zip
```
For users who want to do prediction by using our SAMbinder package:
cd to the “SAMbinder” folder which contains SAMbinder.py. This is the main code which is used to run machine learning and do prediction of SAM interacting residue in a given target protein. It allows users to do prediction using 5 different machine learning models. For more information, run:

```
python3 SAMbinder.py -h
```

Examples for users to do SAM interacting residue prediction.
The input protein sequence for SAMbinder.py should be in fasta format. Please find the example in example folder. The following parameters are required by SAMbinder.py

**COMMAND**
```
python3 sambinder.py -i <input_file> -o <output_file> -m <method> -t <thres>
```
where,
-	<input_file>: Input file having sequence file in FASTA format
-	<output_file>: Output file generated by SAMbinder having prediction result
- <thresh>: User defined threshold score (between 0-9)
- <method>: Machine Learning method and the type of input feature it used
  
1. Binary SVC
2. Binary Random Forest
3. Binary MLP
4. Binary KNN
4. PSSM SVC
5. PSSM Random Forest


For more information type,
```
python3 sambinder.py –h
```

In our package, we have provided 5 different machine learning models which utilizes different features. 
- Method ‘1’ is Support Vector Classifier which utilizes binary profile of the pattern as an input feature.
- Method ‘2’ is Random Forest Classifier which also utilizes binary profile of the pattern as an input feature. 
- Method ‘3’ is Artificial Neural Network model developed using binary profile of the pattern as input feature. 
- Method ‘4’ is K Nearest Neighbor method developed using binary profile of the pattern as input feature.
- Method ‘5’ is Support Vector Classifier which utilizes evolutionary information in the form of PSSM profile as an input feature.
- Method ‘6’ is Random Forest classifier which also utilizes evolutionary information in the form of PSSM profile as an input feature. The PSSM profile is generated using PSI-BLAST by running against the SwissProt database.
