# Smile_Generator
This is a code to systematically generate SMILES of wanted molecules

1. Molecule and its stereoisomer generate.
    ```
    cd run
    ./TaskGenerate.py -n 20
    ./GetMoleculeAndStereoIsomer.py
    ```
2. Similarity analysis based on morgan fingerprint.  
   Get the training set and validation set with given similarity cutoff value.
    ```
    cd scripts
    ./SimilarityAnalysis.py
    ./TrainValidationSetGenerate.py -t FP-similarity --similarity 0.1
    ```
3. Get the training set and validation set using FP-space distance
   ```
   cd run
   ./FPGenerate.py --maxPath 7 --cutoff2 10.0
   cd ../scripts
   ./TrainValidationSetGenerate.py -t FP-distance
   ```
4. Topological fingerprint analysis.
   Get the training set and validation set based on topological fp space distance.
   ```
   cd run
   ./FPGenerate.py
   ./TrainValidationSetGenerate.py -t FP-distance
   ```
